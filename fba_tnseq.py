#!/usr/bin/env python
'''
@name fba_tnseq.py
@author: Juan C. Castro <jcb0@cdc.gov>, Rakin Choudhury <rchoudhury32@gatech.edu>
@version: 1.0.0
@update: 23-Apr-2026
@license: GNU General Public License v3.0

Integrate Tn-seq fitness data into FBA via:
  (A) A modified objective that rewards flux through fitness-important reactions
  (B) Reaction bound constraints scaled by gene fitness
 
Expected Tn-seq input (CSV/TSV): columns `gene_id`, `log2fc` (or `fitness`).
Convention used here:
    fitness = -log2fc
    -> high positive fitness = gene important for growth
    -> high negative fitness = gene loss is beneficial / neutral

'''
import cobra
import numpy as np
import pandas as pd
from cobra.flux_analysis import pfba
 
 
# ---------------------------------------------------------------------------
# 1. Load and prep Tn-seq data
# ---------------------------------------------------------------------------
def load_tnseq(path, gene_col="gene_id", lfc_col="log2fc", fitness_col=None):
    """Load Tn-seq data and produce a `fitness` column (higher = more important)."""
    df = pd.read_csv(path, sep=None, engine="python")
    df = df.rename(columns={gene_col: "gene_id"})
    if fitness_col and fitness_col in df.columns:
        df["fitness"] = df[fitness_col]
    else:
        df["fitness"] = -df[lfc_col]
    return df[["gene_id", "fitness"]]
 
 
# ---------------------------------------------------------------------------
# 2. Map gene-level fitness to reaction-level scores via GPR
# ---------------------------------------------------------------------------
def reaction_fitness_scores(model, gene_fitness, missing_value=0.0):
    """
    Convert gene fitness -> reaction fitness score using GPR logic:
      - 'AND' (complex subunits): minimum fitness across genes (weakest link)
      - 'OR'  (isozymes):         maximum fitness across genes (best option)
    Reactions with no GPR or no matching genes get `missing_value`.
    """
    gf = dict(zip(gene_fitness["gene_id"], gene_fitness["fitness"]))
    scores = {}
 
    for rxn in model.reactions:
        gpr = rxn.gene_reaction_rule
        if not gpr:
            scores[rxn.id] = missing_value
            continue
 
        # Parse GPR as OR of ANDs (the standard form cobrapy uses)
        # e.g. "(g1 and g2) or (g3)"
        or_clauses = [c.strip(" ()") for c in gpr.split(" or ")]
        clause_vals = []
        for clause in or_clauses:
            genes = [g.strip() for g in clause.split(" and ")]
            vals = [gf.get(g, np.nan) for g in genes]
            vals = [v for v in vals if not np.isnan(v)]
            if vals:
                clause_vals.append(min(vals))  # AND -> min
        scores[rxn.id] = max(clause_vals) if clause_vals else missing_value  # OR -> max
 
    return pd.Series(scores, name="fitness_score")
 
 
# ---------------------------------------------------------------------------
# 3. (A) Modified objective: biomass + weighted fitness-flux term
# ---------------------------------------------------------------------------
def set_fitness_objective(model, rxn_scores, biomass_weight=1.0,
                          fitness_weight=0.01, only_positive=True):
    """
    Replace objective with:   biomass_weight * v_biomass
                            + fitness_weight * sum_i (score_i * |v_i|)
 
    `fitness_weight` is small by design — biomass should still dominate.
    Treat it as a tie-breaker that favors flux through important reactions
    among otherwise-equivalent optimal solutions. Tune between 1e-3 and 1e-1.
 
    `only_positive=True` ignores negative fitness scores (reactions whose
    genes are dispensable) rather than actively penalizing them.
    """
    biomass_rxn = _find_biomass(model)
    obj_terms = {biomass_rxn.forward_variable: biomass_weight,
                 biomass_rxn.reverse_variable: -biomass_weight}
 
    for rxn_id, score in rxn_scores.items():
        if only_positive and score <= 0:
            continue
        if not only_positive and score == 0:
            continue
        rxn = model.reactions.get_by_id(rxn_id)
        # Reward |flux|: add score to both forward and reverse variables
        obj_terms[rxn.forward_variable] = (
            obj_terms.get(rxn.forward_variable, 0) + fitness_weight * score
        )
        obj_terms[rxn.reverse_variable] = (
            obj_terms.get(rxn.reverse_variable, 0) + fitness_weight * score
        )
 
    model.objective = model.problem.Objective(
        sum(coef * var for var, coef in obj_terms.items()),
        direction="max",
    )
    return model
 
 
# ---------------------------------------------------------------------------
# 4. (B) Fitness-based bound constraints
# ---------------------------------------------------------------------------
def apply_fitness_bounds(model, rxn_scores, mode="soft",
                         low_pct=10, high_pct=90,
                         lower_scale=0.1, upper_scale=1.0):
    """
    Scale reaction bounds based on fitness percentile.
    - Reactions below the low_pct percentile (low fitness) get tightened
      bounds: multiplied by `lower_scale`.
    - Reactions above the high_pct percentile keep full bounds.
    - Middle tier is scaled linearly between.
 
    `mode='soft'` scales existing bounds (multiplicative).
    `mode='hard'` sets unused low-fitness reactions near zero.
    """
    scored = rxn_scores[rxn_scores != 0]  # skip reactions without data
    if scored.empty:
        print("No scored reactions; skipping bound adjustment.")
        return model
 
    low = np.percentile(scored, low_pct)
    high = np.percentile(scored, high_pct)
    print(f"Bound scaling: low_pct={low:.3f}, high_pct={high:.3f}")
 
    for rxn_id, score in scored.items():
        rxn = model.reactions.get_by_id(rxn_id)
        if score >= high:
            scale = upper_scale
        elif score <= low:
            scale = lower_scale if mode == "soft" else 1e-3
        else:
            # Linear interpolation
            frac = (score - low) / (high - low)
            scale = lower_scale + frac * (upper_scale - lower_scale)
 
        lb, ub = rxn.bounds
        rxn.bounds = (lb * scale if lb < 0 else lb,
                      ub * scale if ub > 0 else ub)
    return model
 
 
# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------
def _find_biomass(model):
    """Find the biomass reaction by common naming patterns."""
    candidates = [r for r in model.reactions
                  if "biomass" in r.id.lower() or "biomass" in r.name.lower()]
    if not candidates:
        raise ValueError("No biomass reaction found. Set model.objective manually.")
    # Prefer the one that's currently the objective, else the first match
    for r in candidates:
        if r.objective_coefficient != 0:
            return r
    return candidates[0]
 
 
def compare_solutions(model, rxn_scores, fitness_weight=0.01):
    """Run baseline FBA, pFBA, and fitness-weighted FBA; compare active reactions."""
    with model as m:
        baseline = m.optimize()
    with model as m:
        parsimonious = pfba(m)
    with model as m:
        set_fitness_objective(m, rxn_scores, fitness_weight=fitness_weight)
        fitness_sol = pfba(m)
 
    def active(fluxes, tol=1e-6):
        return set(fluxes.index[fluxes.abs() > tol])
 
    a_base = active(baseline.fluxes)
    a_pfba = active(parsimonious.fluxes)
    a_fit = active(fitness_sol.fluxes)
 
    print(f"\nActive reactions — baseline FBA: {len(a_base)}")
    print(f"Active reactions — pFBA:        {len(a_pfba)}")
    print(f"Active reactions — fitness FBA: {len(a_fit)}")
    print(f"  Newly active under fitness:   {len(a_fit - a_pfba)}")
    print(f"  Dropped vs. pFBA:             {len(a_pfba - a_fit)}")
    return {"baseline": baseline, "pfba": parsimonious, "fitness": fitness_sol}
 
 
# ---------------------------------------------------------------------------
# Example usage
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    model = cobra.io.read_sbml_model("iPae1146.xml")
    tnseq = load_tnseq("tnseq_fitness.csv", lfc_col="log2fc")
    scores = reaction_fitness_scores(model, tnseq)
 
    print(f"Reactions with fitness data: {(scores != 0).sum()} / {len(scores)}")
    print(f"Score range: {scores.min():.2f} to {scores.max():.2f}")
 
    # --- Approach A: modify objective only ---
    with model as m:
        set_fitness_objective(m, scores, fitness_weight=0.01)
        sol_obj = pfba(m)
        print(f"Biomass under fitness objective: "
              f"{sol_obj.fluxes[_find_biomass(m).id]:.4f}")
 
    # --- Approach B: bound constraints only ---
    with model as m:
        apply_fitness_bounds(m, scores, mode="soft")
        sol_bnd = m.optimize()
        print(f"Biomass under fitness bounds: {sol_bnd.objective_value:.4f}")
 
    # --- Compare all three ---
    compare_solutions(model, scores, fitness_weight=0.01)
