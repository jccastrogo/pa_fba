#!/usr/bin/env python
'''
@name pseudomonas_flux.py
@author: Juan C. Castro <jcb0@cdc.gov>
@version: 1.2.0
@update: 24-Jun-2024
@license: GNU General Public License v3.0
'''

'''1.0 Import modules, define funcions, initialize variables'''
#======================1.1 Import modules=====================#
try :
    import os, sys, subprocess
    import argparse
    import json
    import numpy as np
    import pandas as pd
    import cobra
    from cobra.flux_analysis import mapExpressionToReactions
    import matplotlib.pyplot as plt
except :
    sys.stderr.write('ERROR! Cannot import required modules. pseudomonas_flux.py requores: argparse, subprocess, json and cobra')
#===================1.2 Initialize variables==================#
parser = argparse.ArgumentParser(description="pseudomonas_flux.py: Calculate flux balance analysis for a metabolic model.")
group = parser.add_argument_group('Required arguments') 
group.add_argument('-c', action='store', dest='config_file',
                   required=True, help='A configuraton file in JSON format with options on how to run Flux Balance Analysis.')

#=====================1.3 Define functions====================#
def parseConfigFile(config_file:str):
    '''
    Parse a configuration file to retrieve configuration options of the analysis

    Args:
        config_file (str): Path to the configuration file
    
    Returns:
        dict: A dictionary with the configuration options
    '''
    config_IOS = open(config_file, 'r')
    config_dict = json.load(config_IOS)
    config_IOS.close()
    return(config_dict)

def create_sbml_model(reaction_file, compounds_file, output_file, save_to_file):
    """
    Creates an SBML metabolic model from reaction and compounds files.

    Args:
        reaction_file: Path to the reaction file (e.g., tab-separated format).
        compounds_file: Path to the compounds file (e.g., tab-separated format).
        output_file: Path to save the SBML model.

    Returns:
        None (saves the SBML model to the specified output file).
    """
    # Read reactions and compounds from files (you'll need to parse your specific formats)
    # Example: Read reactions and compounds into appropriate data structures

    # Create a cobra.Model object
    model = cobra.Model()
    
    # Load reactions and compounds from CSV files
    reactions_df = pd.read_csv(reaction_file)
    compounds_df = pd.read_csv(compounds_file)
    
    # Add metabolites to the model
    for index, row in compounds_df.iterrows():
        metabolite = cobra.Metabolite(
            id=row['id'],
            name=row['name'],
            compartment=row['compartment']
        )
        model.add_metabolites([metabolite])     
    # Add reactions to the model
    for index, row in reactions_df.iterrows():
        reaction = cobra.Reaction(id=row['id'])
        reaction.name = row['name']
        reaction.lower_bound = row['lower_bound']
        reaction.upper_bound = row['upper_bound']
        # Add metabolites to the reaction
        metabolites = {}
        for met_id, stoich in zip(row['metabolites'].split(';'), row['stoichiometry'].split(';')):
            metabolites[model.metabolites.get_by_id(met_id)] = float(stoich)
        reaction.add_metabolites(metabolites)
        model.add_reactions([reaction])
    # Save the model to SBML format
    if save_to_file:
        cobra.io.write_sbml_model(model, output_file)    
    return(model)

def optimize_biomass(model: cobra.Model,
                     biomass_rxn: str = None) -> cobra.Solution:
    """Standard FBA: maximize biomass."""
    if biomass_rxn is not None:
        model.objective = model.reactions.get_by_id(biomass_rxn)
    solution = model.optimize()
    _report("biomass", model, solution)
    return solution
 
 
def optimize_pfba(model: cobra.Model,
                  biomass_rxn: str = None,
                  fraction_of_optimum: float = 1.0) -> cobra.Solution:
    """
    Parsimonious FBA (Blank et al. 2005; Lewis et al. 2010).
    Maximize biomass, then minimize sum of |fluxes| subject to biomass
    staying at `fraction_of_optimum` of its maximum.
    """
    if biomass_rxn is not None:
        model.objective = model.reactions.get_by_id(biomass_rxn)
    solution = pfba(model, fraction_of_optimum=fraction_of_optimum)
    _report("pfba", model, solution)
    return solution
 
 
def optimize_min_redox(model: cobra.Model,
                       biomass_rxn: str,
                       redox_metabolite_ids=None,
                       biomass_fraction: float = 1.0) -> cobra.Solution:
    """
    Minimize flux through reactions that produce redox cofactors
    (NADH, NADPH, FADH2 by default), while requiring biomass to stay at
    >= `biomass_fraction` of its unconstrained optimum.
 
    Based on Knorr et al. (2007).
 
    Args:
        model: The constrained model. Must have a biomass reaction.
        biomass_rxn: ID of the biomass reaction (used to set the minimum
                     growth constraint).
        redox_metabolite_ids: Iterable of metabolite IDs for redox cofactors
                              in their reduced form. Defaults to ModelSEED
                              NADH/NADPH/FADH2 in the cytosol.
        biomass_fraction: Required fraction of maximum biomass (0.0 - 1.0).
 
    Returns:
        The cobra Solution.
    """
    if redox_metabolite_ids is None:
        redox_metabolite_ids = DEFAULT_REDOX_METABOLITES
 
    # Use a context manager so our modifications don't persist on the model.
    with model:
        # 1) Find the maximum biomass under current media, set a lower bound.
        biomass = model.reactions.get_by_id(biomass_rxn)
        model.objective = biomass
        max_growth = model.slim_optimize()
        if max_growth is None or max_growth != max_growth:  # NaN check
            raise RuntimeError("Biomass optimization failed; cannot set "
                               "growth floor for min-redox.")
        biomass.lower_bound = biomass_fraction * max_growth
        print(f"min_redox: max biomass = {max_growth:.6f}, "
              f"floor set to {biomass.lower_bound:.6f} "
              f"({biomass_fraction*100:.0f}% of optimum)")
 
        # 2) Find all reactions that *produce* a redox metabolite.
        redox_mets = []
        for mid in redox_metabolite_ids:
            try:
                redox_mets.append(model.metabolites.get_by_id(mid))
            except KeyError:
                print(f"  Warning: redox metabolite '{mid}' not in model; "
                      f"skipping.", file=sys.stderr)
 
        redox_rxns = set()
        for met in redox_mets:
            for rxn in met.reactions:
                coeff = rxn.metabolites[met]
                # A reaction "produces" the metabolite when flux * coeff > 0.
                # Forward (flux > 0) produces if coeff > 0; reverse produces
                # if coeff < 0 and the reaction is reversible. We include the
                # reaction either way and let the solver pick the direction.
                if coeff != 0:
                    redox_rxns.add(rxn.id)
 
        if not redox_rxns:
            raise RuntimeError("No redox-producing reactions found; check "
                               "your --redox_metabolites argument.")
        print(f"min_redox: minimizing across {len(redox_rxns)} "
              f"redox-producing reactions")
 
        # 3) Build a linear objective: minimize sum of (forward - reverse)
        #    fluxes for each redox-producing reaction, signed so that
        #    *production* of the redox cofactor is what gets minimized.
        #    For simplicity and robustness we minimize the sum of absolute
        #    fluxes through these reactions using pFBA-style variable split,
        #    which is what COBRApy's linear objective on reaction forward
        #    variables gives when combined with `direction='min'` after
        #    setting a growth floor.
        obj_terms = {model.reactions.get_by_id(rid): 1.0 for rid in redox_rxns}
        model.objective = obj_terms
        model.objective_direction = "min"
 
        solution = model.optimize()
 
    _report("min_redox", model, solution)
    return solution 

def custom_fba(model, objective_coefficients):
    """
    Perform flux balance analysis (FBA) with a customizable objective function.

    Args:
        model (cobra.Model): The metabolic model (read from SBML or created).
        objective_coefficients (dict): A dictionary mapping reaction IDs to their coefficients.
            Example: {'R1': 1.0, 'R2': -0.5} means maximize R1 flux and minimize R2 flux.

    Returns:
        dict: A dictionary containing the flux values for each reaction.
    """
    # Set the custom objective coefficients
    for reaction_id, coefficient in objective_coefficients.items():
        model.reactions.get_by_id(reaction_id).objective_coefficient = coefficient

    # Optimize the model (FBA)
    solution = model.optimize()

    # Get the fluxes
    fluxes = solution.fluxes

    return fluxes 

def knockout_and_fba(model, gene_ids):
    """
    Perform gene knockout and FBA for each gene in the list.

    Args:
        model (cobra.Model): The metabolic model.
        gene_ids (list): List of gene IDs to knock out.

    Returns:
        dict: A dictionary containing flux values for each reaction after knockout.
    """
    fluxes_after_knockout = {}

    for gene_id in gene_ids:
        # Knock out the gene (set its bounds to zero)
        model.genes.get_by_id(gene_id).knock_out()

        # Calculate fluxes using custom FBA
        fluxes = custom_fba(model, {'biomass_reaction': 1.0})  # Example: maximize biomass production

        # Store fluxes after knockout
        fluxes_after_knockout[gene_id] = fluxes

        # Reset the gene knockout (set bounds back to original)
        model.genes.get_by_id(gene_id).knock_in()

    return fluxes_after_knockout

def read_media_file(media_file_path):
    """
    Read the media file and extract relevant data.

    Args:
        media_file_path (str): Path to the media file (tab-separated format).

    Returns:
        dict: A dictionary containing media data for each compound.
            Example: {'compound1': {'minflux': 0.1, 'maxflux': 10.0, 'concentration': 1.0}, ...}
    """
    media_data = {}
    with open(media_file_path, 'r') as f:
        for line in f:
            compound, _, _, minflux, maxflux, concentration = line.strip().split('\t')
            media_data[compound] = {
                'minflux': float(minflux),
                'maxflux': float(maxflux),
                'concentration': float(concentration)
            }
    return media_data

def fba_with_media(model, media_file_path):
    """
    Perform flux balance analysis (FBA) considering media constraints.

    Args:
        model (cobra.Model): The metabolic model.
        media_file_path (str): Path to the media file (tab-separated format).

    Returns:
        dict: A dictionary containing flux values for each reaction after considering media constraints.
    """
    # Read the media file
    media_data = read_media_file(media_file_path)

    # Set bounds for exchange reactions based on media data
    for compound, data in media_data.items():
        model.reactions.get_by_id(f"{compound}_exchange").lower_bound = data['minflux']
        model.reactions.get_by_id(f"{compound}_exchange").upper_bound = data['maxflux']

    # Optimize the model (FBA)
    solution = model.optimize()

    # Get the fluxes
    fluxes = solution.fluxes

    return fluxes
def read_gene_expression_file(file_path):
    """
    Read gene expression data from a tab-delimited file.

    Args:
        file_path (str): Path to the gene expression file.

    Returns:
        dict: A dictionary mapping gene IDs to expression values.
    """
    gene_expression_data = {}
    try:
        with open(file_path, "r") as file:
            for line in file:
                gene_id, *expression_values = line.strip().split("\t")
                if gene_id:
                    # Assuming the first value after gene ID is the expression value
                    expression_value = float(expression_values[0])
                    gene_expression_data[gene_id] = expression_value
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        # You can handle the exception as needed (e.g., return an empty dictionary)

    return gene_expression_data

def relate_expression_to_fluxes(model, expression_values):
    """
    Relates gene expression values to fluxes in a metabolic model.

    Args:
        model: A cobra.Model object representing the metabolic model.
        expression_values: A dictionary mapping gene IDs to expression values.

    Returns:
        A dictionary mapping reaction IDs to their corresponding fluxes.
    """
    fluxes = {}
    for reaction in model.reactions:
        gpr = reaction.gene_reaction_rule
        if gpr:
            # Calculate the reaction flux based on expression values
            # (You can customize this part based on your specific approach)
            flux = expression_values.get(gpr, 0.0)
            fluxes[reaction.id] = flux
    return fluxes

def map_expression_to_reactions(model, expression_values):
    """
    Maps gene expression values to reactions in a metabolic model using cobrapy's mapExpressionToReactions.

    Args:
        model: A cobra.Model object representing the metabolic model.
        expression_values: A dictionary mapping gene IDs to expression values.

    Returns:
        A dictionary mapping reaction IDs to their corresponding expression values.
    """
    reaction_expression = mapExpressionToReactions(model, expression_values)
    return reaction_expression

def map_genes_to_reactions(model, expression_values, objective_function):
    """
    Relates gene expression values to genes in a metabolic model based on an objective function.

    Args:
        model: A cobra.Model object representing the metabolic model.
        expression_values: A dictionary mapping gene IDs to expression values.
        objective_function: A function that defines the optimization objective (e.g., maximize growth rate).

    Returns:
        A dictionary mapping gene IDs to predicted expression values.
    """
    # Solve the optimization problem to find optimal flux distribution
    model.objective = objective_function
    solution = model.optimize()

    # Extract fluxes and map to gene expression
    fluxes = {reaction.id: solution.fluxes[reaction.id] for reaction in model.reactions}
    predicted_expression = {}
    for reaction in model.reactions:
        gpr = reaction.gene_reaction_rule
        if gpr:
            # Calculate predicted expression based on flux
            gene_expression = np.sum([expression_values.get(gene_id, 0.0) for gene_id in gpr.split(' or ')])
            predicted_expression[gpr] = max(0.0, gene_expression)  # Ensure non-negative values

    return predicted_expression

def create_custom_objective(model, gene_expression_file):
    """
    Create a custom objective function based on gene expression coefficients.

    Args:
        model (cobra.Model): The metabolic model.
        gene_expression_file (str): Path to the gene expression file.

    Returns:
        dict: A dictionary mapping reaction IDs to their coefficients.
    """
    # Read gene expression data from the file (you'll need to implement this part)
    gene_expression_data = read_gene_expression_file(gene_expression_file)

    # Map gene IDs to reaction IDs (you'll need to implement this mapping)
    gene_to_reaction_mapping = map_genes_to_reactions(model, gene_expression_data)

    # Calculate reaction coefficients based on gene expression
    reaction_coefficients = {}
    for gene_id, expression_value in gene_expression_data.items():
        for reaction_id in gene_to_reaction_mapping.get(gene_id, []):
            # Assign a coefficient (e.g., scaled expression value) to the reaction
            reaction_coefficients[reaction_id] = expression_value

    return reaction_coefficients


def plot_observed_vs_predicted(observed_expression, predicted_expression):
    """
    Creates a scatter plot of observed vs. predicted gene expression values.

    Args:
        observed_expression: A dictionary mapping gene IDs to observed expression values.
        predicted_expression: A dictionary mapping gene IDs to predicted expression values.

    Returns:
        None (displays the plot).
    """
    gene_ids = list(observed_expression.keys())
    observed_values = list(observed_expression.values())
    predicted_values = [predicted_expression.get(gene_id, 0.0) for gene_id in gene_ids]

    plt.figure(figsize=(8, 6))
    plt.scatter(observed_values, predicted_values, color='b', alpha=0.7)
    plt.plot([min(observed_values), max(observed_values)], [min(observed_values), max(observed_values)], 'r--')
    plt.xlabel('Observed Expression')
    plt.ylabel('Predicted Expression')
    plt.title('Observed vs. Predicted Gene Expression')
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    #Parse arguments
    args = parser.parse_args()
    #Parse config file
    config_file = args.config_file
    config_dict = parseConfigFile(config_file)
    #Data options
    data_opts = config_dict['data_load']
    rxn_file = data_opts["rxn_file"]
    cpd_file = data_opts["cpd_file"]
    model_file = data_opts["model_file"]
    expr_file = data_opts["expr_file"]
    tn_seq_file = data_opts["tn_seq_file"]
    media_file = data_opts['media_file']
    #Analysis options
    analysis_opts = config_dict["data_analysis"]
    export_sbml = analysis_opts["export_sbml"]
    ko_fba = analysis_opts["ko_fba"]
    epr_to_flux = analysis_opts["epr_to_flux"]
    custom_obj = analysis_opts["custom_obj"]
    predict_expr = analysis_opts["predict_expr"]
    #Output options
    output_opts = config_dict["data_analysis"]
    output_prefix = output_opts["output_prefix"]
        
    if export_sbml:
        output_file = f'{output_prefix}gsm.sbml'
        save_to_file = True
        model = create_sbml_model(rxn_file, cpd_file, output_file, save_to_file)
    else:
        save_to_file = False
        output_file = f'{output_prefix}gsm.sbml'
        model = create_sbml_model(rxn_file, cpd_file, output_file, save_to_file)
    
    if json.loads(ko_fba.lower()):
        gene_ids = pd.read_csv(tn_seq_file)
        knockout_and_fba(model, gene_ids)
            
    media_fluxes = fba_with_media(model, media_file)
    expression_values = read_gene_expression_file(expr_file)
    relate_expression_to_fluxes(model, expression_values)

    tn_seq_opt = create_custom_objective(model, expr_file)
    gene_ids = pd.read_csv(tn_seq_file)
    custom_fba(model, gene_ids)
