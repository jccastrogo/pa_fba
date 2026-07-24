#!/usr/bin/env python
import argparse
import pandas as pd
import cobra
import csv
from sympy .logic.boolalg import to_dnf


def load_xml_model(xml_model: str):
    """
    Loads XML model and validates it. Catches any errors in an error log.

    Args: 
        xml_model: Path to XML file

    Returns:
        XML model as a cobra model object
    """

    model, errors = cobra.io.sbml.validate_sbml_model(xml_model)

    with open('model_validation_errors.tsv', 'w') as error_file:
        for error in errors:
            for message in errors[error]:
                error_file.write(f"{error}\t{message}")

    print(len(model.reactions), "reactions in the model")
    print(len(model.genes), "genes in the model")
    print(len(model.metabolites), "metabolites in the model")

    return model


def incorporate_media(model: cobra.Model, media_file_path: str):
    """
    Incorporates media constraints into the metabolic model by reading a media file and adjusting the bounds of exchange reactions.

    Args:
        model (cobra.Model): The metabolic model.
        media_file_path (str): Path to the media file (tab-separated format).

    Returns:
        model: The metabolic model with updated bounds for exchange reactions.
    """

    # Read the media file and store the data in a dictionary
    media_data = {}
    with open(media_file_path, 'r') as media_file:
        lines = media_file.readlines()
        for line in lines[1:]:  # Skip header line
            compound, _, _, minflux, maxflux, concentration = line.strip().split('\t')
            media_data[compound] = {
                'minflux': float(minflux),
                'maxflux': float(maxflux),
                'concentration': float(concentration)
            }
    # Sets new bounds for exchange reactions based on media data
    for compound, data in media_data.items():
        model.reactions.get_by_id(compound).lower_bound = data['minflux']
        model.reactions.get_by_id(compound).upper_bound = data['maxflux']

    return model


def incorporate_media_alternative(model: cobra.Model, media_file_path: str):
    """
    Alternative (possibly correct?) method to incorporate media constraints into the metabolic model.
    Utilizes media object from cobrapy to replace the original media associated with the model.

    Args:
        model (cobra.Model): The metabolic model.
        media_file_path (str): Path to the media file (tab-separated format).

    Returns:
        model: The metabolic model with updated bounds for exchange reactions.
    """

    # Read the media file and store the data in a dictionary
    scfm2_media_data = {}
    with open(media_file_path, 'r') as media_file:
        lines = media_file.readlines()
        for line in lines[1:]:  # Skip header line
            compound, _, _, _, maxflux, _ = line.strip().split('\t')
            scfm2_media_data[compound] = float(maxflux)
    
    # Replace the original media associated with the model with the new media object
    model.media = scfm2_media_data

    return model


def knockut_gene_essentiality_test(model: cobra.Model):
    """
    Tests gene essentiality by performing single gene knockouts and optimizing for biomass production.

    Args:
        model: A cobra.Model object representing the metabolic model.
    
    Returns:
        A dictionary mapping genes to a true/false indicating essentiality
    """

    gene_essentiality_dict = {}

    gene_ids = [gene.id for gene in model.genes]
    for gene_id in gene_ids:
        knockout_model = model.copy()
        knockout_model.objective = knockout_model.reactions.get_by_id('PAO1_Biomass')
        knockout_model.genes.get_by_id(gene_id).knock_out()
        solution = knockout_model.optimize(objective_sense='maximize')
       
       # Assigns 'True' if the objective value is 0, negative, or extremely close to 0 (less than 1e-6) 
       # False otherwise
        if solution.objective_value <= 0:
            gene_essentiality_dict[gene_id] = True
        elif 0 < solution.objective_value < (1 * 10**-6):
            gene_essentiality_dict[gene_id] = True
        else:
            gene_essentiality_dict[gene_id] = False
    
    return gene_essentiality_dict


def optimize_model_biomass(model: cobra.Model):
    """
    Optimizes metabolic model for biomass production.

    Args:
        model: A cobra.Model object representing the metabolic model.

    Returns:
        solution: The optimization solution object.
    """
    
    # Set the objective of the model to maximize biomass production
    model.objective = model.reactions.get_by_id('PAO1_Biomass')

    # Optimize the model
    biomass_solution = model.optimize(objective_sense='maximize')

    fluxes = [biomass_solution.fluxes[reaction.id] for reaction in model.reactions if biomass_solution.fluxes[reaction.id] != 0.0]
    print("Number of reactions with fluxes for biomass optimization != 0:", len(fluxes))

    return biomass_solution


def optimize_model_min_flux(model: cobra.Model):
    """
    Optimizes metabolic model for minimizing overall flux (parsimonious FBA).

    Args:
        model: A cobra.Model object representing the metabolic model.

    Returns:
        solution: The optimization solution object.
    """

    # optimize the model using parsimonious FBA to minimize overall flux while maintaining the same biomass optimization objectve
    pfba_solution = cobra.flux_analysis.parsimonious.pfba(model)

    fluxes = [pfba_solution.fluxes[reaction.id] for reaction in model.reactions if pfba_solution.fluxes[reaction.id] != 0.0]
    print("Number of reactions with fluxes for parsimonious FBA != 0:", len(fluxes))

    return pfba_solution


def optimize_model_min_redox_flux(model: cobra.Model):
    """
    Optimizes metabolic model by minimizing redox flux (NADH/NADPH/FADH2).
    NADH, NADPH, and FADH2; cpd00004_c, cpd00005_c, and cpd00982_c, respectively.

    Args:
        model: A cobra.Model object representing the metabolic model.
    Returns:
        solution: The optimization solution object.
    """

    # listed below (in a set) are the compound ids for NADH, NADPH, and FADH2
    redox_compounds = set(['cpd00004_c', 'cpd00005_c', 'cpd00982_c'])

    # initialize list entry to store reactions producing redox compounds
    redox_reactions = []

    # identifies reactions that produce the above redox compounds and adds them to the redox_reactions list
    for reaction in model.reactions:
        product_ids = set([metabolite.id for metabolite in reaction.products])
        if len(redox_compounds.intersection(product_ids)) > 0:
            redox_reactions.append(reaction.id)
    
    print("Number of redox reactions:", len(redox_reactions))

    # set the objective to minimize the flux through redox reactions
    model.objective = {model.reactions.get_by_id(rxn_id): 1.0 for rxn_id in redox_reactions}
    model.objective_direction = 'minimize'
    redox_solution = model.optimize()

    fluxes = [redox_solution.fluxes[reaction.id] for reaction in model.reactions if redox_solution.fluxes[reaction.id] != 0.0]
    print("Number of reactions with fluxes for min_redox != 0:", len(fluxes))

    return redox_solution


def map_flux_to_genes(model: cobra.Model, solution: cobra.Solution):
    """
    Maps flux values from an optimization solution to the corresponding genes in the metabolic model.

    Args:
        model: A cobra.Model object representing the metabolic model.
        solution: The optimization solution object containing flux values.
    Returns:
        gene_flux_dict: A dictionary mapping gene IDs to their corresponding flux values.
    """

    # initialize dictionary to store gene fluxes
    gene_flux_dict = {}
    for reaction in model.reactions:
        flux = solution.fluxes[reaction.id]
        # ignores reactions without GPRs or have 'nan' as their GPR
        if reaction.gene_reaction_rule not in ['', 'nan']:
            gpr = reaction.gene_reaction_rule
            # removes parentheses and logical operators from GPR string to allow to easy splitting
            for chr in ['(', ')', ' and', ' or']:
                gpr = gpr.replace(chr, '')
            # adds the absolute value of the flux to the gene in the gene_flux_dict; adds the gene first it if not present
            for gene in set(gpr.split(' ')):
                if gene not in gene_flux_dict:
                    gene_flux_dict[gene] = abs(flux)
                else:
                    gene_flux_dict[gene] += abs(flux)

    with open('genes_found_in_gpr.txt', 'w') as gene_file:
        gene_file.write("Genes\n")
        for gene in gene_flux_dict:
            gene_file.write(f"{gene}\n")

    return gene_flux_dict


def map_flux_to_genes_operator_conscious(model: cobra.Model, solution: cobra.Solution):
    """
    Maps flux values from an optimization solution to the corresponding genes in the metabolic model.

    Pays attention to logical operators in GPRs to avoid overestimating flux values for genes involved in 'or' relationships:
    1. With GPRs with both 'and' and 'or' operators, the flux value is divided by the number of 'or' operators in the GPR and assigned to each of the genes between them
    2. With GPRs with only 'and' operators, the entire flux value is assigned to each gene in the GPR.
    3. With GPRs with only 'or' operators, the flux value is divided by the number of genes in the GPR and assigned to each one.
    4. With GPRs with only one gene (which contain neither operator), the entire flux value is assigned to that gene.

    Args:
        model: A cobra.Model object representing the metabolic model.
        solution: The optimization solution object containing flux values.
    Returns:
        gene_flux_dict: A dictionary mapping gene IDs to their corresponding flux values.
    """

    # initialize dictionary to store gene fluxes
    gene_flux_dict = {}
    
    for reaction in model.reactions:
        flux = solution.fluxes[reaction.id]
        # ignores reactions without GPRs, have 'nan' as their GPR, or includes the following keywords in their GPR (SPONTANEOUS and Unassigned)
        if reaction.gene_reaction_rule not in ['', 'nan', 'SPONTANEOUS', 'Unassigned']:
            gpr = reaction.gene_reaction_rule
            # replaces the 'and' and 'or' operators with '&' and '|' respectively for use with sympy's to_dnf function
            gpr = gpr.replace(' and ', ' & ').replace(' or ', ' | ')

            # uses sympy's to_dnf function to convert to disjunctive normal form (DNF) for easier parsing.
            try:
                gpr_dnf = str(to_dnf(gpr, simplify=True))
            except ValueError: # catches GPRs with too many genes to convert to DNF quickly; these are handled as-is with the rest
                gpr_dnf = str(gpr)
                continue
            
            # flux is divided by the number of 'or' operators in the GPR +1 and assigned to each of the genes between them
            division_count = gpr_dnf.count('|') + 1
            flux /= division_count

            if '|' in gpr_dnf and '&' in gpr_dnf:
                # these will be GPRs with both 'and' and 'or' operators, so the flux value is divided by the number of 'or' operators in the GPR and assigned to each gene between them
                gpr_slices = gpr_dnf.split(' | ')
                for slice in gpr_slices:
                    if '(' and ')' in slice:
                        slice = slice.replace('(', '').replace(')', '')
                    genes = slice.split(' & ')
                    for gene in genes:
                        if gene not in gene_flux_dict:
                            gene_flux_dict[gene] = abs(flux)
                        else:
                            gene_flux_dict[gene] += abs(flux)

            elif '|' in gpr_dnf:
                # these will be GPRs with only 'or' operators, so the flux value is divided by the number of genes in the GPR and assigned to each one
                genes = gpr_dnf.split(' | ')
                for gene in genes:
                    if gene not in gene_flux_dict:
                        gene_flux_dict[gene] = abs(flux)
                    else:
                        gene_flux_dict[gene] += abs(flux)
                        
            elif '&' in gpr_dnf:
                # these will be GPRs with only 'and' operators, so the entire flux value is assigned to each gene in the GPR
                genes = gpr_dnf.split(' & ')
                for gene in genes:
                    if gene not in gene_flux_dict:
                        gene_flux_dict[gene] = abs(flux)
                    else:
                        gene_flux_dict[gene] += abs(flux)
            else:
                # these will be GPRs with only one gene, so the entire flux value is assigned to that gene
                if gpr_dnf not in gene_flux_dict:
                    gene_flux_dict[gpr_dnf] = abs(flux)
                else:
                    gene_flux_dict[gpr_dnf] += abs(flux)

    return gene_flux_dict


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input_XML_model',
                        help='(Req) Path to XML Model for Metabolic Pathways')
    parser.add_argument('-m', '--media_file',
                        help='(Req) Path to media file (tab-separated format with columns: compound, minflux, maxflux, concentration)')

    args = parser.parse_args()

    # load the XML model and validate it, then incorporate media constraints into the model by reading the media file and adjusting the bounds of exchange reactions accordingly
    model = load_xml_model(args.input_XML_model)
    media_model = incorporate_media(model, args.media_file)
    # media_model = incorporate_media_alternative(model, args.media_file)
    
    # make a copy of model for gene_essentialiity
    media_model = model.copy()

    """gene_essentiality = knockut_gene_essentiality_test(media_model)
    print("Gene essentiality test results:")
    for gene in gene_essentiality:
        print(f"{gene}: {gene_essentiality[gene]}")

    with open('gene_essentiality_results.tsv', 'w') as essentiality_file:
        essentiality_file.write("Gene\tBiomass_Objective_Value\n")
        for gene in gene_essentiality:
            essentiality_file.write(f"{gene}\t{gene_essentiality[gene]}\n")
    """

    print(len(media_model.reactions), "reactions in the media model")

    # generate optimization solutions for biomass, minimizing overall flux, and minimizing redox flux, then map the flux values from each to their corresponding genes in the model
    biomass_solution, min_flux_solution, min_redox_solution = optimize_model_biomass(media_model), optimize_model_min_flux(media_model), optimize_model_min_redox_flux(media_model)
    # biomass_fluxes, min_flux_fluxes, min_redox_fluxes = map_flux_to_genes(media_model, biomass_solution), map_flux_to_genes(media_model, min_flux_solution), map_flux_to_genes(media_model, min_redox_solution)
    biomass_fluxes, min_flux_fluxes, min_redox_fluxes = map_flux_to_genes_operator_conscious(media_model, biomass_solution), map_flux_to_genes_operator_conscious(media_model, min_flux_solution), map_flux_to_genes_operator_conscious(media_model, min_redox_solution)

    # load in protein expression data via pandas and filter for only the protein column and the average expression columns for each condition
    protein_expression_data_df = pd.read_csv('/data2/R_Choudhury/pa_fba/individual_scripts/P380_abundance_normalized_averaged_plus_tnseq.csv', sep=',')
    # protein_expression_data_df = protein_expression_data_df[['Accession','AS1_avg','SAS1_avg','ASA1_avg','SASA1_avg']]

    # test; determine the number of genes shared between the model and the protein expression data
    model_genes = set([gene.id for gene in model.genes])
    protein_genes = set(protein_expression_data_df['Accession'])
    shared_genes = model_genes.intersection(protein_genes)

    print("Number of genes in the model:", len(model_genes))
    print("Number of genes in the protein expression data:", len(protein_genes))
    print("Number of shared genes between model and protein expression data:", len(shared_genes))

    # initialize dictionaries to store flux values for each gene for each optimization method
    biomass = {}
    min_flux = {}
    min_redox = {}

    # maps the flux values from each optimization method to the corresponding genes in the protein expression data dataframe and adds new columns for each optimization method's flux values
    for gene in protein_expression_data_df['Accession']:
        biomass[gene] = biomass_fluxes.get(gene, 0.0)
        min_flux[gene] = min_flux_fluxes.get(gene, 0.0)
        min_redox[gene] = min_redox_fluxes.get(gene, 0.0)

    # maps the flux values from each optimization method to the corresponding genes in the protein expression data dataframe and adds new columns for each optimization method's flux values
    protein_expression_data_df['Biomass_Flux'] = protein_expression_data_df['Accession'].map(biomass)
    protein_expression_data_df['Min_Flux'] = protein_expression_data_df['Accession'].map(min_flux)
    protein_expression_data_df['Min_Redox_Flux'] = protein_expression_data_df['Accession'].map(min_redox)
    protein_expression_data_df.to_csv('pao1-proteins01_with_fluxes_all_intersect.csv', index=False)

    # alternatively, filter the protein expression data dataframe to only include rows where the 'Protein' column value is in the set of shared genes
    protein_expression_data_df_filtered = protein_expression_data_df[protein_expression_data_df['Accession'].isin(shared_genes)]
    protein_expression_data_df_filtered.to_csv('pao1-proteins01_with_fluxes_filtered_all_intersect.csv', index=False)


if __name__ == "__main__":
    main()
