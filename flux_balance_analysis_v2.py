#!/usr/bin/env python
import argparse
import pandas as pd
import cobra


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

    return model


def incorporate_media(model: cobra.Model, media_file_path: str):
    """
    Incorporates media constraints into the metabolic model by reading a media file and adjusting the bounds of exchange reactions accordingly.

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
    biomass_solution = model.optimize()

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

    # optimize the model using parsimonious FBA to minimize overall flux while maintaining the same objective value as the biomass optimization
    pfba_solution = cobra.flux_analysis.parsimonious.pfba(model)

    fluxes = [pfba_solution.fluxes[reaction.id] for reaction in model.reactions if pfba_solution.fluxes[reaction.id] != 0.0]
    print("Number of reactions with fluxes for parsimonious FBA != 0:", len(fluxes))

    return pfba_solution


def min_redox_flux(model: cobra.Model):
    """
    Optimizes metabolic model by minimizing redox flux (NADH/NADPH/FADH2).
    NADH, NADPH, and FADH2; cpd00004[c], cpd00005[c], and cpd00982[c], respectively.

    Args:
        model: A cobra.Model object representing the metabolic model.
    Returns:
        solution: The optimization solution object.
    """

    # listed below (in a set) are the compound ids for NADH, NADPH, and FADH2
    redox_compounds = set(['cpd00004[c]', 'cpd00005[c]', 'cpd00982[c]'])

    # initialize list entry to store reactions producing redox compounds
    redox_reactions = []

    # identifies reactions that produce the above redox compounds and adds them to the redox_reactions list
    for reaction in model.reactions:
        product_ids = set([metabolite.id for metabolite in reaction.products])
        if len(redox_compounds.intersection(product_ids)) > 0:
            redox_reactions.append(reaction.id)

    # set the objective to minimize the flux through redox reactions
    model.objective = {model.reactions.get_by_id(rxn_id): 1.0 for rxn_id in redox_reactions}
    model.objective_direction = 'min'
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
        # ignores reactions without GPRs and/or have 'nan' as their GPR
        if reaction.gene_reaction_rule not in ['', 'nan']:
            gpr = reaction.gene_reaction_rule
            # removes parentheses and logical operators from GPR string to allow to easy splitting
            for chr in ['(', ')', ' and', ' or']:
                gpr = gpr.replace(chr, '')
            # splits GPR string by spaces to get individual gene names and adds the absolute flux to each gene in the dictionary
            for gene in set(gpr.split(' ')):
                if gene not in gene_flux_dict:
                    gene_flux_dict[gene] = abs(flux)
                gene_flux_dict[gene] += abs(flux)

    with open('genes_found_in_gpr.txt', 'w') as gene_file:
        gene_file.write("Genes\n")
        for gene in gene_flux_dict:
            gene_file.write(f"{gene}\n")

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

    # generate optimization solutions for biomass, minimizing overall flux, and minimizing redox flux, then map the flux values from each to their corresponding genes in the model
    biomass_solution, min_flux_solution, min_redox_solution = optimize_model_biomass(media_model), optimize_model_min_flux(media_model), min_redox_flux(media_model)
    biomass_fluxes, min_flux_fluxes, min_redox_fluxes = map_flux_to_genes(media_model, biomass_solution), map_flux_to_genes(media_model, min_flux_solution), map_flux_to_genes(media_model, min_redox_solution)

    # load in protein expression data via pandas and filter for only the protein column and the average expression columns for each condition
    protein_expression_data_df = pd.read_csv('/data2/R_Choudhury/pa_fba/example/pao1-proteins01.tsv', sep='\t')
    protein_expression_data_df = protein_expression_data_df[['Protein','Avg_AS1','Avg_SAS1','Avg_ASA1','Avg_SASA1']]

    # test; determine the number of genes shared between the model and the protein expression data
    model_genes = set([gene.id for gene in model.genes])
    protein_genes = set(protein_expression_data_df['Protein'])
    shared_genes = model_genes.intersection(protein_genes)

    print("Number of genes in the model:", len(model_genes))
    print("Number of genes in the protein expression data:", len(protein_genes))
    print("Number of shared genes between model and protein expression data:", len(shared_genes))


    # initialize dictionaries to store flux values for each gene for each optimization method
    biomass = {}
    min_flux = {}
    min_redox = {}

    # maps the flux values from each optimization method to the corresponding genes in the protein expression data dataframe and adds new columns for each optimization method's flux values
    for gene in protein_expression_data_df['Protein']:
        biomass[gene] = biomass_fluxes.get(gene, 0.0)
        min_flux[gene] = min_flux_fluxes.get(gene, 0.0)
        min_redox[gene] = min_redox_fluxes.get(gene, 0.0)

    # maps the flux values from each optimization method to the corresponding genes in the protein expression data dataframe and adds new columns for each optimization method's flux values
    protein_expression_data_df['Biomass_Flux'] = protein_expression_data_df['Protein'].map(biomass)
    protein_expression_data_df['Min_Flux'] = protein_expression_data_df['Protein'].map(min_flux)
    protein_expression_data_df['Min_Redox_Flux'] = protein_expression_data_df['Protein'].map(min_redox)
    protein_expression_data_df.to_csv('pao1-proteins01_with_fluxes.csv', index=False)

    # alternatively, filter the protein expression data dataframe to only include rows where the 'Protein' column value is in the set of shared genes
    protein_expression_data_df_filtered = protein_expression_data_df[protein_expression_data_df['Protein'].isin(shared_genes)]
    protein_expression_data_df_filtered.to_csv('pao1-proteins01_with_fluxes_filtered.csv', index=False)


if __name__ == "__main__":
    main()
