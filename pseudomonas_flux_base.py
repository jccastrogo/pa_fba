#!/usr/bin/env python
import argparse
import pandas as pd
import cobra
import re


def create_sbml_model(compounds_file, reaction_file, output_file, save_to_file):
    """
    Creates an SBML metabolic model from reaction and compounds files.

    Args:
        reaction_file: Path to the reaction file (e.g., tab-separated format).
        compounds_file: Path to the compounds file (e.g., tab-separated format).
        output_file: Path to save the SBML model.
        save_to_file: Boolean value on whether to save the SBML model to file.

    Returns:
        None (saves the SBML model to the specified output file).
    """
    # Read reactions and compounds from files (you'll need to parse your specific formats)
    # Example: Read reactions and compounds into appropriate data structures

    # Create a cobra.Model object
    model = cobra.Model()

    # Load compounds and reactions from the CSV files
    compounds_df = pd.read_csv(compounds_file, sep='\t', encoding='unicode_escape')
    reactions_df = pd.read_csv(reaction_file, sep='\t', encoding='unicode_escape')
    for column in compounds_df.columns:
        compounds_df[column] = compounds_df[column].astype(str).str.rstrip()
    for column in reactions_df.columns:
        reactions_df[column] = reactions_df[column].astype(str).str.rstrip()

    # create entries for the metabolites and then adds them to the model
    for index, row in compounds_df.iterrows():
        metabolite = cobra.Metabolite(
            id=row['Abbreviation'],
            name=row['Name'],
            compartment=row['Abbreviation'][row['Abbreviation'].index('[') + 1:row['Abbreviation'].index(']')],
            formula=str(row['Formula (Charged)']),
            charge=float(row['Charge'])
        )
        model.add_metabolites([metabolite])

    # create entries for reactions (using metabolites by id for linking) and add them to the model
    # parses the csv's reaction column to get appropriate ids and counts of the metabolites per row

    for index, row in reactions_df.iterrows():
        # first checks if the reaction is an exchange, sink, or demand reaction (based on if the Subsystems column is 'Exchange')
        if row['Subsystems'] == 'Exchange':
            # reversible reactions with [c] are demand reactions, irreversible reactions are sink reactions, and reversible reactions without [c] are exchange reactions
            if '->' in row['Reaction']:
                model.add_boundary(model.metabolites.get_by_id(row['Reaction'].split(' ')[0]), type='sink')
            elif '[c]' in row['Reaction'] and '<=>' in row['Reaction']:
                model.add_boundary(model.metabolites.get_by_id(row['Reaction'].split(' ')[0]), type='demand')
            else:
                model.add_boundary(model.metabolites.get_by_id(row['Reaction'].split(' ')[0]), type='exchange')

        else:
            # adds the reaction's id, name, lower_bound, and upper_bound from the csv file
            reaction = cobra.Reaction(id=row['Abbreviation'],
                                      name=row['Name'],
                                      lower_bound=float(row['Lower bound']),
                                      upper_bound=float(row['Upper bound']),
                                      )
            
            # adds gene names to reaction if present in the reactions file
            if pd.notna(row['GPR']) and row['GPR'] not in ['SPONTANEOUS', 'Unassigned']:
                reaction.gene_reaction_rule = row['GPR']

            positive_or_negative = -1.0
            metabolite_dict = {}

            # splits reaction column into individual strings (removes any space characters and '+' from the list)
            reaction_col_list = re.findall(r'\S+', row['Reaction'])
            reaction_col_list = [element for element in reaction_col_list if element != '+']

            """
            Runs through the reaction_col_list in a while loop:
            If the element can be converted into a float, then it will be become a negative (before reaction) or 
            positive (after reaction) float and takes the succeeding element as its key.
            If the element cannot be converted (reagent/product), then the quantity is treated as -1/+1 based on its position.
            """
            i = 0
            while i < len(reaction_col_list):
                if reaction_col_list[i] not in ['->', '<=>']:
                    try:
                        quantity = float(reaction_col_list[i]) * positive_or_negative
                        metabolite_dict[reaction_col_list[i + 1]] = quantity
                        i += 1
                    except ValueError:
                        metabolite_dict[reaction_col_list[i]] = positive_or_negative
                else:
                    # flips the sign for products when the reaction symbol is detected
                    positive_or_negative = 1.0
                i += 1

            # first add reaction to the model    
            model.add_reactions([reaction])
            # adds metabolites to the reaction using the metabolite_dict
            model.reactions.get_by_id(row['Abbreviation']).add_metabolites(metabolite_dict)
            # add objective coefficient to the reaction
            # model.reactions.get_by_id(row['Abbreviation']).objective_coefficient = row['Objective']
            # add reversibility to the reaction
            # check as to why ATPM and rxn00695 cannot have their reversibility set (+1 other)
            model.reactions.get_by_id(row['Abbreviation']).reversibility = True if int(row['Reversible']) == 1 else False
            # add Subsystem to the reaction
            model.reactions.get_by_id(row['Abbreviation']).subsystem = row['Subsystems']

    # Save the model to SBML format
    if save_to_file:
        cobra.io.write_sbml_model(model, output_file)

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

    # initialize dictionary to store gene fluxes
    biomass_gene_flux_dict = {}
    for reaction in model.reactions:
        flux = biomass_solution.fluxes[reaction.id]
        # ignores reactions without GPRs and/or have 'nan' as their GPR
        if reaction.gene_reaction_rule not in ['', 'nan']:
            gpr = reaction.gene_reaction_rule
            # removes parentheses and logical operators from GPR string to allow to easy splitting
            for chr in ['(', ')', ' and', ' or']:
                gpr = gpr.replace(chr, '')
            # splits GPR string by spaces to get individual gene names and adds the absolute flux to each gene in the dictionary
            for gene in gpr.split(' '):
                if gene not in biomass_gene_flux_dict:
                    biomass_gene_flux_dict[gene] = abs(flux)
                biomass_gene_flux_dict[gene] += abs(flux)

    # sets any negative flux values to 0
    for key in biomass_gene_flux_dict:
        if biomass_gene_flux_dict[key] < 0:
            biomass_gene_flux_dict[key] = 0

    return biomass_gene_flux_dict


def optimize_model_min_flux(model: cobra.Model):
    """
    Optimizes metabolic model for minimizing overall flux (parsimonious FBA).

    Args:
        model: A cobra.Model object representing the metabolic model.

    Returns:
        solution: The optimization solution object.
    """

    pfba_solution = cobra.flux_analysis.parsimonious.pfba(model)

    # initialize dictionary to store gene fluxes
    pfba_gene_flux_dict = {}
    for reaction in model.reactions:
        flux = pfba_solution.fluxes[reaction.id]
        # ignores reactions without GPRs and/or have 'nan' as their GPR
        if reaction.gene_reaction_rule not in ['', 'nan']:
            gpr = reaction.gene_reaction_rule
            # removes parentheses and logical operators from GPR string to allow to easy splitting
            for chr in ['(', ')', ' and', ' or']:
                gpr = gpr.replace(chr, '')
            # splits GPR string by spaces to get individual gene names and adds the absolute flux to each gene in the dictionary
            for gene in gpr.split(' '):
                if gene not in pfba_gene_flux_dict:
                    pfba_gene_flux_dict[gene] = abs(flux)
                pfba_gene_flux_dict[gene] += abs(flux)

    # sets any negative flux values to 0
    for key in pfba_gene_flux_dict:
        if pfba_gene_flux_dict[key] < 0:
            pfba_gene_flux_dict[key] = 0
    
    return pfba_gene_flux_dict


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

    # initialize dictionary to store gene fluxes
    redox_gene_flux_dict = {}
    for reaction in model.reactions:
        flux = redox_solution.fluxes[reaction.id]
        # ignores reactions without GPRs and/or have 'nan' as their GPR
        if reaction.gene_reaction_rule not in ['', 'nan']:
            gpr = reaction.gene_reaction_rule
            # removes parentheses and logical operators from GPR string to allow to easy splitting
            for chr in ['(', ')', ' and', ' or']:
                gpr = gpr.replace(chr, '')
            # splits GPR string by spaces to get individual gene names and adds the absolute flux to each gene in the dictionary
            for gene in gpr.split(' '):
                if gene not in redox_gene_flux_dict:
                    redox_gene_flux_dict[gene] = abs(flux)
                redox_gene_flux_dict[gene] += abs(flux)

    # sets any negative flux values to 0
    for key in redox_gene_flux_dict:
        if redox_gene_flux_dict[key] < 0:
            redox_gene_flux_dict[key] = 0
    
    return redox_gene_flux_dict


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--compounds_file',
                        help='(Req) Path to the compounds file (tab-separated format)',
                        metavar='', type=str, required=True)
    parser.add_argument('-r', '--reaction_file',
                        help='(Req) Path to the reaction file (tab-separated format)',
                        metavar='', type=str, required=True)
    parser.add_argument('-o', '--output_file',
                        help='(Opt) Path to save the SBML model.',
                        metavar='', type=str, required=False, default='pseudomonas_model.xml')
    parser.add_argument('-s', '--save_to_file', help="Flag to save the SBML model to file.",
                        required=False, default=False)

    args = parser.parse_args()

    model = create_sbml_model(args.compounds_file, args.reaction_file, args.output_file, args.save_to_file)

    biomass_solution = optimize_model_biomass(model)

    min_flux_solution = optimize_model_min_flux(model)

    min_redox_solution = min_redox_flux(model)

    # load in protein expression data via pandas
    
    protein_expression_data = pd.read_csv('ancestor_pao1_plus_aureus_greater_than_0.csv')
    biomass = {}
    min_flux = {}
    min_redox = {}

    for gene in protein_expression_data['Protein']:
        biomass[gene] = biomass_solution.get(gene, 0.0)
        min_flux[gene] = min_flux_solution.get(gene, 0.0)
        min_redox[gene] = min_redox_solution.get(gene, 0.0)

    protein_expression_data['Biomass_Flux'] = protein_expression_data['Protein'].map(biomass)
    protein_expression_data['Min_Flux'] = protein_expression_data['Protein'].map(min_flux)
    protein_expression_data['Min_Redox_Flux'] = protein_expression_data['Protein'].map(min_redox)

    protein_expression_data.to_csv('pao1_plus_aureus_protein_expression_with_fluxes.csv', index=False)


if __name__ == "__main__":
    main()
