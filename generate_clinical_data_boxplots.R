# Description ------------------------------------------------------------------

# Script is intended to generate the box plots for the transcriptomic data 

# Load Packages ----------------------------------------------------------------

library(dplyr)
library(ggplot2)

# Provide input files ----------------------------------------------------------

# load generated mapping file from Iris's nextflow script that also includes a column mapping the genes
# make said Gene_ID the rownames
input_mapping_csv <- read.csv('GCF_000006765.1_ASM676v1_genomic_merged_with_geneID.csv', header = TRUE)
rownames(input_mapping_csv) <- input_mapping_csv$Gene_ID

# load metadata file that includes the accessions for the transcriptomic samples
# sets rownames as the Accessions
input_metadata_csv <- read.csv('Clinical_Sample_Metadata.csv', header = TRUE)
rownames(input_metadata_csv) <- input_metadata_csv$Accession

#  --------------

# filters only for TPM columns as the selected source of data points, then removes that pieces to line up with the accessions in the metadata file
input_mapping_csv_TPM <- input_mapping_csv[ ,grepl("_TPM", names(input_mapping_csv))]
names(input_mapping_csv_TPM) <- gsub("_TPM", "", names(input_mapping_csv_TPM))

# filters only for the named columns to be for the above data frame to be joined
input_csv_named <- input_csv[,c('Tn.Seq_Biomass_Flux', 'Biomass_Flux', 'Accession')]
rownames(input_csv_named) <- input_csv$Accession

# left join into fluxes
fluxes_csv_plus_TPM <- merge(input_csv_named, input_mapping_csv_TPM, by = "row.names", all.x = TRUE)
rownames(fluxes_csv_plus_TPM) <- fluxes_csv_plus_TPM$Row.names

# removes unecessary rows
fluxes_csv_plus_TPM <- fluxes_csv_plus_TPM %>% select(-Row.names, -Accession)

# generated blank columns for the correlations
input_metadata_csv$Biomass_Cor <- NA
input_metadata_csv$TnSeq_Cor <- NA

# runs through each accession and assigns each row of the input_metadata_csv column with the correlation between the biomass flux or tn-seq flux and the TPM of that accession
for (i in 3:length(colnames(fluxes_csv_plus_TPM))) {
  print(i)
  print(colnames(fluxes_csv_plus_TPM)[i])
  print(cor(fluxes_csv_plus_TPM$Tn.Seq_Biomass_Flux, fluxes_csv_plus_TPM[,i], use = 'complete.obs'))
  input_metadata_csv[colnames(fluxes_csv_plus_TPM)[i], 'TnSeq_Cor'] <- cor(fluxes_csv_plus_TPM$Tn.Seq_Biomass_Flux, fluxes_csv_plus_TPM[, colnames(fluxes_csv_plus_TPM)[i]], use='complete.obs')
  input_metadata_csv[colnames(fluxes_csv_plus_TPM)[i], 'Biomass_Cor'] <- cor(fluxes_csv_plus_TPM$Biomass_Flux, fluxes_csv_plus_TPM[, colnames(fluxes_csv_plus_TPM)[i]], use='complete.obs')
}

# Generate Box Plots ------------------------------------------------------------

# Base Biomass Correlation
ggplot(input_metadata_csv, aes(x=Type, y=Biomass_Cor)) +
  geom_boxplot() +
  geom_jitter() +
  labs(x = "Type", y = "Biomass Correlation", title='Biomass Correlations by Type')

# Biomass Correlation + Tn-Seq
ggplot(input_metadata_csv, aes(x=Type, y=TnSeq_Cor)) +
  geom_boxplot() +
  geom_jitter() +
  labs(x = "Type", y = "Tn-Seq Correlation", title='Tn-Seq Correlations by Type')

