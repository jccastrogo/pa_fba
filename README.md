# Pseudomona aeruginosa Flux Balance Analysis

Use a whole genome metabolic model to perform several analyses

- Perform flux balace analysis using a custom objectuve function
- Recursively knock-out genes in the genomes and perform flux balance analysis
- Perform flux balance analysis in a custom media
- Map genes expression to reactions in the model
- Predict genes expression level

Usage

```bash
src/pseudomonas_flux.py -c config_file.json
```

Where ```config_file.json``` is a comfiguration file in json format that includes the option of which analysis to perform on which datasets.



For additional help use:

```bash
src/pseudomonas_flux.py -h
```



Continued development of this pipeline has moved to [nf_pafba](https://github.com/jccastrogo/nf_pafba) using [![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)# pa_fba
