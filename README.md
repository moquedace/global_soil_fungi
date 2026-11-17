# Global soil fungi

R workflow supporting the article *Soil fungal diversity is greatest where
taxonomic knowledge is poorest, not where ectomycorrhizal lineages are
richest* (Moquedace et al., in preparation).

> **Status:** in preparation. The repository stays private until the article
> is submitted, and the scripts are being adapted from the working code.

- Maps (Zenodo, released at submission): https://doi.org/10.5281/zenodo.24292559
- Google Earth Engine assets: `projects/ee-cassiomoquedace/assets/soil_fungi`

## Structure

```text
global_soil_fungi/
|-- .github/workflows/        # repository validation
|-- assets/                   # banner, map preview, logos
|-- docs/methods_and_code.md  # article -> scripts, data sources, versions
|-- scripts/
|   |-- 00_config.R           # project root and settings
|   |-- functions/            # shared model and figure functions
|   |-- 01_covariates/        # grid, domain, fractions, resampling, gaps, stack
|   |-- 02_responses/         # GloSED download and the six responses
|   |-- 03_modeling/          # extraction, final model, spatial hold-out
|   |-- 04_prediction/        # tiled prediction and mosaics
|   |-- 05_figures_tables/    # figures, tables and numbers of the article
|   `-- 06_distribution/      # COG int16 copy and Earth Engine strips
|-- hpc/                      # generic SLURM job
|-- .zenodo.json
|-- CITATION.cff
|-- DESCRIPTION
|-- LICENSE
`-- README.md
```

## License

Code released under the MIT License (see `LICENSE`).