# README

This contains the R code used to reproduce the analyses in the manuscript for the article `Bayesian Latent Class Regression and Variable Selection with Applications to Sleep Patterns Data`.

## Files
- `run_all.R`: main script to run the full analysis
- `BLCR_master_script.R`: master script that sources functions and runs all simulation study and data analysis scripts
- `01_simulation_study_1.R`: Simulation Study 1
- `01_1_simulation_study_1_sample_sensitivity.R`: sensitivity analysis for Simulation Study 1
- `02_simulation_study_2.R`: Simulation Study 2
- `02_1_simulation_study_2_sample_sensitivity.R`: sensitivity analysis for Simulation Study 2
- `03_cshq_analysis.R`: applied CSHQ analysis

## Requirements
The code requires R and the additional helper scripts in the `R/` folder, along with input data files in the `data/` folder.

## How to run
Run `run_all.R` to execute the full workflow.

## Output files
The full code reproduces all results in the manuscript. To keep the submission size manageable, not all generated output files are included. The scripts will recreate the figures, tables, and saved result objects in the `output/` directory when run. 