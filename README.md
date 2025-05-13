# TADA - Target Aggregate Data Adjustment Method for Transportability Analysis Utilizing Summary-Level Data from the Target Population
R codes to implement the simulation study in "Target Aggregate Data Adjustment Method for Transportability Analysis Utilizing Summary-Level Data from the Target Population".

Here is a brief introduction to the function of each R file. Comprehensive comments provide detailed instructions in each file. Setting the relative working path for storing results and loading data is necessary.

### `TADA-RCS-dataWrangling.R`
Prepares individual-level and aggregate-level datasets for real-case analysis, converting raw SAS files into cleaned CSV formats. This script extracts key covariates, formats censoring indicators, and saves processed IPD and AgD inputs.

### `TADA-realCaseStudy.R`
Implements the TADA method on a real-world survival dataset, incorporating time-varying IPCW and MoM weights to generate adjusted Kaplan-Meier curves and treatment effect estimates.

### `assistFunctions.R`
Defines core simulation functions including `generateTestDataTADA_TTE()` and `generateTestDataTADA_TTE_ABV()`, used to simulate survival data with censoring and transportability scenarios under varying treatment-effect modifiers.

### `TADA-simulation.R`
Conducts the main simulation study evaluating TADA’s performance across 10 scenarios of treatment-covariate interactions, reporting bias, coverage, and confidence interval estimates with and without censoring adjustment.

### `TADA-SA-abnormalValue.R`
Performs sensitivity analysis on TADA by introducing abnormal values in the target aggregate data, examining robustness of estimates under data perturbation including X2 flips and X3 extremes.

### `TADA-SA-sampleSize.R`
Assesses the impact of varying source sample sizes (e.g., n = 300 to 1000) on TADA’s bias, precision, and coverage performance under controlled simulation settings.

### `TADA-SA-truncProp.R`
Explores sensitivity to different final weight truncation thresholds (80%–99%) in TADA, evaluating their effects on estimator stability and coverage across 500 simulation replicates.

We will consistently provide the necessary maintenance for these major R files.
