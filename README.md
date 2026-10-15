# Dengue in Valle del Cauca, Colombia, 2022

Analysis code for *Spatio-temporal clustering of dengue in Valle del Cauca, Colombia*. The R file covers the spatial analysis and the comparisons of maximum temporal windows. The Stata do-file covers hospitalization.

The surveillance data are confidential and are not included. Access to de-identified records from the Secretaría de Salud del Valle del Cauca requires institutional review and a data-use agreement. Input paths are set at the top of each script. Outputs containing case coordinates should remain outside the public repository.

## Spatial analysis in R

`Aim1_spatial_analysis.R` uses R 4.1 or later, the packages `terra`, `readxl` and `lubridate`, and the SaTScan 10.3.2 command-line executable.

The inputs are an Excel case file and the WorldPop 2020 unadjusted population raster for Colombia, at approximately 100 m resolution. The case file has one row per case, with consultation date (`fec_con_`), classification (`tip_cas_`, `ajuste_`) and residential coordinates (`latitude`, `longitude`, WGS84). Laboratory confirmation is coded as 3 in either classification field. Dates can be Excel dates or `YYYY-MM-DD`. Selected cases must have valid consultation dates in 2022. The population raster is a single count layer in geographic coordinates compatible with WGS84.

The script selects cases within the nine study rectangles and crops the population raster to the extent of the selected case points. Counts are grouped by raster cell and consultation date. It writes the `.cas`, `.geo`, `.pop` and parameter files and runs SaTScan for each area with eligible cases. The maximum temporal windows are six months for the primary analysis, and three and one month for the sensitivity analyses. The remaining scan settings are in the parameter block.

Set the case, population, SaTScan and output paths, then run:

```sh
Rscript Aim1_spatial_analysis.R.
```

Alongside the SaTScan files, the script saves case counts in `case_flow.csv` and `cases_by_municipality.csv`, cluster details in `cluster_results.csv`, and summaries in `table1_primary_scan.csv` and `s3_table_temporal_windows.csv`. R and package versions are printed by `sessionInfo()`.

## Hospitalization analysis in Stata

`Aim2_hospitalization.do` uses Stata 18 and a prepared CSV. The classification fields are named `confirm` and `adjust`, corresponding to `tip_cas_` and `ajuste_` in the spatial data. The other fields are `hosp`, `stratum`, `sex`, `ethn`, `insure`, `age`, `severe` and `manage`.

Hospitalization is coded 0/1, SES as strata 1–6, sex as 1 male and 2 female, and ethnicity code 5 identifies Black/Afro-Colombian cases. Insurance codes are 1 uninsured, 2 contributory, 3 subsidized and 4 other. Age is in years. Severity is coded 1 without warning signs, 2 with warning signs and 3 severe dengue. `manage` records the clinical care setting.

The do-file starts with recoding and descriptive comparisons. The primary model includes SES, sex, ethnicity, insurance and age group among laboratory-confirmed cases. It then assesses SES interactions, fits the same model to all notified cases, and fits a severity-adjusted model. Covariates are removed individually from the severity-adjusted model to calculate changes in the SES coefficient. Each model uses complete cases for its variables.

Set the CSV and output paths. The parent of the output folder must already exist. Run in Stata:

```stata
do "Aim2_hospitalization.do"
```

The results, including Hosmer–Lemeshow tests with 10 groups, link tests and VIF, are saved in `hospitalization_analysis.log`. The log also includes the interaction tests and retrospective sample-size calculation. Model coefficients are saved in `primary_coefficients.csv`, `all_notified_coefficients.csv` and `with_severity_coefficients.csv`.
