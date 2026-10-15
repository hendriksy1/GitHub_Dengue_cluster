* Aim 2: hospitalization among dengue cases, Valle del Cauca, 2022 (Stata 18)

version 18.0
clear all
set more off
set varabbrev off

* File paths
local input_csv "path/to/hospitalization_cases.csv"
local output_dir "path/to/private_outputs/hospitalization"
capture mkdir "`output_dir'"

capture log close hospitalization
log using "`output_dir'/hospitalization_analysis.log", name(hospitalization) text replace

* 1. Recoding
* Outcome: hospitalization. Exposure: SES 1-2 versus 3-6
import delimited using "`input_csv'", clear varnames(1) case(lower)
generate byte lab_confirmed = (confirm == 3 | adjust == 3)
recode stratum (1 2 = 1) (3 4 5 6 = 0) (else = .), generate(stratum_grouped)
generate byte ethn_binary = (ethn == 5) if !missing(ethn)
generate age_band70 = floor(age/10)*10
replace age_band70 = 70 if age_band70 >= 70 & !missing(age_band70)
label define ses_lbl 0 "Medium-high (3-6)" 1 "Low (1-2)"
label define ethn_lbl 0 "Other" 1 "Black/Afro-Colombian"
label define insure_lbl 1 "Uninsured" 2 "Contributory" 3 "Subsidized" 4 "Other"
label define sex_lbl 1 "Male" 2 "Female"
label define severe_lbl 1 "Without warning signs" 2 "With warning signs" 3 "Severe dengue"
label values stratum_grouped ses_lbl
label values ethn_binary ethn_lbl
label values insure insure_lbl
label values sex sex_lbl
label values severe severe_lbl
assert inlist(hosp,0,1) if !missing(hosp)
misstable summarize hosp stratum sex ethn insure age severe

* Case counts, confirmation codes and missing SES
count
count if lab_confirmed
tabulate confirm adjust, missing
count if adjust == 6
count if missing(stratum)
count if lab_confirmed == 1 & missing(stratum)
count if lab_confirmed == 1 & age >= 70 & !missing(age)
count if lab_confirmed == 1 & hosp == 1

* 2. S2 Table: cases by laboratory confirmation
bysort lab_confirmed: summarize age, detail
ranksum age, by(lab_confirmed)
foreach v in sex stratum insure severe hosp manage {
    tabulate `v' lab_confirmed, chi2 column
}
tabulate ethn_binary lab_confirmed, exact column
tabulate stratum lab_confirmed, column missing

* 3. Table 2: laboratory-confirmed cases by SES
bysort stratum_grouped: summarize age if lab_confirmed == 1, detail
ranksum age if lab_confirmed == 1, by(stratum_grouped)
foreach v in sex insure severe hosp manage {
    tabulate `v' stratum_grouped if lab_confirmed == 1, chi2 column
}
tabulate ethn_binary stratum_grouped if lab_confirmed == 1, exact column

* 4. Bivariate associations with hospitalization
ranksum age if lab_confirmed == 1, by(hosp)
foreach v in stratum_grouped sex insure severe manage {
    tabulate `v' hosp if lab_confirmed == 1, chi2 column
}
tabulate ethn_binary hosp if lab_confirmed == 1, exact column


* 5. Table 3: primary hospitalization model
* Reference groups: SES 3-6, male, other ethnicity, contributory insurance, age 20-29
local covariates "ib1.sex ib0.ethn_binary ib2.insure ib20.age_band70"
logit hosp ib0.stratum_grouped `covariates' if lab_confirmed == 1, or
assert e(converged) == 1
estimates store primary
generate byte primary_sample = e(sample)
display as text "Table 3 estimation sample: " e(N)
lincom 2.sex, eform
* Hosmer-Lemeshow (10 groups), link test and VIF
estat gof, group(10) table
linktest
quietly regress hosp ib0.stratum_grouped `covariates' if primary_sample
estat vif
estimates restore primary

* Crude SES odds ratio, using the same cases as the primary model
logit hosp ib0.stratum_grouped if primary_sample, or
estimates restore primary

* 6. SES interactions with sex, ethnicity, insurance and age
estimates restore primary
local primary_N = e(N)

* SES by sex
logit hosp ib0.stratum_grouped##ib1.sex ib0.ethn_binary ib2.insure ib20.age_band70 if primary_sample, or
assert e(converged) == 1 & e(N) == `primary_N'
testparm i.stratum_grouped#i.sex

* SES by ethnicity
logit hosp ib0.stratum_grouped##ib0.ethn_binary ib1.sex ib2.insure ib20.age_band70 if primary_sample, or
assert e(converged) == 1 & e(N) == `primary_N'
testparm i.stratum_grouped#i.ethn_binary

* SES by insurance
logit hosp ib0.stratum_grouped##ib2.insure ib1.sex ib0.ethn_binary ib20.age_band70 if primary_sample, or
assert e(converged) == 1 & e(N) == `primary_N'
testparm i.stratum_grouped#i.insure

* SES by age group
logit hosp ib0.stratum_grouped##ib20.age_band70 ib1.sex ib0.ethn_binary ib2.insure if primary_sample, or
assert e(converged) == 1 & e(N) == `primary_N'
testparm i.stratum_grouped#i.age_band70

* 7. Retrospective sample size, with observed SES odds ratio of 1.30
quietly summarize stratum_grouped if primary_sample
local px = r(mean)
quietly summarize hosp if primary_sample & stratum_grouped == 0
local p0 = r(mean)
local p1 = (1.30*`p0')/(1-`p0'+1.30*`p0')
local nratio = `px'/(1-`px')
power twoproportions `p0' `p1', nratio(`nratio') power(0.80) alpha(0.05)

* 8. S1 Table: primary model fitted to all notified cases
logit hosp ib0.stratum_grouped `covariates', or
assert e(converged) == 1
estimates store all_notified
* Hosmer-Lemeshow (10 groups), link test and VIF
generate byte all_notified_sample = e(sample)
estat gof, group(10) table
linktest
quietly regress hosp ib0.stratum_grouped `covariates' if all_notified_sample
estat vif
estimates restore all_notified

* 9. Severity-adjusted model and SES-by-severity interaction
tabulate severe if primary_sample, missing
generate byte severity_sample = primary_sample & inlist(severe,1,2,3)
count if severity_sample
bysort stratum_grouped: tabulate severe hosp if severity_sample, missing
logit hosp ib0.stratum_grouped `covariates' ib1.severe if severity_sample, or
assert e(converged) == 1
estimates store with_severity
local severity_beta = _b[1.stratum_grouped]
* Hosmer-Lemeshow (10 groups), link test and VIF
generate byte with_severity_sample = e(sample)
estat gof, group(10) table
linktest
quietly regress hosp ib0.stratum_grouped `covariates' ib1.severe if with_severity_sample
estat vif
estimates restore with_severity

* All severe-dengue cases in SES 3-6 were hospitalized
* The full interaction is not estimable because of separation
capture noisily logit hosp ib0.stratum_grouped##ib1.severe `covariates' if severity_sample, or

* Remove each covariate and calculate the absolute % change in SES beta
local full_terms "`covariates' ib1.severe"
foreach removed of local full_terms {
    local reduced : list full_terms - removed
    quietly logit hosp ib0.stratum_grouped `reduced' if severity_sample, or
    local change = 100*abs((_b[1.stratum_grouped]-`severity_beta')/`severity_beta')
    display as text "Removed `removed': absolute percent change in SES beta = " %8.2f `change'
}

estimates restore primary

* 10. Coefficient exports: primary, all notified and severity-adjusted models
foreach model in primary all_notified with_severity {
    estimates restore `model'
    tempname coefficient_post
    tempfile coefficient_data
    postfile `coefficient_post' str60 term double beta se odds_ratio lower95 upper95 p_value N using `coefficient_data', replace
    local terms : colnames e(b)
    foreach term of local terms {
        if _se[`term'] > 0 {
            post `coefficient_post' ("`term'") (_b[`term']) (_se[`term']) (exp(_b[`term'])) ///
                (exp(_b[`term']-invnormal(.975)*_se[`term'])) ///
                (exp(_b[`term']+invnormal(.975)*_se[`term'])) ///
                (2*normal(-abs(_b[`term']/_se[`term']))) (e(N))
        }
    }
    postclose `coefficient_post'
    preserve
    use `coefficient_data', clear
    export delimited using "`output_dir'/`model'_coefficients.csv", replace
    restore
}

log close hospitalization
