# Clinical and Economic Impact of Molecular Testing for BRAF Fusion in Pediatric Low-Grade Glioma

**Updated**: October 27 2021 

R code to recreate cost-effectiveness analysis in the 'Clinical and Economic Impact of Molecular Testing for BRAF Fusion in Pediatric Low-Grade Glioma' manuscript. For detailed description of model structure, assumptions and limitations please see published manuscript.

**Authors** : Juan David Rios MA (1), Russanthy Velummailum MPH (1), Julie Bennett MD (2), Liana Nobre MD (2), Derek S. Tsang MD (2,4),  Eric Bouffet MD (2) , Cynthia Hawkins MD (3), Uri Tabori MD (2), Avram Denburg MD (1,2,5)\*, Petros Pechlivanoglou Ph.D (1,5)\* 

**Affiliations:** 

1. Child Health Evaluative Sciences, Peter Gilgan Centre for Research and Learning, The Hospital for Sick Children, Toronto, Canada. 
2. Division of Haematology/Oncology, The Hospital for Sick Children, Toronto, Canada
3. Department of Pathology, Hospital for Sick Children, Toronto, ON, Canada.
4. Radiation Medicine Program, Princess Margaret Cancer Centre, University Health Network, Toronto, Canada.
5. Institute of Health Policy, Management and Evaluation, Dalla Lana School of Public Health, University of Toronto, Toronto, Canada.

\*Co-senior authors

R code by [Juan David Rios](https://github.com/David-rios), Russanthy Velummailum, and [Petros Pechlivanoglou](https://github.com/ppehli).

---

### Running the cost-effectiveness model 

To run the cost-effectiveness model go to `Cost-effectiveness_basecase.R` and run the file line by line. Afterwards ` Manuscript results.Rmd` will summarize the results. 

Running the model for 100,000 individuals, 1,000 times for the radiation benefit and no radiation benefit scenarios required 100 days of computational time on a high-performance cluster. If running locally we recommend using a smaller number of individuals, or the deterministic scenario available in `Cost-effectiveness_scenario.R`.

The analysis is structured using the recommendations of 'A Need for Change! A Coding Framework for Improving Transparency in Decision Modeling' (https://doi.org/10.1007/s40273-019-00837-x). 



---

### Modifying the cost-effectiveness model 

* To modify costs: modify the function `gen_costs()` in `Costs.R`.
* To modify utilities: modify the function `gen_Effects()` in `Effects.R`.
* To modify incidence of radiation related AE, modify the `EstProbs.AE()` function in `Fun_surv_est.R`.
* To modify PLGG related transitions, generate a dataframe with the same structure as the files in `02_data/deter` or `02_data/flexsurv`.

To run the cost-effectiveness model go to `Cost-effectiveness_basecase.R` and run the file line by line. Afterwards ` Manuscript results.Rmd` will summarize the results. 

---

### Project structure

This project contains the following files

```
01_data_raw: not included in this repository. 

02_data: data
  |- Utilities_CochranSE.csv: Utility estimates and coefficients 
  |- SRS_cancer.RDS: Standardized risk ratio for Cancer
  |- AE-digitized.RDS: Pseudo IPD from Digitized KM curve for adverse events incidence
  |- AE-atrisk.RDS: Numbers at risk from adverse events
  |- Rate-cancer.RDS: rate of cancer per 100000 by age 
  |- Rate-stroke.RDS: rate of stroke per 100000 by age
  |- Rate-cardio.RDS: rate of cardiovascular disease per 100000 by age
  |- life-exp-smn.RDS: rate of all cause mortality per 100000 by age excluding SN
  |- life-exp-minus.RDS: rate of all cause mortality per 100000 by age excluding SN
  |- CardioMortPar.RDS: Mortality after Cardio AE
  |- cherlow_dig.RDS: Pseudo IPD from Digitized KM curve from Cherlow et al.
  |- IPD_ae.RDS: Pseudo IPD data from Effinger et al KM 
  |- syn_data.RDS: A synthetic version of IPD (100,000 fake individuals with corelation structure of age, sex, and fusion status held the same). 
    |- bootstrap: a folder with 1000 RDS files, each RDS specifying to the flexsurv fit of the transition probabilities. Flexsurv objects have been de-removed to only have summary functions and variance and covariance matrices. 
       |- boostrapped_1_TRUE.RDS: 1st Bootstrap flexsurv fit 
       |- ...
       |- bootstrapped_1000_TRUE.RDS: 1000th Bootstrapped flexsurv fit
    |deter: a folder with 1 value for the flexsurv objects fit without bootstrapping
       
03_ancilliary-functions: Set of functions that run the economic model 
 |-Microsim.R: Core of the microsimulation model
 |- Fun_surv_est.R: functions to estimate PLGG and treatment related AE  probailities 
 |- microsim-helpers.R: Helpers for Probs and Microsimulation functions
 |- ae_functions.R: Generate probabilities for non treatment related AEs
 |- Probs.R: Assigns probabilities of individuals transitioning 
 |- Effects.R generates HrQoL estimates
 |- Costs.R generates cost estimates
 |- Radiation-Benefit.R: Generates probabilities according to benefit of Radiation therapy
 |- summary.R generates summaries of model results 
 
04_analysis: 
  |- Cost-effectiveness_basecase.R: Generates model results for base case results in manuscript. Wrapper to Microsim.R
  |- Cost-effectiveness_scenario.R: Generates model results for base case results in manuscript. Wrapper to Microsim.R
    |- output: stores output from Cost-effectiveness_basecase.R
      |- res_rad   
      |- res_model  
      |- res_OS  
      |- res_tmat  
      |- res_CI   
05_report 
  |-  Manuscript results.Rmd: Generates manuscript results from markdown file

```
---

### Required R packages 

The R code requires that you have the following packages installed.

```
here
tidyverse
pryr
formattable
knitr
scales 
survminer
survival
flexsurv
tidyr 
dfoptim
stringr 
matrixStats
```
