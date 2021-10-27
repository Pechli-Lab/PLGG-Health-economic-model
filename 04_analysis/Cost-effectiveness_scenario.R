### R Script for running cost-effectiveness analysis 
require(dplyr)
require(here)
require(pryr)


# location of files, and output 
loc_here <- paste0(here::here(),"/")
loc_here_out <- paste0(here::here(),"/04_analysis/outputscenario/")

# create folders in the location of the output
dir.create(loc_here_out)
dir.create(paste0(loc_here_out,"res_rad"))
dir.create(paste0(loc_here_out,"res_model"))
dir.create(paste0(loc_here_out,"res_OS"))
dir.create(paste0(loc_here_out,"res_tmat"))
dir.create(paste0(loc_here_out,"res_CI"))

## Required functions  ----
### Microsimulation.R Core of the microsimulation model 
source(paste0(loc_here,"03_ancilliary-functions/MicrosimScenario.R"))

## Sourcing functions ----
source(paste0(loc_here, "03_ancilliary-functions/Fun_surv_est.R"))
# Fun_surv_est.R functions to estimate PLGG and treatment related AE probailities 
# EstProbs.AE: estimated AE probabilities
# EstProbs.plgg: estimates PLGG probs
# FunFlexFit: helper
# FunFlexFit_AE: helper

source(paste0(loc_here,"03_ancilliary-functions/microsim-helpers.R"))
# Helpers for Probs and Microsimulation functions
#Fun.GenI: Bootstrap wrapper to generate sample population
#fun_index: given a time from an event and returns the index of how many cycles since 0
#Probs.plggFused:  Probs helper
#Probs.plggWT:  Probs helper
#Probs.trtNoRad: Probs helper
#Probs.trtRad: Probs helper
#samplev2: used to estimate probabilities
#SurvProbFun: flexsurv objects to probabilities

source(paste0(loc_here,"03_ancilliary-functions/ae_functions.R"))
# AE non treatment related function
#AE_outside_pr: Generates matrix of probability of having various AE
#MortCardio_est: Probability of death related to Cardiovascular event
#PrDeath_SM: Probability of death realted to SMN event
#PrDeath_Stroke: Probability of death related to stroke AE
#RSR_gen: ratios of RSR ratios 

source(paste0(loc_here,"03_ancilliary-functions/Probs.R"))
# Probs: workhorse of natural history PLGG

source(paste0(loc_here,"03_ancilliary-functions/Costs.R"))
# gen_costs: generates costs
# Costs_AE: assining AE costs
# Costs_plgg: assinging PLGG related costs

source(paste0(loc_here,"03_ancilliary-functions/Effects.R"))
# gen_costs: generates costs
# Effects: assining utilities costs
# gen_Effects: generating effects
# Effects.R

source(paste0(loc_here,"03_ancilliary-functions/Radiation-Benefit.R"))
# gen_costs: generates costs
# Effects: assining utilities costs
# gen_Effects: generating effects
# Effects.R

source(paste0(loc_here, "03_ancilliary-functions/summary.R"))
# to summarize model ouptut
source(paste0(loc_here, "03_ancilliary-functions/gen_synpop.R"))
# to generate fake data sets (no correlation structure)



RSR_mat <- RSR_gen(cycles_year = 12,
                   SRS_small = readRDS(paste0(loc_here,"02_data/SRS_cancer.RDS")))


# Non-treatment related AE
AE_outside_mat1  <-  AE_outside_pr(
  Rate_cancer.RDS = readRDS(paste0(loc_here,"02_data/Rate-cancer.RDS")),
  Rate_cardio.RDS = readRDS(paste0(loc_here,"02_data/Rate-cardio.RDS")),
  Rate_stroke.RDS = readRDS(paste0(loc_here,"02_data/Rate-stroke.RDS")),
  life_exp_minus.RDS = readRDS(paste0(loc_here,"02_data/life-exp-minus.RDS")),
  life_exp_smn.RDS =  readRDS(paste0(loc_here,"02_data/life-exp-smn.RDS")),
  cycles_inyear = 12)

# Probaiblity of death aftern cardiovascular AE
PrMortCardio_mat1 <- MortCardio_est(est_par =readRDS(paste0(loc_here,"02_data/CardioMortPar.RDS")), 
                                    cycleinyear = 12)

# Costs inputs
cost_input.l <- gen_costs()

Utils <- read.csv(file = paste0(loc_here,"02_data/Utilities_CochranSE.csv"),row.names = 1)

pars <- betaPar(Utils$Wtd_Mean,Utils$St_Error)

df_beta_Utils <- data.frame( alpha = pars$a,  beta = pars$b,
                             mean = pars$a / (pars$a + pars$b),
                             row.names = rownames(Utils))

dig_radbenefit1  <- readRDS(paste0(loc_here, "02_data/cherlow_dig.RDS"))

k1 <- fun_genRadBenefit(dig_radbenefit =dig_radbenefit1, deter =T)




sim_num <- 1
  # Actual Simulation loop
  set.seed(sim_num*300) 
  
  ## load input data
  df_ae_ipd <- readRDS(paste0(loc_here,"02_data/IPD_ae.RDS"))
  
  ## Estimating model inputs ----
  # AE treatment related probabilities
  Estimated_ae1 <- EstProbs.AE(df_ae_ipd,returnAll = F ,deter = T)
  
  # PLGG related probabilities
  # Estimated_plgg1  <-EstProbs.plgg(Surv1$SurvEvent,returnAll = F)
  # Estimating radiated paitents benefit
  
  # here stuck at generating flexsurv ojbects will be done after model run 
  
   rad_benefit_flex <- F
  
  Estimated_plgg1 <- readRDS(
    paste0(loc_here,"02_data/Deter/deter_1_T.RDS"))
  Util_i <- gen_Effects(deter = T,  beta_Utils =df_beta_Utils)
  
  
  n.y_g <- 80              # time horizon, 70 years
  cyc.t_g <- 1/12           # cycle length 
  seed_n.g <- 1
  df_char_g <- gen_synpop(n1 = 10000)
  
  df_char_g$Radiation <- 0
  n.i_g <- dim(df_char_g)[1]         # number of individuals        
  
  pr_vect <- seq(0.1,1, by=  0.1)
  res_list <- vector(mode ="list", length = length(pr_vect)*2)
  index1 <- 1
   for( j in seq_along(pr_vect)){
    
     rad_benefit_flex<-    fun_genRadBenefit(dig_radbenefit =dig_radbenefit1, deter  = T)
     rad_benefit_flex <- FALSE
     
     qr1 <- Microsimulation_scenario(n.i = n.i_g,n.y = n.y_g,cyc.t = 1/12,monitor = T,seed_n = sim_num*200,
                                     df_char = df_char_g,Estimated_plgg = Estimated_plgg1,Estimated_ae = Estimated_ae1,
                                     AE_outside_mat = AE_outside_mat1,PrMortCardio_mat = PrMortCardio_mat1,
                                     SRS_mat1 = RSR_mat,cost_input = cost_input.l, rad_benefit  =  rad_benefit_flex,
                                     util_input = Util_i, uindx = 1, 
                                     d.c = 0.00246627, d.e = 0.00246627,
                                     sim_numi = 1, rri = 1,loc_out1 = loc_here_out,prob_scenario = pr_vect[j])
     res_list[[index1]] <- qr1$res_model
     index1 <- index1 + 1
     
     qr1 <- Microsimulation_scenario(n.i = n.i_g,n.y = n.y_g,cyc.t = 1/12,monitor = T,seed_n = sim_num*200,
                                     df_char = df_char_g,Estimated_plgg = Estimated_plgg1,Estimated_ae = Estimated_ae1,
                                     AE_outside_mat = AE_outside_mat1,PrMortCardio_mat = PrMortCardio_mat1,
                                     SRS_mat1 = RSR_mat,cost_input = cost_input.l, rad_benefit  =  rad_benefit_flex,
                                     util_input = Util_i, uindx = 1, 
                                     d.c = 0.00246627, d.e = 0.00246627,
                                     sim_numi = 1, rri = 0,loc_out1 = loc_here_out,prob_scenario = pr_vect[j])
     res_list[[index1]] <- qr1$res_model
     res_list[[index1]]$pr_v <- pr_vect[j]
     
     index1 <- index1 + 1 
     
     rm(qr1)
     
     ## assuming benefit of radiation
     rad_benefit_flex<-    fun_genRadBenefit(dig_radbenefit =dig_radbenefit1, deter =F)
     rad_benefit_flex$rr <- 1
     
     qr2 <-Microsimulation_scenario(n.i = n.i_g,n.y = n.y_g,cyc.t = 1/12,monitor = T,seed_n = sim_num*200,
                                    df_char = df_char_g,Estimated_plgg = Estimated_plgg1,Estimated_ae = Estimated_ae1,
                                    AE_outside_mat = AE_outside_mat1,PrMortCardio_mat = PrMortCardio_mat1,
                                    SRS_mat1 = RSR_mat,cost_input = cost_input.l, rad_benefit  =  rad_benefit_flex,
                                    util_input = Util_i, uindx = 1, 
                                    d.c = 0.00246627, d.e = 0.00246627,
                                    sim_numi = 1, rri = 1,loc_out1 = loc_here_out,prob_scenario = pr_vect[j])
     
     res_list[[index1]] <- qr2$res_model
     res_list[[index1]]$pr_v <- pr_vect[j]
     index1 <- index1 + 1 
   }
  

for(j in seq_along(res_list)){
  
  res_list[[i]] <-  res_list[[i]]  %>% 
    filter(subset %in% c("all"), Type %in% c("Total","LY","QALY"), discount == 1) %>% 
    transmute(Value, Variable, intervention, rr, model = pr_v) %>% 
    spread(intervention, Value) %>% 
    mutate(delta =Norad_Fused-Rad_Fused )
}
  
  
  res_list <-bind_rows(res_list)
  
  res_list %>% 
    mutate(model = str_remove_all(model, ".RDS"),
           k1 = str_remove_all(model, pattern = "[:alpha:]"),
           k2 = str_remove_all(model, pattern = "[:digit:]|\\.")
    ) %>% 
    select(-model, -rr) %>% 
    gather(Norad_Fused:delta, key= "int", value = "val") %>% 
    pivot_wider(names_from = c(Variable,int),values_from = val) %>% 
    mutate(k1= as.numeric(k1)) %>% 
    arrange(k2, k1)
