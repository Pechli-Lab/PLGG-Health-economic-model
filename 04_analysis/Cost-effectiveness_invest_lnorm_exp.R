### R Script for running cost-effectiveness analysis 
require(dplyr)
require(here)
require(pryr)
library(darthtools)

scen <- "1d"

# location of files, and output 
loc_here <- paste0(here::here(),"/")
loc_here_out <- paste0(here::here(),"/04_analysis/output/")

# create folders in the location of the output
dir.create(loc_here_out)
dir.create(paste0(loc_here_out,"res_rad"))
dir.create(paste0(loc_here_out,"res_model"))
dir.create(paste0(loc_here_out,"res_OS"))
dir.create(paste0(loc_here_out,"res_tmat"))
dir.create(paste0(loc_here_out,"res_CI"))
dir.create(paste0(loc_here_out,"trace"))
dir.create(paste0(loc_here_out,"AE_track"))


## Required functions  ----
### Microsimulation.R Core of the microsimulation model 
source(paste0(loc_here,"03_ancilliary-functions/Microsim.R"))

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

### Survival curves
## Targeted
# Investigator assessment
#load(paste0(here::here("04_analysis","control scenarios","Erics firstline RCT"),"/fitted_PFS_targeted_firstline_investigator_exp_100.RData")) # pre-prog to prog1
load(paste0(here::here("04_analysis","control scenarios","Erics firstline RCT"),"/fitted_PFS_targeted_firstline_investigator_lnorm_100.RData")) 

# Independent reviewer
#load(paste0(here::here("04_analysis","control scenarios","Erics firstline RCT"),"/fitted_PFS_targeted_firstline_exp_100.RData")) # pre-prog to prog1
#load(paste0(here::here("04_analysis","control scenarios","Erics firstline RCT"),"/fitted_PFS_targeted_firstline_lnorm_100.RData")) 

## Control 
# Investigator assessment
load(paste0(here::here("04_analysis","control scenarios","Erics firstline RCT"),"/fitted_PFS_control_firstline_investigator_exp_100.RData")) # pre-prog to prog1
#load(paste0(here::here("04_analysis","control scenarios","Erics firstline RCT"),"/fitted_PFS_control_firstline_investigator_lnorm_100.RData")) # pre-prog to prog1

# Independent reviewer
#load(paste0(here::here("04_analysis","control scenarios","Erics firstline RCT"),"/fitted_PFS_control_firstline_lnorm_100.RData")) # pre-prog to prog1

# SickKids (control arm, match with only independent reviewer targeted)
#load(paste0(here::here("04_analysis","control scenarios","Bryans SickKids data firstline"),"/fitted_PFS_control_sk_firstline_exp_100.RData")) # pre-prog to prog 1     


# Prog 1 to prog 2
load(paste0(here::here("04_analysis","control scenarios","Bryans SickKids data"),"/fitted_PFS_control_sk_exp_100.RData")) # prog1 to prog2

# plgg OS adjusted
load(paste0(here::here("literature"),"/fitted_plgg_os_adjust_exp_100.RData")) # prog1 to death


## Global Parameters
N_sim <- 50
torun1 <- 1:N_sim

start.time <- Sys.time()
for(sim_num1 in torun1){
  
  sim_num <- sim_num1
  # Actual Simulation loop
  set.seed(sim_num*300) 
  
  ## load input data
  df_ae_ipd <- readRDS(paste0(loc_here,"02_data/IPD_ae.RDS"))
  
  ## Estimating model inputs ----
  # AE treatment related probabilities
  Estimated_ae1 <- EstProbs.AE(df_ae_ipd,returnAll = F ,deter = F)
  
  # PLGG related probabilities
  # Estimated_plgg1  <-EstProbs.plgg(Surv1$SurvEvent,returnAll = F)
  # Estimating radiated paitents benefit
  
  # here stuck at generating flexsurv ojbects will be done after model run 
  rad_benefit_flex    <- fun_genRadBenefit(dig_radbenefit = dig_radbenefit1, deter = F)
  rad_benefit_flex$rr <- 1
  #rad_benefit_flex <- F
  
  Estimated_plgg1 <- readRDS(
    paste0(loc_here,"02_data/bootstrap/bootstrapped_",sim_num1,"_TRUE.RDS"))
  Util_i <- gen_Effects(deter = F,  beta_Utils =df_beta_Utils)
  
  # Bootstrap sampling with replacement from Liana's paper's PFS swimmer plot data
  Estimated_plgg1$flexobj[[2]][1][[1]]  <- fitted_PFS_targeted[[sim_num1]]  # preprog to prog1, targeted
  Estimated_plgg1 <- rbind(Estimated_plgg1, Estimated_plgg1[1,])
  Estimated_plgg1$flexobj[[11]][1][[1]] <- fitted_PFS_chemo[[sim_num1]]     # preprog to prog1, control
  Estimated_plgg1$flexobj[[10]][1][[1]] <- fitted_PFS_chemo_sk[[sim_num1]]  # prog1 to prog2, both arms
  Estimated_plgg1$flexobj[[4]][1][[1]]  <- fitted_os_adjusted[[sim_num1]]   # prog1 to death, both arms
  
  n.y_g     <- 80             # time horizon, 80 years
  cyc.t_g   <- 1/12           # cycle length 
  seed_n.g  <- 1
  df_char_g <- gen_synpop(n1 = 10000)
  
  df_char_g$Radiation <- 0
  n.i_g <- dim(df_char_g)[1]         # number of individuals        
  
  ## assuming no benefit to radiation 
  qr1 <- Microsimulation(n.i              = n.i_g,
                         n.y              = n.y_g,
                         cyc.t            = 1/12,
                         monitor          = T,
                         seed_n           = sim_num*200,
                         df_char          = df_char_g,
                         Estimated_plgg   = Estimated_plgg1,
                         Estimated_ae     = Estimated_ae1,
                         AE_outside_mat   = AE_outside_mat1,
                         PrMortCardio_mat = PrMortCardio_mat1,
                         SRS_mat1         = RSR_mat,
                         cost_input       = cost_input.l, 
                         rad_benefit      = rad_benefit_flex,
                         util_input       = Util_i, 
                         uindx            = 1, 
                         d.c              = 0.015/12, # 1.5% annual discount rate, to monthly
                         d.e              = 0.015/12,
                         sim_numi         = sim_num, 
                         rri              = 0,
                         loc_out1         = loc_here_out,
                         trt_duration     = 2,
                         combo            = TRUE)
  
  # # # Testing
  # n.i              = n.i_g;
  # n.y              = n.y_g;
  # cyc.t            = 1/12;
  # monitor          = T;
  # seed_n           = sim_num*200;
  # df_char          = df_char_g;
  # Estimated_plgg   = Estimated_plgg1;
  # Estimated_ae     = Estimated_ae1;
  # AE_outside_mat   = AE_outside_mat1;
  # PrMortCardio_mat = PrMortCardio_mat1;
  # SRS_mat1         = RSR_mat;
  # cost_input       = cost_input.l;
  # rad_benefit      = rad_benefit_flex;
  # util_input       = Util_i
  # uindx            = 1;
  # d.c              = 0.015/12; # 1.5% annual discount rate; to monthly
  # d.e              = 0.015/12;
  # sim_numi         = sim_num;
  # rri              = 0;
  # loc_out1         = loc_here_out;
  # trt_duration     = 1;
  # combo            = T;
  
  message("finished assuming no radiation benefit")
  
  for( i in seq_along(qr1)){
    #saveRDS(qr1[[i]], paste0(loc_here_out, names(qr1)[i], "/", "model_", sim_num, "_RR_0.RDS" ))
    saveRDS(qr1[[i]], paste0(loc_here_out, "scenarios/", scen, "/basecase lnorm exp lifetime investigator combo/", names(qr1)[i], "/", "model_", sim_num, "_RR_0.RDS" ))
  }

  rm(qr1)
  ## assuming benefit of radiation
  # rad_benefit_flex<-    fun_genRadBenefit(dig_radbenefit =dig_radbenefit1, deter =F)
  # rad_benefit_flex$rr <- 1
  # 
  # qr2 <- Microsimulation(n.i = n.i_g,n.y = n.y_g,cyc.t = 1/12,monitor = F,seed_n = sim_num*200,
  #                        df_char = df_char_g,Estimated_plgg = Estimated_plgg1,Estimated_ae = Estimated_ae1,
  #                        AE_outside_mat = AE_outside_mat1,PrMortCardio_mat = PrMortCardio_mat1,
  #                        SRS_mat1 = RSR_mat,cost_input = cost_input.l, rad_benefit  =  rad_benefit_flex,
  #                        util_input = Util_i, uindx = 1, 
  #                        d.c = 0.00246627, d.e = 0.00246627,
  #                        sim_numi = sim_num, rri = 1,loc_out1 = loc_here_out )
  # 
  # for( i in seq_along(qr2)){
  #   saveRDS(qr2[[i]], paste0(loc_here_out, names(qr2)[i], "/", "model_", sim_num, "_RR_0.RDS" ))
  # }
  # rm(qr2)
  # message(paste0( "run:",  floor(sim_num / 100), "percent:"  ,round(sim_num1/100, digits = 2)))
  
}
Time.elapsed <- Sys.time()-start.time





