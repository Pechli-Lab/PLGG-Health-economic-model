## Costs function
# Create generator for Costs

gen_costs <- function(){
  # input:
  # output: a list with Inpatient_Costs.l , Outpatient_costs.l, AE_costs.l, Lookback_costs.l
  # see 06_report, Tables_out for source of values
  inflat <- 151.2/133.4 # StatsCan CPI 2022 vs. 2018
  discount <- 0 # x% discount, e.g. 20% discount = 0.2
  
  Inpatient_Costs.l  <- list(
    pre_prog = c(14955.52, 357.26)*inflat, # Cycle 1; cycle 2 until patient is 18 years old (Supp table 1)
    prog_1   = c(3876.69, 806.4)  *inflat, # Cycle 1; cycle 2 until patient is 18 years old (Supp table 1)
    prog_2   = c(7358.42, 1367.02)*inflat) # Cycle 1; cycle 2 until patient is 18 years old (Supp table 1) # Cycle 1; cycle 2 until patient is 18 years old (Supp table 1)
  
  Outpatient_costs.l <- list(
    pre_prog = c(3823.31)*inflat,  # Cycles 1-17 (Supp table 1)
    prog_1   = c(3786.90)*inflat,  # Cycles 1 to 13 (0-52 weeks after progression) (Supp table 1)
    prog_2   = c(3786.90)*inflat)  # Cycles 1 to 13 (0-52 weeks after progression) (Supp table 1)
  
  AE_cost.l <- list(Neurologic = c(30500,7400,1500)*inflat, # Cycle 1-2 post-injury; Cycles 3 to 11 after injury; Cycles 12 to (5 cycles before death) 
                    Auditory   = c(8.33)*inflat,            # Cycle 1 to death
                    Visual     = c(17.95,8.98,5.98)*inflat, # Cycles 1 to 12; Cycles 13 to 24; Cycles 25 to death
                    Stroke     = c(21735, 4933)*inflat,     # Cycle 1; Cycles 2 to 20
                    Cardiac    = c(20565,2997)*inflat,      # Cycle 1; Cycles 2 to 21
                    SMN        = c(4566.84,664.5)*inflat)   # Cycles 1 - 6; Cycles 7- (12 cycles before death)
  
  Lookback_cost.l  <- list(death_plgg  = 25960.21*inflat,   # 12 Cycles before death 
                           death_SN    = 5494.93*inflat,    # 12 cycles before death (Terminal)
                           Dx_SN       = 2842.71*inflat)    # Pre-diagnosis (3 cycles)
  
  General_pop_costs.l <- list(five        = 2383.36/12*inflat, # Supp table 1
                              ten         = 1412.68/12*inflat, 
                              fifteen     = 1486.44/12*inflat, 
                              twenty      = 1846.25/12*inflat, 
                              twentyfive  = 1955.36/12*inflat, 
                              thirty      = 2349.73/12*inflat, 
                              thirtyfive  = 2618.95/12*inflat, 
                              fourty      = 2543.75/12*inflat, 
                              fourtyfive  = 2518.44/12*inflat, 
                              fifty       = 2819.49/12*inflat, 
                              fiftyfive   = 3281.00/12*inflat, 
                              sixty       = 3948.83/12*inflat, 
                              sixtyfive   = 4885.19/12*inflat, 
                              seventy     = 6504.60/12*inflat, 
                              seventyfive = 8438.69/12*inflat, 
                              eighty      = 11313.55/12*inflat, 
                              eightyfive  = 15676.13/12*inflat, 
                              ninty       = 23267.29/12*inflat, 
                              death       = 29464.58/12*inflat)
  
  Costs_testing.l = list(FISH       = 0, # One-time cost (Supp table 1)
                         Nanostring = 0) # One-time cost (Supp table 1)
  
  Costs_trt.l = list(Targeted = 3958.42*inflat * (1-discount), # debraf
                     tram = 76.96 * 30.44 * 1.3269 * (151.2/136.0)) # monthly
  
  # 76.96/day from Petros on teams, 1.3269 is average USD to CAD exchange rate in 2019 (Bank of Canada)
  # 30.44 days a month
  # StatsCan CPI 2022 vs. 2019 
  
  return(list(Inpatient_Costs.l   = Inpatient_Costs.l,
              Outpatient_costs.l  = Outpatient_costs.l,
              AE_cost.l           = AE_cost.l,
              Lookback_cost.l     = Lookback_cost.l,
              General_pop_costs.l = General_pop_costs.l,
              Costs_testing.l     = Costs_testing.l,
              Costs_trt.l         = Costs_trt.l))  
}




# Inpatient_costs.l





Costs_plgg <- function(Inpatient_Costs,Outpatient_costs, m_Dur,curent_age_v ){
  # input: 
  #   Inpatient_costs.l: list of inpatient costs for pre_prog, prog1, prog2 and death
  #   M_dur: m.Dur 
  #   cur_age_v: current age vector

# Creating vector i will be returning
ret_costs_out <-ret_costs_in <- vector(mode = "numeric", length = dim(m_Dur)[1])

# Create index for which state they are in 
# pre_prog = 1, prog1 = 2, prog2 = 3, death 
cur_state  <- rowSums(m_Dur[,c("pre_prog","prog1","prog2")] > 0)
cur_state[rowSums(m_Dur[,c("death_plgg","death_outside","death_stroke","death_SN","death_cardio")]) > 0] <- 4    

# Finding how long the've been in the most upto date state
loc_i<- 1:dim(m_Dur)[1] +  (cur_state-1)*dim(m_Dur)[1]
t_state <-as.vector(m_Dur)[loc_i]
first_cycle <-abs(t_state - 1/12)  < (1/24)
age_log <-  curent_age_v <= 18 

# Inpatient Costs
## Assinging based on what state, whether under 18, and if first cycle
ret_costs_in[cur_state == 1 & age_log & first_cycle]  <- Inpatient_Costs$pre_prog[1]
ret_costs_in[cur_state == 1 & age_log & !first_cycle] <- Inpatient_Costs$pre_prog[2]
ret_costs_in[cur_state == 2 & age_log & first_cycle]  <- Inpatient_Costs$prog_1[1]
ret_costs_in[cur_state == 2 & age_log & !first_cycle] <- Inpatient_Costs$prog_1[2]
ret_costs_in[cur_state == 3 & age_log & first_cycle]  <- Inpatient_Costs$prog_2[1]
ret_costs_in[cur_state == 3 & age_log & !first_cycle] <- Inpatient_Costs$prog_2[2]

# Outpatinet Costs
# Each only occur for a sepcific cycle length
ret_costs_out[cur_state == 1 & age_log & m_Dur[,"pre_prog"] <= (17.01*1/12)]    <- Outpatient_costs$pre_prog
ret_costs_out[cur_state == 2 & age_log & m_Dur[,"prog1"]    <= (13.01*1/12)]    <- Outpatient_costs$prog_1
ret_costs_out[cur_state == 3 & age_log & m_Dur[,"prog2"]    <= (13.01*1/12)]    <- Outpatient_costs$prog_2

return(ret_costs_out + ret_costs_in)

}



Costs_AE <- function(m_Dur,AE_cost){
  # Each AE has differnetlogicals 
  # Neurological


# is alive 
alive_indicator <- rowSums(m_Dur[,c("death_plgg","death_outside","death_stroke","death_SN","death_cardio")]) == 0  
  
vect_ae_costs_Neurologic    <- vector(mode="numeric", length = dim(m_Dur)[1])    
vect_ae_costs_Auditory    <- vector(mode="numeric", length = dim(m_Dur)[1])    
vect_ae_costs_Visual   <- vector(mode="numeric", length = dim(m_Dur)[1])    
vect_ae_costs_Stroke   <- vector(mode="numeric", length = dim(m_Dur)[1])    
vect_ae_costs_Cardiac   <- vector(mode="numeric", length = dim(m_Dur)[1])    
vect_ae_costs_SMN   <- vector(mode="numeric", length = dim(m_Dur)[1])    
  
# Neurologic , cycles 1-2, cycles 3-11, cycles 12 
has_neurologic_cycle0102  <- m_Dur[,c("neurologic")] > 0 & m_Dur[,c("neurologic")] < (1/12)*2.01     
has_neurologic_cycle311  <- m_Dur[,c("neurologic")] > (1/12)*2.01 & m_Dur[,c("neurologic")] < (1/12)*11.01  
has_neurologic_cycle12  <- m_Dur[,c("neurologic")] > (1/12)*11.01  

vect_ae_costs_Neurologic[has_neurologic_cycle0102 & alive_indicator] <- AE_cost$Neurologic[1]
vect_ae_costs_Neurologic[has_neurologic_cycle311 & alive_indicator] <- AE_cost$Neurologic[2]
vect_ae_costs_Neurologic[has_neurologic_cycle12 & alive_indicator] <- AE_cost$Neurologic[3]

# Auditory
has_auditory <- m_Dur[,c("auditory")] > 0
vect_ae_costs_Auditory[has_auditory & alive_indicator] <- AE_cost$Auditory

# Visual
has_visual_cycle112 <- m_Dur[,c("visual")] > 0 & m_Dur[,c("visual")] < (1/12)*12.01  
has_visual_cycle1324 <- m_Dur[,c("visual")] > (1/12)*12.01   & m_Dur[,c("visual")] < (1/12)*24.01  
has_visual_cycle25 <- m_Dur[,c("visual")] > (1/12)*24.01  

vect_ae_costs_Visual[has_visual_cycle112 & alive_indicator] <- AE_cost$Visual[1]
vect_ae_costs_Visual[has_visual_cycle1324 & alive_indicator] <- AE_cost$Visual[2]
vect_ae_costs_Visual[has_visual_cycle25 & alive_indicator] <- AE_cost$Visual[3]

# Stroke
time_stroke  <- rowSums(m_Dur[,c("stroke","stroke_gen")])
has_stroke_1 <- time_stroke > 0 & time_stroke <  (1/12)*1.01
has_stroke_20 <- time_stroke > (1/12)*1.01 & time_stroke <  (1/12)*20.01

vect_ae_costs_Stroke[has_stroke_1 & alive_indicator] <- AE_cost$Stroke[1]
vect_ae_costs_Stroke[has_stroke_20 & alive_indicator] <- AE_cost$Stroke[2]

# Cardiac
time_cardiac <-rowSums(m_Dur[,c("cardiovascular_gen","cardiovascular")])
has_cardiac_1 <- time_cardiac > 0 & time_cardiac < (1/12)*1.01
has_cardiac_21 <- time_cardiac >  (1/12)*1.01 & time_cardiac < (1/12)*20.01

vect_ae_costs_Cardiac[has_cardiac_1 & alive_indicator] <- AE_cost$Cardiac[1]
vect_ae_costs_Cardiac[has_cardiac_21& alive_indicator] <- AE_cost$Cardiac[2]

# SMN
time_smn <- rowSums(m_Dur[,c("SN_gen","SN")])

has_smn6 <- time_smn >0 & time_smn < (1/12)*6.01
has_smn7p <- time_smn > (1/12)*6.01

vect_ae_costs_SMN[has_smn6 & alive_indicator] <- AE_cost$SMN[1]
vect_ae_costs_SMN[has_smn7p &  alive_indicator] <- AE_cost$SMN[2]


ret_vect <- vect_ae_costs_Neurologic + vect_ae_costs_Auditory + vect_ae_costs_Visual + vect_ae_costs_Stroke + vect_ae_costs_Cardiac + vect_ae_costs_SMN
return(ret_vect)

}



