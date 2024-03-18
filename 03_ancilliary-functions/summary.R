### Summary function

# function of the following inputs
sum_function <- function(sim_num, intervention1, rr1, d_rate,mM, mDur, mVis,
                         costsplgg, costsAE, costsgen, coststrt, costsrad, util, vRAD, vFusion  ){
  # What simulation number,
  # Parameter on radiation benefit
  # Intervention
  #sim_num = 1
  #intervention1 = "Intervention"
  #rr1 = 0
  # d_rate <- 0.00246627
  v_dwc <- 1 / (1 + d_rate) ^ (0:(dim(mM)[2]-1)) # calculate the cost discount weight based on the discount rate d.c 
  ## Requires
  require(dplyr)
  require(stringr)
  
  # output
  
  res_ly <- vector(mode = "numeric", length = 3)
  type   <- c("all","fused","fused_rad")
  
  ## Matrix of model inputs
  # LY 
  #   All
  res_ly[[1]] <- mean(rowSums(util == 0)*cyc.t)
  #   Fused
  res_ly[[2]] <- mean(rowSums(util[vFusion == 1,] == 0)*cyc.t)
  #   Fused and Prog 1 
  res_ly[[3]] <- mean(rowSums(util[vFusion == 1 & mVis[,2],] == 0)*cyc.t)
  
  df_LY       <- tibble(Value = res_ly,  Variable = "LY", Type = "LY" ,  subset = type , discount = 0)
  
  res_qaly    <- vector(mode = "numeric", length = 6)
  
  # QALY 
  #   All
  res_qaly[[1]] <- mean(rowSums(util*cyc.t*v_dwc))
  #   Fused
  res_qaly[[2]] <- mean(rowSums(util[vFusion == 1,]*cyc.t*v_dwc))
  #   Fused and Prog 1 
  res_qaly[[3]] <- mean(rowSums(util[vFusion == 1 & mVis[,2],]*cyc.t*v_dwc))
  
  # QALY - undiscounted 
  #   All
  res_qaly[[4]] <- mean(rowSums(util*cyc.t))
  #   Fused
  res_qaly[[5]] <- mean(rowSums(util[vFusion == 1,]*cyc.t))
  #   Fused and Prog 1 
  res_qaly[[6]] <- mean(rowSums(util[vFusion == 1 & mVis[,2],]*cyc.t))
  
  df_QALY <- tibble(Value    = res_qaly, 
                    Variable = "QALY", 
                    Type     = "QALY",
                    subset   = rep(type, 2), 
                    discount = c(rep(1, times = 3), rep(0, times = 3)))
  
  ### Costs -----
  costs_plgg <- vector(mode = "numeric", length = 3)
  
  # Costs PLGG 
  #   All 
  costs_plgg[[1]] <- mean(rowSums(costsplgg*v_dwc))
  #   Fused
  costs_plgg[[2]] <- mean(rowSums(costsplgg[vFusion == 1,]*v_dwc))
  #   Fused and Prog 1 
  costs_plgg[[3]] <- mean(rowSums(costsplgg[vFusion == 1 & mVis[,2],]*v_dwc))
  
  
  df_costs_plgg <- tibble(Value    = costs_plgg, 
                          Variable = "Cost", 
                          Type     = "PLGG",  
                          subset   = unlist(type), 
                          discount = 1)
  
  # Costs AE 
  costs_ae <- vector(mode = "numeric", length = 3)
  #   All 
  costs_ae[[1]] <-mean(rowSums(costsAE*v_dwc))
  #   Fused
  costs_ae[[2]] <- mean(rowSums(costsAE[vFusion == 1,]*v_dwc))
  #   Fused and Prog 1 
  costs_ae[[3]] <- mean(rowSums(costsAE[vFusion == 1 & mVis[,2],]*v_dwc))
  
  df_costs_ae   <- tibble(Value    = costs_ae, 
                          Variable = "Cost", 
                          Type     = "AE",  
                          subset   = unlist(type), 
                          discount = 1)
  
  # Costs General
  costs_gen <- vector(mode = "numeric", length = 3)
  #   All 
  costs_gen[[1]] <-mean(rowSums(costsgen*v_dwc))
  #   Fused
  costs_gen[[2]] <- mean(rowSums(costsgen[vFusion == 1,]*v_dwc))
  #   Fused and Prog 1 
  costs_gen[[3]] <- mean(rowSums(costsgen[vFusion == 1 & mVis[,2],]*v_dwc))
  
  df_costs_gen   <- tibble(Value    = costs_gen, 
                          Variable = "Cost", 
                          Type     = "gen",  
                          subset   = unlist(type), 
                          discount = 1)
  
  # Costs treatment 
  costs_trt <- vector(mode = "numeric", length = 3)
  #   All 
  costs_trt[[1]] <-mean(rowSums(coststrt*v_dwc))
  #   Fused
  costs_trt[[2]] <- mean(rowSums(coststrt[vFusion == 1,]*v_dwc))
  #   Fused and Prog 1 
  costs_trt[[3]] <- mean(rowSums(coststrt[vFusion == 1 & mVis[,2],]*v_dwc))
  
  df_costs_trt   <- tibble(Value    = costs_trt, 
                           Variable = "Cost", 
                           Type     = "Treatment",  
                           subset   = unlist(type), 
                           discount = 1)
  
  # Costs radiation 
  costs_rad <- vector(mode = "numeric", length = 3)
  #   All 
  costs_rad[[1]] <-mean(rowSums(costsrad*v_dwc))
  #   Fused
  costs_rad[[2]] <- mean(rowSums(costsrad[vFusion == 1,]*v_dwc))
  #   Fused and Prog 1 
  costs_rad[[3]] <- mean(rowSums(costsrad[vFusion == 1 & mVis[,2],]*v_dwc))
  
  df_costs_rad   <- tibble(Value    = costs_rad, 
                           Variable = "Cost", 
                           Type     = "Radiation",  
                           subset   = unlist(type), 
                           discount = 1)
  
  # Costs Total
  #   All 
  costs_total <- costs_ae + costs_gen + costs_plgg + costs_trt + costs_rad
  
  df_costs_total <- tibble(Value    = costs_total, 
                           Variable = "Cost", 
                           Type     = "Total",  
                           subset   = unlist(type),
                           discount = 1)
  
  # Costs Discounted
  
  costs_plgg_undis <- vector(mode = "numeric", length = 3)
  
  #   All 
  costs_plgg_undis[[1]] <-mean(rowSums(costsplgg))
  #   Fused
  costs_plgg_undis[[2]] <- mean(rowSums(costsplgg[vFusion == 1,]))
  #   Fused and Prog 1 
  costs_plgg_undis[[3]] <- mean(rowSums(costsplgg[vFusion == 1 & mVis[,2],]))
  
  
  df_costs_plgg_undis <- tibble(Value    = costs_plgg_undis, 
                                Variable = "Cost", 
                                Type     = "PLGG" ,  
                                subset   = unlist(type), 
                                discount = 0)
  
  costs_ae_undis<- vector(mode = "numeric", length = 3)
  
  # Costs AE 
  #   All 
  costs_ae_undis[[1]] <-mean(rowSums(costsAE))
  #   Fused
  costs_ae_undis[[2]] <- mean(rowSums(costsAE[vFusion == 1,]))
  #   Fused and Prog 1 
  costs_ae_undis[[3]] <- mean(rowSums(costsAE[vFusion == 1 & mVis[,2],]))
  
  df_costs_ae_undis <- tibble(Value    = costs_ae_undis, 
                              Variable = "Cost", 
                              Type     = "AE",  
                              subset   = unlist(type), 
                              discount = 0)
  
  costs_gen_undis<- vector(mode = "numeric", length = 3)
  
  # Costs general 
  #   All 
  costs_gen_undis[[1]] <-mean(rowSums(costsgen))
  #   Fused
  costs_gen_undis[[2]] <- mean(rowSums(costsgen[vFusion == 1,]))
  #   Fused and Prog 1 
  costs_gen_undis[[3]] <- mean(rowSums(costsgen[vFusion == 1 & mVis[,2],]))
  
  df_costs_gen_undis <- tibble(Value    = costs_gen_undis, 
                              Variable = "Cost", 
                              Type     = "gen",  
                              subset   = unlist(type), 
                              discount = 0)
  
  costs_trt_undis<- vector(mode = "numeric", length = 3)
  
  # Costs trt 
  #   All 
  costs_trt_undis[[1]] <-mean(rowSums(coststrt))
  #   Fused
  costs_trt_undis[[2]] <- mean(rowSums(coststrt[vFusion == 1,]))
  #   Fused and Prog 1 
  costs_trt_undis[[3]] <- mean(rowSums(coststrt[vFusion == 1 & mVis[,2],]))
  
  df_costs_trt_undis <- tibble(Value    = costs_trt_undis, 
                               Variable = "Cost", 
                               Type     = "Treatment",  
                               subset   = unlist(type), 
                               discount = 0)
  
  costs_rad_undis <- vector(mode = "numeric", length = 3)
  
  # Costs rad 
  #   All 
  costs_rad_undis[[1]] <-mean(rowSums(costsrad))
  #   Fused
  costs_rad_undis[[2]] <- mean(rowSums(costsrad[vFusion == 1,]))
  #   Fused and Prog 1 
  costs_rad_undis[[3]] <- mean(rowSums(costsrad[vFusion == 1 & mVis[,2],]))
  
  df_costs_rad_undis <- tibble(Value    = costs_rad_undis, 
                               Variable = "Cost", 
                               Type     = "Radiation",  
                               subset   = unlist(type), 
                               discount = 0)
  
  costs_total_undis <- costs_ae_undis + costs_gen_undis + costs_plgg_undis + costs_trt_undis + costs_rad_undis
  # Costs Total
  #   All 
  df_costs_total_undis <- tibble(Value    = costs_total_undis, 
                                 Variable = "Cost", 
                                 Type     = "Total" ,  
                                 subset   = unlist(type), 
                                 discount = 0)
  
  df_Costs <- rbind(df_costs_ae,
                    df_costs_gen,
                    df_costs_plgg,
                    df_costs_trt, 
                    df_costs_rad,
                    df_costs_total,
                    df_costs_ae_undis, 
                    df_costs_gen_undis,
                    df_costs_plgg_undis,
                    df_costs_trt_undis,
                    df_costs_rad_undis,
                    df_costs_total_undis)
  
  ### Model results
  res_model  <- rbind(df_Costs, df_QALY, df_LY)
  
  ## OS -----
  
  # OS all
  OS_all        <- tibble( OS    = colMeans(util != 0),
                           cycle = seq(0, n.y_g, by = cyc.t),
                           type  = "OS")
  # OS fused
  OS_fused      <- tibble( OS    = colMeans(util[vFusion == 1,] != 0),
                           cycle = seq(0, n.y_g, by = cyc.t),
                           type  = "Fused")
  # OS Fused and Progressed
  OS_fused_prog <- tibble( OS    = colMeans(util[vFusion == 1 & mVis[,2],] != 0),
                           cycle = seq(0, n.y_g, by = cyc.t),
                           type  = "FusedProg")
  
  res_OS <-  rbind(OS_all,OS_fused, OS_fused_prog)
  

## Time in State figure

pos_states <- colnames(mVis)

pre_prog           <- colMeans(mM == "pre_prog")
prog1              <- colMeans(mM == "prog1")
prog2              <- colMeans(mM == "prog2")
death_plgg         <- colMeans(mM == "death_plgg")
neurologic         <- colMeans(mM == "neurologic")
auditory           <- colMeans(mM == "auditory")
visual             <- colMeans(mM == "visual")
stroke             <- colMeans(mM == "stroke")
cardiovascular     <- colMeans(mM == "cardiovascular")
SN                 <- colMeans(mM == "SN")
stroke_gen         <- colMeans(mM == "stroke_gen")
cardiovascular_gen <- colMeans(mM == "cardiovascular_gen")
SN_gen             <- colMeans(mM == "SN_gen")
death_outside      <- colMeans(mM == "death_outside")
death_stroke       <- colMeans(mM == "death_stroke")
death_SN           <- colMeans(mM == "death_SN")
death_cardio       <- colMeans(mM == "death_cardio")

res_tmat <- tibble( per_state = c( pre_prog, 
                                   prog1, 
                                   prog2, 
                                   death_plgg, 
                                   neurologic, 
                                   auditory, 
                                   visual, 
                                   stroke, 
                                   cardiovascular, 
                                   SN, 
                                   stroke_gen, 
                                   cardiovascular_gen, 
                                   SN_gen, 
                                   death_outside, 
                                   death_stroke,
                                   death_SN,
                                   death_cardio),      
                    cycle = rep( seq(0, n.y_g, by = cyc.t), times = length(pos_states)),
                    state = unlist(lapply(pos_states, function(x){rep(x, times = dim(mM)[2])})))



## Cumulative incidence accounting for competing risk
# ci_states <-pos_states[!str_detect(pos_states, "death")]
ci_states <- pos_states
res_vect  <- vector(mode = "list", length = length(ci_states))
time_vect <- vector(mode = "list", length = length(ci_states))
per_alive <- colMeans(util != 0)
time_v    <- seq(0, n.y_g, by = cyc.t)
tol       <- cyc.t/10
mDur1     <- (n.y_g+cyc.t)-mDur 

for(j in seq_along(res_vect)){
  k1 <- vector( mode =  "numeric", length = length(time_v))  
  
  for(t1 in seq_along(k1)){
    
    if(str_detect(ci_states[j], "death")){
      k1[t1] <- mean(abs(mDur1[,ci_states[j]] - time_v[t1]) < tol)
      
      
    } else{
      k1[t1] <- mean(abs(mDur1[,ci_states[j]] - time_v[t1]) < tol)/per_alive[t1]
    }
  }
  res_vect[[j]]  <- cumsum(k1)  
  time_vect[[j]] <- time_v
  
}

res_CI <-tibble( 
  ci    = unlist(res_vect),
  time  = unlist(time_vect),
  state = unlist(lapply(ci_states, function(x){ rep(x, times = length(time_v)) })))

## Number of events
res_rad <- tibble(rad          = sum(vRAD),
                  Type         = "Number Radiated",
                  sim          = sim_num, 
                  intervention = intervention1, 
                  rr           = rr1) 

res_CI    <- mutate(res_CI,    sim = sim_num, intervention = intervention1, rr = rr1)
res_model <- mutate(res_model, sim = sim_num, intervention = intervention1, rr = rr1)
res_OS    <- mutate(res_OS,    sim = sim_num, intervention = intervention1, rr = rr1)
res_tmat  <- mutate(res_tmat,  sim = sim_num, intervention = intervention1, rr = rr1)

return(list(res_rad   = res_rad, 
            res_CI    = res_CI,
            res_model = res_model,
            res_OS    = res_OS,
            res_tmat  = res_tmat,
            util      = util
            
))
}
