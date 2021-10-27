# MicrosimScenario.R
Microsimulation_scenario <- function(n.i, n.y, cyc.t, monitor,seed_n, df_char,Estimated_plgg,
                                     Estimated_ae,AE_outside_mat,PrMortCardio_mat, SRS_mat1, cost_input,
                                     util_input,uindx, rad_benefit, d.c,d.e,sim_numi, rri, loc_out1, prob_scenario ){
  # input:  
  # n.i: number of individuals
  # n.y: time horizon
  # cyc.t: cycle length
  # monitor: Logical indicating whether to monitor 
  # seed_n: seed number
  # df_char: a data frame of characteristics including a column Radiation
  # Estimated_plgg: output from EstProbs.plgg
  # Estimated_ae: output from EstProbs.ae
  # AE_outside_mat: output from AE_outside_pr
  # PrMortCardio_mat: output from MortCardio_est
  # SRS_mat1: output from SRS_gen
  # output: 
  # list of model run outcomes with the same parameters but varying intervention
  
  # This function runs the microsimulation and is structured in a similar manner as the 
  # DARTH Microsimualtion format. It consists of X sections:
  # 1. Model inputs: brining into function evironment model inputs
  # 2. Initializing: creating matrices to store simulation results
  # 3. Initial cycle: information on status at cycle 1 and costs and utilities associated 
  # 4. Cycle 2- End: Modelling natural history, costs, utilities, and radiation decision
  # 5. Summarizing: Summarizing model results 
  # Sections 2-5, are done for intervention and control  
  
  
  #1. Model inputs -----
  # Due to scoping issues have to create some variables as global variables.  
  n.i <<- n.i      # number of individuals
  cyc.t <<- cyc.t  # cycle length (1/12)
  n.y <<- n.y      # max cycle
  v.n  <<-  c("pre_prog", "prog1","prog2","death_plgg",  # PLGG related
              "neurologic","auditory","visual","stroke", # Treatment related AE
              "cardiovascular","SN", 
              "stroke_gen", "cardiovascular_gen","SN_gen", # General population risk of AE
              "death_outside","death_stroke","death_SN","death_cardio") # Mortality events
  
  n.t <<-  n.y/cyc.t      # number of cycles
  n.s <<- length(v.n)     # the number of states
  v.M_1 <- rep("pre_prog", n.i)  # Everyone begins in the itial pre-progression state(1)
  v.time <- seq(from = cyc.t, to = n.y , by = cyc.t) # time vector
  
  v.dwc <- 1 / (1 + d.c) ^ (0:n.t)   # calculate the cost discount weight based on the discount rate d.c 
  v.dwe <- 1 / (1 + d.e) ^ (0:n.t)   # calculate the QALY discount weight based on the discount rate d.e
  
  # Creating costs inputs into global environment
  for(i in seq_along(cost_input)){ 
    assign(names(cost_input)[i], value = cost_input[[i]])
  }
  
  # This for loop is for interventions Rad_fused (fused patients are radiated)
  # Norad_Fused (fused patients are not radiated)
  
  
  
  for( intervention in c("Rad_Fused","Norad_Fused")){  
    
    # Assinging a logical to whether fused patients can be radiated  
    if(intervention == "Rad_Fused"){
      FusedNoRad <- FALSE
    } else {
      FusedNoRad <- TRUE
    }
    
    # Reseeting radiation account
    
    df_char$Radiation <- 0 
    
    # 2. Initializing Values -----       
    #   m.Utilities (utilities) and m.Costs.plgg, m.Costs.ae
    m.Utilities <- m.Costs.ae <-  m.Costs.plgg <- matrix(nrow = n.i, 
                                                         ncol = n.t +1,
                                                         dimnames = list(paste("ind", 1:n.i, sep = " "), 
                                                                         paste( 0:n.t)),data = 0)
    # m.M (state tracking)
    m.M <- matrix(nrow = n.i, 
                  ncol = n.t +1,
                  dimnames = list(paste("ind", 1:n.i, sep = " "), 
                                  paste(0:n.t)),data = "")
    
    # Creating a m.Dur counts how long you've been in a state
    # m.Vis counts if have been in a state
    
    m.Dur <- matrix(nrow = n.i,
                    ncol = n.s, 
                    dimnames = list(paste("ind", 1:n.i, sep = " "), 
                                    v.n), data = 0)
    m.Vis <- matrix(nrow = n.i,
                    ncol = n.s,
                    dimnames = list(paste("ind", 1:n.i, sep = " "), v.n), data = F) 
    
    
    # Creating curent_age vector (Aware of spelling mistake but propogated through sub functions)
    # which is updated with the current age
    curent_age <- df_char$AgeDx - df_char$AgeDx %% cyc.t
    
    
    # 3. Initial cycle ----
    
    # Start everyone into the pre_progression state that they've been there for one cycle
    # and have visited pre_progression
    m.M[, 1] <- v.M_1
    m.Dur[,1] <- cyc.t
    m.Vis[,1] <-  T
    
    ## First cylce costs
    m.Costs.plgg[,1]  <- Costs_plgg(Inpatient_Costs = Inpatient_Costs.l,
                                    Outpatient_costs = Outpatient_costs.l,m_Dur = m.Dur,curent_age_v = curent_age) 
    
    # Utilities first cycle
    m.Utilities[,1] <- Effects(lst.Util = util_input,cur_age_v = curent_age,m_V = m.Vis)
    
    
    # If testing (Noread_fused) don't radiate 
    if(intervention == "Norad_Fused"){
      m.Costs.plgg[,1]   <- m.Costs.plgg[,1]   +  Costs_testing.l$Nanostring
    }
    
    # Update age 
    curent_age <- curent_age + cyc.t
    set.seed(seed_n)
    
    # 4. Cycle 2 - End ----
    for (t in 1:n.t) {  # for t in 1 number of cycles
      #if(t == 100){debug(Probs)}
      # if(any(m.Vis[,2] & df_char$Radiation == 1)){debug(Probs)}
      # Create m.p: transition matrices for next period  
      m.p <- Probs_scenario(M_it = m.M[,t], D_m = m.Dur, df.i = df_char, Vis_m = m.Vis ,CycleLength = cyc.t,
                            GlobT = v.time[t], EstPLGG = Estimated_plgg ,EstAE = Estimated_ae,
                            AE_outside_mat1 = AE_outside_mat,PrMortCardio_mat1 = PrMortCardio_mat,rad_benefit1 = rad_benefit,
                            SRS_mat = SRS_mat1,curent_age_v = curent_age,  prob_scenario1= prob_scenario)
      
      # if(t == 100){undebug(Probs)}
      
      
      # Tests to ensure m.p is working correctly
      if(any(is.nan(m.p)) | any(is.na(m.p))){ stop("NAN or NA in M.P")}
      if(any((m.p) < 0  )){ stop("less than zero")}
      
      # Use samplev to sample from the transition matrices and assign the next state
      m.M[, t + 1] <- samplev2( prob = m.p, m = 1)  
      # Tests to ensure that samplev is working 
      if(any(m.M[, t + 1] == "cardiovascular_gen" & m.M[, t ] == "cardiovascular")){ 
        stop("returning to state")}
      
      # updating if i've visited a state the Vis and Dur matrices
      for(j in v.n){
        m.Vis[ m.M[, t + 1] == j ,j] <- T
        m.Dur[ m.Vis[ , j ],j] <- m.Dur[  m.Vis[ , j ] , j ] + cyc.t
      }
      
      # Creating radiation decision for those that progress
      RadDec.patients <-m.M[, t + 1] == "prog1" & m.M[, t] != "prog1"
      
      if(any(RadDec.patients)){
        if(FusedNoRad){ # Which intervention 
          
          AgeLog   <-  curent_age[RadDec.patients]   > 8
          FusionLog <- df_char[RadDec.patients,"Fusion"] == 0     
          df_char$Radiation[RadDec.patients][AgeLog & FusionLog] <-  1 
          
        } else if (!FusedNoRad){
          
          AgeLog   <- curent_age[RadDec.patients] > 8
          df_char$Radiation[RadDec.patients][AgeLog] <-  1    
          
        }
        
        
        ## Costs treatment        
        m.Costs.plgg[, t +1 ] <- Costs_plgg(Inpatient_Costs = Inpatient_Costs.l,
                                            Outpatient_costs = Outpatient_costs.l,
                                            m_Dur = m.Dur,curent_age_v = curent_age)
        
        
        m.Costs.ae[, t +1 ] <- Costs_AE(AE_cost = AE_cost.l,m_Dur = m.Dur)
        
        
        ##  Costs look back
        # if dying from plgg related $ 25,960.21 for 12 cycles prior to
        
        death_plgg_lookback  <- m.M[ , t +1  ] == "death_plgg" & m.M[, t ] != "death_plgg"
        look_back_plgg <- seq( to = t +1 , by = 1, length.out = 12)[seq( to = t +1, by = 1, length.out = 12) > 0] 
        m.Costs.plgg[ death_plgg_lookback , look_back_plgg ] <- Lookback_cost.l$death_plgg
        
        
        # if dying SMN 12 cycles prior to assigning cost of death_SMN 
        death_sn_lookback  <- m.M[ , t +1  ] == "death_SN" & m.M[, t ] != "death_SN"
        
        look_back_sn <- seq( to = t +1 , by = 1, length.out = 12)[seq( to =  t +1, by = 1, length.out = 12) > 0]
        for(j in which(death_sn_lookback)){
          look_back_sn <- look_back_sn[look_back_sn*cyc.t <= sum(m.Dur[j,c("SN_gen","SN")])] 
          m.Costs.ae[ j , look_back_sn ] <- m.Costs.ae[ j , look_back_sn ]  + Lookback_cost.l$death_SN - AE_cost.l$SMN[2]
          
        }
        
        # If diagnosed with SMN prior three periods assinging 
        new_mn <- m.M[ , t +1  ] %in% c("SN","SN_gen") & rowSums(m.Dur[,c("SN","SN_gen")] < cyc.t +0.01) 
        look_back_newmn <- seq( to = t , by = 1, length.out = 3)[seq( to = t , by = 1, length.out = 3) > 0]
        m.Costs.ae[ new_mn , look_back_newmn ] <- m.Costs.ae[ new_mn , look_back_newmn ]  +  Lookback_cost.l$Dx_SN
      }
      
      # General population costs
      # Assinged when alive and no costs
      alive_indc <-!(m.M[, t+ 1] %in% c("death_plgg","death_SN","death_outside","death_stroke","death_SN","death_cardio"))   
      nocost <-(m.Costs.ae[,t +1 ]  + m.Costs.plgg[,t +1 ]) == 0
      
      m.Costs.ae[alive_indc & nocost & curent_age < 5,t+1] <- General_pop_costs.l$five
      m.Costs.ae[alive_indc & nocost & (curent_age >= 5 & curent_age < 10),t+1] <- General_pop_costs.l$ten
      m.Costs.ae[alive_indc & nocost & (curent_age >= 10 & curent_age < 15),t+1] <- General_pop_costs.l$fifteen
      m.Costs.ae[alive_indc & nocost & (curent_age >= 15 & curent_age < 20),t+1] <- General_pop_costs.l$twenty
      m.Costs.ae[alive_indc & nocost & (curent_age >= 20 & curent_age < 25),t+1] <- General_pop_costs.l$twentyfive
      m.Costs.ae[alive_indc & nocost & (curent_age >= 25 & curent_age < 30),t+1] <- General_pop_costs.l$thirty
      m.Costs.ae[alive_indc & nocost & (curent_age >= 30 & curent_age < 35),t+1] <- General_pop_costs.l$thirtyfive
      m.Costs.ae[alive_indc & nocost & (curent_age >= 35 & curent_age < 40),t+1] <- General_pop_costs.l$fourty
      m.Costs.ae[alive_indc & nocost & (curent_age >= 40 & curent_age < 45),t+1] <- General_pop_costs.l$fourtyfive
      m.Costs.ae[alive_indc & nocost & (curent_age >= 45 & curent_age < 50),t+1] <- General_pop_costs.l$fifty
      m.Costs.ae[alive_indc & nocost & (curent_age >= 50 & curent_age < 55),t+1] <- General_pop_costs.l$fiftyfive
      m.Costs.ae[alive_indc & nocost & (curent_age >= 55 & curent_age < 60),t+1] <- General_pop_costs.l$sixty
      m.Costs.ae[alive_indc & nocost & (curent_age >= 60 & curent_age < 65),t+1] <- General_pop_costs.l$sixtyfive
      m.Costs.ae[alive_indc & nocost & (curent_age >= 65 & curent_age < 70),t+1] <- General_pop_costs.l$seventy
      m.Costs.ae[alive_indc & nocost & (curent_age >= 70 & curent_age < 75),t+1] <- General_pop_costs.l$seventyfive
      m.Costs.ae[alive_indc & nocost & (curent_age >= 75 & curent_age < 80),t+1] <- General_pop_costs.l$eighty
      m.Costs.ae[alive_indc & nocost & (curent_age >= 80 & curent_age < 85),t+1] <- General_pop_costs.l$eightyfive
      m.Costs.ae[alive_indc & nocost & (curent_age >= 85 & curent_age < 90),t+1] <- General_pop_costs.l$ninty
      m.Costs.ae[alive_indc & nocost & (curent_age >= 90 ),t+1] <- General_pop_costs.l$death
      
      
      
      # Utilities
      m.Utilities[,t +1 ] <-  Effects(lst.Util = util_input,cur_age_v = curent_age,m_V = m.Vis)
      
      if(t %% (n.t/10) == 0 & monitor) {
        message(paste("\r", round(t/(n.t), digits = 2),"%", "sim: ",sim_numi,print(pryr::mem_used()) ))}
      #message(t)
      
      
      # Update age 
      curent_age <- curent_age + cyc.t
    } # time loop ends  
    
    
    # Creating an event matrix 
    
    #trying to create event
    # globT <- n.t*cyc.t 
    # 
    # ind_v <- rep(1:n.i, times = length(v.n))
    # event1  <-unlist(lapply(v.n,   FUN = function(x){rep(x, times = n.i ) })) 
    # m.DurV <- globT + cyc.t -  as.vector(m.Dur)   
    # 
    # event_df  <-data.frame(id = ind_v, time = m.DurV, state = event1)
    
    
    t1 <- sum_function(sim_num = sim_numi,
                       intervention1 = intervention,
                       rr1 = rri,d_rate = d.c,mM = m.M,mDur = m.Dur,mVis = m.Vis,
                       costsplgg = m.Costs.plgg,costsAE = m.Costs.ae,util = m.Utilities,
                       vRAD = df_char$Radiation, vFusion = df_char$Fusion )
    
    
    # 
    assign(intervention, t1)
    
  } # intervention loop ends
  
  
  res_rad <- rbind(Rad_Fused$res_rad,Norad_Fused$res_rad)
  res_model <- rbind(Rad_Fused$res_model,Norad_Fused$res_model)
  res_OS <- rbind(Rad_Fused$res_OS,Norad_Fused$res_OS)
  res_tmat <- rbind(Rad_Fused$res_tmat,Norad_Fused$res_tmat)
  res_CI <- rbind(Rad_Fused$res_CI,Norad_Fused$res_CI)
  return(list(res_rad = res_rad,
              res_model = res_model,
              res_OS = res_OS,
              res_tmat = res_tmat,
              res_CI = res_CI
  ))
  
  # return(paste0("completed",sim_numi))
  
}
