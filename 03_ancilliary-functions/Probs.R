# Probs.R

Probs <- function( M_it, D_m,  df.i, Vis_m, CycleLength, GlobT, EstPLGG, EstAE,
                   AE_outside_mat1,PrMortCardio_mat1,SRS_mat, curent_age_v, rad_benefit1, inter, tt, trt_dur){
  # input: 
  # M_it: a numeric vector indicating current state 
  # D_m: matrix (numeric) indicating the duration of stay in each state
  # df.i: data.frame of cohort characteristics
  # Vis_m: a matrix indicating whether a patient has visited a specific state
  # CycleLength: numeric length of cycle
  # GlobT: global time
  # EstPLGG: output from EstProbs.plgg
  # EstAE: output from EstProvs.ae
  # AE_outside_mat1: output from AE_outside_pr
  # PrMortCardio_mat: a matrix indicating probaility of Cardiovascular related mortality
  # SRS_mat: SRS_mat
  # output  
  #  m.p.it: a matrix of dim number of individuals,number of states  
  #           where the i,jth element is the probability
  #           of transitioning from the current state to state j for individual i
  
  
  # First create empty vector that i'll be adding to
  m.p.it <- matrix(0, n.i, n.s)   
  colnames(m.p.it) <- v.n
  
  # Identify which cases (if any ) are fused
  FusedLog <- df.i$Fusion == 1
  anyFused <- any(FusedLog)
  # Do the same for Radiation
  RadIndc <- df.i$Radiation == 1
  
  
  # Identifying which PLGG related state each patient is in not counting death
  # Can do this sequentially since patients go from pre-prog to prog 1 than prog 2
  plgg.it <- rowSums(D_m[,c("pre_prog","prog1","prog2")] > 0 )
  plgg.char  <- c("pre_prog","prog1","prog2")[plgg.it]
  
  
  # Identify time in state for the most recent progression. 
  # This requiers some interesting subsetting. 
  # Basically when we collpase D_m as a vector it does so by rows which we can
  # correct index values using plgg.timein. This make sense if you do a simple test
  plgg.timeinstate <- as.vector(D_m[,c(1,2,3)])[(plgg.it-1)*n.i + 1:n.i]
  
  # Get errors if attempt to do this with empty vectors
  # so must have atleast some fused patients
  
  
  if(anyFused){
    # Use probs functions created earlier to estimate probabilities 
    # m.p.it[FusedLog, 1:4] <- Probs.plggFused(t_instate   = plgg.timeinstate[FusedLog] +CycleLength , 
    #                                          Cycle1      =  CycleLength,
    #                                          StateIndx   = plgg.char[FusedLog],
    #                                          FlexSurvRes = EstPLGG,
    #                                          inter       = inter,
    #                                          tt          = tt,
    #                                          trt_dur     = trt_dur)
    # 
    # # Similar for WT patients
    # m.p.it[!FusedLog,c(1,2,3,4)] <- Probs.plggWT(t_instate = plgg.timeinstate[!FusedLog] + CycleLength, 
    #                                              Cycle1 = CycleLength,
    #                                              StateIndx = plgg.char[!FusedLog],
    #                                              FlexSurvRes = EstPLGG)
    
  } else { # No fused
    m.p.it[,c(1,2,3,4)] <- Probs.plggWT(t_instate   = plgg.timeinstate+ CycleLength, 
                                        Cycle1      = CycleLength,
                                        StateIndx   = plgg.char,
                                        FlexSurvRes = EstPLGG,
                                        inter       = inter,
                                        tt          = tt,
                                        trt_dur     = trt_dur)
  }
  
  ## Radiation Benefit
  if(!(is.data.frame(rad_benefit1))){ # baseline
    # do nothing
    
  } else if(is.data.frame(rad_benefit1)){
    
    # if any have been radiated and currently in prog2
    if(any( RadIndc & plgg.char == "prog1")){
      
      m.p.it[ (RadIndc  & plgg.char == "prog1"),"prog2"] <- SurvProbFun(object = rad_benefit1$flexobj[[1]][[1]], 
                                                                               plgg.timeinstate[RadIndc & plgg.char == "prog1"] + CycleLength,
                                                                               
                                                                               cycle1 = CycleLength) * rad_benefit1$rr
      
    }
}
  
  
  
# Radiation related 
v.RadPr <- Probs.trtRad(t  = GlobT+ CycleLength,cycle1 = CycleLength,flexAE_Rad = EstAE)
v.NoRadPr <- Probs.trtNoRad(t  = GlobT+ CycleLength,cycle1 = CycleLength,flexAE_NoRad =  EstAE)
  
  # Assinging based on radiation related indicator created earlier have to 
  # turn into matrixces sometimes freaks out if trying to create an empty matrix
  # nrow = 0 but data > 0 so suppressingWarnings
  
  
m.p.it[RadIndc,c("neurologic","auditory","visual","stroke","cardiovascular","SN")] <- suppressWarnings(matrix(ncol = length(v.RadPr), 
                                                                                                                nrow = sum(RadIndc), 
                                                                                                                data = v.RadPr,byrow = T))
m.p.it[!RadIndc,c("neurologic","auditory","visual","stroke","cardiovascular","SN")] <- suppressWarnings(matrix(ncol = length(v.NoRadPr), 
                                                                                                                 nrow = sum(!RadIndc), 
                                                                                                                 data = v.NoRadPr,byrow = T))
  
  #indexing AE_outside_mat1 to get probability of storke, cardiovascular, and SN 
  m.p.it[,c("stroke_gen","cardiovascular_gen","SN_gen","death_outside")] <- AE_outside_mat1[fun_index(cur_time = curent_age_v,cycles_year = 1/cyc.t),c("Pr_stroke_vect","Pr_cardio_vect","Pr_cancer_vect","Pr_death_vect")]  
  
  
  # Mortality estimates
  # Non-plgg, non-cancer,stroke or cardiovascular risk of mortality.
  
  # Probability of death due to Cardiovascular event, 
  
  HasCardioVascularEvent <- (rowSums(Vis_m[, c("cardiovascular","cardiovascular_gen") ]) > 0)
  if(any(HasCardioVascularEvent)){
    time_sincecardiodx <-rowSums(D_m[HasCardioVascularEvent, c("cardiovascular","cardiovascular_gen"),drop = FALSE])
    
    m.p.it[HasCardioVascularEvent, c("death_cardio")]  <- PrMortCardio_mat1[fun_index(time_sincecardiodx, cycles_year = 1/cyc.t),2]
  }
  
  # Stroke within the last year
  Stroke_mort  <- (rowSums(Vis_m[,c("stroke","stroke_gen")]) >0  ) & rowSums(Vis_m[,c("stroke","stroke_gen"), drop = FALSE])  < 1 + cyc.t
  
  if(any(Stroke_mort)){
    age_strok1 <- curent_age_v[Stroke_mort]  - rowSums(Vis_m[Stroke_mort,c("stroke","stroke_gen"),drop = FALSE])
    m.p.it[Stroke_mort, "death_stroke"] <- PrDeath_Stroke(AgeStroke = age_strok1,Gender = df.i$Sex[Stroke_mort],cycles_year = 1/cyc.t)
  }
  
  
  # Death SMN
  
  HasSMN <- (rowSums(Vis_m[, c("SN","SN_gen"),drop = FALSE ]) > 0)
  if(any(HasSMN)){
    pr_mortnonSMN <- AE_outside_mat1[ fun_index(cur_time = curent_age_v[HasSMN],cycles_year = 1/cyc.t  ), "Pr_death_smn"]
    time_dx_ae <- rowSums(Vis_m[HasSMN, c("SN","SN_gen") ,drop = FALSE ])
    time_dx_ae[time_dx_ae>10] <- 10
    
    m.p.it[HasSMN, "death_SN"]  <-   PrDeath_SM( CurAge = curent_age_v[HasSMN], VectSRSdif = SRS_mat[fun_index(time_dx_ae,cycles_year = 1/cyc.t),"SRS_dif"],MortRiskAge = pr_mortnonSMN
                                                 
                                                 
    )
  }
  
  
  for(j in c("neurologic","auditory","visual")){
    
    m.p.it[Vis_m[,j],j]  <- 0
  }
  
  for(j in c("stroke", "cardiovascular","SN", "stroke_gen","cardiovascular_gen","SN_gen")){
    if(j %in%  c("stroke","stroke_gen")){
      m.p.it[(Vis_m[,"stroke"] |Vis_m[,"stroke_gen"]) ,c("stroke","stroke_gen")]  <- 0
      
    } else if (j %in% c("cardiovascular","cardiovascular_gen")){
      m.p.it[(Vis_m[,"cardiovascular"] |Vis_m[,"cardiovascular_gen"]) , c("cardiovascular","cardiovascular_gen")]  <- 0
      
    } else if(j %in% c("SN","SN_gen")){
      m.p.it[(Vis_m[,"SN"] |Vis_m[,"SN_gen"]) ,c("SN","SN_gen")]  <- 0
    }
  }
  
  
  
  # Calculate the probability of reamining in current state
  
  
  # Calculating the probability of remaining in state as 
  # 1 - sum(probability i leave to another state) 
  
  for(i in v.n ){
    if(any(M_it == i)){
      m.p.it[M_it == i,i] <-1-rowSums(m.p.it[M_it == i, c(v.n != i),drop = FALSE ])
    }
  }
  
  # Ensuringthat if death state you can't leave.
  
  
  m.p.it[M_it == "death_plgg",] <- 0
  m.p.it[M_it == "death_plgg","death_plgg"] <-1   
  
  m.p.it[M_it == "death_outside",] <- 0
  m.p.it[M_it == "death_outside","death_outside"] <-1   
  
  
  m.p.it[M_it == "death_stroke",] <- 0
  m.p.it[M_it == "death_stroke","death_stroke"] <-1   
  
  m.p.it[M_it == "death_SN",] <- 0
  m.p.it[M_it == "death_SN","death_SN"] <-1  
  
  m.p.it[M_it == "death_cardio",] <- 0
  m.p.it[M_it == "death_cardio","death_cardio"] <-1  
  
  # Returning Probs
  return(m.p.it)
  
}



Probs_scenario <- function( M_it, D_m,  df.i, Vis_m, CycleLength, GlobT, EstPLGG, EstAE,
                   AE_outside_mat1,PrMortCardio_mat1,SRS_mat, curent_age_v,rad_benefit1, prob_scenario1){
  # input: 
  # M_it: a numeric vector indicating current state 
  # D_m: matrix (numeric) indicating the duration of stay in each state
  # df.i: data.frame of cohort characteristics
  # Vis_m: a matrix indicating whether a patient has visited a specific state
  # CycleLength: numeric length of cycle
  # GlobT: global time
  # EstPLGG: output from EstProbs.plgg
  # EstAE: output from EstProvs.ae
  # AE_outside_mat1: output from AE_outside_pr
  # PrMortCardio_mat: a matrix indicating probaility of Cardiovascular related mortality
  # SRS_mat: SRS_mat
  # output  
  #  m.p.it: a matrix of dim number of individuals,number of states  
  #           where the i,jth element is the probability
  #           of transitioning from the current state to state j for individual i
  
  
  # First create empty vector that i'll be adding to
  m.p.it <- matrix(0, n.i, n.s)   
  colnames(m.p.it) <- v.n
  
  # Identify which cases (if any ) are fused
  FusedLog <- df.i$Fusion == 1
  anyFused <- any(FusedLog)
  # Do the same for Radiation
  RadIndc <- df.i$Radiation == 1
  
  
  # Identifying which PLGG related state each patient is in not counting death
  # Can do this sequentially since patients go from pre-prog to prog 1 than prog 2
  plgg.it <- rowSums(D_m[,c("pre_prog","prog1","prog2")] > 0 )
  plgg.char  <- c("pre_prog","prog1","prog2")[plgg.it]
  
  
  # Identify time in state for the most recent progression. 
  # This requiers some interesting subsetting. 
  # Basically when we collpase D_m as a vector it does so by rows which we can
  # correct index values using plgg.timein. This make sense if you do a simple test
  plgg.timeinstate <- as.vector(D_m[,c(1,2,3)])[(plgg.it-1)*n.i + 1:n.i]
  
  # Get errors if attempt to do this with empty vectors
  # so must have atleast some fused patients
  
  
  if(anyFused){
    # Use probs functions created earlier to estimate probabilities 
    m.p.it[FusedLog, 1:4] <- Probs.plggFused(t_instate = plgg.timeinstate[FusedLog] +CycleLength , 
                                             Cycle1 =  CycleLength,
                                             StateIndx = plgg.char[FusedLog],FlexSurvRes = EstPLGG)
    
    # Similar for WT patients
    m.p.it[!FusedLog,c(1,2,3,4)] <- Probs.plggWT(t_instate = plgg.timeinstate[!FusedLog] + CycleLength, 
                                                 Cycle1 = CycleLength,
                                                 StateIndx = plgg.char[!FusedLog],
                                                 FlexSurvRes = EstPLGG)
    
  } else { # No fused
    m.p.it[,c(1,2,3,4)] <- Probs.plggWT(t_instate = plgg.timeinstate+ CycleLength, 
                                        Cycle1 =CycleLength,
                                        StateIndx = plgg.char,
                                        FlexSurvRes = EstPLGG)
  }
  
  ## Radiation Benefit
  if(!(is.data.frame(rad_benefit1))){ # baseline
    # do nothing
    
  } else if(is.data.frame(rad_benefit1)){
    
    # if any have been radiated and currently in prog2
    if(any( RadIndc & plgg.char == "prog1")){
      
      m.p.it[ (RadIndc  & plgg.char == "prog1"),"prog2"] <- SurvProbFun(object = rad_benefit1$flexobj[[1]][[1]], 
                                                                        plgg.timeinstate[RadIndc & plgg.char == "prog1"] + CycleLength,
                                                                        
                                                                        cycle1 = CycleLength) * rad_benefit1$rr
      
    }
  }
  
  
  
  # Radiation related 
  v.RadPr <- Probs.trtRadScenario(t  = GlobT+ CycleLength,cycle1 = CycleLength,flexAE_Rad = EstAE, prob_delta_scenario = prob_scenario1 )
  v.NoRadPr <- Probs.trtNoRad(t  = GlobT+ CycleLength,cycle1 = CycleLength,flexAE_NoRad =  EstAE)
  
  # Assinging based on radiation related indicator created earlier have to 
  # turn into matrixces sometimes freaks out if trying to create an empty matrix
  # nrow = 0 but data > 0 so suppressingWarnings
  
  
  m.p.it[RadIndc,c("neurologic","auditory","visual","stroke","cardiovascular","SN")] <- suppressWarnings(matrix(ncol = length(v.RadPr), 
                                                                                                                nrow = sum(RadIndc), 
                                                                                                                data = v.RadPr,byrow = T))
  m.p.it[!RadIndc,c("neurologic","auditory","visual","stroke","cardiovascular","SN")] <- suppressWarnings(matrix(ncol = length(v.NoRadPr), 
                                                                                                                 nrow = sum(!RadIndc), 
                                                                                                                 data = v.NoRadPr,byrow = T))
  
  #indexing AE_outside_mat1 to get probability of storke, cardiovascular, and SN 
  m.p.it[,c("stroke_gen","cardiovascular_gen","SN_gen","death_outside")] <- AE_outside_mat1[fun_index(cur_time = curent_age_v,cycles_year = 1/cyc.t),c("Pr_stroke_vect","Pr_cardio_vect","Pr_cancer_vect","Pr_death_vect")]  
  
  
  # Mortality estimates
  # Non-plgg, non-cancer,stroke or cardiovascular risk of mortality.
  
  # Probability of death due to Cardiovascular event, 
  
  HasCardioVascularEvent <- (rowSums(Vis_m[, c("cardiovascular","cardiovascular_gen") ]) > 0)
  if(any(HasCardioVascularEvent)){
    time_sincecardiodx <-rowSums(D_m[HasCardioVascularEvent, c("cardiovascular","cardiovascular_gen"),drop = FALSE])
    
    m.p.it[HasCardioVascularEvent, c("death_cardio")]  <- PrMortCardio_mat1[fun_index(time_sincecardiodx, cycles_year = 1/cyc.t),2]
  }
  
  # Stroke within the last year
  Stroke_mort  <- (rowSums(Vis_m[,c("stroke","stroke_gen")]) >0  ) & rowSums(Vis_m[,c("stroke","stroke_gen"), drop = FALSE])  < 1 + cyc.t
  
  if(any(Stroke_mort)){
    age_strok1 <- curent_age_v[Stroke_mort]  - rowSums(Vis_m[Stroke_mort,c("stroke","stroke_gen"),drop = FALSE])
    m.p.it[Stroke_mort, "death_stroke"] <- PrDeath_Stroke(AgeStroke = age_strok1,Gender = df.i$Sex[Stroke_mort],cycles_year = 1/cyc.t)
  }
  
  
  # Death SMN
  
  HasSMN <- (rowSums(Vis_m[, c("SN","SN_gen"),drop = FALSE ]) > 0)
  if(any(HasSMN)){
    pr_mortnonSMN <- AE_outside_mat1[ fun_index(cur_time = curent_age_v[HasSMN],cycles_year = 1/cyc.t  ), "Pr_death_smn"]
    time_dx_ae <- rowSums(Vis_m[HasSMN, c("SN","SN_gen") ,drop = FALSE ])
    time_dx_ae[time_dx_ae>10] <- 10
    
    m.p.it[HasSMN, "death_SN"]  <-   PrDeath_SM( CurAge = curent_age_v[HasSMN], VectSRSdif = SRS_mat[fun_index(time_dx_ae,cycles_year = 1/cyc.t),"SRS_dif"],MortRiskAge = pr_mortnonSMN
                                                 
                                                 
    )
  }
  
  
  for(j in c("neurologic","auditory","visual")){
    
    m.p.it[Vis_m[,j],j]  <- 0
  }
  
  for(j in c("stroke", "cardiovascular","SN", "stroke_gen","cardiovascular_gen","SN_gen")){
    if(j %in%  c("stroke","stroke_gen")){
      m.p.it[(Vis_m[,"stroke"] |Vis_m[,"stroke_gen"]) ,c("stroke","stroke_gen")]  <- 0
      
    } else if (j %in% c("cardiovascular","cardiovascular_gen")){
      m.p.it[(Vis_m[,"cardiovascular"] |Vis_m[,"cardiovascular_gen"]) , c("cardiovascular","cardiovascular_gen")]  <- 0
      
    } else if(j %in% c("SN","SN_gen")){
      m.p.it[(Vis_m[,"SN"] |Vis_m[,"SN_gen"]) ,c("SN","SN_gen")]  <- 0
    }
  }
  
  
  
  # Calculate the probability of reamining in current state
  
  
  # Calculating the probability of remaining in state as 
  # 1 - sum(probability i leave to another state) this is wonky 
  # for the death state but I fix that in the next line.
  for(i in v.n ){
    if(any(M_it == i)){
      m.p.it[M_it == i,i] <-1-rowSums(m.p.it[M_it == i, c(v.n != i),drop = FALSE ])
    }
  }
  
  # Ensuringthat if death state you can't leave.
  
  
  m.p.it[M_it == "death_plgg",] <- 0
  m.p.it[M_it == "death_plgg","death_plgg"] <-1   
  
  m.p.it[M_it == "death_outside",] <- 0
  m.p.it[M_it == "death_outside","death_outside"] <-1   
  
  
  m.p.it[M_it == "death_stroke",] <- 0
  m.p.it[M_it == "death_stroke","death_stroke"] <-1   
  
  m.p.it[M_it == "death_SN",] <- 0
  m.p.it[M_it == "death_SN","death_SN"] <-1  
  
  m.p.it[M_it == "death_cardio",] <- 0
  m.p.it[M_it == "death_cardio","death_cardio"] <-1  
  
  # Returning Probs
  return(m.p.it)
  
}
