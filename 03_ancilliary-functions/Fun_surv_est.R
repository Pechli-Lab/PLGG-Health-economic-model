# Fun_surv_est.Rmd 
## Est.probs.plgg


# FunFlexFit is used by PLGG probs to estimate survival probabilities.

FunFlexFit <- function(tr1, dist1, fuse1, Survdat1){
  # FunFlex Fits estimates flexsurv survival curves
  # Inputs
  #   tr1: transitions
  #   dist1: distribution to be tested
  #   fuse1: whether fused or not fused
  #   SurvDat1: survival data to be fitted
  # Returns
  #  either a flexsurv object or an error message

  ret <- suppressWarnings(try( flexsurvreg(Surv(Time,   status ) ~ 1 , data = subset(Survdat1, tr == tr1 & Fusion == fuse1), dist = dist1), silent = T))
  return(ret)
}

# EstProbs.plgg #####
# Estimates transition probabilites
EstProbs.plgg <- function(SurvDat, returnAll = FALSE){
  ## Input:
  # SurvDat: Survival data of the structure (Time, Status, Fusion)
  # returnALL: logical indicating if to return all fits (TRUE) or only best (FALSE)
  # Output: 
  # BestFitsAIC a data frame with the best fitting curve (according to AIC)
  # for the 5 transitions and fusion and non-fused status
  # transitions: from: c(1,1,2,2,3). to c(2,4,3,4,4) 
  
  # Need to have these packages loaded
  require(dfoptim , quietly = T)
  require(flexsurv, quietly = T)
  require(dplyr,    quietly = T)
  require(tidyr,    quietly = T)
  
  # Get the names of the survival distributions
  pos.dist <- names(flexsurv::flexsurv.dists)
  # Remove repeated ones
  pos.dist <- pos.dist[!(pos.dist %in% c("genf.orig","gengamma","gengamma.orig",
                                         "weibullPH","exponential","lognormal","genf"))]
  
  # create a data frame to store results 
  res.df <- SurvDat %>% 
    dplyr::count(tr,Fusion) %>% select(-n) %>% 
    mutate(dist = list(pos.dist)) %>% unnest(dist) 
  res.df$flexobj <- ""
  res.df$AIC <- ""
  
  
  # For each transition,fusion status, and distribution attempt to fit a survcurv
  # If returns an error make it NA
  # If not save the flexsurv object as well as AIC
  for(i in 1:dim(res.df)[1]){
    fit1   <- FunFlexFit(tr1 = res.df$tr[i], fuse1 = res.df$Fusion[i],
                         dist1 = res.df$dist[i], Survdat1 = SurvDat)
    if(typeof(fit1) != "character"){
      res.df$flexobj[[i]] <- list( fit1)
      res.df$AIC[[i]] <- fit1$AIC  
    } else{
      res.df$flexobj[[i]] <- NA
      res.df$AIC[[i]]  <- NA
    }
  }
  
  # Select the best fits for each transfusion and fusion status by selecting
  # the distribution with the lowest AIC.
  if(returnAll == FALSE){
    
    BestFitsAIC<- res.df %>% group_by(tr, Fusion) %>%  
      mutate(n = min(as.numeric(AIC), na.rm = TRUE)) %>% 
      filter(AIC == n) %>% select(-n,-AIC)
    # return that data frame which includes, fusion,tr, distribution, flexsurvobj
    return(BestFitsAIC) }
  if(returnAll == TRUE){
    return(res.df)
  }
}

# DeidentifiedSurv <- readRDS("~/OneDrive - SickKids/PLGG_Economic_Evaluation/PLGG-R_2019-07-10/DeidentifiedSurv.RDS")
# R1 <-EstProbs.plgg(filter(DeidentifiedSurv, Time > 0),returnAll = TRUE)


## Estimating AE events

FunFlexFit_AE <- function(dist1,split_i, df_IPD1){
  # FunFlex Fits estimates flexsurv survival curves for AE
  # Inputs
  #  split_i: split estimating adverse events
  #  df_IPD  (here::here("02_data","IPD_ae.RDS")) long IPD from digitized  
  # Returns
  #  either a flexsurv object or an error message

  if(any(df_IPD1$event ==3)){
    ret <- suppressWarnings(try( flexsurvreg(Surv(time = Time0, time2 =  Time1,event =  event, type = "interval") ~ 1 , data = subset(df_IPD1, split1 == split_i), dist = dist1), silent = T))
    return(ret)
  } else{
    ret <- suppressWarnings(try( flexsurvreg(Surv(time =  time2,event =  event) ~ 1 , data = subset(df_IPD1, split1 == split_i), dist = dist1), silent = T))
  }
}


EstProbs.AE <- function(df_IPD, returnAll = FALSE, deter = T ){
  ## Input:
  # df_IPD: digitized ipatiend data
  # returnALL: logical indicating if to return all fits (TRUE) or only best (FALSE)
  # deter: a logical inicating if deterministic or bootstrapped
  # Output: 
  # BestFitsAIC a data frame with the best fitting curve (according to AIC)
  # for the 12 AE and radiation no radiation combination

## Bootstrapping procedure

  
  
if(deter == F){
  
  # create a list with the names of each item in df_IPD$split1
  # 
bs_list <- vector(mode = "list",length = length(unique(df_IPD$split1)))  
names(bs_list) <- unique(df_IPD$split1)

for( i in seq_along(bs_list)){
  tempdf <- df_IPD[df_IPD$split1 == names(bs_list)[i],]
  rows_tosample <- seq(from = 1, to = dim(tempdf)[1])
  tempdf <- tempdf[sample(x = rows_tosample , size = length(rows_tosample),replace = T), ]
  bs_list[[i]] <- tempdf
  
  
  
}

df_IPD <- do.call(what = dplyr::bind_rows, bs_list)
}  

  # Need to have these packages loaded
  require(dfoptim , quietly = T)
  require(flexsurv, quietly = T)
  require(dplyr,    quietly = T)
  require(tidyr,    quietly = T)
  
  # Get the names of the survival distributions
  pos.dist <- names(flexsurv::flexsurv.dists)
  # Remove repeated ones
  pos.dist <- pos.dist[!(pos.dist %in% c("genf.orig","gengamma","gengamma.orig",
                                         "weibullPH","exponential","lognormal","genf"))]
  
  # create a data frame to store results 
  
  res.df <- tibble(split1 = unique(df_IPD$split1),
                   dist = list(pos.dist),
                   flexobj = "",
                   AIC = "") %>% unnest(dist) 
  
  
  # For each split1, and distribution attempt to fit a survcurv
  # If returns an error make it NA
  # If not save the flexsurv object as well as AIC
  for(i in 1:dim(res.df)[1]){
    fit1   <- FunFlexFit_AE(dist1 = res.df$dist[i],split_i = res.df$split1[i],df_IPD1 = df_IPD)
    if(typeof(fit1) != "character"){
      res.df$flexobj[[i]] <- list( fit1)
      res.df$AIC[[i]] <- fit1$AIC  
    } else{
      res.df$flexobj[[i]] <- NA
      res.df$AIC[[i]]  <- NA
    }
  }
  
  # Select the best fits for each transfusion and fusion status by selecting
  # the distribution with the lowest AIC.
  if(returnAll == FALSE){
    
    BestFitsAIC<- res.df %>% group_by(split1) %>%  
      mutate(n = min(as.numeric(AIC), na.rm = TRUE)) %>% 
      filter(AIC == n) %>% select(-n,-AIC)
    # return that data frame which includes, fusion,tr, distribution, flexsurvobj
    return(BestFitsAIC) }
  if(returnAll == TRUE){
    return(res.df)
  }
}

# df_IPC <- readRDS( here::here("02_data","IPD_ae.RDS"))
# R2 <-EstProbs.AE(df_IPD,returnAll = T)

