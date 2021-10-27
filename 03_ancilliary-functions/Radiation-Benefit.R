# Radiation Benefit

# This R script creates a function to generate deterministic and bootstrapped radiation fitting
# fitsdistributions

FunFlexFit1 <- function(tr1, dist1,  Survdat1){
  # FunFlex Fits estimates flexsurv survival curves
  # Inputs
  #   tr1: transitions
  #   dist1: distribution to be tested
  #   fuse1: whether fused or not fused
  #   SurvDat1: survival data to be fitted
  # Returns
  #  either a flexsurv object or an error message

  ret <- suppressWarnings(try( flexsurvreg(Surv(time,   status ) ~ 1 , data = Survdat1 , dist = dist1), silent = T))
  return(ret)
}



## Read data
### Input to move over
# dig_radbenefit <- readRDS(here::here("02_data","cherlow_dig.RDS"))
# allows for a decrease in radiation benefit


fun_genRadBenefit <- function(dig_radbenefit, deter){
  # digitized_radbenefit: digitized data from cherlow
  # deter : a logical indicating if deterministic
  # requires: flexsurv
require(flexsurv, quietly = T)
require(dplyr)
pos.dist <- names(flexsurv::flexsurv.dists)
  # Remove repeated ones
pos.dist <- pos.dist[!(pos.dist %in% c("genf.orig","gengamma","gengamma.orig",
                                         "weibullPH","exponential","lognormal","genf"))]


# Storing data  
res1  <- data.frame(post.dist = pos.dist, tr  = 2)
res1$AIC <- Inf
res1$flexobj <-""
res1$post.dist <- as.character(res1$post.dist)
# if determiistic no bootstrapping
if(deter){
  
  for(j in seq_along(res1$post.dist)){
    
    fit1 <- FunFlexFit1(tr1 = res1$tr[j], dist1 = res1$post.dist[j], Survdat1 = dig_radbenefit)
    
    if(typeof(fit1) != "character"){
      res1$flexobj[[j]] <- list( fit1)
      res1$AIC[[j]] <- fit1$AIC  
    } else{
      res1$flexobj[[j]] <- NA
      res1$AIC[[j]]  <- NA
    }
    
} # for each distribution loop

  
  
  
  } # deterministic loop   

if(!deter){
  
rows_tosample <-  sample(x = seq(1,dim(dig_radbenefit)[1]), size = length(seq(1,dim(dig_radbenefit)[1])), replace = T)
dig_radbenefit <- dig_radbenefit[rows_tosample,]


  for(j in seq_along(res1$post.dist)){
    
    fit1 <- FunFlexFit1(tr1 = res1$tr[j], dist1 = res1$post.dist[j], Survdat1 = dig_radbenefit)
    
    if(typeof(fit1) != "character"){
      res1$flexobj[[j]] <- list( fit1)
      res1$AIC[[j]] <- fit1$AIC  
    } else{
      res1$flexobj[[j]] <- NA
      res1$AIC[[j]]  <- NA
    }
  
  } 
  
 } # bootstrapped loop

BestFitsAIC<- res1 %>% group_by(tr) %>%  
  mutate(n = min(as.numeric(AIC), na.rm = TRUE)) %>% 
  filter(AIC == n) %>% select(-n,-AIC)
return(BestFitsAIC)
# Function end
}
