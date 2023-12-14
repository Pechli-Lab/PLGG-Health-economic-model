# Made a small modification where the simulation runs using numeric state indicators not colnames.

# Replaced lev <- dimnames(m.p) to lev <- 1:length(dimnmaes) 
# Function to run the vectorized version of microsimulation 

samplev2 <- function (probs, m) {
  d <- dim(probs)
  n <- d[1]
  k <- d[2]
  lev <- dimnames(probs)[[2]]
  if (!length(lev)) 
    lev <- 1:k
  ran <- matrix(lev[1], ncol = m, nrow = n)
  U <- t(probs)
  for(i in 2:k) {
    U[i, ] <- U[i, ] + U[i - 1, ]
  }
  if (any((U[k, ] - 1) > 1e-05))
    stop("error in multinom: probabilities do not sum to 1")
  
  for (j in 1:m) {
    un <- rep(runif(n), rep(k, n))
    ran[, j] <- lev[1 + colSums(un > U)]
  }
  ran
}

# Fun.GenI generates individuals. Need the file location of the analysis set and 
# the number to generate. NumToGen == n.i


## tempdat$Sex = 0 
## tempdat$Sex[tempdat$Gender == male]  =1 
## although I am not sure why you set NAs equal to zero

Fun.GenI<- function(loc.analysis.set,NumToGen, Seed_n1){
  # input: 
  set.seed(Seed_n1)
  
  tempdat   <- readRDS(loc.analysis.set)$df.diagnosis
  tempdat <- tempdat[,c("Fusion","AgeDx", "Gender")]
  tempdat$Sex <- 0
  tempdat$Sex[tempdat$Gender == "male"] <- 1
  tempdat <- tempdat[,c("Fusion","AgeDx","Sex")]
  tempdat <- tempdat[sample(x = 1:dim(tempdat)[1], size = NumToGen, replace = T),]
  return(tempdat)
}

# SurvProbFun ####
# This is the guts of the flexsurv summary function made to create probabilities
# that are vectorized. Not commented if interested run once with debug
SurvProbFun <- function(object, t , cycle1){
  start  = t - cycle1 
  x <- object
  dat <- x$data
  fn <-  function(t,start,...) {
    s_t <-(1 - x$dfns$p(start,...))
    s_t1 <- (1 - x$dfns$p(t,...))
    ret <- 1- s_t1/s_t
    
    ret[t<start] <- 1 # prob[t<start] was previously 0
    ret[s_t- s_t1 < .Machine$double.eps  ] <- 0
    ret
  }
  fncall <- list(t, start)
  beta <- if (x$ncovs == 0) {0}
  X <- as.matrix(0, nrow = 1, ncol = max(x$ncoveffs, 1))
  
  dlist <- x$dlist
  ret <- vector(nrow(X), mode = "list")
  
  basepars.mat <- flexsurv:::add.covs(x, x$res.t[dlist$pars, "est"], 
                                      beta, X[1, , drop = FALSE], transform = FALSE)
  basepars <- as.list(as.data.frame(basepars.mat))
  fncall[dlist$pars] <- basepars
  y <- do.call(fn, fncall)  
  return(y)
}

# Probs.plggFused and Probs.plggWT ########
# Similarly to treatment related adverse events creating a 
# helper function this time based on BRAF value
# They both follow a similar structure

## Flexsurv res is structured by state and fusion status
# row 1 tr:1 & fusion:0
# row 2 tr:1 & fusion:1
# row 3 tr:2 & fusion:0
# row 4 tr:2 & fusion:1
# row 5 tr:3 & fusion:0
# row 6 tr:3 & fusion:1
# row 7 tr:4 & fusion:0
# row 8 tr:4 & fusion:1
# row 9 tr:5 & fusion:0
# row 10 tr:5 & fusion:1

# rows 1,3,5,7,9 are associated with Fusion == 0

Probs.plggWT <- function(t_instate,StateIndx, Cycle1, FlexSurvRes, inter, tt, trt_dur){
  # Inputs
  # t_instate: time in current state, a vector
  # StateIndx: what state are they in
  # Cycle1: cycle length
  # Outputs probability natrix
  # Probability matrix which i'll return
  
  if (inter == "Targeted") {
    arm <- 2  
  } else if (inter == "SoC") {
    arm <- 11 
  }
  if (trt_dur == 2 & tt > 24) {arm <- 11} # treatment effect lasts for 2 years, goes back to SoC 
  
  Prob.mat <-  matrix(nrow = length(t_instate), ncol = 4,data = 0)
  
  # If any state is equal to the first state uses the SurvProbFun 
  # to estimate transition probabilities
  # only for those whose state is equal to 1. And so on.

  if(any(StateIndx ==  "pre_prog")){
    state1.log <- StateIndx ==  "pre_prog"
    Prob.mat[state1.log,2] <- SurvProbFun(object = FlexSurvRes$flexobj[[arm]][[1]], t = t_instate[state1.log], cycle1 = Cycle1)
    Prob.mat[state1.log,4] <- SurvProbFun(object = FlexSurvRes$flexobj[[3]][[1]],   t = t_instate[state1.log], cycle1 = Cycle1)
  }
  
  if(any(StateIndx ==   "prog1")){
    state2.log <- StateIndx ==   "prog1"
    Prob.mat[state2.log,3] <- SurvProbFun(object = FlexSurvRes$flexobj[[5]][[1]], t = t_instate[state2.log], cycle1 = Cycle1)
    Prob.mat[state2.log,4] <- SurvProbFun(object = FlexSurvRes$flexobj[[7]][[1]], t = t_instate[state2.log], cycle1 = Cycle1)
  }
  
  if(any(StateIndx ==   "prog2")){
    state3.log <- StateIndx ==   "prog2"
    Prob.mat[state3.log,4] <- SurvProbFun(object = FlexSurvRes$flexobj[[9]][[1]], t = t_instate[state3.log], cycle1 = Cycle1)
  }
  
  return(Prob.mat)
}
# Estimated_plgg1$flexobj[[5]][[1]] -> prog1to2
# prog1to2$res[,1] -> res2
# plot(0:10,plnorm(0:10, meanlog=res2[1], sdlog=res2[2], lower.tail = F), ylim=c(0,1), type="l")

# This is exact same structure as above only differencce is indexing in FlexSurvRes
# Probs.plggFused <- function(t_instate,StateIndx, Cycle1,FlexSurvRes ){
#   Prob.mat <-  matrix(nrow = length(t_instate), ncol = 4,data = 0)
#   if(any(StateIndx == "pre_prog")){
#     state1.log <- StateIndx == "pre_prog"
#     Prob.mat[state1.log,2] <- SurvProbFun(object = FlexSurvRes$flexobj[[2]][[1]], t = t_instate[state1.log], cycle1 = Cycle1)
#     Prob.mat[state1.log,4] <- SurvProbFun(object = FlexSurvRes$flexobj[[4]][[1]], t = t_instate[state1.log], cycle1 = Cycle1)
#   }
#   if(any(StateIndx == "prog1")){
#     state2.log <- StateIndx == "prog1"
#     Prob.mat[state2.log,3] <- SurvProbFun(object = FlexSurvRes$flexobj[[6]][[1]], t = t_instate[state2.log], cycle1 = Cycle1)
#     Prob.mat[state2.log,4] <- SurvProbFun(object = FlexSurvRes$flexobj[[8]][[1]], t = t_instate[state2.log], cycle1 = Cycle1)
#   }
#   if(any(StateIndx == "prog2")){
#     state3.log <- StateIndx == "prog2"
#     Prob.mat[state3.log,4] <- SurvProbFun(object = FlexSurvRes$flexobj[[10]][[1]], t = t_instate[state3.log], cycle1 = Cycle1)
#   }
#   
#   return(Prob.mat)
# }

## fun_index is a function that takes an input as time from an event and returns the index 
# of what time that would be assocaited with if the data frame started from 0 and went up by 
# cycles in a year 

# Creating `PrDeath_SM` which is the probability of death for those who have a secondary malignancy.


fun_index <- function(cur_time, cycles_year){
  # input: curr_time (either a digit or a vector) indicating the current time
  #        cylces_year: how many cycles in a year we want to be indexing from
  index_1 <- floor((cur_time*cycles_year)+0.0001) + 1
  return(index_1)
}


Probs.trtRad <- function(t, cycle1,flexAE_Rad){
  #Inputs
  # t: current time in model (since model starts at diagnosis)
  # cycle1: length of cycle
  # flexAE_Rad a list of flexsurv object (with the ith object being associated with the AE)
  # Outputs:
  # res.vect: A vector of the probabilities with the ith probably being the probaiblity of AE occuring at btween t-cycle1
  
  flexAE_Rad <-  flexAE_Rad[1:dim(flexAE_Rad)[1] %% 2 == 1,]
  
  res.vect <- vector(mode = "numeric", length = dim(flexAE_Rad)[1])
  
  for( i in 1:length(flexAE_Rad$flexobj)){
    res.vect[i] <- SurvProbFun(flexAE_Rad$flexobj[[i]][[1]], t, cycle1)
  }
  
  
  return(res.vect) 
}




Probs.trtNoRad <- function(t, cycle1,flexAE_NoRad){
  #Inputs
  # t: current time in model (since model starts at diagnosis)
  # cycle1: length of cycle
  # flexAE_NoRad a list of flexsurv object (with the ith object being associated with the AE)
  # Outputs:
  # res.vect: A vector of the probabilities with the ith probably being the probaiblity of AE occuring at btween t-cycle1
  flexAE_NoRad <-  flexAE_NoRad[1:dim(flexAE_NoRad)[1] %% 2 == 0,]
  res.vect <- vector(mode = "numeric", length = dim(flexAE_NoRad)[1])
  for( i in 1:length(flexAE_NoRad$flexobj)){
    res.vect[i] <- SurvProbFun(flexAE_NoRad$flexobj[[i]][[1]], t, cycle1)
  }
  
  
  return(res.vect) 
  
  
}

# test that it works

# Probs.trtNoRad(1,0.1,BestAE_NoRad)
# Probs.trtRad(1,0.1,BestAE_Rad)

## New AE function 

Probs.trtRadScenario <- function(t, cycle1,flexAE_Rad, prob_delta_scenario){
  #Inputs
  # t: current time in model (since model starts at diagnosis)
  # cycle1: length of cycle
  # flexAE_Rad a list of flexsurv object (with the ith object being associated with the AE)
  # Outputs:
  # res.vect: A vector of the probabilities with the ith probably being the probaiblity of AE occuring at btween t-cycle1
  flexAE_NoRad1 <-  flexAE_Rad[1:dim(flexAE_Rad)[1] %% 2 == 0,]
  flexAE_Rad1 <-    flexAE_Rad[1:dim(flexAE_Rad)[1] %% 2 == 1,]
  
  res.vect <- vector(mode = "numeric", length = dim(flexAE_Rad1)[1])
  
  for( i in 1:length(flexAE_Rad1$flexobj)){
    res.vect[i] <- (SurvProbFun(flexAE_Rad1$flexobj[[i]][[1]], t, cycle1)-SurvProbFun(flexAE_NoRad1$flexobj[[i]][[1]], t, cycle1))*prob_delta_scenario  + SurvProbFun(flexAE_NoRad1$flexobj[[i]][[1]], t, cycle1)
  }
  
  
  return(res.vect) 
}
