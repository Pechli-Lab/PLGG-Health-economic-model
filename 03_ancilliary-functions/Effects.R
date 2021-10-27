# Effects R consists of three functions 
# gen_effects that generates model input
# betaPar, getting beta parameters from the mean and SE of utilities sourche
# Effects that assigns costs

# get the beta parameters from the mean and se
require(dplyr)
betaPar=function(m,s){
  a <-  m * ((m * (1 - m)/s ^ 2) - 1)
  b <-  (1 - m) * ((m * (1 - m)/s ^ 2) - 1)
  list(a = a, b = b)
}


genBeta <- function(n, betaPar)
{
  res <- rbeta(n = n, betaPar$a, betaPar$b)
  return(res)
}





gen_Effects <- function(deter = TRUE, 
                        beta_Utils  ){
  # Inputs:
  #   deter:  Logical if deterministic or not 
  #   list_beta_input: beta parameter inputs
  # output: 
  #   utilities list to be used by gen_effects

# mean of a beta distribution is alpha/(alpha + beta)
# creating deterministic mean. 

if(deter){ # deterministic
  
  util_list <- list(
    pre_prog        = c(beta_Utils["Children_PLGG", "mean"], beta_Utils["Adults_PLGG", "mean"]),
    prog1           = c(beta_Utils["Children_PLGG", "mean"], beta_Utils["Adults_PLGG", "mean"]),
    prog2           = c(beta_Utils["Children_PLGG", "mean"], beta_Utils["Adults_PLGG", "mean"]),
    Death           = 0,
    auditory        = c(beta_Utils["Children_auditory", "mean"], beta_Utils["Adults_auditory", "mean"]),
    cardiovascular  =   beta_Utils["heart",       "mean"],
    neurologic      =   beta_Utils["neurologic",  "mean"],
    SN              =   beta_Utils["cancer_ever", "mean"],
    stroke          =   beta_Utils["stroke",      "mean"],
    visual          =   beta_Utils["visual",      "mean"],
    
    Nonpre_prog        = c(beta_Utils["NonChildren_PLGG", "mean"], beta_Utils["NonAdults_PLGG", "mean"]),
    Nonprog1           = c(beta_Utils["NonChildren_PLGG", "mean"], beta_Utils["NonAdults_PLGG", "mean"]),
    Nonprog2           = c(beta_Utils["NonChildren_PLGG", "mean"], beta_Utils["NonAdults_PLGG", "mean"]),
    Nonauditory        = c(beta_Utils["NonChildren_auditory", "mean"], beta_Utils["NonAdults_auditory", "mean"]),
    Noncardiovascular  =   beta_Utils["Nonheart",       "mean"],
    Nonneurologic      =   beta_Utils["Nonneurologic",  "mean"],
    NonSN              =   beta_Utils["Noncancer_ever", "mean"],
    Nonstroke          =   beta_Utils["Nonstroke",      "mean"],
    Nonvisual          =   beta_Utils["Nonvisual",      "mean"],
    Non_all            =   beta_Utils["Noncancer_noAE", "mean"])
  
  
  
  util_list_ret <- list( 
    pre_prog        = util_list$pre_prog       /util_list$Nonpre_prog,
    prog1           = util_list$prog1          /util_list$Nonprog1,
    prog2           = util_list$prog2          /util_list$Nonprog2,
    Death           = util_list$Death,
    auditory        = util_list$auditory       /util_list$Nonauditory,
    cardiovascular  = util_list$cardiovascular /util_list$Noncardiovascular,
    neurologic      = util_list$neurologic     /util_list$Nonneurologic,
    SN              = util_list$SN             /util_list$NonSN,
    stroke          = util_list$stroke         /util_list$Nonstroke,
    visual          = util_list$visual         /util_list$Nonvisual,
    Non_all         = util_list$Non_all
   
  )
  
   
 } else if(!deter){

# if not determenistic creating random draw from beta distribution  
plgg_rel  <-  c(genBeta(n = 1, beta_Utils["Children_PLGG", c("alpha","beta")]),
                genBeta(n = 1, beta_Utils["Adults_PLGG",   c("alpha","beta")]))
Nonplgg_rel  <-  c(genBeta(n = 1, beta_Utils["NonChildren_PLGG", c("alpha","beta")]),
                genBeta(n = 1, beta_Utils["NonAdults_PLGG",   c("alpha","beta")]))

util_list <-   list(
    pre_prog        = plgg_rel,
    prog1           = plgg_rel,
    prog2           = plgg_rel,
    Death           = 0,
    auditory        = c(genBeta(n = 1, beta_Utils["Children_auditory", c("alpha","beta")]),
                        genBeta(n = 1, beta_Utils["Adults_auditory"   , c("alpha","beta")])),
    cardiovascular  =   genBeta(n = 1, beta_Utils["heart",             c("alpha","beta")]),
    neurologic      =   genBeta(n = 1, beta_Utils["neurologic"       , c("alpha","beta")]),
    SN              =   genBeta(n = 1, beta_Utils["cancer_ever"      , c("alpha","beta")]),
    stroke          =   genBeta(n = 1, beta_Utils["stroke"           , c("alpha","beta")]),
    visual          =   genBeta(n = 1, beta_Utils["visual"           , c("alpha","beta")]),
    
    Nonpre_prog        = Nonplgg_rel,
    Nonprog1           = Nonplgg_rel,
    Nonprog2           = Nonplgg_rel,
    Nonauditory        = c(genBeta(n = 1, beta_Utils["NonChildren_auditory", c("alpha","beta")]),
                           genBeta(n = 1, beta_Utils["NonAdults_auditory"  , c("alpha","beta")])),
    Noncardiovascular  =   genBeta(n = 1, beta_Utils["Nonheart",             c("alpha","beta")]),
    Nonneurologic      =   genBeta(n = 1, beta_Utils["Nonneurologic"       , c("alpha","beta")]),
    NonSN              =   genBeta(n = 1, beta_Utils["Noncancer_ever"      , c("alpha","beta")]),
    Nonstroke          =   genBeta(n = 1, beta_Utils["Nonstroke"           , c("alpha","beta")]),
    Nonvisual          =   genBeta(n = 1, beta_Utils["Nonvisual"           , c("alpha","beta")]),
    NonAll             =   genBeta(n = 1, beta_Utils["Noncancer_noAE"      , c("alpha","beta")]))


util_list_ret <- list( 
  pre_prog        = util_list$pre_prog       /util_list$Nonpre_prog,
  prog1           = util_list$prog1          /util_list$Nonprog1,
  prog2           = util_list$prog2          /util_list$Nonprog2,
  Death           = util_list$Death,
  auditory        = util_list$auditory       /util_list$Nonauditory,
  cardiovascular  = util_list$cardiovascular /util_list$Noncardiovascular,
  neurologic      = util_list$neurologic     /util_list$Nonneurologic,
  SN              = util_list$SN             /util_list$NonSN,
  stroke          = util_list$stroke         /util_list$Nonstroke,
  visual          = util_list$visual         /util_list$Nonvisual,
  Non_all         = util_list$NonAll
  
)


}

  return(util_list_ret)
}
  






Effects <- function(lst.Util, cur_age_v, m_V){
  # inputs:
  #   lst.Util a list of length n.v +1 with varying length of items
  #   m_V: matrix indicating whether an individual has visited a state
  #   cur_age_v: a vector indicating current age of PLGG patient
  #   
  #   Output
  #   Util.Res: A vector of the current utility in each state.
  
  # Going to store utilities in m.Util
  m.Util <- matrix(nrow = n.i, ncol = n.s + 1, data = 1 ,dimnames = list(NULL, c(v.n,"gen_pop")))
  
  # First identifying which state they are in
  PLGG.relstate <- c("pre_prog", "prog1","prog2")[rowSums(m_V[,c("pre_prog", "prog1","prog2")])]
  DeathIdc <- rowSums(m_V[,c("death_plgg","death_outside","death_stroke","death_SN","death_cardio")] > 0) == 1
  
  
  ## PLGG Related Utilities
  ## Currently don't have varying utilities varying by age
  m.Util[m_V[,"pre_prog"] & cur_age_v <= 18 ,"pre_prog"] <- lst.Util$pre_prog[1]
  m.Util[m_V[,"pre_prog"] & cur_age_v > 18 ,"pre_prog"] <- lst.Util$pre_prog[2]
  
  # Utility when dead. Again don't have mulitple disutilieis so assinging 0 to m.Util
  m.Util[DeathIdc,"death_plgg"] <- 0
  
  ## Treatment related Adverse events
  # Neurologic
  m.Util[(m_V[,c("neurologic")] & cur_age_v > 18) ,"neurologic"] <- lst.Util$neurologic
  
  # Auditory
  m.Util[(m_V[,"auditory"] & cur_age_v <= 18) , "auditory" ] <- lst.Util$auditory[1]
  m.Util[(m_V[,"auditory"] & cur_age_v > 18) , "auditory" ]  <- lst.Util$auditory[2]
  
  # Cardiovascular only adults
  m.Util[((m_V[,"cardiovascular"] |m_V[,"cardiovascular_gen"])& cur_age_v > 18) , "cardiovascular"] <- lst.Util$cardiovascular 

  # SN
  m.Util[(m_V[,"SN"] | m_V[,"SN_gen"]),"SN"] <- lst.Util$SN[1]
  
  # Stroke
  m.Util[m_V[,"stroke"] | m_V[,"stroke_gen"],"stroke"] <- lst.Util$stroke[1]
    
  # Visual
  m.Util[m_V[,c("visual")],"visual"] <- lst.Util$visual[1]
  m.Util[, "gen_pop"] <- lst.Util$Non_all[1]
  
  
  
  Util.ret <- matrixStats::rowProds(m.Util)
  return(Util.ret)
  
}




