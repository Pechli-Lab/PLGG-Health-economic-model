gen_synpop <- function( n1=  100000, deter = T){
  #l1: list of mean ,sd agedx, and neseted glms fit to fusion and sex
    if(deter = T){
  AgeDx = round(runif(n = n1, min = 0,max = 18))
  coef1 <- c(0.03758116, -0.006658198) # intercenpt + AgeDx, binomial 
  o1 <- exp(cbind(1, AgeDx) %*% coef1)
  Sex = rbinom(n1,1, prob = o1/(1+o1))
  coef2 <- c(0.1107713, -0.05921511, -0.2586516) # interept + Agedx + Sex
  
  o2 <- exp(cbind(1, AgeDx,Sex) %*% coef2)
  Fusion =  rbinom(n1,1, prob = o2/(1+o2))
  
  
  ret<-  data.frame(AgeDx = AgeDx,
             Sex = Sex,
             Fusion = Fusion)

  } else{
    coef1<- mvtnorm::rmvnorm(n = 1,mean = c(0.03758116, -0.006658198), sigma =
                               matrix(c(0.05869754, -0.004873916, -0.004873916, 0.0004982868), nrow = 2) )
    coef2<- mvtnorm::rmvnorm(n = 1,mean = c(0.1107713, -0.05921511, -0.2586516), sigma =
                               matrix(c(0.07417177, -0.005295126, -0.02501394,
                                        -0.005295126, 0.0005518633, 0.0001718476,
                                        -0.02501394, 0.0001718476, 0.04914024), nrow = 3) )
    
    AgeDx = round(runif(n = n1, min = 0,max = 18))
    o1 <- exp(cbind(1, AgeDx) %*% t(coef1))
    Sex = rbinom(n1,1, prob = o1/(1+o1))
    o2 <- exp(cbind(1, AgeDx,Sex) %*% t(coef2))
    Fusion =  rbinom(n1,1, prob = o2/(1+o2))
    
    
    ret<-  data.frame(AgeDx = AgeDx,
                      Sex = Sex,
                      Fusion = Fusion)
    

    
  }
    return(ret)
}


