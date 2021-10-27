gen_synpop <- function(l1, n1=  100000){
  #l1: list of mean ,sd agedx, and neseted glms fit to fusion and sex
  
ret<-  data.frame(AgeDx = round(runif(n = n1, min = 0,max = 18)),
             Sex = as.numeric(runif(n= n1, min= 0,max =1) < 0.4931129),
             Fusion = as.numeric(runif(n= n1, min = 0, max =1 ) < 0.3581267)
             )
    return(ret)
}
