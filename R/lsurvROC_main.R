
lsurvROC_main <- function(dat.long, dat.short, time.window,
                     cutoff.type.basis, sens.type.basis,
                     covariate1, covariate2, tau, nResap = 200, nknot = NULL){
  n <-  nrow(dat.short)
  # dat.short$delta.censor <- 1- dat.short$delta
  # dat.long$W <- dat.long$vtime + time.window
  # dat.long$logt = log(dat.long$vtime + 0.1)
  # dat.long$sqrtt = sqrt(dat.long$vtime + 0.1)
  # dat.long$sqrtt_inv = sqrt(1/(dat.long$vtime + 0.1))
  
  resap.weight.mat <- matrix(rexp(n*nResap, 1), nrow =n, ncol = nResap) #generate weight matrix for pertubation
  #fit models 
  res <- lsurvROC_tmp(1, dat.long, dat.short, time.window, 
                      cutoff.type.basis, sens.type.basis, 
                      covariate1,covariate2, tau, nknot)
  
  #perturbation
  resap.res <- apply(resap.weight.mat, 2, lsurvROC_tmp, dat.long, dat.short, 
                     time.window, cutoff.type.basis, sens.type.basis, covariate1,
                     covariate2, tau, nknot)
  
  return(list(model.results = res, resap.results = resap.res))
}

