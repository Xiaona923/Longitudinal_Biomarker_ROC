
lsurvROC.tmp <- function(model.weight, dat.long, dat.short, time.window, 
                           cutoff.type.basis, sens.type.basis, 
                           covariate1 = 1, covariate2 = 1, tau, nknot = NULL){
  #function to fit threshold model and sensitivity model as well as output the predicted values
  # model.weight: without pertubation = 1, with perturbation draw from exp(1)
  # dat.long: long data, each patients can have multiple row representing different visit time and biomarker measurements
  # dat.short: short data, each paitents only have one recode with id, event time, delta, baseline covariates
  # time.window: pre-specified time window
  # cutoff.type.basis: type of basis function for cutoff (threshold) model, could be FP, linear, constant, intercept
  # sens.type.basis: type of basis function for sensitivity model, could be FP, linear, constant, intercept
  # covariate1: a vector of covariates name, baseline covariates adjusted in the cutoff (threshold) model
  # covariate2: a vector of covariates name, baseline covariates adjusted in the sensitivity model
  # tau: a single value or set of tau values from (0, 1)
  # newdata: data.frame for data prediction. default is NULL.
  
  processed.dat <- process_dat(model.weight, dat.long, dat.short, time.window)
  cutoff.data <- processed.dat$cutoff.data
  sensitivity.data <- processed.dat$sensitivity.data
  #cutoff of specificity  
  spec.fit0 = fit_spec(cutoff.data, covariate1, type.basis = cutoff.type.basis, tau, nknot = nknot)
  spec.fit = spec.fit0$model
  #prepare data
  design.mat <- model.matrix(formula(spec.fit), data = sensitivity.data)
  est.cutoff <- design.mat %*% spec.fit$coefficients
  est.outcome <- apply(sensitivity.data$Xt >= est.cutoff, 2, as.numeric)
  #estimated sensitivity
  sens.fit <- apply(est.outcome, 2, fit_sens, data = sensitivity.data, covariates = covariate2, 
                    type.basis = sens.type.basis, nknot = nknot)
  
  return(list(cutoff.model = spec.fit0,
              sensitivity.model = sens.fit))
    
  
}



lsurvROC <- function(dat.long, dat.short, time.window,
                     cutoff.type.basis, sens.type.basis,
                     covariate1, covariate2, tau, nResap = 200, nknot = NULL){
  #function to conduct a full data analysis with variance estimation (perturbation)
  #nResap: number of perturbation, default = 200, the larger the number, the better precision of variance estimation
  n <-  nrow(dat.short)
  dat.short$delta.censor <- 1- dat.short$delta
  dat.long$W <- dat.long$vtime + time.window
  dat.long$logt = log(dat.long$vtime + 0.1)
  dat.long$sqrtt = sqrt(dat.long$vtime + 0.1)
  dat.long$sqrtt_inv = sqrt(1/(dat.long$vtime + 0.1))
  
  resap.weight.mat <- matrix(rexp(n*nResap, 1), nrow =n, ncol = nResap) #generate weight matrix for pertubation
  #fit models 
  res <- lsurvROC.tmp(1, dat.long, dat.short, time.window, 
                      cutoff.type.basis, sens.type.basis, 
                      covariate1,covariate2, tau, nknot)
  
  #perturbation
  resap.res <- apply(resap.weight.mat, 2, lsurvROC.tmp, dat.long, dat.short, 
                     time.window, cutoff.type.basis, sens.type.basis, covariate1,
                     covariate2, tau, nknot)
  
  return(list(model.results = res, resap.results = resap.res))
}

