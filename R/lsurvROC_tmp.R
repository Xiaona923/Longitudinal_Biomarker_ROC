lsurvROC_tmp <- function(model.weight, dat.long, dat.short, time.window, 
                         cutoff.type.basis, sens.type.basis, 
                         covariate1 = 1, covariate2 = 1, tau, nknot = NULL){
  
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