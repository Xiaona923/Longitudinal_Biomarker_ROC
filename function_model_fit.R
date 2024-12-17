
process_dat <- function(my.weights, dat.long, dat.short, time.window){
  #function to process data for modeling
  #my.weights: without pertubation = 1, with perturbation draw from exp(1)
  #dat.long: long data, each patients can have multiple row representing different visit time and biomarker measurements
  #dat.short: short data, each paitents only have one recode with id, event time, delta, baseline covariates
  #time.window: pre-specified time window
  
  dat.short$wt_rsap <- my.weights
  km.fit <- survfit(Surv(Y, delta.censor) ~ 1, data = dat.short, weights = wt_rsap) 
  survest <- stepfun(km.fit$time, c(1, km.fit$surv)) 
  dat.short$G_Y <- pmax(survest(dat.short$Y), 0.05)
  dat.long$G_s <- pmax(survest(dat.long$vtime), 0.05)
  
  #merge two data set
  dat.full <- left_join(dat.long, dat.short, by = c("obs_id" = "obs_id")) #create formatted time for model fitting 
  
  #return processed data
  ## specificity
  dat.spec <- dat.full %>% filter(W <= Y)
  
  ## sensitivity
  #filter data - sensitivity
  dat.sens <- dat.full %>% filter(Y <= W) %>%
    mutate(mod.weight = wt_rsap * delta/(G_Y/G_s))
  return(list(cutoff.data = dat.spec, sensitivity.data = dat.sens))
  
}



fit.spec <- function(data, covariates, type.basis, tau){
  #function to estimate targeted threshold under a give tau(s)
  #data: cutoff.data from process_dat()
  #covariates: a vector of covariates name, baseline covariates adjusted in the model
  #type.basis: could be FP, linear, constant, intercept
  #tau: numbers, a single value or set of tau values from (0, 1)
  
  part2 = paste("(", paste(covariates, collapse = "+"), ")")
  
  if(type.basis == "FP"){
    part1 = "(1  + logt + sqrtt +sqrtt_inv) *"
    form <- as.formula(paste("Xt ~", part1 , part2))
  }else if(type.basis == "linear"){
    part1 = "(1 + vtime) *"
    form <- as.formula(paste("Xt ~", part1, part2))
  }else if(type.basis == "constant"){
    form <- as.formula(paste("Xt ~", part2))
  }else if(type.basis == "intercept"){
    form <- as.formula(paste("Xt ~", 1))
  }
  mod <- rq(form, data = data, tau = tau, weights = wt_rsap)
  
  return(mod)
}


pred.spec <- function(model, newdata, type.basis){
  #function to calculate fitted threshold value
  #model: results from fit.spec()
  #newdata: newdata for prediction
  #type.basis: type of basis used in fitting model
  if(type.basis == "FP"){
    newdata$logt = log(newdata[,"vtime"]+0.1)
    newdata$sqrtt = sqrt(newdata[,"vtime"]+0.1)
    newdata$sqrtt_inv = 1/(sqrt(newdata[,"vtime"]+0.1))
  }
  predicted <- predict(model, newdata = newdata)
  return(predicted)
}


fit.sens <- function(outcome_t, data, covariates, type.basis){
  #function to estimate sensitivity level under a give tau(s)
  #outcome_t: a vector with all 0, 1 values, 1 = death, 0 = alive
  #data: sensitivity.data from process_dat()
  #covariates: a vector of covariates name, baseline covariates adjusted in the model
  #type.basis: could be FP, linear, constant, intercept
  
  part2 = paste("(", paste(covariates, collapse = "+"), ")")
  
  if(type.basis == "FP"){
    part1 = "(1  + logt + sqrtt +sqrtt_inv) *"
    form <- as.formula(paste("outcome_t ~", part1, part2))
  }else if(type.basis == "linear"){
    part1 = "(1 + vtime) *"
    form <- as.formula(paste("outcome_t ~", part1, part2))
  }else if(type.basis == "constant"){
    #covariates a vector of covariates name
    form <- as.formula(paste("outcome_t ~", part2))
  }else if(type.basis == "intercept"){
    form <- as.formula(paste("outcome_t ~", 1))
  }
  
  mod <- glm(form, data = data, weights = mod.weight, family = binomial(link = "logit"))
  
  return(mod)
}


pred.sens <- function(model, newdata, type.basis){
  #function to calculate fitted threshold value
  #model: results from fit.sens()
  #newdata: newdata for prediction
  #type.basis: type of basis used in fitting model
  if(type.basis == "FP"){
    newdata$logt = log(newdata[,"vtime"]+0.1)
    newdata$sqrtt = sqrt(newdata[,"vtime"]+0.1)
    newdata$sqrtt_inv = 1/(sqrt(newdata[,"vtime"]+0.1))
  }
  predicted <- predict(model, newdata = newdata, type="response")
  return(predicted)
}


main <- function(model.weight, dat.long, dat.short, time.window, 
                 cutoff.type.basis, sens.type.basis, 
                 covariate1, covariate2, tau, newdata = NULL){
  #function to fit threshold model and sensitivity model as well as output the predicted values
  # model.weight: without pertubation = 1, with perturbation draw from exp(1)
  # dat.long: long data, each patients can have multiple row representing different visit time and biomarker measurements
  # dat.short: short data, each paitents only have one recode with id, event time, delta, baseline covariates
  # time.window: pre-specified time window
  # cutoff.type.basis: type of basis function for cutoff (threshold) model, could be FP, linear, constant, intercept
  # sens.type.basis: type of basis function for sensitivity model, could be FP, linear, constant, intercept
  # covariate1: a vector of covariates name, baseline covariates adjusted in the cutoff (threshold) model
  # covariate2: a vector of covariates name, baseline covariates adjusted in the sensitivity model
  # tau: numbers, a single value or set of tau values from (0, 1)
  # newdata: data.frame for data prediction. default is NULL.
  
  processed.dat <- process_dat(model.weight, dat.long, dat.short, time.window)
  cutoff.data <- processed.dat$cutoff.data
  sensitivity.data <- processed.dat$sensitivity.data
  #cutoff of specificity  
  spec.fit = fit.spec(cutoff.data, covariate1, type.basis = cutoff.type.basis, tau)
  #prepare data
  design.mat <- model.matrix(formula(spec.fit), data = sensitivity.data)
  est.cutoff <- design.mat %*% spec.fit$coefficients
  est.outcome <- apply(sensitivity.data$Xt >= est.cutoff, 2, as.numeric)
  #estimated sensitivity
  sens.fit <- apply(est.outcome, 2, fit.sens, data = sensitivity.data, covariates = covariate2, type.basis = sens.type.basis)
  
  if(length(newdata) > 0){
    cutoff.pred = pred.spec(model = spec.fit, newdata = newdata, type.basis = cutoff.type.basis)
    sensitivity.pred = pred.sens(model = sens.fit, newdata = newdata, type.basis = sens.type.basis)
    return(list(cutoff.model = spec.fit, cutoff.pred = cutoff.pred, 
                sensitivity.model = sens.fit, sens.pred = sensitivity.pred))
  }
  return(list(cutoff.model = spec.fit, sensitivity.model = sens.fit))
}


analysis_model_fit <- function(dat.long, dat.short, time.window, cutoff.type.basis, sens.type.basis, covariate1, covariate2, tau, nResap){
  #function to conduct a full data analysis with variance estimation (perturbation)
  #nResap: number of perturbation, the larger the number, the better precision of variance estimation
  n <-  nrow(dat.short)
  dat.short$delta.censor <- 1- dat.short$delta
  dat.long$W <- dat.long$vtime + time.window
  dat.long$logt = log(dat.long$vtime + 0.1)
  dat.long$sqrtt = sqrt(dat.long$vtime + 0.1)
  dat.long$sqrtt_inv = sqrt(1/(dat.long$vtime + 0.1))
  
  resap.weight.mat <- matrix(rexp(n*nResap, 1), nrow =n, ncol = nResap) #generate weight matrix for pertubation
  #fit models 
  res <- main(1, dat.long, dat.short, time.window, cutoff.type.basis, sens.type.basis, covariate1,
              covariate2, tau)
  
  #perturbation
  resap.res <- apply(resap.weight.mat, 2, main, dat.long, dat.short, time.window, cutoff.type.basis, sens.type.basis, covariate1,
                     covariate2, tau)
  
  return(list(model.results = res, resap.results = resap.res))
}


