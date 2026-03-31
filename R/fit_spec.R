fit_spec <- function(data, covariates, type.basis, tau, nknot = NULL){
  #function to estimate targeted threshold under a give tau(s) (1-specificity)
  #data: cutoff.data from process_dat()
  #covariates: a vector of covariates name, baseline covariates adjusted in the model
  #type.basis: could be FP, linear, constant, intercept, BS
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
  }else if(type.basis == "BS"){
    knot <- quantile(unique(data$vtime), ((1:nknot) /(nknot + 1)))
    part1 = "(1 + bs(vtime, knots = knot)) *"
    form <- as.formula(paste("Xt ~", part1 , part2))
  }
  mod <- rq(form, data = data, tau = tau, weights = mod.weight)
  
  return(list(model = mod, 
              basis = type.basis,
              tau = tau))
}
