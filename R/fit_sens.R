fit_sens <- function(outcome_t, data, covariates, type.basis, nknot = NULL){

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
  }else if(type.basis == "BS"){
    knot <- quantile(unique(data$vtime), ((1:nknot) /(nknot + 1)))
    part1 = "(1 + bs(vtime, knots = knot)) *"
    form <- as.formula(paste("outcome_t ~", part1 , part2))
  }
  
  suppressWarnings(mod <- glm(form, data = data, weights = mod.weight, family = binomial(link = "logit")))
  
  return(list(model = mod, 
              basis = type.basis))
}
