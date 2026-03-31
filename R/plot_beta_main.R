
plot_beta_main <- function(coefficients, beta.Resap, basis,
                      covariates,
                      visit.time,
                      titles,
                      reverse = 0, nknot = NULL, tol = 100){
  
  coefficients <- as.data.frame(coefficients)
  if(basis == "FP"){
    time.mat = data.frame(1, 
                          logt = log(visit.time +0.1),
                          sqrtt = sqrt(visit.time + 0.1),
                          sqrtt_inv =sqrt(1/(visit.time + 0.1)))
  }else if(basis == "BS"){
    knot <- quantile(visit.time,  ((1:nknot) /(nknot + 1)))
    bs.func <- bs(visit.time, knots = knot)
    time.mat <- data.frame(1, bs.func)
  }else if(basis == "linear"){
    time.mat <- data.frame(1, vtime = visit.time)
  }
  
  if(reverse == 1){
    beta_all <- -1 * rowMeans(coefficients)
    beta.Resap <- -1 * beta.Resap
  }else{
    beta_all <- rowMeans(coefficients)
  }
  
  check.ntau <- dim(coefficients)[2]
  coef.names = rownames(coefficients)
  if(check.ntau > 1){
    beta.Resap1 <- do.call(cbind, lapply(beta.Resap, rowMeans, na.rm = T))
  }else{
    beta.Resap1 <- do.call(cbind, beta.Resap)
  }
  beta.Resap.check <- apply(beta.Resap1, 2, function(x){max(abs(x)) < tol}) #delete extreme values
  beta.Resap2 <- beta.Resap1[,beta.Resap.check]
  
  intercept_pattern <- paste0(paste(covariates, "$", sep = ""), collapse = "|")
  intercept_index <- which(!grepl(intercept_pattern, coef.names))
  intercept.res <- calculate_beta_t(beta = beta_all, beta.Resap = beta.Resap2, index = intercept_index, time.mat)
  plot_beta(beta.t = intercept.res$beta.t, beta.t.CI = intercept.res$beta.t.CI, visit.time, my.title = titles[1])
  
  for (k in 1:length(covariates)) {
    tmp.pattern <- paste(covariates[k], "$", sep = "")
    tmp.index <- which(grepl(tmp.pattern, coef.names))
    tmp.res <- calculate_beta_t(beta_all, beta.Resap2, tmp.index, time.mat)
    plot_beta(beta.t = tmp.res$beta.t, beta.t.CI = tmp.res$beta.t.CI, visit.time, my.title = titles[k+1])
  }
  
}

