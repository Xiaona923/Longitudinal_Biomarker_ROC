calculate_beta_t <- function(beta, beta.Resap, index, time.mat){
  beta.t <- as.matrix(time.mat) %*% as.matrix(beta[index])
  beta.t.Resap <- as.matrix(time.mat) %*% as.matrix(beta.Resap[index, ])
  beta_CI <- t(apply(beta.t.Resap, 1, quantile, c(0.025, 0.975), na.rm = TRUE))
  return(list(beta.t = beta.t, beta.t.CI = beta_CI))
}
