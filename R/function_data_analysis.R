
#function for data analysis - fit cutoff and sensitivity models
#dat.long: long data with id, visit time and biomarker values
#dat.short: short data with id, baseline covariates, event time and event indicator
#cutoff.type.basis: basis used for fitting biomarker cutoff
#sens.type.basis: basis for fitting senstivity
#covariate1: covariates used for biomarker cutoff model
#covariate2: covariates used for sensitivity model
#tau: vector of tau that we are interested in
#W: time-window
#nResap: the number of perturbation

plot_beta <- function(beta.t, beta.t.CI, visit.time, my.title = NULL){
  #need to make sure the sequence of beta and time.mat need to be consistent
  x <- visit.time
  #prepare data to plot
  min.y <- floor(min(beta.t.CI[, 1])) 
  max.y <- ceiling(max(beta.t.CI[, 2]))
  
  plot.dat <- unique(data.frame(x, beta.t, beta.t.CI))
  reorder <- order(plot.dat$x)
  plot(x = plot.dat$x[reorder], y = plot.dat$beta.t[reorder], 
       ylim = c(min.y, max.y), 
       type = "l",
       lwd = 1,
       main = my.title,
       font.main = 2,
       yaxt='n',xaxt='n',
       xlab = "s",
       ylab = ifelse(grepl("gamma",my.title), expression(paste(gamma[tau], " (s)")), expression(paste(beta[tau], " (s)"))))
  axis(1, cex.axis=1, tck=-0.02, lwd = 0.8)
  axis(2, cex.axis=1, tck=-0.02, lwd = 0.8)
  
  lines(x = plot.dat$x[reorder], plot.dat[,3][reorder], col = "black", lwd = 0.8, lty = 3)
  lines(x = plot.dat$x[reorder], plot.dat[,4][reorder], col = "black", lwd = 0.8, lty = 3)
  abline(h = 0, col = "#4292C6", lty = 2, lwd = 1)
}


calculate_beta_t <- function(beta, beta.Resap, index, time.mat){
  beta.t <- as.matrix(time.mat) %*% as.matrix(beta[index])
  beta.t.Resap <- as.matrix(time.mat) %*% as.matrix(beta.Resap[index, ])
  beta_CI <- t(apply(beta.t.Resap, 1, quantile, c(0.025, 0.975), na.rm = TRUE))
  return(list(beta.t = beta.t, beta.t.CI = beta_CI))
}


plot_main <- function(coefficients, beta.Resap, covariates, visit.time, titles, reverse = 0){
  coefficients <- as.data.frame(coefficients)
  time.mat = data.frame(1, 
                        logt = log(visit.time +0.1),
                        sqrtt = sqrt(visit.time + 0.1),
                        sqrtt_inv =sqrt(1/(visit.time + 0.1)))
  
  if(reverse == 1){
    beta_all <- -1 * rowMeans(coefficients)
    beta.Resap <- -1 * beta.Resap
  }else{
    beta_all <- rowMeans(coefficients)
  }
  
  check.ntau <- dim(coefficients)[2]
  coef.names = rownames(coefficients)
  if(check.ntau > 1){
    beta.Resap1 <- do.call(cbind, lapply(beta.Resap, rowMeans))
  }else{
    beta.Resap1 <- do.call(cbind, beta.Resap)
  }
  beta.Resap.check <- apply(beta.Resap1, 2, function(x){max(abs(x)) < 100}) #delete extreme values
  beta.Resap2 <- beta.Resap1[,beta.Resap.check]
  
  
  intercept_index <- which(coef.names %in% c("(Intercept)","logt","sqrtt","sqrtt_inv"))
  intercept.res <- calculate_beta_t(beta = beta_all, beta.Resap = beta.Resap2, index = intercept_index, time.mat)
  plot_beta(beta.t = intercept.res$beta.t, beta.t.CI = intercept.res$beta.t.CI, visit.time, my.title = titles[1])
  
  for (k in 1:length(covariates)) {
    tmp.pattern <- paste(covariates[k], "$", sep = "")
    tmp.index <- which(grepl(tmp.pattern, coef.names))
    tmp.res <- calculate_beta_t(beta_all, beta.Resap2, tmp.index, time.mat)
    plot_beta(beta.t = tmp.res$beta.t, beta.t.CI = tmp.res$beta.t.CI, visit.time, my.title = titles[k+1])
  }
  
}


analysis_main <- function(dat.long, dat.short, cutoff.type.basis, sens.type.basis, covariate1 = 1, covariate2 = 1, tau, W, nResap){
  model_results <- fit.2model.main(dat.long, dat.short, time.window, cutoff.type.basis, sens.type.basis, covariate1,
                                       covariate2, tau, nResap)
  vtime = sort(unique(dat.long$vtime))
  
  #extract coefficients
  cutoff_coef = coefficients(model_results$model.results$cutoff.model)
  sens_coef = as.data.frame(lapply(model_results$model.results$sensitivity.model, coefficients))
  
  cutoff_resap = lapply(model_results$resap.results, 
                        function(x){cutoff_resap = coefficients(x$cutoff.model)})
  sens_resap = lapply(model_results$resap.results, 
                      function(x){sens_resap = as.data.frame(lapply(x$sensitivity.model, coefficients))})
  sens_converge <- lapply(model_results$resap.results, 
                          function(x){sens_conv = as.data.frame(lapply(x$sensitivity.model,
                                                                       function(y){y$converged}))})
  
  sens_resap_clean <- mapply(function(x, y){Z = x[unlist(y)]}, sens_resap, sens_converge, SIMPLIFY = FALSE)
  sens_resap_clean2 <- lapply(sens_resap_clean, 
                              function(x){})
  
  
  
  #plot cutoff model time-dependent coef
  ncovari1 <- length(covariate1)
  cutoff_plot <- function(){
    tau_values <- unique(range(tau))
    plot_main(coefficients = cutoff_coef, beta.Resap = cutoff_resap, 
              covariates = covariate1, visit.time = vtime, 
              titles = unlist(sapply(0:ncovari1, function(i) bquote(beta[.(paste(i)) ~ ","~ tau ~ "= ["~.(paste(tau_values, collapse = ", ")) ~ "]"] ~ "(s)"))),
              reverse = 0)}
  
  #plot sens model time-dependent coef
  ncovari2 <- length(covariate2)
  sens_plot <- function(){ 
    tau_values <- unique(range(tau))
    plot_main(coefficients = sens_coef, beta.Resap = sens_resap_clean, 
             covariates = covariate2, visit.time = vtime, 
             titles = unlist(sapply(0:ncovari2, function(i) bquote(gamma[.(paste(i)) ~ ","~ tau ~ "= ["~.(paste(tau_values, collapse = ", ")) ~ "]"] ~ "(s)"))),
             reverse = 0)}
  return(list(model_results = model_results,
              cutoff_plots = cutoff_plot,
              sens_plots = sens_plot))
}



