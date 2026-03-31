
#function for data analysis - fit cutoff and sensitivity models
#dat.long: long data with id, visit time and biomarker values
#dat.short: short data with id, baseline covariates, event time and event indicator
#cutoff.type.basis: basis used for fitting biomarker cutoff
#sens.type.basis: basis for fitting senstivity
#covariate1: covariates used for biomarker cutoff model
#covariate2: covariates used for sensitivity model
#tau: vector of tau that we are interested in
#time-window
#nResap: the number of perturbation

plot_beta <- function(beta.t, beta.t.CI, visit.time, my.title = NULL){
  #the sequence of beta and time.mat need to be consistent
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


plot_main <- function(coefficients, beta.Resap, basis, 
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




plot_coef_lsurvROC <- function(dat.long, dat.short,
                               cutoff.type.basis, sens.type.basis, 
                               covariate1 = 1, covariate2 = 1, 
                               tau, time.window, nResap, 
                               show_plots = TRUE, nknot = NULL, model = NULL, tol = 1e3){
  
  
  if(is.null(model)){
    model_results <- lsurvROC(dat.long, 
                              dat.short, 
                              time.window,
                              cutoff.type.basis,
                              sens.type.basis, 
                              covariate1, 
                              covariate2, 
                              tau, 
                              nResap,
                              nknot)
  }else{
    model_results <- model
  }
  
  vtime = sort(unique(dat.long$vtime))
  
  #extract coefficients
  cutoff_coef = coefficients(model_results$model.results$cutoff.model$model)
  sens_coef = as.data.frame(lapply(model_results$model.results$sensitivity.model,
                                   function(x) coefficients(x$model)))
  #check convergence and extreme values
  sens_converge = sapply(model_results$model.results$sensitivity.model, 
                         function(x){y = x$model; y$converged & !any(y$coefficients >= tol)})

  if (any(!sens_converge)) {
    warning(
      paste0(sum(!sens_converge),
             " sensitivity model(s) did not converge; delete tau = ",
             paste0(tau[which(!sens_converge)],collapse = ", "))
    )
  }
  
  cutoff_resap = lapply(model_results$resap.results, 
                        function(x){
                          cutoff_resap = coefficients(x$cutoff.model$model)
                          })
  sens_resap = lapply(model_results$resap.results, 
                      function(x){sens_resap = as.data.frame(
                        lapply(x$sensitivity.model, function(x) coefficients(x$model))
                                                             )})
  sens_resap_converge <- lapply(model_results$resap.results, 
                                function(x){sens_conv = as.data.frame(
                                  lapply(x$sensitivity.model, function(y){y$converged})
                                  )})
  
  #retain results only if the model has converged
  
  sens_resap_clean <- mapply(function(x, y){
                            extreme_val <- apply(x, 2, function(q) any(abs(q) > tol)) 
                            x[ , !y | extreme_val] <- NA
                            Z = x[, sens_converge]
                              }, 
                       sens_resap, sens_resap_converge, SIMPLIFY = FALSE)

  
 
  #plot cutoff model time-dependent coef
  ncovari1 <- length(covariate1)
  cutoff_plot <- function(){
    tau_values <- unique(range(tau))
    plot_main(coefficients = cutoff_coef, 
              beta.Resap = cutoff_resap, 
              basis = cutoff.type.basis,
              covariates = covariate1,
              visit.time = vtime, 
              titles = unlist(sapply(0:ncovari1, function(i) bquote(beta[.(paste(i)) ~ ","~ tau ~ "= ["~.(paste(tau_values, collapse = ", ")) ~ "]"] ~ "(s)"))),
              reverse = 0, 
              nknot, 
              tol)}
  
  #plot sens model time-dependent coef
  ncovari2 <- length(covariate2)
  sens_plot <- function(){ 
    tau_values <- unique(range(tau[sens_converge]))
    plot_main(coefficients = sens_coef[sens_converge], 
              beta.Resap = sens_resap_clean, 
              basis = sens.type.basis,
              covariates = covariate2,
              visit.time = vtime, 
              titles = unlist(sapply(0:ncovari2, function(i) bquote(gamma[.(paste(i)) ~ ","~ tau ~ "= ["~.(paste(tau_values, collapse = ", ")) ~ "]"] ~ "(s)"))),
              reverse = 0, 
              nknot, tol)}
  
  if(isTRUE(show_plots)){
    print(cutoff_plot()) 
    print(sens_plot())
    
  }
  
  invisible((list(cutoff_plots = cutoff_plot,
                  sens_plots = sens_plot)))
}


