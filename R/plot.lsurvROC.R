#' @title Plot method for lsurvROC objects
#' @description Custom plot for lsurvROC class.
#' @param x a lsurvROC object
#' @param ROC a logical value that indicates whether plot the ROC curve
#' @method plot lsurvROC
#' @export
plot.lsurvROC <- function(x, ROC = FALSE, ...) {
  model_results <- x$model
  vtime <- x$param$times
  tau <- x$param$tau
  cutoff.type.basis <- x$param$cutoff.type.basis
  sens.type.basis <- x$param$sens.type.basis
  covariate1 <- x$param$covariate1
  covariate2 <- x$param$covariate2
  tol <- x$param$tol
  nknot <- x$param$nknot
  
  if(ROC == TRUE){
    data <- x$ROC
    f <- approxfun(data$FalsePos, data$new_meas, method = "constant", f = 0)
    par(pty="s")
    curve(f(x), from = 0, to = 1, ylim = c(0, 1), lwd = 1.5, 
          xlab = "",
          ylab = "", 
          axes=FALSE,
          frame=TRUE,
          main = paste0(sprintf("AUC = %.3f", x$AUC$AUC)), 
          ...)
    title(xlab = "1 - Specificity", line = 3)
    title(ylab = "Sensitivity", line = 3)
    axis(1, cex.axis=1, tck=-0.02, lwd = 0.8)
    axis(2, cex.axis=1, tck=-0.02, lwd = 0.8)
    
    
  }
  
  
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
  tau_values <- unique(range(tau))
  plot_beta_main(coefficients = cutoff_coef,
                 beta.Resap = cutoff_resap,
                 basis = cutoff.type.basis,
                 covariates = covariate1,
                 visit.time = vtime,
                 titles = unlist(sapply(0:ncovari1, function(i) bquote(beta[.(paste(i)) ~ ","~ tau ~ "= ["~.(paste(tau_values, collapse = ", ")) ~ "]"] ~ "(s)"))),
                 reverse = 0,
                 nknot,
                 tol = tol)
  
  #plot sens model time-dependent coef
  ncovari2 <- length(covariate2)
  tau_values <- unique(range(tau[sens_converge]))
  plot_beta_main(coefficients = sens_coef[sens_converge],
                 beta.Resap = sens_resap_clean,
                 basis = sens.type.basis,
                 covariates = covariate2,
                 visit.time = vtime,
                 titles = unlist(sapply(0:ncovari2, function(i) bquote(gamma[.(paste(i)) ~ ","~ tau ~ "= ["~.(paste(tau_values, collapse = ", ")) ~ "]"] ~ "(s)"))),
                 reverse = 0,
                 nknot, tol = tol)

  
  
}