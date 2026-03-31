#' lsurvROC
#'
#' It runs flexible regression models to estimate biomarker thresholds at selected quantile(s) and the corresponding sensitivities.
#' It returns model results
#' and plots the time-varying coefficients with 95% confidence interval.
#'
#' @param dat.long a data frame with long data format. It should contain subject id (id), longitudinal marker (Xt), measurement time (vtime)
#' @param dat.short a data frame. It should contain unique subject id (id), observed event time (Y), event indicator (delta, 1 = observed event, 0 = censored), covariates
#' @param cutoff.type.basis basis function to estimate marker threshold ("FP" = fractional polynomial, "bs" = "B splines", "linear" = linear basis)
#' @param sens.type.basis basis function to estimate sensitivity levels ("FP" = fractional polynomial, "bs" = "B splines", "linear" = linear basis)
#' @param covariate1 a vector of name of covariates to estimate marker thresholds at level tau(s)
#' @param covariate2 a vector of name of covariates to estimate sensitivity levels
#' @param tau the quantile(s) of marker threshold to be estimated
#' @param time.window a numeric value. This value will be used to specify the future time point to perform the evaluation, meaning vtime + time.window 
#' @param nResap a numeric value to perform perturbation resampling. The default value = 200. Larger value can generate more accurate ASE value
#' @param show_plots if TRUE then the time-varying coefficients plot will be displayed
#' @param nknot if basis function is B splines, need to specify the number of knot
#' @param tol a constant used for removing unconverged results
#' @examples
#'   data(example_data)
#'   head(example_data$data.long)
#'   head(example_data$data.short)
#'   
#'   res <- lsurvROC(dat.long = example_data$data.long, 
#'     dat.short = example_data$data.short,
#'     cutoff.type.basis = "FP",
#'     sens.type.basis = "FP", 
#'     covariate1 = c("Z", "Zcont"), 
#'     covariate2 = c("Z", "Zcont"), 
#'     tau = c(0.7, 0.8, 0.9), 
#'     time.window = 1, 
#'     nResap = 50, 
#'     show_plots = TRUE
#'     )
#' @import ggplot2
#' @import dplyr
#' @importFrom quantreg rq
#' @import survival
#' @import splines
#' @export
lsurvROC <- function(dat.long, dat.short,
                               cutoff.type.basis, sens.type.basis, 
                               covariate1 = 1, covariate2 = 1, 
                               tau, time.window, nResap = 200, 
                               newdata = NULL,
                               nknot = NULL,
                               tol = 1e3
                     ){
 #fit models
  model_results <- lsurvROC_main(dat.long, 
                            dat.short, 
                            time.window,
                            cutoff.type.basis,
                            sens.type.basis, 
                            covariate1, 
                            covariate2, 
                            tau, 
                            nResap,
                            nknot)
 #get all visit time 
  vtime = sort(unique(dat.long$vtime))
  
 #get ROC curve if newdata == TRUE
  if(!is.null(newdata)){
    #get the main ROC curve
    my.ROC <- get_ROC(model = model_results$model.results$sensitivity.model,
                      basis = sens.type.basis, 
                      my.newdat = newdata,
                      tau = tau, 
                      tol = tol)
    rownames(my.ROC$ROC) <- NULL
    
    #get perturbed ROC curves
    ROC.resap <- lapply(model_results$resap.results, 
                        function(x){get_ROC(x$sensitivity.model,
                                            basis = sens.type.basis, 
                                            my.newdat = newdata, 
                                            tau = tau)})
    
    AUC.sd = sd(unlist(lapply(ROC.resap, "[", "AUC")))
    
    output <- list(model = model_results,
                   ROC = my.ROC$ROC,
                   AUC = list(AUC = my.ROC$AUC, sd = AUC.sd),
                   param = list(times = vtime,
                                taus = tau,
                                cutoff.type.basis = cutoff.type.basis,
                                sens.type.basis = sens.type.basis, 
                                covariate1 = covariate1, 
                                covariate2 = covariate2,
                                nknot = nknot, 
                                tol = tol))
    
  }else{
    output <- list(model = model_results,
                   param = list(times = vtime,
                                taus = tau,
                                cutoff.type.basis = cutoff.type.basis,
                                sens.type.basis = sens.type.basis, 
                                covariate1 = covariate1, 
                                covariate2 = covariate2,
                                nknot = nknot, 
                                tol = tol))
                   
  }
  

                
  
  class(output) <- "lsurvROC"
  return(output)
  
}

#' @keywords package function 
