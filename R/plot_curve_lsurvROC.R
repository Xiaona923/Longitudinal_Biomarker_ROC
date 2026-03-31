plot_curve_lsurvROC <- function(model, my.newdat, tau, basis, tol = 1e3, add = FALSE,
                                col = "black", lty = 1, main = NULL){

  model_results <- model
  
  #get the main ROC curve
  my.ROC <- get_ROC(model = model_results$model.results$sensitivity.model,
                  basis, 
                  my.newdat,
                  tau, 
                  tol)
  rownames(my.ROC$ROC) <- NULL
  
  #get perturbed ROC curves
  ROC.resap <- lapply(model_results$resap.results, 
                      function(x){get_ROC(x$sensitivity.model, basis, my.newdat, tau)})
  
  AUC.sd = sd(unlist(lapply(ROC.resap, "[", "AUC")))
  
  plot_ROC(my.ROC$ROC, my.add = add, my.col = col, my.lty = lty, my.main = main)
  
  return(list(ROC.results = my.ROC, 
              AUC.sd =AUC.sd,
              ROC.resap = ROC.resap))
}
