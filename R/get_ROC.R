get_ROC <- function(model, basis, my.newdat, tau, tol = 1e3){
  #function to get an ROC curve
  predicted.sens <- unlist(lapply(model, 
                                  function(x){pred_sens(x$model, my.newdat, basis)}))
  model_sens_clean <- unlist(lapply(model, function(y){ 
    x = y$model
    x$converged  & !any(x$coefficients >= tol)
    }))
  
  origin.res <- data.frame(FalsePos = 1 - tau, pred.sens = predicted.sens) %>%
    filter(model_sens_clean == TRUE) %>% 
    arrange(FalsePos)
  
  mono.roc <-  monotone_ROC(FalsePos = origin.res$FalsePos, 
                            ROC.original = origin.res$pred.sens)
  
  auc.val <- get_AUC(mono.roc$FalsePos,
                     mono.roc$new_meas)
  return(list(ROC = mono.roc, AUC = auc.val))
}
