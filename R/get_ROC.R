
get_ROC <- function(model, method = "both", cutoff.type.basis = NULL, sens.type.basis = NULL, my.newdat, tau, tol = 1e3){
  #function to get monotoned ROC curve
  #get original value for threshold & sensitivity
  predicted.cutoff <- pred_spec(model$cutoff.model$model, my.newdat, cutoff.type.basis)
  predicted.sens <- unlist(lapply(model$sensitivity.model, 
                                  function(x){pred_sens(x$model, my.newdat, sens.type.basis)}))
  #only keep converged results
  model_sens_clean <- unlist(lapply(model$sensitivity.model, function(y){ 
    x = y$model
    x$converged  & !any(x$coefficients >= tol)
  }))
  
  origin.res <- data.frame(tau = tau, FalsePos = 1 - tau, pred.sens = predicted.sens) %>%
    filter(model_sens_clean == TRUE) %>% 
    arrange(FalsePos)
  
  if(method == "sensitivity"){#only monotone sensitivity
    mono.roc <-  monotone_value(x = origin.res$FalsePos, 
                                y.origin = origin.res$pred.sens, ROC = TRUE)
    auc.val <- get_AUC(mono.roc$x,
                       mono.roc$y.mono)
    colnames(mono.roc) <- c("FalsePos", "TruePos")
    thres.out <- data.frame(tau = tau, threshold = as.vector(predicted.cutoff)) %>%
      filter(as.character(tau) %in% as.character(1-mono.roc$FalsePos))
    
    rownames(mono.roc) = 1:nrow(mono.roc)
    rownames(thres.out) = 1:nrow(thres.out)
    return(list(ROC = mono.roc, AUC = auc.val, threshold = thres.out))
  }else if(method == "threshold"){#only monotone threshold
    mono.threshold <-  monotone_value(x = tau, 
                                      y.origin = predicted.cutoff, ROC = FALSE)
    colnames(mono.threshold) <- c("tau", "threshold")
    roc.out <- origin.res %>% 
      filter(as.character(tau) %in% as.character(mono.threshold$tau)) %>%
      dplyr::select(FalsePos, pred.sens) %>% rename(TruePos = pred.sens) %>%
      dplyr::add_row(FalsePos = 0, TruePos = 0) %>%
      dplyr::add_row(FalsePos = 1, TruePos = 1) %>%
      arrange(FalsePos)
      
    auc.val <- get_AUC(roc.out$FalsePos, roc.out$TruePos)
    rownames(roc.out) = 1:nrow(roc.out)
    rownames(mono.threshold) = 1:nrow(mono.threshold)
    return(list(ROC = roc.out, AUC = auc.val, threshold = mono.threshold))
    
  }else if(method == "both"){#monotone both
    #monotone threshold 
    mono.threshold <-  monotone_value(x = tau, y.origin = predicted.cutoff, ROC = FALSE)
    origin.res.sub <- origin.res %>% 
      filter(as.character(tau) %in% as.character(mono.threshold$x))
    
    mono.roc <-  monotone_value(x = origin.res.sub$FalsePos, 
                                y.origin = origin.res.sub$pred.sens)
    
    auc.val <- get_AUC(mono.roc$x, mono.roc$y.mono)
    thres.out <- mono.threshold %>% filter(as.character(x) %in% as.character(1-mono.roc$x))
    
    colnames(mono.roc) <- c("FalsePos", "TruePos")
    colnames(thres.out) <- c("tau", "threshold")
    rownames(mono.roc) = 1:nrow(mono.roc)
    rownames(thres.out) = 1:nrow(thres.out)
    return(list(ROC = mono.roc, AUC = auc.val, threshold = thres.out))
  }
}

# get_ROC <- function(model, basis, my.newdat, tau, tol = 1e3){
#   #function to get an ROC curve
#   predicted.sens <- unlist(lapply(model, 
#                                   function(x){pred_sens(x$model, my.newdat, basis)}))
#   model_sens_clean <- unlist(lapply(model, function(y){ 
#     x = y$model
#     x$converged  & !any(x$coefficients >= tol)
#     }))
#   
#   origin.res <- data.frame(FalsePos = 1 - tau, pred.sens = predicted.sens) %>%
#     filter(model_sens_clean == TRUE) %>% 
#     arrange(FalsePos)
#   
#   mono.roc <-  monotone_ROC(FalsePos = origin.res$FalsePos, 
#                             ROC.original = origin.res$pred.sens)
#   
#   auc.val <- get_AUC(mono.roc$FalsePos,
#                      mono.roc$new_meas)
#   return(list(ROC = mono.roc, AUC = auc.val))
# }
