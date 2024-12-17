library(MESS)
source("Monotone_ROC.R")
get_AUC<-function(Abs,Ord){
  nobs<-length(Abs)
  dAbs<-Abs[-1]-Abs[-nobs]
  mil<-(Ord[-nobs]+Ord[-1])/2
  area<-sum(dAbs*mil)
  return(area)
}


get_ROC <- function(sensitivity.model, type.basis, my.newdat, tau){
  #function to get an ROC curve
  predicted.sens <- unlist(lapply(sensitivity.model, pred.sens, my.newdat, type.basis))
  model_sens_clean <- unlist(lapply(sensitivity.model, function(x) x$converged))
  
  origin.res <- data.frame(FalsePos = 1 - tau, pred.sens = predicted.sens) %>%
    filter(model_sens_clean == TRUE) %>% 
    arrange(FalsePos)
  
  mono.roc <-  monotone_ROC(FalsePos = origin.res$FalsePos, ROC.original = origin.res$pred.sens)
  auc.val <- get_AUC(mono.roc$FalsePos, mono.roc$new_meas)
  return(list(ROC = mono.roc, AUC = auc.val))
}

ROC.main <- function(my.newdat, dat.long, dat.short, tau, time.window, cutoff.type.basis, sens.type.basis, covariate1, covariate2, nResap){
  #fit model with perturbations
  model_results <- analysis_model_fit(dat.long, dat.short, time.window, cutoff.type.basis, sens.type.basis, covariate1,
                                       covariate2, tau, nResap)
  #get the main ROC curve
  ROC1 <- get_ROC(sensitivity.model = model_results$model.results$sensitivity.model, sens.type.basis, my.newdat, tau)
  #get perturbed ROC curves
  ROC.resap <- lapply(model_results$resap.results, 
                      function(x){get_ROC(x$sensitivity.model, sens.type.basis, my.newdat, tau)})
  AUC.sd = sd(unlist(lapply(ROC.resap, function(x){x$AUC})))
  

  return(list(ROC.results = ROC1, AUC.sd =AUC.sd, ROC.resap = ROC.resap))
}

plot_ROC <- function(data, my.add, my.col, my.lty, my.main =""){
  f <- approxfun(data$FalsePos, data$new_meas, method = "constant", f = 0)
  par(pty="s")
  curve(f(x), from = 0, to = 1, ylim = c(0, 1), lwd = 1.5, lty = my.lty,
        xlab = "",
        ylab = "", main = my.main,
        col = my.col,
        add = my.add,
        axes=FALSE,
        frame=TRUE)
  title(xlab = "1 - Specificity", line = 3)
  title(ylab = "Sensitivity", line = 3)
  axis(1, cex.axis=1, tck=-0.02, lwd = 0.8)
  axis(2, cex.axis=1, tck=-0.02, lwd = 0.8)
}


library(MESS)

get_res_tab <- function(result_list, cutoff, s_range, pkg_auc = NULL){
  s = unlist(lapply(strsplit(names(result_list), "_"), "[", 1))
  gender = unlist(lapply(strsplit(names(result_list), "_"), "[", 2))
  auc.res = sprintf("%.3f", unlist(lapply(result_list, "[", 2)))
  sd = sprintf("%.3f", unlist(lapply(lapply(result_list, "[", 3), function(x){sd(unlist(x))})))
  

  if(length(pkg_auc) == 0){res.tab <- data.frame(s, gender, cutoff = cutoff, s_range = s_range, auc.res, sd)}
  else{ res.tab <- data.frame(s, gender, cutoff = cutoff, s_range = s_range, auc.res, sd, pkg_auc)}
  return(res.tab)
}

