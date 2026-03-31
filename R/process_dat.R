process_dat <- function(my.weights, dat.long, dat.short, time.window){
  #function to process data for modeling
  #my.weights: without pertubation = 1, with perturbation draw from exp(1)
  #dat.long: long data, each patients can have multiple row representing different visit time and biomarker measurements
  #dat.short: short data, each paitents only have one recode with id, event time, delta, baseline covariates
  #time.window: pre-specified time window
  #id: the column of subject id
  
  dat.short$delta.censor <- 1- dat.short$delta
  dat.long$W <- dat.long$vtime + time.window
  dat.long$logt = log(dat.long$vtime + 0.1)
  dat.long$sqrtt = sqrt(dat.long$vtime + 0.1)
  dat.long$sqrtt_inv = sqrt(1/(dat.long$vtime + 0.1))
  
  dat.short$wt_rsap <- my.weights
  km.fit <- survfit(Surv(Y, delta.censor) ~ 1, data = dat.short, weights = wt_rsap) 
  survest <- stepfun(km.fit$time, c(1, km.fit$surv)) 
  dat.short$G_Y <- pmax(survest(dat.short$Y), 0.05)
  dat.long$G_s <- pmax(survest(dat.long$vtime), 0.05)
  
  #merge two data set
  dat.full <- left_join(dat.long, dat.short, by = c("id" = "id")) #create formatted time for model fitting 
  
  #specificity
  dat.spec <- dat.full %>% filter(W <= Y)  %>%
    mutate(mod.weight = wt_rsap)
  #sensitivity
  dat.sens <- dat.full %>% filter(Y <= W) %>%
    mutate(mod.weight = wt_rsap * delta/(G_Y/G_s))
  
  return(list(cutoff.data = dat.spec, 
              sensitivity.data = dat.sens))
  
}
