monotone_value <- function(x, y.origin, ROC = TRUE){
  startTau = 0.5
  leftRec <- data.frame(x = rep(NA, length(x)),
                        y.mono = rep(NA, length(x)))
  rightRec <- data.frame(x = rep(NA, length(x)),
                         y.mono = rep(NA, length(x)))
  
  myi = 1
  
  y.orig = y.origin
  
  #find the closest point from startTau
  leftRec$x[myi] = x[which.min(abs(startTau-x))]
  leftRec$y.mono[myi] = y.orig[which.min(abs(startTau-x))]
  new_idx = 1
  
  while (!is.na(new_idx)) {
    tmp <- which(y.orig <= leftRec$y.mono[myi] & x < leftRec$x[myi])
    len = length(tmp)
    if (len >= 1) {
      new_idx <- max(tmp)
      myi <- myi+1
      leftRec$x[myi] <- x[new_idx]
      leftRec$y.mono[myi] <- y.orig[new_idx]
    } else {
      new_idx = NA
    }
  }
  
  leftRec <- na.omit(leftRec)
  
  myi = 1
  rightRec$x[myi] = startTau
  rightRec$y.mono[myi] = y.orig[which.min(abs(startTau-x))]
  
  new_idx = 1
  while (!is.na(new_idx)) {
    tmp <- which(y.orig >= rightRec$y.mono[myi] & x > rightRec$x[myi])
    len = length(tmp)
    if (len >= 1) {
      new_idx <- min(tmp)
      myi <- myi+1
      rightRec$x[myi] <- x[new_idx]
      rightRec$y.mono[myi] <- y.orig[new_idx]
    } else {
      new_idx = NA
    }
  }
  
  rightRec <- na.omit(rightRec)
  
  
  monoRes <- rbind(leftRec, rightRec[-1,]) %>%
    arrange(x)
  
  if(ROC){
    if(min(monoRes$x) ==0){
      monoRes$y.mono[monoRes$x == 0] = 0
      mono_roc2 <- monoRes
    }else{
      mono_roc2 <- monoRes %>% add_row(x = 0, y.mono = 0)
    }
    
    if(max(mono_roc2$x) == 1){
      mono_roc2$y.mono[mono_roc2$x == 1] = 1
      mono_roc3 <- mono_roc2
    }else{
      mono_roc3 <- mono_roc2 %>% add_row(x = 1, y.mono = 1) %>% arrange(x)
    }
    return(mono_roc3)
  }else{
    return(monoRes)
  }
  
}

