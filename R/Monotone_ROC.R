monotone_ROC <- function(FalsePos, ROC.original){
  startTau = 0.5
  leftRec <- data.frame(FalsePos = rep(NA, length(FalsePos)),
                        new_meas = rep(NA, length(FalsePos)))
  rightRec <- data.frame(FalsePos = rep(NA, length(FalsePos)),
                         new_meas = rep(NA, length(FalsePos)))
  
  myi = 1
  
  sens.orig = ROC.original
  
  #find the closest point from startTau
  leftRec$FalsePos[myi] = FalsePos[which.min(abs(startTau-FalsePos))]
  leftRec$new_meas[myi] = sens.orig[which.min(abs(startTau-FalsePos))]
  new_idx = 1
  
  while (!is.na(new_idx)) {
    tmp <- which(sens.orig <= leftRec$new_meas[myi] & FalsePos < leftRec$FalsePos[myi])
    len = length(tmp)
    if (len >= 1) {
      new_idx <- max(tmp)
      myi <- myi+1
      leftRec$FalsePos[myi] <- FalsePos[new_idx]
      leftRec$new_meas[myi] <- sens.orig[new_idx]
    } else {
      new_idx = NA
    }
  }
  
  leftRec <- na.omit(leftRec)
  
  myi = 1
  rightRec$FalsePos[myi] = startTau
  rightRec$new_meas[myi] = sens.orig[which.min(abs(startTau-FalsePos))]
  
  new_idx = 1
  while (!is.na(new_idx)) {
    tmp <- which(sens.orig >= rightRec$new_meas[myi] & FalsePos > rightRec$FalsePos[myi])
    len = length(tmp)
    if (len >= 1) {
      new_idx <- min(tmp)
      myi <- myi+1
      rightRec$FalsePos[myi] <- FalsePos[new_idx]
      rightRec$new_meas[myi] <- sens.orig[new_idx]
    } else {
      new_idx = NA
    }
  }
  
  rightRec <- na.omit(rightRec)
  
  
  monoRes <- rbind(leftRec, rightRec[-1,]) %>%
    arrange(FalsePos)
  
  if(min(monoRes$FalsePos) ==0){
    monoRes$new_meas[monoRes$FalsePos == 0] = 0
    mono_roc2 <- monoRes
  }else{
    mono_roc2 <- monoRes %>% add_row(FalsePos = 0, new_meas = 0)
  }
  
  if(max(mono_roc2$FalsePos) == 1){
    mono_roc2$new_meas[mono_roc2$FalsePos == 1] = 1
    mono_roc3 <- mono_roc2
  }else{
    mono_roc3 <- mono_roc2 %>% add_row(FalsePos = 1, new_meas = 1) %>% arrange(FalsePos)
  }
  return(mono_roc3)
}

