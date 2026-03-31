pred_spec <- function(model, newdata, type.basis){
  #function to calculate fitted threshold value
  #model: results from fit.spec()
  #newdata: newdata for prediction
  #type.basis: type of basis used in fitting model
  if(type.basis == "FP"){
    newdata$logt = log(newdata[,"vtime"]+0.1)
    newdata$sqrtt = sqrt(newdata[,"vtime"]+0.1)
    newdata$sqrtt_inv = 1/(sqrt(newdata[,"vtime"]+0.1))
  }
  predicted <- predict(model, newdata = newdata)
  return(predicted)
}

