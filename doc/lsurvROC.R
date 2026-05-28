## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----install, eval=FALSE------------------------------------------------------
# library(remotes)
# remotes::install_github("Xiaona923/Longitudinal_Biomarker_ROC", force = TRUE)
# library(lsurvROC)

## ----load, include=FALSE------------------------------------------------------
library(lsurvROC)

## -----------------------------------------------------------------------------
data("example_data")

## ----message=FALSE, warning=FALSE---------------------------------------------
res = lsurvROC(dat.long = example_data$data.long, 
              dat.short = example_data$data.short,
              cutoff.type.basis = "FP",
              sens.type.basis = "FP", 
              covariate1 = c("Z", "Zcont"), 
              covariate2 = c("Z", "Zcont"), 
              tau = c(0.7, 0.8, 0.9), 
              time.window = 1, 
              nResap = 50,
              newdata = NULL
              )

par(mfrow=c(1,3))
plot(res)


## ----message=FALSE, warning=FALSE---------------------------------------------
res = lsurvROC(dat.long = example_data$data.long, 
              dat.short = example_data$data.short,
              cutoff.type.basis = "FP",
              sens.type.basis = "FP", 
              covariate1 = c("Z", "Zcont"), 
              covariate2 = c("Z", "Zcont"), 
              tau = seq(0.1, 0.9, 0.05), 
              time.window = 1, 
              nResap = 50,
              newdata = data.frame(vtime = 0.5, Z = 1, Zcont = 0.25)
              )

plot(res, ROC = TRUE)


