    remotes::install_github("Xiaona923/Longitudinal_Biomarker_ROC", force = TRUE)

    ## 
    ## ── R CMD build ─────────────────────────────────────────────────────────────────
    ##      checking for file ‘/private/var/folders/xf/l79029316mxg61sybqflzb_40000gn/T/RtmpTb4ZsZ/remotes96d66e377a38/Xiaona923-Longitudinal_Biomarker_ROC-735d48c/DESCRIPTION’ ...  ✔  checking for file ‘/private/var/folders/xf/l79029316mxg61sybqflzb_40000gn/T/RtmpTb4ZsZ/remotes96d66e377a38/Xiaona923-Longitudinal_Biomarker_ROC-735d48c/DESCRIPTION’
    ##   ─  preparing ‘lsurvROC’:
    ##      checking DESCRIPTION meta-information ...  ✔  checking DESCRIPTION meta-information
    ##   ─  checking for LF line-endings in source and make files and shell scripts
    ##   ─  checking for empty or unneeded directories
    ##   ─  building ‘lsurvROC_0.0.1.tar.gz’
    ##      
    ## 

    library(lsurvROC)

## Introduction

We propose a flexible regression framework to evaluate the prognostic
capability of longitudinal biomarker with survival outcome accounting
for the covariate and biomaker measurement time. Here is an example
using a simulated dataset to estimate the biomarker threshold at a
specified specificity level along with the corresponding sensitivity,
and to construct a covariate- and measurement time-specific ROC curve.

## Read in example data

    data("example_data")

## Estimate time-varying coefficients

-   dat.long: long format data, each subject may have multiple biomarker
    records
-   dat.short: short format data, each subject only has one record,
    including observed event time, event indicator, covariates (eg. Z,
    Zcont in this example)
-   cutoff.type.basis: type of basis function to estimate biomarker
    threshold, (FP: Fractional polynomial, BS: B Splines, linear: linear
    basis)
-   sens.type.basis: type of basis function to estimate sensitivity
    level, (FP: Fractional polynomial, BS: B Splines, linear: linear
    basis)
-   covariate1: a vector of covariate names to estimate the biomarker
    threshold
-   covariate2: a vector of covariate names to estimate the sensitivity
    level
-   tau: target specificity level, could be a single value or a set of
    values between 0 and 1
-   time.window: a time window for biomarker evaluation
-   nResap: number of perturbation resampling

<!-- -->

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

![](README_files/figure-markdown_strict/unnamed-chunk-3-1.png)![](README_files/figure-markdown_strict/unnamed-chunk-3-2.png)

## Conditional ROC curve

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

![](README_files/figure-markdown_strict/unnamed-chunk-4-1.png)![](README_files/figure-markdown_strict/unnamed-chunk-4-2.png)![](README_files/figure-markdown_strict/unnamed-chunk-4-3.png)![](README_files/figure-markdown_strict/unnamed-chunk-4-4.png)![](README_files/figure-markdown_strict/unnamed-chunk-4-5.png)![](README_files/figure-markdown_strict/unnamed-chunk-4-6.png)![](README_files/figure-markdown_strict/unnamed-chunk-4-7.png)
