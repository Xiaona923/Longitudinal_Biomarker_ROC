# Longitudinal_Biomarker_ROC
Covariate- and Measurement Time-Specific Evaluation of ROC curve with Survival Outcome for Longitudinal Biomarkers

#--------------CODE-----------------#

**function_model_fit.R: contains functions to estimate biomarker threshold and corresponding sensitivity levels**
  - process_dat(): process short and long data for modeling
  - fit.spec(): estimate biomarker threshold conditioned on covariates and tau
  - pred.spec(): get the fitted values at different biomarker measurement time
  - fit.sens(): estimate corresponding sensitivity levels 
  - pred.sens(): get the fitted values
  - analysis_model_fit(): warp-up function for the model fitting
    
**function_data_analysis.R: contains functions to summarize the model results**
  - plot_beta(): intermediate function to plot the time-varying beta and 95% CI
  - calculate_beta_t(): intermediate function for plotting
  - plot_main(): function to plot the results
  - analysis_main(): KEY FUNCTION, warp-up all functions to fit the models and plot the results  

**monotone_ROC.R**
 function to monotone the ROC curve
 
**Get_ROC.R**

  - get_ROC()
  - ROC.main()
  - plot_ROC()

#--------------data-----------------#
**simulation data under regular visits**
- reg_data_sim_long.csv: long format data, each subject has multiple records
- reg_data_sim_short.csv:  short format data, each subject has 1 records
