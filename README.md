# Bayesian-Analysis-of-Survival-Data-with-Missing-Censoring-Indicators
Code for paper with the same name

Code is organized in two folders, "Simulations" and "ApplyToOPPERA", each with subfolders based on either the current Bayesian method or previous methods. 

The "Simulations" folder contains code to run the simulations. Together, these files create Table 1.
The "Simulations_Bayes" subfolder contains code to run the simulations on the Bayesian method. The "sims_FB.R" file runs the simulations based on the model, programmed in "model-Final_Biometrics.R". 
The "Simulations_run_other_methods" subfolder contains code to run the previous methods (no missing data, complete case, treat all as censored, treat all as failures, multiple imputation (Brownstein Cai, Slade and Bair 2015 in Statistics in Medicine), and the bootstrapping based method of Cook and Kosorok (2004). The "Table1_Biom.R" file runs the simulations for these methods and calls the other two code files ("marfuncs.R" and "MI_binary.R"), which contain code for the functions used.

The "ApplyToOPPERA" folder contains code to apply the methods to the OPPERA study. Together, these files create Table 2. 
The "OPPERA_run_Bayes" subfolder contains code to apply the Bayesian method to the OPPERA study, i.e. the rightmost column in Table 2. This code can be run using "code_biometrics_JAGS_OPPERA.R"
The "OPPERA_run_other_methods" subfolder contains code to apply the multiple imputation method and ad-hoc method to treat all missing data as censored. The files "do_mi_clin_rev.R", "do_mi_psych_rev.R", "do_mi_QST_rev.R" run these methods on the clinical, psych, and QST data, respectively. Necessary functions are included in "MI_binary.R", "read_data.R", and "read_triggers_biom_nb.R".
