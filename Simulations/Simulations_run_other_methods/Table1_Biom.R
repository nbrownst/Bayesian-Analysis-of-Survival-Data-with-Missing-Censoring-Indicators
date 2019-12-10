#library(survival)
#library(boot)
#library(doMC)
#registerDoMC()

setwd("M:\\lab\\Lab_Brownstein\\Research\\MissingCensoringIndicators\\NB_code\\Sims\\try_SIM_code_windows")


source("marfuncs_rev.R")

set.seed(12321)
#options(error=traceback)
#syntax is as follows: dosims(truebeta,mcar,misspec)


#run with 250
dosims(-0.5,0,0,nsim=250,nx=1000)
dosims(-1.5,0,0,nsim=250,nx=1000)
dosims(-3,0,0,nsim=250,nx=1000)
