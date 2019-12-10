#library(survival)
#library(boot)
#library(doMC)
#registerDoMC()

setwd("M:\\lab\\Lab_Brownstein\\Research\\MissingCensoringIndicators\\NB_code\\Sims\\try_SIM_code_windows")


source("marfuncs_rev.R")

set.seed(12321)
#options(error=traceback)
#syntax is as follows: dosims(truebeta,mcar,misspec)

#test
#dosims(-0.5,0,0,nsim=2)


#run with 250
dosims(-0.5,0,0,nsim=250,nx=1000)
dosims(-1.5,0,0,nsim=250,nx=1000)
dosims(-3,0,0,nsim=250,nx=1000)

######not needed
##change sample size
#set.seed(12321)
#dosims(-0.5,0,0,nsim=250,nx=2000)
#dosims(-1.5,0,0,nsim=250,nx=2000)
#dosims(-3,0,0,nsim=250,nx=2000)

#set.seed(12321)

#dosims(-0.5,0,0,nsim=250,nx=5000)
#dosims(-1.5,0,0,nsim=250,nx=5000)
#dosims(-3,0,0,nsim=250,nx=5000)
