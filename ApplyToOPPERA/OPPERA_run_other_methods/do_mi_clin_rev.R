n.mi <- 10
source("MI_binary.R")
library(survival)
#previous version using Eric's linux machine
#library(doMC)

#registerDoMC(cores=multicore:::detectCores())

#try this for windows from https://stackoverflow.com/questions/23926334/how-do-i-parallelize-in-r-on-windows-example
# process in parallel
library(doParallel) 
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)
#there will be an off key at the end

#to get recode from mi.binary to work
#library(dplyr)
library(car)


#need to setwd or source won't work
setwd("M:\\lab\\Lab_Brownstein\\Research\\MissingCensoringIndicators\\NB_code\\OPPERA\\oppera analysis NB")


source("read_triggers_biom_nb.R")

clin.all <- read.csv("Clinical2.csv")
clin.all <- merge(trigger.default, clin.all, all=TRUE)

clin.all[is.na(clin.all[,10]),10] <- 0
clin.all[clin.all[,1]==2703,10] <- 1

clin.all[clin.all[,5]==99,5] <- 4
clin.all[clin.all[,10]==1,9] <- 0

clin.x <- clin.all[,11:121]
clin.x <- clin.x[,-(39:43)]

out1 <- matrix(nrow=ncol(clin.x), ncol=4)

for (i in 1:ncol(clin.x)) {
  if (length(unique(clin.x[!is.na(clin.x[,i]),i]))>2) {
    cur.x <- scale(clin.x[,i])
  }
  else {
    cur.x <- clin.x[,i]
  }
  cur.cox <- summary(coxph(Surv(clin.all[,8], clin.all[,9])~cur.x+
                           factor(clin.all[,4])+factor(clin.all[,5])+
                           factor(clin.all[,6])+clin.all[,7]))
  if (cur.cox$coefficients[1,1]<0) {
    cur.x <- -1*cur.x
    cur.cox <- summary(coxph(Surv(clin.all[,8], clin.all[,9])~cur.x+
                             factor(clin.all[,4])+factor(clin.all[,5])+
                             factor(clin.all[,6])+clin.all[,7]))
  }
  out1[i,1] <- cur.cox$conf.int[1,1]
  out1[i,2] <- cur.cox$conf.int[1,3]
  out1[i,3] <- cur.cox$conf.int[1,4]
  out1[i,4] <- cur.cox$coefficients[1,5]
}

out1 <- as.data.frame(out1)
names(out1) <- c("HR", "L95CI", "U95CI", "P")
write.csv(out1, "Clinical_NoMI.csv", row.names=names(clin.x))

out2 <- matrix(nrow=ncol(clin.x), ncol=5)

for (i in 1:ncol(clin.x)) {
  if (length(unique(clin.x[!is.na(clin.x[,i]),i]))>2) {
    cur.x <- scale(clin.x[,i])
  }
  else {
    cur.x <- clin.x[,i]
  }
  cur.y <- clin.all[,8]
  mi.out <- foreach(j=1:n.mi, .combine=cbind, .packages = c("car", "arm", "mi", "survival")) %dopar% {
    ic.star <- clin.all[,9]
    no.rdc <- is.na(clin.all[,2])&!is.na(clin.all[,3])
    ic.star[no.rdc] <- mi.binary(tmd~siteid+b8facesymptomcount+timetoqhu,
              data=trigger.imp)#@random #removed 5/31/2019
    cur.cox <- summary(coxph(Surv(cur.y, ic.star)~cur.x+
                             factor(clin.all[,4])+factor(clin.all[,5])+
                             factor(clin.all[,6])+clin.all[,7]))
    if (!is.na(cur.cox$coefficients[1,1])) {
      if (cur.cox$coefficients[1,1]<0) {
        cur.x <- -1*cur.x
        cur.cox <- summary(coxph(Surv(cur.y, ic.star)~cur.x+
                                 factor(clin.all[,4])+factor(clin.all[,5])+
                                 factor(clin.all[,6])+clin.all[,7]))
      }
      c(cur.cox$coefficients[1,1], cur.cox$coefficients[1,3])
    }
  }
  
  
  
#    if (!is.null(mi.out)) {
  out2[i,1] <- exp(mean(mi.out[1,]))
  cur.ubar <- mean(mi.out[2,]^2)
  cur.b <- var(mi.out[1,])
  cur.se <- sqrt(cur.ubar + (1+1/n.mi)*cur.b)
  cur.df <- (n.mi-1)*(1 + n.mi*cur.ubar/((n.mi+1)*cur.b))^2
  cur.r <- (1+1/n.mi)*cur.b/cur.ubar
  out2[i,2] <- exp(mean(mi.out[1,]) - qt(0.975, cur.df)*cur.se)
  out2[i,3] <- exp(mean(mi.out[1,]) + qt(0.975, cur.df)*cur.se)
  out2[i,4] <- 2*pt(-1*abs(mean(mi.out[1,]))/cur.se, cur.df)
  out2[i,5] <- (cur.r+2/(cur.df+3))/(cur.r+1)
#    }
}

# turn parallel processing off and run sequentially again:
registerDoSEQ()


out2 <- as.data.frame(out2)
names(out2) <- c("HR", "L95CI", "U95CI", "P")
write.csv(out2, "Clinical_MI.csv", row.names=names(clin.x))

rm(i, n.mi, cur.x, cur.cox, mi.out, cur.ubar, cur.b, cur.se, cur.df)
