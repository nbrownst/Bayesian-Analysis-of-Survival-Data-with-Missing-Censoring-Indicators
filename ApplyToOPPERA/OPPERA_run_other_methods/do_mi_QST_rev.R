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
setwd("M:\\lab\\Lab_Brownstein\\Research\\MissingCensoringIndicators\\OPPERA\\oppera analysis NB")

#source("read_triggers.R")
source("read_triggers_biom_nb.R")



qst.all <- read.csv("QST_data.csv")
qst.all <- merge(trigger.default, qst.all, all=TRUE)

qst.all[qst.all[,5]==99,5] <- 4
qst.all[qst.all[,15]==1,14] <- 0

qst.x <- qst.all[,16:48]

out1 <- matrix(nrow=ncol(qst.x), ncol=4)

for (i in 1:ncol(qst.x)) {
    cur.x <- scale(qst.x[,i])
    cur.cox <- summary(coxph(Surv(qst.all[,11], qst.all[,14])~cur.x+
                             factor(qst.all[,4])+factor(qst.all[,5])+
                             factor(qst.all[,6])+qst.all[,7]))
    if (cur.cox$coefficients[1,1]<0) {
        cur.x <- -1*cur.x
        cur.cox <- summary(coxph(Surv(qst.all[,11], qst.all[,14])~cur.x+
                                 factor(qst.all[,4])+factor(qst.all[,5])+
                                 factor(qst.all[,6])+qst.all[,7]))
    }
    out1[i,1] <- cur.cox$conf.int[1,1]
    out1[i,2] <- cur.cox$conf.int[1,3]
    out1[i,3] <- cur.cox$conf.int[1,4]
    out1[i,4] <- cur.cox$coefficients[1,5]
}

out1 <- as.data.frame(out1)
names(out1) <- c("HR", "L95CI", "U95CI", "P")
write.csv(out1, "QST_NoMI.csv", row.names=names(qst.x))

out2 <- matrix(nrow=ncol(qst.x), ncol=5)

for (i in 1:ncol(qst.x)) {
    cur.x <- scale(qst.x[,i])
    cur.y <- qst.all[,11]
    mi.out <- foreach(j=1:n.mi, .combine=cbind, .packages = c("car", "arm", "mi", "survival")) %dopar% {
      ic.star <- qst.all[,14]
      no.rdc <- is.na(qst.all[,2])&!is.na(qst.all[,3])
      ic.star[no.rdc] <- mi.binary(tmd~siteid+b8facesymptomcount+timetoqhu,
              data=trigger.imp)#@random #removed 5/31/2019
      cur.cox <- summary(coxph(Surv(cur.y, ic.star)~cur.x+
                               factor(qst.all[,4])+factor(qst.all[,5])+
                               factor(qst.all[,6])+qst.all[,7]))
      if (cur.cox$coefficients[1,1]<0) {
        cur.x <- -1*cur.x
        cur.cox <- summary(coxph(Surv(cur.y, ic.star)~cur.x+
                                 factor(qst.all[,4])+factor(qst.all[,5])+
                                 factor(qst.all[,6])+qst.all[,7]))
      }
      c(cur.cox$coefficients[1,1], cur.cox$coefficients[1,3])
    }
    out2[i,1] <- exp(mean(mi.out[1,]))
    cur.ubar <- mean(mi.out[2,]^2)
    cur.b <- var(mi.out[1,])
    cur.se <- sqrt(cur.ubar + (1+1/n.mi)*cur.b)
    cur.df <- (n.mi-1)*(1 + n.mi*cur.ubar/((n.mi+1)*cur.b))^2
    ur.r <- (1+1/n.mi)*cur.b/cur.ubar
    out2[i,2] <- exp(mean(mi.out[1,]) - qt(0.975, cur.df)*cur.se)
    out2[i,3] <- exp(mean(mi.out[1,]) + qt(0.975, cur.df)*cur.se)
    out2[i,4] <- 2*pt(-1*abs(mean(mi.out[1,]))/cur.se, cur.df)
    out2[i,5] <- (cur.r+2/(cur.df+3))/(cur.r+1)
}

# turn parallel processing off and run sequentially again:
registerDoSEQ()

out2 <- as.data.frame(out2)
names(out2) <- c("HR", "L95CI", "U95CI", "P")
write.csv(out2, "QST_MI.csv", row.names=names(qst.x))

rm(i, n.mi, cur.x, cur.cox, mi.out, cur.ubar, cur.b, cur.se, cur.df)
