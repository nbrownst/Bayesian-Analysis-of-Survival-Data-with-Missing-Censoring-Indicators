library(mi)
library(car)
library(car)
library(arm)

source("read_data.R")

triggers.all <- read.csv("alltriggers.csv", na.strings="L")

#qhudate2 <- rep("01/01/42", nrow(triggers.all))
#qhudate2 <- as.Date(qhudate2, "%m/%d/%y")
#for (i in 1:nrow(triggers.all)) {
#  qhudate2[i] <- as.Date(toString(triggers.all$qhudatefilledout[i]),
#                         format="%d%b%Y")
#}
#triggers.all <- cbind(triggers.all, qhudate2)

	qhudate.r=rep(NA,nrow(triggers.all))
	enroldate.r=rep(NA,nrow(triggers.all))
	for (i in seq(nrow(triggers.all))){
	qhudate.r[i]=
	as.Date(toString(triggers.all$qhudatefilledout[i]), format="%d%b%Y")
	enroldate.r[i]=
	as.Date(toString(triggers.all$enrolmentdate[i]), format="%d%b%Y")
	timetoqhu=(qhudate.r-enroldate.r)/365.25
	}
	triggers.all=cbind(triggers.all,qhudate.r,enroldate.r,timetoqhu)




triggers$QHU_BodyAches[is.na(triggers$QHU_BodyAches)] <- 0
triggers.all$QHU_BodyAches[is.na(triggers.all$QHU_BodyAches)] <- 0
triggers$QHU_InjuryEventMouthOpen[is.na(triggers$QHU_InjuryEventMouthOpen)] <- 0
triggers.all$QHU_InjuryEventMouthOpen[is.na(triggers.all$QHU_InjuryEventMouthOpen)] <- 0
triggers$QHU_VisitProvider[is.na(triggers$QHU_VisitProvider)] <- 0
triggers.all$QHU_VisitProvider[is.na(triggers.all$QHU_VisitProvider)] <- 0

#source("qhu_model.R")

triggers.all$siteid <- factor(triggers.all$siteid)

#x.mat <- model.matrix(~siteid+QHU_BodyAches+b8facesymptomcount+
#                      QHU_InjuryEventMouthOpen+QHU_VisitProvider,
#                      data=triggers.all)
#x.mat <- model.matrix(~siteid+b8facesymptomcount+timetoqhu, data=triggers.all)
#logodds.hat <- x.mat %*% fixef(trigger.glmm3)
#prob.hat <- exp(logodds.hat) / (1+exp(logodds.hat))

#triggers.all <- cbind(triggers.all, prob.hat)

#extract respid, tmd, qhudate.r, prob.hat (also previously included 132,133: what are these?)
trigger.imp <- triggers.all[,c(1,11,2,113,134,132)]
trigger.imp <- trigger.imp[order(trigger.imp[,1], trigger.imp[,6],
                                     decreasing=TRUE),]
trigger.imp <- trigger.imp[!duplicated(trigger.imp[,1]),]
trigger.imp <- trigger.imp[,-6]
trigger.imp <- trigger.imp[order(trigger.imp[,1]),]

trigger.default <- trigger.imp
trigger.default$tmd.imp <- trigger.default$tmd
trigger.default$tmd.imp[is.na(trigger.default$tmd)] <-
    mi.binary(tmd~siteid+b8facesymptomcount+timetoqhu,
              data=trigger.imp)#@random #had to take this out, it wasn't needed (check MI_binary.R), still works and is random

trigger.default <- trigger.default[,c(1,2,6)]

rm(i)
