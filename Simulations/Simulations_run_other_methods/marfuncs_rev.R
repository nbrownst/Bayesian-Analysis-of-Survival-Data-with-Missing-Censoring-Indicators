library(survival)
library(boot)
library(arm)
library(car)
library(mi)
#library(foreach) probably not needed for windows


#parallel processing added later so we can call functions

#old: for unix
#library(doMC)
#library(mi)
##registerDoMC(cores=multicore:::detectCores())
#registerDoMC(cores=1)
##library(multicore)

#install.packages("doMC",repos="http://lib.stat.cmu.edu/R/CRAN",lib="~/Rlibs")


setwd("M:\\lab\\Lab_Brownstein\\Research\\MissingCensoringIndicators\\NB_code\\Sims\\try_SIM_code_windows")

#NB added 6/7/19
source("MI_binary.R")


#function to create the data
create.data <- function(beta,mcar,samp.size){
	#create data
	x <- rnorm(samp.size, 2)
	y.true <- -1*log(runif(samp.size)) / exp(beta*x)
	c.time <- rexp(samp.size, 1/5)
	ic.ideal <- 1*(y.true<c.time)
	y <- y.true
	y[ic.ideal==0] <- c.time[ic.ideal==0]

	#create indicators
	ic <- ic.ideal

	x2 <- rnorm(samp.size, ic.ideal, 0.3) #censoring indicator is now related to a covariate
	#create Q variable
	#Q=1 means last QHU is positive, 0 means last QHU is negative (censored)
	Q=rep(0,length(x2))
	Q[x2>0.5]=1


#	x3 <- rnorm(samp.size,ic.ideal/10,1)
	#in practice we don't observe all censoring indicators
	#if Q=0 then no positive QHU so they are censored automatically
	ic[Q==0]=0

        if(mcar==1){
		#missing completely at random
		ic[Q==1&(runif(samp.size)<0.4)] <- NA
	}
	else if(mcar==0){
		#missing at random
		#ic[(Q==1)&(runif(samp.size)<exp(-0.25*x)/(1+exp(-0.25*x)))] <- NA
		#ic[(Q==1)&(runif(samp.size)<(exp(-0.5*x+0.2*y)/(1+exp(-0.5*x+0.2*y))))] <- NA
		ic[(Q==1)&(runif(samp.size)<(exp(-0.2-0.3*x+0.1*y)/(1+exp(-0.2-0.3*x+0.1*y))))] <- NA
	}
        else if(mcar==2){
		#missing at random
		#ic[(Q==1)&(runif(samp.size)<exp(-0.25*x)/(1+exp(-0.25*x)))] <- NA
		#ic[(Q==1)&(runif(samp.size)<(exp(-0.5*x+0.2*y)/(1+exp(-0.5*x+0.2*y))))] <- NA
		ic[(Q==1)&(runif(samp.size)<(exp(-0.2-0.3*x+0.1*y)/(1+exp(-0.2-0.3*x+0.1*y))))] <- NA
                ic[Q==0&(runif(samp.size)<0.4)] <- NA
	}
        else if(mcar==-1){
		#missing not at random
		ic[(Q==1)&
		  (((ic.ideal==0)&(runif(1000)<0.3))
		  |((ic.ideal==1)&(runif(1000)<0.5)))] <- NA
	}
        else if(mcar==-2){
		#missing not at random
		ic[(Q==1)&
		  (((ic.ideal==0)&(runif(1000)<0.2))
		  |((ic.ideal==1)&(runif(1000)<0.6)))] <- NA
	}

#	indata <- data.frame(x,x2,y,ic,ic.ideal,Q,x3)
        indata <- data.frame(x,x2,y,ic,ic.ideal,Q)
	#beta=coxmodels(indata)
	return(indata)
}

#Once the data is created, here is the function definition for the EM algorithm
#to get predicted probabilities
#however we no longer iterate, just do it once
#it takes in the covariates, the survival time, the censoring indicator for the raw data (some missing)
#and the censoring indicator for the complete data
#also notes if model is misspecified
#(1=extra covariate, -1=left one out, 0=neither/the model is correct)
#returns predicted probabilities

EM.algorithm <- function(indata,misspecified){
        indata2 <- indata[indata$Q==1|is.na(indata$ic),]
#	ones=rep(1,length(x))

	#fit a logistic regression model predicting ic based on x.qhu
	#note: doesn't include people who didn't come in, b/c ic will be missing
	#for them, so it's only Q=1 people
	if(misspecified==0){
#		logreg <- summary(glm(ic ~ x+x2+x3+y, family=binomial, subset=Q==1))
#                logreg <- glm(ic ~ x+x2+x3+y, family=binomial, subset=Q==1)
                cur.mi <- mi.binary(ic~x+x2+y, data=indata2)#@random #removed 6/7/2019 for windows
		#this uses the correct model
#		designmat=cbind(ones,x,x2,x3,y)
	}
	else if(misspecified==1){
		indata2$x3 <- rnorm(nrow(indata2))
#		logreg <- summary(glm(ic ~ x+x2+y+x3+x4, family=binomial,subset=Q==1))
#                logreg <- glm(ic ~ x+x2+y+x3+x4, family=binomial,subset=Q==1)
                cur.mi <- mi.binary(ic~x+x2+x3+y, data=indata2)#@random #removed 6/7/2019 for windows
		#this puts in an extra covariate (x3)
#		designmat=cbind(ones,x,x2,y,x3,x4)
	}
	else if(misspecified==-1){
#	logreg <- summary(glm(ic ~ x+x2+y, family=binomial, subset=Q==1)) #leaves out x2
#        logreg <- glm(ic ~ x+x2+y, family=binomial, subset=Q==1) #leaves out x2
        cur.mi <- mi.binary(ic~x+y, data=indata2)#@random #removed 6/7/2019 for windows
		#this leaves out a covariate (x3)
#		designmat=cbind(ones,x,x2,y)
	}
#	coef <- logreg$coefficients
#	parms <- coef[,1]

	#use the model fit to get predicted probabilities
#	xbeta=designmat%*%parms
#	pred.prob=exp(xbeta)/(1+exp(xbeta))
#        cur.pred <- predict(logreg, newdata=indata, se.fit=TRUE)
#	return(list(logodds=cur.pred$fit, se.lo=cur.pred$se.fit))
        return(cur.mi)
}

EM.algorithm2 <- function(indata,misspecified){
        x <- indata$x
        x2 <- indata$x2
#        x3 <- indata$x3
        y <- indata$y
        ic <- indata$ic
        Q <- indata$Q
	ones=rep(1,length(x))

	#fit a logistic regression model predicting ic based on x.qhu
	#note: doesn't include people who didn't come in, b/c ic will be missing
	#for them, so it's only Q=1 people
	if(misspecified==0){
#		logreg <- summary(glm(ic ~ x+x2+x3+y, family=binomial, subset=Q==1))
                logreg <- glm(ic ~ x+x2+y, family=binomial, subset=Q==1)
		#this uses the correct model
#		designmat=cbind(ones,x,x2,x3,y)
	}
	else if(misspecified==1){
		x3 <- rnorm(length(x))
#		logreg <- summary(glm(ic ~ x+x2+y+x3+x4, family=binomial,subset=Q==1))
                logreg <- glm(ic ~ x+x2+y+x3, family=binomial,subset=Q==1)
		#this puts in an extra covariate (x3)
#		designmat=cbind(ones,x,x2,y,x3,x4)
	}
	else if(misspecified==-1){
#	logreg <- summary(glm(ic ~ x+x2+y, family=binomial, subset=Q==1)) #leaves out x2
        logreg <- glm(ic ~ x+y, family=binomial, subset=Q==1) #leaves out x2
		#this leaves out a covariate (x3)
#		designmat=cbind(ones,x,x2,y)
	}
#	coef <- logreg$coefficients
#	parms <- coef[,1]

	#use the model fit to get predicted probabilities
#	xbeta=designmat%*%parms
#	pred.prob=exp(xbeta)/(1+exp(xbeta))
        cur.pred <- predict(logreg, newdata=indata)
	return(exp(cur.pred)/(1+exp(cur.pred)))
}


#takes in the data, uses it to find the predicted probabilities
#the return is the vector of beta parameter estimates from the cox models
#for each data analysis method
coxmodels <- function(indata,misspecified0) {

	#get predicted probabilities
	pred.probs = EM.algorithm(indata,misspecified0)
	ic.star=indata$ic
	ic.star[is.na(indata$ic)] =
	rbinom(sum(is.na(indata$ic)),1,pred.probs[is.na(indata$ic)])

	#now use regular cox model on the completed data
	completed=summary(coxph(Surv(indata$y,ic.star) ~ indata$x))

        return(completed$coef[1])
}

# fitting the Cox model using the Cook and Kosorok approach (I think)
coxmodels2 <- function(indata,misspecified0) {

	#get predicted probabilities
	pred.probs = EM.algorithm2(indata,misspecified0)
#        logodds <- EM.algorithm2(indata,misspecified0)$logodds
#        pred.probs <- exp(logodds)/(1+exp(logodds))
#        pred.probs[logodds>700] <- 1

        newx <- c(indata$x, indata$x[is.na(indata$ic)])
        newy <- c(indata$y, indata$y[is.na(indata$ic)])
        newic <- c(indata$ic, rep(0, sum(is.na(indata$ic))))
        newic[is.na(newic)] <- 1
        pred.probs[!is.na(indata$ic)] <- 1
        wts <- c(pred.probs, 1-pred.probs[is.na(indata$ic)])
		#print(sum(wts==0))
		#print(indata[wts==0,])

	#drop things that have wts=0
		newxd=newx[wts!=0]
		newyd=newy[wts!=0]
		newicd=newic[wts!=0]
		wtsd=wts[wts!=0]

	#now use regular cox model on the weighted data
	completed=summary(coxph(Surv(newyd, newicd)~newxd, weights=wtsd))

        return(completed$coef[1])
}

#function that considers the certain set of indices so we can bootstrap
coxmodels.boot <- function(indata,i,mis) {
  indata.star <- indata[i,]
  return(coxmodels(indata.star,mis))
}

coxmodels.boot2 <- function(indata, i,mis) {
  indata.star <- indata[i,]
  return(coxmodels2(indata.star,mis))
}

#one run of bootstrapping with 1000 replicates
one.boot <- function(indata,truebeta,missp)
{
	cox.boot <- boot(indata,coxmodels.boot,R=1000,mis=missp)
#                         parallel="multicore",
#                         ncpus=multicore:::detectCores())
	completed.ci <-boot.ci(cox.boot)
	beta=completed.ci$t0
	#conf.int=completed.ci$bca[4:5]
	#cover=as.numeric((conf.int[2]>truebeta)&(conf.int[1]<truebeta))
	#cover.bca=as.numeric((completed.ci$bca[5]>truebeta)&(completed.ci$bca[4]<truebeta))
	#width.bca=completed.ci$bca[5]-completed.ci$bca[4]

	cover.pct=as.numeric((completed.ci$percent[5]>truebeta)&(completed.ci$percent[4]<truebeta))
	width.pct=completed.ci$percent[5]-completed.ci$percent[4]

	#get rid of bca stuff #cover.bca,width.bca,
	bootinfo=cbind(beta,completed.ci$percent[4],completed.ci$percent[5],
			   cover.pct,width.pct)
	#previously it was this:
	#bootinfo=cbind(beta,completed.ci$bca[4],completed.ci$bca[5],
	#cover.bca,width.bca,completed.ci$percent[4],completed.ci$percent[5],
	#cover.pct,width.pct)

	return(bootinfo)
}

one.boot2 <- function(indata,truebeta,missp)
{
	cox.boot <- boot(indata,coxmodels.boot2,R=1000,mis=missp)
#                         parallel="multicore",
#                         ncpus=multicore:::detectCores())
	#print(142) #don't need: removed 6/7/2019
	#save.image()
	completed.ci <-boot.ci(cox.boot,type="perc")

	beta=completed.ci$t0
	#conf.int=completed.ci$bca[4:5]
	#cover=as.numeric((conf.int[2]>truebeta)&(conf.int[1]<truebeta))
	#get rid of bca stuff #
	#cover.bca=as.numeric((completed.ci$bca[5]>truebeta)&(completed.ci$bca[4]<truebeta))
	#width.bca=completed.ci$bca[5]-completed.ci$bca[4]

	cover.pct=as.numeric((completed.ci$percent[5]>truebeta)&(completed.ci$percent[4]<truebeta))
	width.pct=completed.ci$percent[5]-completed.ci$percent[4]


	bootinfo=cbind(beta,completed.ci$percent[4],completed.ci$percent[5],
			   cover.pct,width.pct)
	#previously it was this!
	#bootinfo=cbind(beta,completed.ci$bca[4],completed.ci$bca[5],
	#cover.bca,width.bca,completed.ci$percent[4],completed.ci$percent[5],
	#cover.pct,width.pct)
	#print(bootinfo)
	return(bootinfo)
}

#try this for windows from https://stackoverflow.com/questions/23926334/how-do-i-parallelize-in-r-on-windows-example
# process in parallel
library(doParallel) 
cl <- makeCluster(detectCores(), type='PSOCK')
clusterExport(cl,c("EM.algorithm")) #need this to pass functions
registerDoParallel(cl)


#do the whole procedure many times
#external function to do the sims
#takes in an indicator of whether the data is mcar or not
#and an indicator of model misspecification
#nsim=number of simulation runs
#nx=sample size
dosims=function(truebeta,mcar,misspec,nsim=1000,nx=1000){

	#create matrices to store the data
	sims = seq(nsim)
	betas=matrix(NA,nrow=length(sims),ncol=7)
	lcl=matrix(NA,nrow=length(sims),ncol=7)
	ucl=matrix(NA,nrow=length(sims),ncol=7)
	cover=matrix(NA,nrow=length(sims),ncol=7)
	width=matrix(NA,nrow=length(sims),ncol=7)
        rtime=matrix(NA,nrow=length(sims),ncol=7)
        miri <- rep(NA, length(sims))

        
  #mse: 6/10/2019
  mse=matrix(NA,nrow=length(sims),ncol=7)
  
	#do each scenario length(sims) number of times (usually 1000)
	for (i in sims){
	        print(i)
		#create data
		input.data <- create.data(truebeta,mcar,nx)
	        x <- input.data$x
	        x2 <- input.data$x2
	        y <- input.data$y
	        ic <- input.data$ic
	        ic.ideal <- input.data$ic.ideal
		  Q <- input.data$Q

		#check to make sure there are enough controls and enough events
		#in the full data
		print(table(ic,ic.ideal)) #maybe remove this later
		flag1 = (sum(ic.ideal==1)==0)+(sum(ic.ideal==0)==0)
		flag2 = (sum(ic==1&!is.na(ic))==0)+(sum(ic==0&!is.na(ic))==0)
		print(flag1) #maybe remove this later
		print(flag2) #maybe remove this later
		#print(sum(is.na(ic.ideal)))
		#print(length(ic))
		#print(sum(ic[!is.na(ic)]))
		#print(sum(is.na(ic)))
		#print(-9)

		#if enough cases/controls then proceed
		if(flag1==0){
			#print(5) #not needed, remove 6/7/2019
			#need compare to complete data to get the true value of beta
                        startt <- Sys.time()
                        complete.data <- summary(coxph(Surv(y,ic.ideal)~x))$coef
			beta.complete=complete.data[1]
			betas[i,1]=beta.complete
			lcl[i,1]=complete.data[1]-qnorm(.975)*complete.data[3]
			ucl[i,1]=complete.data[1]+qnorm(.975)*complete.data[3]
			cover[i,1]=as.numeric((ucl[i,1]>truebeta)&(lcl[i,1]<truebeta))
			width[i,1]=ucl[i,1]-lcl[i,1]
                        endt <- Sys.time()
                        rtime[i,1]=endt-startt
      
                        
			#check to make sure still enough cases/controls even with MCIs
			#if so proceed
			if(flag2==0){
			#compare to complete case
                        startt <- Sys.time()
			complete.case <- summary(coxph(Surv(y,ic)~x))$coef
			betas[i,2]=complete.case[1]
			lcl[i,2]=complete.case[1]-qnorm(.975)*complete.case[3]
			ucl[i,2]=complete.case[1]+qnorm(.975)*complete.case[3]
			cover[i,2]=as.numeric((ucl[i,2]>truebeta)&(lcl[i,2]<truebeta))
			width[i,2]=ucl[i,2]-lcl[i,2]
                        endt <- Sys.time()
                        rtime[i,2]=endt-startt

			#compare to treating missing as all censored or as all cases
                        startt <- Sys.time()
                        ic.bad <- ic
			ic.bad[is.na(ic)] <- 0 #for comparison treat missing as censored
			bad1=summary(coxph(Surv(y,ic.bad) ~ x))$coef
			betas[i,3]=bad1[1]
			lcl[i,3]=bad1[1]-qnorm(0.975)*bad1[3]
			ucl[i,3]=bad1[1]+qnorm(0.975)*bad1[3]
			cover[i,3]=as.numeric((ucl[i,3]>truebeta)&(lcl[i,3]<truebeta))
			width[i,3]=ucl[i,3]-lcl[i,3]
                        endt <- Sys.time()
                        rtime[i,3]=endt-startt

                        startt <- Sys.time()
                        ic.bad2 <- ic
			ic.bad2[is.na(ic)] <- 1 #for comparison treat missing as case
			bad2=summary(coxph(Surv(y,ic.bad2) ~ x))$coefficients

			betas[i,4]=bad2[1]
			lcl[i,4]=bad2[1]-qnorm(0.975)*bad2[3]
			ucl[i,4]=bad2[1]+qnorm(0.975)*bad2[3]
			cover[i,4]=as.numeric((ucl[i,4]>truebeta)&(lcl[i,4]<truebeta))
			width[i,4]=ucl[i,4]-lcl[i,4]
                        endt <- Sys.time()
                        rtime[i,4]=endt-startt

      		#do the Cook & Kosorok bootstrapping
                        startt <- Sys.time()
			#print(505) #not needed, remove 6/7/2019
			boot.out=one.boot2(input.data,truebeta,mis=misspec)
			#print(506) #not needed, remove 6/7/2019
			#betas[i,5]=boot.out[1]
			#lcl[i,5]=boot.out[2]
			#ucl[i,5]=boot.out[3]
			#cover[i,5]=boot.out[4]
			#width[i,5]=boot.out[5]

			#add percentile intervals
			lcl[i,6]=boot.out[2] #6
			ucl[i,6]=boot.out[3] #7
			cover[i,6]=boot.out[4] #8
			width[i,6]=boot.out[5] #previously 9 check this!
			betas[i,6]=boot.out[1] #betas[i,5]
                        endt <- Sys.time()
                        rtime[i,6]=endt-startt

			##print(42)
			#compare to multiple imputation
                        startt <- Sys.time()
			n.mi=10
			mi.out=matrix(NA,nrow=n.mi,ncol=2)
			#now do multiple imputation
#	      	  pred.probs = EM.algorithm(input.data,misspec)
#                        em.out = EM.algorithm(input.data,misspec)
			mi.out <-foreach(j=1:n.mi, .combine=rbind, .packages = c("car", "arm", "mi", "survival"),
			                 .export=c("EM.algorithm")) %dopar% { 
			  #old: foreach(j=seq(n.mi),.combine=rbind) %dopar% {
#                                cur.lo <- rnorm(length(em.out$logodds),
#                                                em.out$logodds,
#                                                em.out$se.lo)
#                                pred.probs <- exp(cur.lo)/(1+exp(cur.lo))
#                                pred.probs[cur.lo>700] <- 1
			   #need to define mi.binary here
         source("MI_binary.R", local = TRUE)
			                   
				ic.star=input.data$ic
				ic.star[is.na(input.data$ic)] =
                                    EM.algorithm(input.data,misspec)
	               cur.cox <- summary(coxph(Surv(input.data$y,ic.star) ~ input.data$x))
      	          c(cur.cox$coefficients[1,1], cur.cox$coefficients[1,3])
			} #closes loop for the imputations
		  } #closes loop for flag2
	     } #closes loop containing flag1

		#get necessary info
	      mi.avg <- mean(mi.out[,1])
	      mi.avgvar <- mean(mi.out[,2]^2)
	      mi.var <- var(mi.out[,1])
	      mi.se <- sqrt(mi.avgvar + (1+1/n.mi)*mi.var)
	      mi.df <- (n.mi-1)*(1 + n.mi*mi.avgvar/((n.mi+1)*mi.var))^2
              cur.r <- (1+1/n.mi)*mi.var/mi.avgvar
              miri[i] <- (cur.r+2/(mi.df+3))/(cur.r+1)
	      betas[i,7] = mi.avg
	      lcl[i,7] <- mi.avg - qt(0.975, mi.df)*mi.se
	      ucl[i,7] <- mi.avg + qt(0.975, mi.df)*mi.se
	      #pval <- 2*pt(-1*abs(mi.avg)/mi.se, mi.df)
		width[i,7] = ucl[i,7]-lcl[i,7]
		cover[i,7]=as.numeric((ucl[i,7]>truebeta)&(lcl[i,7]<truebeta))
                endt <- Sys.time()
                rtime[i,7]=endt-startt
                
                
#    ###mse=
#    mse[i]=(betahat-trubeta)^2 #or bias^2
                            
	} #closes loop for the nsim simulations

	bias=betas-truebeta
	bias2=bias^2

	#define mean and variance function that gets rid of NA values
	mean.na=function(x){
		return(mean(x,na.rm=TRUE))
	}
	var.na=function(x){
		return(var(x,na.rm=TRUE))
	}
	sum.not.na=function(x){
		return(sum(!is.na(x)))
	}
	
	#from SIM
	avg.bias=round(apply(bias,2,mean.na),4)
	se.bias=round(sqrt(apply(bias,2,var.na)/sum.not.na(bias)),4)
	avg.width=round(apply(width,2,mean.na),4)
	se.width=round(sqrt(apply(width,2,var.na)/sum.not.na(width)),4)
	cover.prob=round(apply(cover,2,mean.na),3)
	cover.mcerror=round(rep(sqrt(0.95*0.05/nsim),sum.not.na(cover.prob)),3)
        mir <- round(rep(mean(miri), sum.not.na(cover.prob)),3)
        avg.rtime=round(apply(rtime,2,mean.na),3)

        
  #mse overall
  bias.2.vec=apply(bias2,2,mean.na)
  mse=round(sqrt(bias.2.vec),4)
    #round(sqrt(mean(bias.2.vec,na.rm=TRUE)),4) #mse from Mauricio
              
	print(bias)
	print(bias2)
	print(bias.2.vec)
	print(mse)
	print(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime))



	#test if bias is different from zero if there are enough observations
	#ttest=function(x){return(t.test(x,mu=0)$p.value)}
	#testbias=rep(NA,7)
	#if(sum.not.na(bias)>30){
	#	testbias=apply(bias,2,ttest)
	#}

	#test if coverage is 0.95
	#testcover=rep(NA,7)
	#testprop=function(s){prop.test(x=sum(s),n=length(sims),p=0.95)$p.value}
	#if(sum.not.na(testprop)>30){
	#	testcover=apply(cover,2,testprop)
	#}

	#test if lengths are different
	#tests=matrix(NA,nrow=5,ncol=2)
	#boot vs. ideal
	#if(sum.not.na(width)>30){
	#	tests[1,1]=t.test(width[,6],width[,1])$p.value
	#}
	#boot vs. complete case
	#if(sum.not.na(width)>30){
	#	tests[2,1]=t.test(width[,6],width[,2])$p.value
	#}
	#boot vs. treat all as censored
	#if(sum.not.na(width)>30){
	#	tests[3,1]=t.test(width[,6],width[,3])$p.value
	#}
	#boot vs. treat all as failures
	#if(sum.not.na(width)>30){
	#	tests[4,1]=t.test(width[,6],width[,4])$p.value
	#}
	#mi vs. ideal
	#if(sum.not.na(width)>30){
	#	tests[1,2]=t.test(width[,7],width[,1])$p.value
	#}
	#mi vs. complete case
	#if(sum.not.na(width)>30){
	#	tests[2,2]=t.test(width[,7],width[,2])$p.value
	#}
	##mi vs. treat all as censored
	#if(sum.not.na(width)>30){
	#	tests[3,2]=t.test(width[,7],width[,3])$p.value
	#}
	##mi vs. treat all as failures
	#if(sum.not.na(width)>30){			tests[4,2]=t.test(width[,7],width[,4])$p.value
	#}
	##mi vs. boot
	#if(sum.not.na(width)>30){
	#	tests[5,2]=t.test(width[,6],width[,7])$p.value
	#}
	#print(tests)
	#write the tables to the appropriate filenames
	#mcar: probably not needed for Biom
	if(mcar==1 & misspec==0 & truebeta==-1.5){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mcarbet1pt5.txt", sep = "&")
	}
	else if(mcar==1 & misspec==0 & truebeta==-0.5){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mcarbetpt5.txt", sep = "&")
	}
	else if(mcar==1 & misspec==0 & truebeta==-3){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mcarbet3.txt", sep = "&")
	}
	else if(mcar==1 & misspec==1 & truebeta==-3){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mcarbet3toobig.txt", sep = "&")
	}
	else if(mcar==1 & misspec==1 & truebeta==-0.5){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mcarbetpt5toobig.txt", sep = "&")
	}
	else if(mcar==1 & misspec==1 & truebeta==-1.5){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mcarbet1pt5toobig.txt", sep = "&")
	}
	else if(mcar==1 & misspec==-1 & truebeta==-3){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mcarbet3toosmall.txt", sep = "&")
	}
	else if(mcar==1 & misspec==-1 & truebeta==-0.5){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mcarbetpt5toosmall.txt", sep = "&")
	}
	else if(mcar==1 & misspec==-1 & truebeta==-1.5){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mcarbet1pt5toosmall.txt", sep = "&")
	}

        #mcar including Q=0: not for BIom
        else if(mcar==2 & misspec==0 & truebeta==-1.5){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mcar2bet1pt5.txt", sep = "&")
	}
	else if(mcar==2 & misspec==0 & truebeta==-0.5){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mcar2betpt5.txt", sep = "&")
	}
	else if(mcar==2 & misspec==0 & truebeta==-3){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mcar2bet3.txt", sep = "&")
	}
	else if(mcar==2 & misspec==1 & truebeta==-3){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mcar2bet3toobig.txt", sep = "&")
	}
	else if(mcar==2 & misspec==1 & truebeta==-0.5){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mcar2betpt5toobig.txt", sep = "&")
	}
	else if(mcar==2 & misspec==1 & truebeta==-1.5){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mcar2bet1pt5toobig.txt", sep = "&")
	}
	else if(mcar==2 & misspec==-1 & truebeta==-3){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mcar2bet3toosmall.txt", sep = "&")
	}
	else if(mcar==2 & misspec==-1 & truebeta==-0.5){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mcar2betpt5toosmall.txt", sep = "&")
	}
	else if(mcar==2 & misspec==-1 & truebeta==-1.5){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mcar2bet1pt5toosmall.txt", sep = "&")
	}

        #mnar scenario 1: probably not for Biom
	else if(mcar==-1 & misspec==0 & truebeta==-1.5){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mnarbet1pt5.txt", sep = "&")
	}
	else if(mcar==-1 & misspec==0 & truebeta==-0.5){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mnarbetpt5.txt", sep = "&")
	}
	else if(mcar==-1 & misspec==0 & truebeta==-3){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mnarbet3.txt", sep = "&")
	}
	else if(mcar==-1 & misspec==1 & truebeta==-3){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mnarbet3toobig.txt", sep = "&")
	}
	else if(mcar==-1 & misspec==1 & truebeta==-0.5){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mnarbetpt5toobig.txt", sep = "&")
	}
	else if(mcar==-1 & misspec==1 & truebeta==-1.5){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mnarbet1pt5toobig.txt", sep = "&")
	}
	else if(mcar==-1 & misspec==-1 & truebeta==-3){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mnarbet3toosmall.txt", sep = "&")
	}
	else if(mcar==-1 & misspec==-1 & truebeta==-0.5){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mnarbetpt5toosmall.txt", sep = "&")
	}
	else if(mcar==-1 & misspec==-1 & truebeta==-1.5){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mnarbet1pt5toosmall.txt", sep = "&")
	}

        #mnar scenario 2: probably not for biom
	else if(mcar==-2 & misspec==0 & truebeta==-1.5){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mnar2bet1pt5.txt", sep = "&")
	}
	else if(mcar==-2 & misspec==0 & truebeta==-0.5){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mnar2betpt5.txt", sep = "&")
	}
	else if(mcar==-2 & misspec==0 & truebeta==-3){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mnar2bet3.txt", sep = "&")
	}
	else if(mcar==-2 & misspec==1 & truebeta==-3){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mnar2bet3toobig.txt", sep = "&")
	}
	else if(mcar==-2 & misspec==1 & truebeta==-0.5){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mnar2betpt5toobig.txt", sep = "&")
	}
	else if(mcar==-2 & misspec==1 & truebeta==-1.5){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mnar2bet1pt5toobig.txt", sep = "&")
	}
	else if(mcar==-2 & misspec==-1 & truebeta==-3){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mnar2bet3toosmall.txt", sep = "&")
	}
	else if(mcar==-2 & misspec==-1 & truebeta==-0.5){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mnar2betpt5toosmall.txt", sep = "&")
	}
	else if(mcar==-2 & misspec==-1 & truebeta==-1.5){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "mnar2bet1pt5toosmall.txt", sep = "&")
	}

#####################this is what we are looking at: beginning#################
	#mar not mcar: for Biom
	else if(mcar==0 & misspec==0 & truebeta==-1.5){
#		write.table(cbind(avg.bias,se.bias,
#		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
	  write.table(cbind(avg.bias,se.bias,mse,avg.width,se.width,cover.prob,cover.mcerror),
	  file = "marnewbet1pt5.txt", sep = "&")
	}
	else if(mcar==0 & misspec==0 & truebeta==-0.5){
#		write.table(cbind(avg.bias,se.bias,
#		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
	  write.table(cbind(avg.bias,se.bias,mse,avg.width,se.width,cover.prob,cover.mcerror),
	  file = "marnewbetpt5.txt", sep = "&")
	}
	else if(mcar==0 & misspec==0 & truebeta==-3){
#		write.table(cbind(avg.bias,se.bias,
#		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
	  write.table(cbind(avg.bias,se.bias,mse,avg.width,se.width,cover.prob,cover.mcerror),
	  file = "marnewbet3.txt", sep = "&")
	}
	#####################this is what we are looking at: end#################

		#misspecified: probably don't need these for Biom
	else if(mcar==0 & misspec==1 & truebeta==-3){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "marnewbet3toobig.txt", sep = "&")
	}
	else if(mcar==0 & misspec==1 & truebeta==-0.5){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "marnewbetpt5toobig.txt", sep = "&")
	}
	else if(mcar==0 & misspec==1 & truebeta==-1.5){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "marnewbet1pt5toobig.txt", sep = "&")
	}
	else if(mcar==0 & misspec==-1 & truebeta==-3){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "marnewbet3toosmall.txt", sep = "&")
	}
	else if(mcar==0 & misspec==-1 & truebeta==-0.5){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "marnewbetpt5toosmall.txt", sep = "&")
	}
	else if(mcar==0 & misspec==-1 & truebeta==-1.5){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "marnewbet1pt5toosmall.txt", sep = "&")
	}
	#add new possibilities: probably not needed for Biom
	else if(mcar==0 & misspec==0 & truebeta==-10){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "marnewbet10.txt", sep = "&")
	}
	else if(mcar==0 & misspec==0 & truebeta==-5){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "marnewbet5.txt", sep = "&")
	}
	else if(mcar==0 & misspec==0 & truebeta==-4){
		write.table(cbind(avg.bias,se.bias,mse,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),
		file = "marnewbet4.txt", sep = "&")
	}
	else if(mcar==0 & misspec==0 & truebeta==-3.5){
		write.table(cbind(avg.bias,se.bias,mse,#testbias,
		avg.width,se.width,cover.prob,cover.mcerror,mir,avg.rtime),#,testcover),
		file = "marnewbet3pt5.txt", sep = "&")
	}
	#print warnings
	warnings()
	#save the results
	#return(cbind(avg.bias,se.bias,testbias,
	#	avg.width,se.width,cover.prob,testcover))
	#return(cbind(avg.bias,se.bias,
	#            avg.width,se.width,cover.prob))
	
	return(cbind(avg.bias,se.bias,mse,#mse.bias,
		avg.width,cover.prob))
	save.image()
	
	#turn off multicore for windows 6/7/2019
	stopCluster(cl)
	registerDoSEQ()
}
