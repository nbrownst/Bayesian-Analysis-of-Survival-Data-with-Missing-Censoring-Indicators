###Library

library(rjags)

library(runjags)

library(survival)

library(coda)

#Create tau

CreateBin2 <- function(v, min.bin.width, min.bin.count) {
	N <- length(v)
	max.v <- max(v) 
	v.sort <- sort(v)
	counts <- n.used <- times <- bin.min <- 0
	bin.max <- min.bin.width
	while(bin.max < max.v) {
		bin.count <- sum(bin.min < v & v <= bin.max)
		n.used <- sum(v < bin.max)
		if (bin.count < min.bin.count) {
			n.add <- min.bin.count - bin.count
			bin.max <- v.sort[min(N, n.used + n.add)]
		}
		bin.count <- sum(bin.min < v & v <= bin.max)
		counts <- c(counts, bin.count)
		times <- c(times, bin.max)
		bin.min <- bin.max
		bin.max <- bin.max + min.bin.width
	}
	times[length(times)] <- max.v + 0.001
	if((bin.max -min.bin.width) < max.v) {max.v <- max.v + 0.001;times <- c(times, max.v)}
	return(times)
}
    
############Function to create the data###################    

create.data <- function(beta,mcar,samp.size){

#create data
	z <- rnorm(samp.size, mean=2,sd=1) 
	y.true <- -1*log(runif(samp.size)) / exp(beta*z)
	c.time <- rexp(samp.size, 1/5)
	ic.ideal <- 1*(y.true<c.time)
	y <- y.true
	y[ic.ideal==0] <- c.time[ic.ideal==0]
  
	m<-rep(99,samp.size)

	V <- rep(0, samp.size)
	for(i in 1:samp.size){
		V[i] <- min(y.true[i], c.time[i])
	}
	ic <- ic.ideal 
  
	x <- rnorm(samp.size, ic.ideal, 0.3)
	Q<-rep(0,length(x))
	Q[x>0.5]<-1
  
	ic[Q==0]<-0
  
	if(mcar==1){
		#missing completely at random
		ic[Q==1&(runif(samp.size)<0.4)] <- NA
	}
	else if(mcar==0){
		#missing at random
		#ic[(Q==1)&(runif(samp.size)<exp(-0.25*z)/(1+exp(-0.25*z)))] <- NA
		#ic[(Q==1)&(runif(samp.size)<(exp(-0.5*z+0.2*y)/(1+exp(-0.5*z+0.2*y))))] <- NA
		ic[(Q==1)&(runif(samp.size)<(exp(-0.2-0.3*z+0.1*y)/(1+exp(-0.2-0.3*z+0.1*y))))] <- NA
	}
	else if(mcar==2){
		#missing at random
		#ic[(Q==1)&(runif(samp.size)<exp(-0.25*z)/(1+exp(-0.25*z)))] <- NA
		#ic[(Q==1)&(runif(samp.size)<(exp(-0.5*z+0.2*y)/(1+exp(-0.5*z+0.2*y))))] <- NA
		ic[(Q==1)&(runif(samp.size)<(exp(-0.2-0.3*z+0.1*y)/(1+exp(-0.2-0.3*z+0.1*y))))] <- NA
		ic[Q==0&(runif(samp.size)<0.4)] <- NA
	}
	else if(mcar==-1){
		#missing not at random
		ic[(Q==1)&(((ic.ideal==0)&(runif(samp.size)<0.3))|((ic.ideal==1)&(runif(samp.size)<0.5)))] <- NA
	}
	else if(mcar==-2){
		#missing not at random
		ic[(Q==1)&(((ic.ideal==0)&(runif(samp.size)<0.2))|((ic.ideal==1)&(runif(samp.size)<0.6)))] <- NA
	}
  
	indata <- data.frame(z,x,y,ic,ic.ideal,Q,V, y.true, c.time,m, row.names=NULL)
	return(indata)
}


#Functions

mean.na=function(x){
	return(mean(x,na.rm=TRUE))
}

var.na=function(x){
	return(var(x,na.rm=TRUE))
}

sum.not.na=function(x){
	return(sum(!is.na(x)))
}

sims <- 250
chains <- 1
burnin <- 5000
iter <- 8000
sampsize<-250
truebeta<--3
textfile<-'marbet3.txt'

#############Functions to run the simulations###############################
doFBsims<-function(sims,chains,burnin,iter,sampsize,truebeta,textfile){

  #NB added to store these as vectors
	bias.vec<-rep(NA,length=sims)
	bias.2.vec<-rep(NA,length=sims)
	width.vec<-rep(NA,length=sims)
	coverage.vec<-rep(NA,length=sims)
  
  
	output <- matrix(NA, ncol=4, nrow=sims) 
	startt <- Sys.time()
  
	results<-matrix(NA,nrow=sims,ncol=4) #4 columns b/c of mean, sd, 2.5 and 97.5 percentiles
	for(i in 1:sims){
    
		###########
		#this creates the data for each simulation with one observation per individual 
		#(doesn't create time intervals)
		#############    
    
		init.data <- create.data(beta=truebeta, mcar=0, samp.size=sampsize)
		data<- init.data[order(init.data$V),c("z", "V","ic", "ic.ideal","m", "Q", "y", "c.time","x")]
		colnames(data) <- c("Z", "V", "delta", "del.true","m", "Q", "f.time", "c.time", "X")
    
		#print(-1)
		#print(truebeta)
		complete.data <- summary(coxph(Surv(y,ic.ideal)~z,data=init.data))$coef
		print(complete.data)
		#print(init.data[1:100,])
		#print(dim(init.data))
		#print(names(init.data))
		#print(summary(init.data))
  

		times <- CreateBin2(data$V, 0.33, 15)
		T <- length(times)-1
		eps <- 0.0001
		if(max(times) < max(data$V)){times[T+1] <- max(data$V)}

		N <- dim(data)[1]
		for(j in 1:N){
			if(is.na(data$delta[j])){data$zeta[j] <- 0}
			else(data$zeta[j] <- 1)
			data$e[j] <- data$Q[j]*data$delta[j]
			if(is.na(data$e[j])){data$e[j] <- 0}
			if(is.na(data$delta[j])){data$delta[j]<-9}##9 is missing for delta
		}

		for(j in 1:N){
			aux<-0
			for(k in 2:length(times)){
				if(data$V[j]<times[k]){aux<-aux+1}
			}
			data$m[j]<-length(times)-aux
		}


		D1<-D2<-D3<-D4<-rep(0,N)

		for (j in 1:N){

			if (data$Q[j]==0){D1[j]<-0}

			if (data$Q[j]==1){D1[j]<-0}
	
			if (data$Q[j]==1 & data$zeta[j]==1 & data$delta[j]==1){D2[j]<-1;D3[j]<-0}

			if (data$Q[j]==1 & data$zeta[j]==1 & data$delta[j]==0){D2[j]<-0;D3[j]<-1}

			if (data$Q[j]==1 & data$zeta[j]==0 & data$delta[j]==9){D4[j]<-1}
		}

		##Auxiliary variable for Jags

		Ind<-NULL

		for (j in 1:N){

			if (data$Q[j]==0){Ind[j]<-1}
	
			if (data$Q[j]==1 & data$zeta[j]==1 & data$delta[j]==1){Ind[j]<-2}

			if (data$Q[j]==1 & data$zeta[j]==1 & data$delta[j]==0){Ind[j]<-3}

			if (data$Q[j]==1 & data$zeta[j]==0 & data$delta[j]==9){Ind[j]<-4}

		}

		ind1<-ind2<-ind3<-ind4<-rep(0,N)

		ind1[Ind==1]<-1
		ind2[Ind==2]<-1
		ind2[Ind==3]<-1
		ind3[Ind==2]<-1
		ind3[Ind==3]<-1
		ind4[Ind==4]<-1


		new.data<-data.frame(Ind,ind1,ind2,ind3,ind4,D1,D2,D3,D4,data)
		
		datas <- list ("m"=new.data$m,"V"=new.data$V,"D1"=new.data$D1,"D2"=new.data$D2,"D3"=new.data$D3,
		"D4"=new.data$D4,"Z"=new.data$Z,"N"=N,"ind1"=new.data$ind1,"ind2"=new.data$ind2,"ind3"=new.data$ind3,"ind4"=new.data$ind4,"T"=T+1)  

		inits <- function(){list(beta=rnorm(1),lambda=rgamma(T+1,1,1),h=rgamma(1,1,1))}

		param<-c("beta","lambda","h")

		model.biometrics <- "model-Final_Biometrics.R"


		set.seed(1241)

		lineout<-run.jags(model=model.biometrics,monitor= param,data= datas,
		inits= inits,n.chains=1,sample=8000,adapt=1000,burnin=4000,thin=5,modules='runjags',jags.refresh=40)

		#median.V <- 1.52
		#dH0.beta <- 1
		#dH0.alpha <- 1
		#linedata <- list(data$Z, data$X, data$V, data$zeta, data$e, times, N=250, T, eps, dH0.beta, dH0.alpha, data$Q, truebeta)
		#names(linedata) <- c("Z", "X", "V", "zeta", "e", "times", "N", "T", "eps", "dH0.beta", "dH0.alpha", "Q", "truebeta")
    
		#lineinits <- function(){list(beta= truebeta, dH0=rep(1,T))}
    
		#debugging check
		#print("linedata")
		#print(linedata)
        
		#lineout <- bugs(data=linedata, inits=lineinits, parameters.to.save = c("beta"),
		#model.file=modelfile, n.chains=chains, n.burnin=burnin, n.iter=iter, digits=5)
		print(0)
		print(lineout)
		print(0.5)
		print(lineout$summaries)
		#print(0.75)
		#print(lineout$summary[,"Rhat"])
		#print(all(lineout$summary[,"Rhat"]<1.1))
    
		#out.coda <- read.bugs(lineout)
		#print(0.8)
		#print(out.coda)
		#print(0.9)
		#print(xyplot(out.coda))
    
    
		print(1)
		print(lineout$summaries[1,c(4, 5, 1,3)], 1, 4)
		print(2)
		print(lineout$summaries[1,c(4, 5, 1,3)])
		print(3)
		print(lineout$summaries)
		print(4)
		print(matrix(lineout$summaries[1,c(4, 5, 1,3)], 1, 4))
		print(5)
		print(matrix(lineout$summaries[1,]))
		dim(t(matrix(lineout$summaries[1,])))
		print(5.5)
		#results[i,] <- t(matrix(lineout$summaries[1,])) 
		results[i,] <- matrix(lineout$summaries[1,c(4, 5, 1,3)], 1, 4)
		print(6)
    
		######################################################################################
		#retain results for each Bayesian simulation run
  		######################################################################################
    
		bias.vec[i] <- results[i,1] - truebeta
		bias.2.vec[i] <- (results[i,1] - truebeta)^2 
		#MSE <- sqrt(bias.2) 
		se.beta <- mean(results[,2], na.rm=TRUE)
		width.vec[i] <- results[i,4] - results[i,3] 
		coverage.vec[i] <- ifelse(results[i,3] < truebeta && results[i,4] > truebeta, 1, 0) 
    
		#debugging
		print(7)
		print(dim(output))
		print(bias.vec[i])
		print(bias.2.vec[i])
		#print(se.beta)
		print(width.vec[i])
		print(coverage.vec[i])
		#print(8)
		#output[i,]  <- cbind(bias, MSE, se.beta, width, coverage)
		#print(9)
		#print(output[i,])
		print(i)

		write.table(results,file = "results_bet3.txt", sep = " ",row.names=FALSE, col.names=FALSE)

	}
  
	output  <- cbind(bias.vec, bias.2.vec, width.vec, coverage.vec)
	#print(9)
	#print(output)
  
  
	#overall results based on all simulations
	#bias <- mean(results[,1] - truebeta, na.rm=TRUE)
	#bias.2 <- mean((results[,1] - truebeta)^2, na.rm=TRUE) #used to be sum 6/17/19
	#MSE <- sqrt(bias.2) #not divided by i 6/17/19
	#se.beta <- mean(results[,2], na.rm=TRUE)
	#width <- results[i,3] - results[i,7] #needed to add an i: 06/14/2019
	#coverage <- ifelse(results[,3] < truebeta && results[,7] > truebeta, 1, 0) #< and > were switched: 06/17/2019

  
	avg.bias<-round(mean(bias.vec, na.rm=TRUE),4)#round(apply(bias.vec,2,mean.na),4)
	mse<-round(sqrt(mean(bias.2.vec,na.rm=TRUE)),4) ##THIS IS THE WAY TO COMPUTE MSE
	se.bias<-round((sqrt(var.na(bias.vec))/sum.not.na(bias.vec)),4)
	avg.width<-round(mean(width.vec,na.rm=TRUE),4)
	se.width<-round(sqrt(var.na(width.vec)/sum.not.na(width.vec)),4)
	cover.prob<-round(mean(coverage.vec,na.rm=TRUE),3)
	cover.mcerror<-round(rep(sqrt(0.95*0.05/sims),sum.not.na(cover.prob)),3)

	######################################################################################
	#Puts the results into one table (1 line) and prints that to a file (named with textfile)
	######################################################################################
	table1results<-c(avg.bias,se.bias,se.beta,mse,avg.width,se.width,cover.prob)
	write.table(cbind(avg.bias,se.bias,se.beta,mse,avg.width,se.width,cover.prob),file = textfile, sep = "&")

	print(width.vec)
	print(avg.width)
	print(se.width)
	print(table1results)
  
	endt <- Sys.time()
	rtime1 <- endt-startt
	print(rtime1)
	#sqrt of average of 250 (bias^2) values  
}

######################################################################################
#These run the simulations, but the first few lines are just for testing purposes (small runs: change as needed)
######################################################################################

#testing
doFBsims(sims=2,chains=1,burnin=4000,iter=8000,sampsize=250,truebeta=-0.5,textfile="test_marbet0pt5.txt")
doFBsims(sims=2,chains=1,burnin=5000,iter=8000,sampsize= 250,truebeta=-1.5,textfile="test_marbet1pt5.txt")
doFBsims(sims=2,chains=1,burnin=5000,iter=8000,sampsize= 250,truebeta=-3,textfile="test_marbet3.txt")


######################################################################################
#These set the seed and run the simulations
######################################################################################


#real
set.seed(54321)
doFBsims(sims=250,chains=1,burnin=4000,iter=8000,sampsize=250,truebeta=-0.5,textfile="marbet0pt5.txt")
doFBsims(sims=250,chains=1,burnin=4000,iter=8000,sampsize=250,truebeta=-1.5,textfile="marbet1pt5.txt")
doFBsims(sims=250,chains=1,burnin=4000,iter=8000,sampsize=250,truebeta=-3,textfile="marbet3.txt")

