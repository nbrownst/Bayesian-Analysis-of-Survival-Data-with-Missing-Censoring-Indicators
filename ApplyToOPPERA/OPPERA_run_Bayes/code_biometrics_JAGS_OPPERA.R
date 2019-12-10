###Reading data

#data.all <- read.csv("OPPERA_clinical_data.csv", na.strings="L")

#data.all <- read.csv("OPPERA_psych_data.csv", na.strings="L")

#data.all <- read.csv("OPPERA_qst_data.csv", na.strings="L")

D1<-rep(0,dim(data.all)[1])

D2<-rep(0,dim(data.all)[1])

D3<-rep(0,dim(data.all)[1])

D4<-rep(0,dim(data.all)[1])

for (i in 1:dim(data.all)[1]){

	if (data.all$Q[i]==0){D1[i]<-0}

	if (data.all$Q[i]==1){D1[i]<-0}
	
	if (data.all$Q[i]==1 & data.all$zeta[i]==1 & data.all$delta[i]==1){D2[i]<-1;D3[i]<-0}

	if (data.all$Q[i]==1 & data.all$zeta[i]==1 & data.all$delta[i]==0){D2[i]<-0;D3[i]<-1}

	if (data.all$Q[i]==1 & data.all$zeta[i]==0 & data.all$delta[i]==9){D4[i]<-1}
}

##Auxiliary variable for Jags

Ind<-NULL

for (i in 1:dim(data.all)[1]){

	if (data.all$Q[i]==0){Ind[i]<-1}
	
	if (data.all$Q[i]==1 & data.all$zeta[i]==1 & data.all$delta[i]==1){Ind[i]<-2}

	if (data.all$Q[i]==1 & data.all$zeta[i]==1 & data.all$delta[i]==0){Ind[i]<-3}

	if (data.all$Q[i]==1 & data.all$zeta[i]==0 & data.all$delta[i]==9){Ind[i]<-4}

}

ind1<-rep(0,dim(data.all)[1])
ind2<-rep(0,dim(data.all)[1])
ind3<-rep(0,dim(data.all)[1])
ind4<-rep(0,dim(data.all)[1])

ind1[Ind==1]<-1
ind2[Ind==2]<-1
ind2[Ind==3]<-1
ind3[Ind==2]<-1
ind3[Ind==3]<-1
ind4[Ind==4]<-1

###W

new.data<-data.frame(Ind,ind1,ind2,ind3,ind4,D1,D2,D3,D4,data.all)

write.csv(new.data,"OPPERA_Biometrics_Approach2_clin_vs2.csv")

#write.csv(new.data,"OPPERA_Biometrics_Approach2_psych_vs2.csv")

#write.csv(new.data,"OPPERA_Biometrics_Approach2_qst.csv")

###RJags code

install.packages("rjags",repos="http://dirichlet.mat.puc.cl")

install.packages("runjags",repos="http://dirichlet.mat.puc.cl")

library(rjags)

library(runjags)

rm(list=ls())

model.biometrics<-'model{

##Priors
	for (k in 1:15){
		lambda[k]~dgamma(1,1)
	 }
	h~dgamma(1,1) #h_j=h for all j.

	for (k in 1:5){
		beta[k]~dnorm(0,1.0E-6)
	}

	eta~dnorm(0,1.0E-6)

##Transformations

	cslambda[1]<-0
	for (j in 1:14) {
		cslambda[j+1] <- cslambda[j] + lambda[j]
	}
	hazard<-exp(eta)

##Likelihood
	phi<- -log(1)+0.0001
	for(i in 1:n){
		##1:Q_i=0:
		D_{1i}=0,D_{1i}~Poi(theta_i+mu_{im_i})
		mu_i[i]<-lambda[m[i]]*(V[i]-0.33*(m[i]-1))##mu_i
		theta_i[i]<-ifelse(m[i]>1,0.33*cslambda[(m[i])],0)##theta_i
		mu11[i]<-mu[i]*(mu_i[i]+theta_i[i])
		theta_iplusGamma_i[i]<-ifelse(m[i]==1,lambda[1]*V[i]*mu[i]+h*V[i],0.33*cslambda[(m[i])]*mu[i]+h*(m[i]-1))##theta_i+gamma_i
		mu1[i]<-ind1[i]*mu11[i]+(1-ind1[i])*theta_iplusGamma_i[i]
		D1[i]~dpois(mu1[i])

		##2:Q_i=1,zeta_i=1,Delta_i=1,D_{2i}=1,D_{3i}=0,D_{2i}~Poi(mu_{im_i}),D_{3i}~Poi(G_{im_i})
		##3:Q_i=1,zeta_i=1,Delta_i=0,D_{2i}=0,D_{3i}=1,D_{2i}~Poi(mu_{im_i}),D_{3i}~Poi(G_{im_i})
		mu21[i]<-mu[i]*mu_i[i]
		mu2[i]<-ind2[i]*mu21[i]+(1-ind2[i])*phi
		D2[i]~dpois(mu2[i])
		mu31[i]<-h*(V[i]-0.33*(m[i]-1))
		mu3[i]<-ind3[i]*mu31[i]+(1-ind3[i])*phi
		D3[i]~dpois(mu3[i])

		##4:Q_i=1,zeta_i=0,Delta_i=9,D_{4i}=1,D_{4i}~Poi(mu_{im_i}+G_{im_i})
		mu41[i]<-mu21[i]+mu31[i]
		mu4[i]<-ind4[i]*mu41[i]+(1-ind4[i])*phi
		D4[i]~dpois(mu4[i])   
	}

	for (i in 1:n){
		mu[i] <-exp(beta[1]*z1[i]+beta[2]*z2[i]+beta[3]*z3[i]+beta[4]*z4[i]+beta[5]*z5[i]+eta*y_1[i]) 
	}

}
'

###Data and initial values for JAGS

data.biom<-read.csv("OPPERA_Biometrics_Approach2_clin_vs2.csv", na.strings="L")

n<-dim(data.biom)[1]

attach(data.biom)

datas <- list ("m"=m,"V"=V,"D1"=D1,"D2"=D2,"D3"=D3,"D4"=D4,"z1"=z1,"z2"=z2,"z3"=z3,"z4"=z4,"z5"=z5,"y_1"=clin_var1,"n"=n,"ind1"=ind1,"ind2"=ind2,"ind3"=ind3,"ind4"=ind4)  

inits <- function(){list(beta=rnorm(5),eta=rnorm(1),lambda=rgamma(15,1,1),h=rgamma(1,1,1))}

param<-c("beta","eta","hazard","lambda","h")

set.seed(1241)

model.out<-run.jags(model=model.biometrics,monitor= param,data= datas,inits= inits,n.chains=2,sample=10000,adapt=2000,burnin=5000,thin=10,modules='runjags',jags.refresh=40)

