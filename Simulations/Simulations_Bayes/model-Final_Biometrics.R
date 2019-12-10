model{

	##Priors

	for (k in 1:T){

		lambda[k]~dgamma(1,1)
	}	

	h~dgamma(1,1)

	beta~dnorm(0,1.0E-6)


	##Transformations

	cslambda[1]<-0
	
	for (j in 1:(T-1)) {

    		cslambda[j+1] <- cslambda[j] + lambda[j]
	}

	##Likelihood

	phi<- -log(1)+0.0001

	for(i in 1:N){

		##1:Q_i=0:D_{1i}=0,D_{1i}~Poi(theta_i+mu_{im_i})##

		mu_i[i]<-lambda[m[i]]*(V[i]-0.33*(m[i]-1))##mu_i
    
		theta_i[i]<-ifelse(m[i]>1,0.33*cslambda[(m[i])],0)##theta_i

		mu11[i]<-mu[i]*(mu_i[i]+theta_i[i])

		theta_iplusGamma_i[i]<-ifelse(m[i]==1,lambda[1]*V[i]*mu[i]+h*V[i],0.33*cslambda[(m[i])]*mu[i]+h*(m[i]-1))##theta_i+gamma_i

		mu1[i]<-ind1[i]*mu11[i]+(1-ind1[i])*theta_iplusGamma_i[i]

		D1[i]~dpois(mu1[i])

		##2:Q_i=1,zeta_i=1,Delta_i=1,D_{2i}=1,D_{3i}=0,D_{2i}~Poi(mu_{im_i}),D_{3i}~Poi(G_{im_i})##

		##3:Q_i=1,zeta_i=1,Delta_i=0,D_{2i}=0,D_{3i}=1,D_{2i}~Poi(mu_{im_i}),D_{3i}~Poi(G_{im_i})##

		mu21[i]<-mu[i]*mu_i[i]

		mu2[i]<-ind2[i]*mu21[i]+(1-ind2[i])*phi

		D2[i]~dpois(mu2[i])

		mu31[i]<-h*(V[i]-0.33*(m[i]-1))

		mu3[i]<-ind3[i]*mu31[i]+(1-ind3[i])*phi

		D3[i]~dpois(mu3[i])

		##4:Q_i=1,zeta_i=0,Delta_i=9,D_{4i}~Poi(mu_{im_i}+G_{im_i})##

		mu41[i]<-mu21[i]+mu31[i]

		mu4[i]<-ind4[i]*mu41[i]+(1-ind4[i])*phi

		D4[i]~dpois(mu4[i])
	}

	for (i in 1:N){

    		mu[i]<-exp(beta*Z[i])

	}
}

