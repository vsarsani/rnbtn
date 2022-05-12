# TODO: Add comment
# 
# Author: Brian
###############################################################################

# Infer the z value's distribution using EM with finite mixture model
# Input: z value-z[vector], number of non-null components-K[scalar,4]
#        Maximum iteration number-maxIter[scalar,1000]
#        Tolerance-tol[scalar, 1e-8], Show iteration info-info[bool, T]
# Output: The proportion of null components-pi0[scale]
#         The proportion of non-null components-Pi1[vector]
#         The estimated variance of non-null components-sigma2[vector]
#         The probability of being null component-h0[vector]
#		  The probability of being each of the non-null components-h[list]
#         Iteration number-iter[scalar]
#         Expectation of complete log likelihood-Qval[scalar]
snpEM<-function(z, K=2, maxIter=1000, tol=1e-4, beta0=length(z)/5, info=TRUE){
	m<-length(z)
	#initialization
	pi0_0<-0.9
	Pi1_0<-rep((1-pi0_0)/K,K)
	sigma2_0<-rgamma(K,1)
	
	pi0_t<-pi0_0
	Pi1_t<-Pi1_0
	sigma2_t<-sigma2_0
	h<-list()
	h0<-0
	tmpH0<-0
	tmpH<-list()
	nanVal<-function(x) ifelse(is.nan(x),0,x)
	for(iter in 1:maxIter){
		#E step
		for(i in 0:K){
			if(i==0) tmpH0<-pi0_t*dnorm(z)
			else tmpH[[i]]<-Pi1_t[i]*dnorm(z,sd=sqrt(1+sigma2_t[i]))
		}
		norH<-tmpH0+Reduce('+',tmpH)
		h0<-nanVal(tmpH0/norH)
		h<-lapply(tmpH,FUN=function(x) return(nanVal(x/norH)))
		
		#M step
		pi0_t1<-(sum(h0)+beta0)/(m+beta0)
		Pi1_t1<-sapply(h,FUN=sum)/(m+beta0)
		sigma2_t1<-sapply(h, FUN=function(x){
			if(sum(x)==0) return(0)
			else return(max(sum(x*z^2)/sum(x)-1,0))
		} )
		
		if( (abs(nanVal((pi0_t1-pi0_t)/pi0_t))<tol) && sqrt(nanVal(sum((Pi1_t1-Pi1_t)^2)/sum(Pi1_t^2)))<tol && 
				Reduce('+',Map(function(x,y){return(ifelse(sum(y^2)==0, 0, sqrt(sum((x-y)^2)/sum(y^2))))},sigma2_t1,sigma2_t))<tol) break
		else{
			pi0_t<-pi0_t1
			Pi1_t<-Pi1_t1
			sigma2_t<-sigma2_t1
		}
	}
	
#	Qval<-ifelse(tmpH0==0,0,sum(h0*log(tmpH0)))+sum(Reduce('+',Map('*',h,lapply(tmpH,function(x){return(ifelse(x==0,0,log(x)))}))))
	logF<-function(x) ifelse(x==0,0,log(x))
	Qval<-sum(h0*logF(tmpH0))+sum(Reduce('+',Map('*',h,lapply(tmpH,logF))))
	if(info){
		cat('pi0:',pi0_t,'\n')
		cat('Pi1:\n')
		print(Pi1_t)
		cat('sigma^2:\n')
		print(sigma2_t)
		cat('Iteration:',iter,'\nLog-likelihood:',Qval,'\n')
	}
	return(list(pi0=pi0_t,Pi1=Pi1_t,sigma2=sigma2_t,h0=h0,h=h,iter=iter,Qval=Qval))
}

#Bayesian Predictive Power (Probability under H1)
predPower<-function(lower, Pi1, sigma2){	
	prob<-0
	K<-length(Pi1)
	lower<-abs(lower)
	for(i in 1:K)
		prob<-prob+2*Pi1[i]*(1-pnorm(lower,sd=sqrt(1+sigma2[i])))
	prob<-prob/sum(Pi1)
	return(prob)
}

#Calculate Bayes Fdr
calBayesFdr<-function(t, Pi1, sigma2){
	t<-abs(t)
	F0<-2*(1-pnorm(t))
	F1<-predPower(t, Pi1, sigma2)
	pi0<-1-sum(Pi1)
	return(pi0*F0/(pi0*F0+(1-pi0)*F1))
}

#Bisection to find positive x such that fun(x)=y, fun is a monotonic increasing function
incSearch<-function(fun,y,init=1000,tol=0.01){
	x<-init
	lowSign<-FALSE
	highSign<-FALSE
	repeat{
		currY<-fun(x)
		if(currY<y){
			lowX<-x
			x<-2*x
			lowSign<-TRUE
		}else if(currY>y){
			highX<-x
			x<-x/2
			highSign<-TRUE
		}else return(x)
		if(lowSign && highSign) break
	}
	while(highX-lowX>tol){
		x<-(lowX+highX)/2
		currY<-fun(x)
		if(currY<y){
			lowX<-x
		}else if(currY>y){
			highX<-x
		}else return(x)
	}
	return(x)
}

#Bayes Fdr Control
BayesFdr<-function(z,q,K=2, beta0=length(z)/5){
	EMresult<-snpEM(z, K=K, beta0=beta0, info=TRUE)
	thld<-incSearch(function(x) return(1-calBayesFdr(x,EMresult$Pi1,EMresult$sigma2)),1-q,init=2)
	rejectIdx<-(1:length(z))[abs(z)>=thld]
	return(list(thld=thld,idx=rejectIdx))
}
