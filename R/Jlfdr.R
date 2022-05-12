# TODO: Joint Local False Discovery Rate (Jlfdr) control in joint analysis of two GWAS 
# 
# Author: wjiangaa
###############################################################################

#install.packages('mvtnorm',dependencies=T)
#install.packages("lib/mvtnorm_1.0-2.zip") #For windows
library(mvtnorm)

# Infer the z value pairs' distribution using EM with finite mixture model
# Input: z value pairs-z1, z2[vector], number of non-null components-K[scalar,4]
#        Maximum iteration number-maxIter[scalar,1000]
#        Tolerance-tol[scalar, 1e-8], Show iteration info-info[bool, T]
# Output: The proportion of null components-pi0[scale]
#         The proportion of non-null components-Pi1[vector]
#         The estimated covariance matrix of non-null components-Sigma[list]
#         The probability of being null component-h0[vector]
#		  The probability of being each of the non-null components-h[list]
#         Iteration number-iter[scalar]
#         Expectation of complete log likelihood-Qval[scalar]
snpEM2<-function(z1,z2, K=2, maxIter=1000, tol=1e-4, beta0=length(z1)/5, info=TRUE){
	m<-length(z1)
	#initialization
	pi0_0<-0.9
	Pi1_0<-rep((1-pi0_0)/K,K)
	sigma2<-rgamma(K,1)
	Sigma_0<-list()
	for(i in 1:K){
		Sigma_0[[i]]<-matrix(rep(sigma2[i],4),2,2)
	}
	
	pi0_t<-pi0_0
	Pi1_t<-Pi1_0
	Sigma_t<-Sigma_0
	h<-list()
	h0<-0
	tmpH0<-0
	tmpH<-list()
	nanVal<-function(x) ifelse(is.nan(x),0,x)
	for(iter in 1:maxIter){
		#E step
		for(i in 0:K){
			if(i==0) tmpH0<-pi0_t*dmvnorm(cbind(z1,z2),sigma=diag(2))
			else tmpH[[i]]<-Pi1_t[i]*dmvnorm(cbind(z1,z2),sigma=diag(2)+Sigma_t[[i]])
		}
		norH<-tmpH0+Reduce('+',tmpH)
		h0<-nanVal(tmpH0/norH)
		h<-lapply(tmpH,FUN=function(x) return(nanVal(x/norH)))
		
		#M step
		pi0_t1<-(sum(h0)+beta0)/(m+beta0)
		Pi1_t1<-sapply(h,FUN=sum)/(m+beta0)
		z11<-z1^2; z12<-z1*z2; z22<-z2^2
		Sigma_t1<-lapply(h, FUN=function(x) {
					sumH<-sum(x)
					if(sumH==0) return(matrix(0,2,2))
					Sigma11<-max(sum(x*z11)/sumH-1,0)
					Sigma12<-max(sum(x*z12)/sumH,0)
					Sigma22<-max(sum(x*z22)/sumH-1,0)
					return(matrix(c(Sigma11,Sigma12,Sigma12,Sigma22),2,2))
				})
		if(abs(nanVal((pi0_t1-pi0_t)/pi0_t))<tol && sqrt(nanVal(sum((Pi1_t1-Pi1_t)^2)/sum(Pi1_t^2)))<tol && 
				Reduce('+',Map(function(x,y){return(ifelse(sum(y^2)==0, 0, sqrt(sum((x-y)^2)/sum(y^2))))},Sigma_t1,Sigma_t))<tol) break
		else{
			pi0_t<-pi0_t1
			Pi1_t<-Pi1_t1
			Sigma_t<-Sigma_t1
		}
	}
	logF<-function(x) ifelse(x==0,0,log(x))
	Qval<-sum(h0*logF(tmpH0))+sum(Reduce('+',Map('*',h,lapply(tmpH,logF))))
	if(info){
		cat('pi0:',pi0_t,'\n')
		cat('Pi1:\n')
		print(Pi1_t)
		cat('Sigma:\n')
		print(Sigma_t)
		cat('Iteration:',iter,'\nLog-likelihood:',Qval,'\n')
	}
	return(list(pi0=pi0_t,Pi1=Pi1_t,Sigma=Sigma_t,h0=h0,h=h,iter=iter,Qval=Qval))
}
	
#Mapping original z value pair to |z1|-sign(z1)z2 plane (Right Half Plane)
map2RHP<-function(z1,z2){
	newZ1<-abs(z1)
	newZ2<-ifelse(z1>=0,z2,-z2)
	return(list(z1=newZ1,z2=newZ2))
}

#Mapping original z value pair to sign(z2)z1-sign(z1)z2 plane (Cross Signed Plane)
map2CSP<-function(z1,z2){
	newZ1<-ifelse(z2>=0,z1,-z1)
	newZ2<-ifelse(z1>=0,z2,-z2)
	return(list(z1=newZ1,z2=newZ2))
}

#Bayesian Predictive Power (Probability under H1) with rectangle rejection region
#BPowerRect<-function(lower, SE1, SE2, sigma02, upper=c(Inf, Inf)){
BPowerRect<-function(lower, Pi1, Sigma,  upper=c(Inf, Inf)){	
	prob<-0
#	m<-length(SE1)
	K<-length(Pi1)
	if(lower[1]<0){
		upper<-map2RHP(lower[1],lower[2])
		upper<-c(upper$z1,upper$z2)
		lower<-map2RHP(upper[1],upper[2])
		lower<-c(lower$z1,lower$z2)
	}else{
		lower<-map2RHP(lower[1],lower[2])
		lower<-c(lower$z1,lower$z2)
		upper<-map2RHP(upper[1],upper[2])
		upper<-c(upper$z1,upper$z2)
	}
	for(i in 1:K)
		prob<-prob+2*Pi1[i]*pmvnorm(lower,upper,sigma=diag(2)+Sigma[[i]])
	prob<-prob/sum(Pi1)

#	for(i in 1:m){
#		Sigma<-matrix(c(1+sigma02/SE1[i]^2,sigma02/(SE1[i]*SE2[i]),sigma02/(SE1[i]*SE2[i]),1+sigma02/SE2[i]^2),2,2)
		#only consistent sign of z values is considered
#		prob<-prob+2*pmvnorm(lower,upper,sigma=Sigma)
#	}
#	prob<-prob/m
	
	return(prob[[1]])
}

#adjust local Jlfdr with monotonic in cross signed plane
adjustJlfdr2<-function(z1, z2, Jlfdr){
	newZ<-map2CSP(z1,z2)
	z1<-newZ$z1; z2<-newZ$z2
	
	#fdr must monotonic in each coordinates of sign(z2)z1-sign(z1)z2 plane
	orderIdx<-order(z1,z2,decreasing=T)
	for(i in 2:length(z1)){
		tmpZ2<-z2[orderIdx[i]]
		tmpIdx<-orderIdx[1:(i-1)][z2[orderIdx[1:(i-1)]]>=tmpZ2]
		if(length(tmpIdx)!=0){
			pos1<-tail(tmpIdx,1) #left most and lowest
			pos2<-tmpIdx[tail(which(z2[tmpIdx]==min(z2[tmpIdx])),1)] #lowest and left-most
			maxOfdr<-max(Jlfdr[pos1],Jlfdr[pos2])
			if(Jlfdr[orderIdx[i]]<maxOfdr) Jlfdr[orderIdx[i]]<-maxOfdr
		}
	}
	return(Jlfdr)
}

#Calculate Jlfdr 
#calJlfdr<-function(z1, z2, SE1, SE2, pi0, sigma02){
calJlfdr2<-function(z1, z2, Pi1, Sigma){
	f0<-dmvnorm(cbind(z1,z2))
#   m<-length(z1)
#   f1<-rep(0,m)
#	for(i in 1:m){
#		Sigma<-matrix(c(1+sigma02/SE1[i]^2,sigma02/(SE1[i]*SE2[i]),sigma02/(SE1[i]*SE2[i]),1+sigma02/SE2[i]^2),2,2)
#		f1<-f1+dmvnorm(cbind(z1,z2),sigma=Sigma)
#	}
#	f1<-f1/m
	K<-length(Sigma)
	pi0<-1-sum(Pi1)
	f1<-0
	for(i in 1:K)
		f1<-f1+Pi1[i]*dmvnorm(cbind(z1,z2),sigma=diag(2)+Sigma[[i]])
	f1<-f1/sum(Pi1)
	Jlfdr<-pi0*f0/(pi0*f0+(1-pi0)*f1)
	
#	Jlfdr<-adjustOfdr2(z1,z2,Jlfdr)
	return(Jlfdr)
}


#Predictive power (Probability under H1) with a determined region in the absolute value scale
# Region is a matrix recording the keypoints of z value pair with decreasing order
#predPowerRegion<-function(region, SE1, SE2, sigma02){
BPowerRegion<-function(region, Pi1, Sigma){
	k<-nrow(region)
#	prob<-BPowerRect(region[1,],SE1,SE2,sigma02)
	prob<-BPowerRect(region[1,],Pi1, Sigma)
	for(i in 2:k)
		prob<-prob+BPowerRect(region[i,],Pi1, Sigma, upper=c(region[i-1,1], Inf))
#		prob<-prob+BPowerRect(region[i,],SE1,SE2,sigma02,upper=c(region[i-1,1],Inf))
	return(prob[[1]])
}

#fdr to Fdr mapping
fdr2Fdr<-function(fdr){
	fdrVal<-sort(unique(fdr))
	FdrVal<-rep(0,length(fdrVal))
	sortfdr<-sort(fdr)
	m<-length(fdr)
	sumfdr<-sortfdr[1]
	FdrVal[1]<-sortfdr[1]
	j<-1
	for(i in 2:m){
		if(sortfdr[i]>sortfdr[i-1]){ 
			j<-j+1
			FdrVal[j]<-sumfdr/(i-1)
		}
		sumfdr<-sumfdr+sortfdr[i]
	}
	j<-j+1
	FdrVal[j]<-sumfdr/i
	return(list(fdrVal=fdrVal,FdrVal=FdrVal))
}

#Find fdr threshold such that corresponding Fdr<q
fdrThld<-function(fdrVal, FdrVal, q){
	pos<-tail((1:length(FdrVal))[FdrVal<=q],1)
	return(fdrVal[pos])
}

#Find points (z value pairs) in the rejection region based on fdr threshold
rejectPointsfdr<-function(fdr, fdrThld){
	reject<-(fdr<=fdrThld)
	return((1:length(fdr))[reject])
}

#Find the region of the points
points2region<-function(points){
#	orderIdx<-order(points[,1],points[,2],decreasing=T)
#	findKeyPts<-function(L,x){
#		if(points[x,1]<=points[L[[length(L)]],1] && points[x,2]<=points[L[[length(L)]],2]) L[[length(L)]]<-x
#		else L[[length(L)+1]]<-x
#		return(L)
#	}
#	keyPts<-Reduce(findKeyPts,orderIdx)
#	keyPts<-Reduce(findKeyPts,keyPts)
#	return(points[keyPts,])
	
	##Similar to Graham scan in finding convex hull
	orderIdx<-order(points[,1],points[,2],decreasing=T)
	tmpOrder<-orderIdx
	total<-length(orderIdx)
	remainNum<-length(orderIdx)
	keyPts<-NULL
	while(TRUE){
		lowpos<-which.min(points[tmpOrder,2])
		keyPts<-c(keyPts,total-remainNum+lowpos)
		remainNum<-remainNum-lowpos
		if(remainNum==0) break
		tmpOrder<-tmpOrder[-(1:lowpos)]
	}
	return(points[orderIdx[keyPts],])
}

plotRegion<-function(region,xlab='',ylab='',main='',...){
	keyNum<-nrow(region)
	points<-matrix(0,nrow=2*keyNum+1,ncol=2)
	params<-list(...)
	if(is.null(params$xlim)) points[1,1]<-region[1,1]*1.1
	else points[1,1]<-params$xlim[2]
	points[1,2]<-region[1,2]
	points[2,]<-region[1,]
	for(i in 2:keyNum){
		points[2*i-1,1]<-region[i-1,1]
		points[2*i-1,2]<-region[i,2]
		points[2*i,]<-region[i,]
	}
	points[2*keyNum+1,1]<-region[keyNum,1]
	if(is.null(params$ylim)) points[2*keyNum+1,2]<-region[keyNum,2]*1.1
	else points[2*keyNum+1,2]<-params$ylim[2]
	plot(points,xlab=xlab,ylab=ylab,main=main,type='l',col='blue',lwd=2,...)
#	points(region,pch=4,col='red',lwd=2)
}

#Fdr control via controlling Jlfdr
FdrControl2<-function(z1, z2, K=2, q=0.05, beta0=length(z1)/5, plot=T, output=T, dir='output'){
	EMresult<-snpEM2(z1,z2, K,beta0=beta0)
	Jlfdr<-EMresult$h0
#	Jlfdr<-adjustJlfdr2(z1,z2,EMresult$h0)
	
#	sampleNum<-1e5
#	componentIdx<-sample(0:K,sampleNum,replace=T,prob=c(EMresult$pi0,EMresult$Pi1))
#	randZ<-matrix(0,nrow=sampleNum,ncol=2)
#	randZ[componentIdx==0,]<-rmvnorm(sum(componentIdx==0),mean=rep(0,2))
#	for(i in 1:K){
#		randZ[componentIdx==i,]<-rmvnorm(sum(componentIdx==i),sigma=diag(2)+EMresult$Sigma[[i]])
#	}
#	randJlfdr<-calJlfdr2(randZ[,1],randZ[,2], EMresult$Pi1, EMresult$Sigma)
#	Jlfdr2FdrMap<-fdr2Fdr(randJlfdr)
#	JlfdrThldVal<-fdrThld(Jlfdr2FdrMap$fdrVal,Jlfdr2FdrMap$FdrVal,q)
#	tmpPts<-map2CSP(randZ[randJlfdr<=JlfdrThldVal,1],randZ[randJlfdr<=JlfdrThldVal,2])
#	region<-points2region(cbind(tmpPts$z1,tmpPts$z2))
	Jlfdr2FdrMap<-fdr2Fdr(Jlfdr)
	JlfdrThldVal<-fdrThld(Jlfdr2FdrMap$fdrVal,Jlfdr2FdrMap$FdrVal,q)
	rejected<-rejectPointsfdr(Jlfdr,JlfdrThldVal)
	rejectedPts<-map2RHP(z1[rejected],z2[rejected])
#	region<-points2region(cbind(rejectPtsByfdr$z1,rejectPtsByfdr$z2))
#	powerByfdr<-predPowerRegion(region, EMresult$Pi1, EMresult$Sigma)
	
	if(plot){
		if(missing(dir)) dev.new()
		else pdf(paste(dir,'/rejectedPts.pdf',sep=''))
		maxX<-max(abs(z1))*1.1; maxY<-max(sign(z1)*z2)*1.1
		minY<-min(sign(z1)*z2)*1.1
#		plotRegion(region,xlab=expression(paste('|',z^(1),'|')),ylab=expression(paste('sgn(',z^(1),') ',z^(2))),				
#				main='Rejected Points',xlim=c(0,maxX),ylim=c(0,maxY))
		plot(rejectedPts$z1,rejectedPts$z2,col='orange',pch=1,lwd=1.5,xlab=expression(paste('|',z^(1),'|')),ylab=expression(paste('sgn(',z^(1),') ',z^(2))),				
				main='Rejected Points',xlim=c(0,maxX),ylim=c(minY,maxY))
		if(!missing(dir)) dev.off()
	}
	if(output){
		pval1<-2*(1-pnorm(abs(z1)))
#		pval2<-ifelse(z1!=0,1-pnorm(sign(z1)*z2),2*(1-pnorm(abs(z2))))
		pval2<-2*(1-pnorm(abs(z2)))
		dir.create(dir,showWarnings=F)
		write.table(cbind(rejected, z1[rejected],z2[rejected],pval1[rejected],pval2[rejected],Jlfdr[rejected]),
				paste(dir,'rejected.txt',sep='/'), col.names=c('Index','Z1','Z2','P1','P2','Jlfdr'),row.names=F,sep='\t',quote=F)
		write.table(cbind(JlfdrThldVal),
				paste(dir,'threshold.txt',sep='/'), 
				col.names=c('Jlfdr-Thld'),row.names=F,sep='\t',quote=F)
	}
	
	return(list(Pi1=EMresult$Pi1,Sigma=EMresult$Sigma,rejected=rejected,JlfdrThld=JlfdrThldVal, Jlfdr=Jlfdr))
}
