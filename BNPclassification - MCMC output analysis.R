##################################################
### Read in MCMC output from either the full BNP or stage 2 model

# load('BNPclassification - full BNP model.RData')
# load('BNPclassification - full BNP model_method1.RData')

library(coda)
#	grDevices::devAskNewPage(ask=TRUE)




#################################
## Cluster Membership Characterization
## Full BNP model

clust.pairs <- array(NA, c(n.pat, n.pat, niter))
for( i in 1:n.pat ){ for( j in i:n.pat ){
	clust.pairs[i,j,] <- clust.pairs[j,i,] <- keep.clust[i,]==keep.clust[j,]
}}
clust.hat <- apply(clust.pairs,1:2,mean)
image( (1:n.pat), (1:n.pat), clust.hat, zlim=c(0,1), 
	col=gray((32:0)/32),
	axes=T, xlab='', ylab='', main='Pairwise cluster probs')
box()

## clust.import selects that clusters that have an average of at least one 
## observation for viewing full BNP output.
clust.import <- sapply(1:H.max, function(s) 
	mean(apply(keep.clust==s,2,sum))>.99 )
# to show all clusters use
# clust.import <- rep( TRUE, H.max)


# Number of non-empty clusters
# blue is rolling average
clust.nonemp <- apply(keep.clust,2,function(s) length(unique(s)) )
clust.med <- apply(keep.clust,2,function(s) sum(table(s)>5) )
plot( clust.nonemp, type='l', main='Number of non-empty clusters' )
points( sapply(1:niter, function(s)
	mean( clust.nonemp[max(1,s-100):min(niter,s+100)]) ), type='l', 
	col='blue', lwd=2)
acf(clust.nonemp, main='Number of non-empty clusters' )
effectiveSize(clust.nonemp)
geweke.diag(clust.nonemp)

mean(clust.nonemp)
quantile(clust.nonemp, c(.025,.975))



plot( clust.med, type='l', main='Number of med clusters (n_h>5)' )
points( sapply(1:niter, function(s)
	mean( clust.med[max(1,s-100):min(niter,s+100)]) ), type='l', 
	col='blue', lwd=2)
acf( clust.med, main='Number of med clusters (n_h>5)' )
effectiveSize(clust.med)
geweke.diag(clust.med)


# Number of patients per cluster
for( k in (1:H.max)[clust.import] ){
	plot(1:niter, apply(keep.clust==k,2,sum), type='l',  
		main=paste('Cluster',k))
}



#################################
## Cluster Membership Characterization
## Stage 2 Model

optimal.method
	#1=Dahl(MCMC); 2=HCL Average, h=.75; 3=HCL Avg, k=median number of clusters
	#4=HCL Ward, k=median; 6=Avg Silhouette, 7=Ward Silhouette
	#8=Avg Gamma, 9=Avg Tau, 10=Ward Gamma, 11=Ward Tau
	#12=VI, 13=Dahl(all)
H.max
apply(cluster.ind,2,sum)
clust.hat <- cluster.ind %*% t(cluster.ind)
image( (1:n.pat), (1:n.pat), clust.hat, zlim=c(0,1), 
	col=gray((32:0)/32),
	axes=T, xlab='', ylab='', main='Pairwise cluster probs')
box()
keep.clust <- array(apply(cluster.ind,1,which.max),c(n.pat,niter))
clust.import <- rep( TRUE, H.max)

# estimated of cluster probabilities
apply( cluster.ind, 2, sum)/n.pat

apply( keep.clust.prob, 1, mean)
apply( keep.clust.prob, 1, quantile, c(.025,.975))
apply( keep.clust.prob, 1, effectiveSize)

for( k in (1:H.max)[clust.import]){
	plot(1:niter, keep.clust.prob[k,], type='l', ylim=c(0,1),
		main=paste('Membership Prob Cluster',k))
}


#################################
## MCMC posterior samples
## Note: In full BNP many are subject to label-switching, and these plots may not be informative.



## Dirichlet concentration parameter
## (full BNP only)
mean(keep.alpha)
quantile(keep.alpha, c(.025,.975))
effectiveSize(keep.alpha)
geweke.diag(keep.alpha)
plot(keep.alpha, type='l', main='alpha')
acf(keep.alpha, main='alpha')




## Disease probability
phi.hat <- apply(keep.phi, 1, mean)
phi.hat[clust.import]
apply(keep.phi, 1, quantile, c(.025,.975))[,clust.import]
effectiveSize(t(keep.phi))[clust.import]
geweke.diag(t(keep.phi))$z[clust.import]

for( k in (1:H.max)[clust.import] ){
	plot(1:niter, keep.phi[k,], type='l', 
		main=paste('Disease Prob Cluster',k))
	for( iter in (1:niter)[apply(keep.clust==k,2,sum)==0] ){
		polygon( c(iter-.75,iter-.75,iter+.75,iter+.75), 
			c(-.1,1.1,1.1,-0.1), col='white',lty=0)
	}
}



# sigma
1/mean(1/keep.sigma2)
quantile(keep.sigma2,c(.025,.975))
effectiveSize(keep.sigma2)
geweke.diag(keep.sigma2)
plot(keep.sigma2, type='l', main='sigma2')
acf(keep.sigma2, main='sigma2')
# MH acceptance rate for sigma2:
mean(keep.MHcounter[2,-1])/keep.prop



# gamma
1/mean(1/keep.gamma2)
quantile(keep.gamma2,c(.025,.975))
effectiveSize(keep.gamma2)
geweke.diag(keep.gamma2)
plot(keep.gamma2, type='l', main='gamma2')
acf(keep.gamma2, main='gamma2')
# MH acceptance rate for gamma2:
mean(keep.MHcounter[3,-1])/keep.prop


# rho
mean(keep.rho)
quantile(keep.rho, c(.025,.975))
effectiveSize(keep.rho)
geweke.diag(keep.rho)
plot(keep.rho, type='l', main='rho')
acf(keep.rho, main='rho')
# MH acceptance rate for rho:
mean(keep.MHcounter[4,-1])/keep.prop



# Longitudinal theta parameters
# MH acceptance rate for an average cluster, per iteration
mean(keep.MHcounter[1,-1])/keep.prop
plot(apply(keep.adpMCMC,2,mean), type='l', main='Adaptive MCMC parameter')


theta.hat <- apply(keep.theta, 1:2, mean)
cbind(theta.hat,phi.hat)[ clust.import, ]
apply(keep.theta,1:2,quantile, c(.025,.975))[,clust.import,]
apply(keep.theta,1:2, function(s) effectiveSize(s) )[clust.import,]
apply(keep.theta,1:2, function(s) geweke.diag(s)$z )[clust.import,]



for( k in (1:H.max)[clust.import] ){ for( j in 1:3 ){
	plot(1:niter, keep.theta[k,j,], type='l', 
		main=paste('Theta Paramter',j, 'for Cluster',k))
	for( iter in (1:niter)[apply(keep.clust==k,2,sum)==0] ){
		polygon( c(iter-.75,iter-.75,iter+.75,iter+.75), c(-20,20,20,-20),
			col='white',lty=0)
	}
}}


# Plot of average trend 
color.set <- c('black','red','blue','orange','green','yellow','grey','purple')
plot(NULL, xlim=c(10,80), ylim=c(1,5), ylab='')
counter <- 1
for( k in (1:H.max)[clust.import] ){
	points(10:80, apply(apply(keep.theta[k,,], 2, function(s) 
		theta.function(10:80, s)),1,mean), type='l',
		col=color.set[counter], lwd=2)
	text(15+counter*5,1.2, col=color.set[counter],
		round(mean(keep.phi[k,]),2))
	text(15+counter*5,1, col=color.set[counter],
		round(mean(apply(keep.clust==k,2,sum))))
	counter <- counter+1
}




# Theta.star
apply(keep.th.star,1,mean)	
apply(keep.th.star,1,quantile,c(.025,.975))	
effectiveSize(t(keep.th.star))
geweke.diag(t(keep.th.star))
for( k in 1:3 ){
	plot( keep.th.star[k,], type='l', main=paste('Theta Star', k))
	acf( keep.th.star[k,], main=paste('Theta Star', k))
}


# Cov matrix of theta
apply(keep.th.cov,1:2,mean)
th.cov.hat <- solve(array(apply(apply(keep.th.cov,3,solve),1,mean), c(3,3)))
th.cov.hat
apply(keep.th.cov,1:2,effectiveSize)
apply(keep.th.cov,1:2,function(s) geweke.diag(s)$z)
for( j in 1:3 ){ for( k in j:3 ){
	plot( keep.th.cov[j,k,], type='l', main=paste('Theta Cov',j,k))
	acf( keep.th.cov[j,k,], main=paste('Theta Cov', j,k))
}}


#################################
## Compute log-likelihood function by iteration
attach(data)
library(mvtnorm)
loglik.marg <- rep(0, niter)
loglik.clust <- rep(0, niter)
for( iter in 1:niter ){
for( i in (1:n.pat) ){
	set <- (i==ID)
	loglik.i <- rep(0,dim(keep.theta)[1])
	Var.temp <- keep.sigma2[iter]*keep.rho[iter]^( abs(array(Times[set], 
		c(length(Times[set]),length(Times[set])))-
		t(array(Times[set], c(length(Times[set]),
		length(Times[set])))))/corr.unit) + keep.gamma2[iter]
	for( k in 1:(dim(keep.theta)[1]) ){
		loglik.i[k] <- loglik.i[k] + dmvnorm(logHCG[set],
			theta.function(Times[set], keep.theta[k,,iter]),
			Var.temp, log=T)
	}
	loglik.i <- loglik.i + dbinom(Disease[i],1,keep.phi[,iter],log=T)
	loglik.clust[iter] <- loglik.clust[iter] + 
		loglik.i[ keep.clust[i,iter] ]
	loglik.i <- loglik.i + log(keep.clust.prob[,iter])
	loglik.marg[iter] <- loglik.marg[iter] + max(loglik.i) + 
		log(sum(exp(loglik.i - max(loglik.i))))
	}
}
plot(loglik.marg, type='l', main='Log-likelihood (marginal over cluster)')
acf(loglik.marg)

plot(loglik.clust, type='l', main='Log-likelihood (conditional on current cluster)')
acf(loglik.clust)

effectiveSize(loglik.marg) 
effectiveSize(loglik.clust) 

geweke.diag(loglik.marg)
geweke.diag(loglik.clust)




#################################
## Within Sample prediction accuracy


plot(apply(keep.test,1,mean), pch=21-6*(Disease==1), main='Training Data',
	ylab='Predicted Probability')
legend('topleft',pch=c(21,15),c('Healthy','Diseased'))

# Mean Square Predictive Error
mean((Disease - apply(keep.test,1,mean))^2)

# Misclassification Rate
mean(abs(Disease - apply(keep.test,1,mean))>.5)





