##################################################
### Read in MCMC output from two component model

# load('BNPclassification - two component model.RData')

library(coda)
#	grDevices::devAskNewPage(ask=TRUE)




#################################
## Posterior estimates

# MH acceptance rates
apply(keep.MHcounter[,-1],1,mean)/keep.prop
	#want values to typically be between 25-45%

# Traceplots for adaptive MH parameters
plot(keep.adpMCMC[1,],type='l')
plot(keep.adpMCMC[2,],type='l')
plot(keep.adpMCMC[3,],type='l')
plot(keep.adpMCMC[4,],type='l')



# Disease probability
mean(keep.phi)
quantile(keep.phi,c(.025,.975))
plot(keep.phi,type='l')
geweke.diag(keep.phi)


# Longitudinal theta parameters
# row 1 for healthy
# row 2 for disease
apply(keep.theta,1:2,mean)
apply(keep.theta,1:2,quantile,prob=c(.025,.975))
apply(keep.theta,1:2, function(s) geweke.diag(s)$z)
plot(keep.theta[1,1,],type='l')
plot(keep.theta[1,2,],type='l')
plot(keep.theta[1,3,],type='l')
plot(keep.theta[2,1,],type='l')
plot(keep.theta[2,2,],type='l')
plot(keep.theta[2,3,],type='l')



# Estimated sig2/gamma2/rho for healthy, diseased patients
1/apply(1/keep.sigma2,1,mean)
1/apply(1/keep.gamma2,1,mean)
apply(keep.rho,1,mean)

apply(keep.sigma2,1,quantile,prob=c(.025,.975))
apply(keep.gamma2,1,quantile,prob=c(.025,.975))
apply(keep.rho,1,quantile,prob=c(.025,.975))

apply(keep.sigma2,1, function(s) geweke.diag(s)$z)
apply(keep.gamma2,1, function(s) geweke.diag(s)$z)
apply(keep.rho,1, function(s) geweke.diag(s)$z)

plot(1/keep.sigma2[1,],type='l')
plot(1/keep.sigma2[2,],type='l')
plot(1/keep.gamma2[1,],type='l')
plot(1/keep.gamma2[2,],type='l')
plot(keep.rho[1,],type='l')
plot(keep.rho[2,],type='l')



# Longitudinal Trend by Group
# Black=healthy; Blue=disease
plot(10:80,theta.function(10:80, apply(keep.theta,1:2,mean)[1,]), type='l')
points(10:80,theta.function(10:80, apply(keep.theta,1:2,mean)[2,]), 
	type='l',col='blue')





#################################
## Compute log-likelihood function by iteration
attach(data)
cluster.ind <- array(FALSE, c(n.pat,2))
cluster.ind[,1] <- as.logical(!Disease)
cluster.ind[,2] <- as.logical(Disease)

logHCG.clust <- Times.clust <- ID.clust <- Corr.matrix <-
	as.list(rep(NA,ncol(cluster.ind)))
for( k in 1:ncol(cluster.ind) ){
	logHCG.clust[[k]] <- logHCG[ ID %in% (1:n.pat)[cluster.ind[,k]] ]
	Times.clust[[k]] <- Times[ ID %in% (1:n.pat)[cluster.ind[,k]] ]
	ID.clust[[k]] <- ID[ ID %in% (1:n.pat)[cluster.ind[,k]] ]
}
library(mvtnorm)


loglik <- dbinom(sum(Disease, na.rm=T), sum(!is.na(Disease)), keep.phi, log=T)
for( k in 1:ncol(cluster.ind) ){
	logHCG.temp <- logHCG[ ID %in% (1:n.pat)[cluster.ind[,k]] ]
	Times.temp <- Times[ ID %in% (1:n.pat)[cluster.ind[,k]] ]
	for( iter in 1:niter ){
		Corr.temp <- keep.rho[k,iter]^( abs(array(Times.clust[[k]], 
			c(length(Times.clust[[k]]),length(Times.clust[[k]])))-
			t(array(Times.clust[[k]], c(length(Times.clust[[k]]),
			length(Times.clust[[k]])))))/7) * (1*(array(ID.clust[[k]], 
			c(length(ID.clust[[k]]),length(ID.clust[[k]])))==
			t(array(ID.clust[[k]], c(length(ID.clust[[k]]),
			length(ID.clust[[k]]))))))
		loglik[iter] <- loglik[iter] + dmvnorm( logHCG.temp, 
			theta.function(Times.temp,keep.theta[k,,iter]),
			keep.sigma2[k,iter]*Corr.temp + 
			keep.gamma2[k,iter]*(Corr.temp!=0), log=T)
	}
}
plot(loglik, type='l')
acf(loglik)
effectiveSize(loglik) 
geweke.diag(loglik)




mean(keep.alpha)
quantile(keep.alpha, c(.025,.975))
effectiveSize(keep.alpha)
geweke.diag(keep.alpha)
plot(keep.alpha, type='l', main='alpha')
acf(keep.alpha, main='alpha')



#################################
## Within Sample prediction accuracy


plot(apply(keep.test,1,mean), pch=21-6*(Disease==1), main='Training Data',
	ylab='Predicted Probability')
legend('topleft',pch=c(21,15),c('Healthy','Diseased'))

# Mean Square Predictive Error
mean((Disease - apply(keep.test,1,mean))^2)

# Misclassification Rate
mean(abs(Disease - apply(keep.test,1,mean))>.5)



