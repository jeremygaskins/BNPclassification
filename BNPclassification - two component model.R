burn.in <- 10000		#number of burn-in iterations
niter <- 2000 		#number of iterations to store
keep.prop <- 5		#proportion of post-burn-in iterations to store
				#total number of iterations is (burn.in + niter*keep.prop)
MCMC.seed <- 5111984
file.out <- 'BNPclassification - two component model.RData'



##################################################
### Simulation Specification Choices

library(mvtnorm); library(MCMCpack)

theta.function <- function(time,theta){
	theta[1]/(1+exp(-time*theta[2]-theta[3]))
}
# theta.function is the form of the longitudinal mean function

tr <- function(A) sum(diag(A))
rinvgamma <- function (n, shape = 1, scale = 1) 
{
    x <- rgamma(n = n, shape = shape, rate = scale)
    x[which(x < 1e-300)] <- 1e-300
    return(1/x)
}

# Parameters for the Adaptive MH step for theta
MH.cov <- array(diag(.01,3),c(3,3,2))
MH.mean <- array(0,c(3,2))	
MH.scale <- rep(1,2) 
MH.ss <- function(s){ (max(burn.in,iter))^(-.5) }
alpha.star <- 0.234

# MH proposal parameters for the covariance terms
MH.sig.sd <- 0.3
MH.rho.sd <- 0.25
MH.gam.sd <- 1.2


##################################################
### Prior Choices

phi.alpha <- .5
phi.beta <- .5
theta.mean <- rep(1,3)
theta.cov <- diag(100, 3)
sigma2.shape <- 0.1
sigma2.scale <- 0.1
gamma2.shape <- 0.1
gamma2.scale <- 0.1

corr.unit <- 7
# rho corresponds to the correlation between observations corr.unit apart



##################################################
### Load Data

load('BNCLD_sim1_1.RData')
# The RData file contains 3 objects defining the dataset.
# data is 3-column data.frame.
###  ID represents the ID number corresponding to which patient that measurement corresponds.
###  NOTE: IDs must be consecutive.
###  Times gives the observation time for the measurement.
###  logHCG gives the longitudinal response value for that patient ID and measurement time.
# Disease gives the binary indicators of classification/disease status.
# n.pat gives the total number of patients/observational units.

attach(data)




##################################################
### Initialize MCMC parameters

phi <- .5
theta <- array(rnorm(6), c(2,3))
sigma2 <- c(1,1)
gamma2 <- c(1,1)
rho <- c(.5,.5)

cluster.ind <- array(FALSE, c(n.pat,2))
cluster.ind[,1] <- as.logical(!Disease)
cluster.ind[,2] <- as.logical(Disease)

logHCG.clust <- Times.clust <- fhat.clust <- ID.clust <- Corr.matrix <-
	as.list(rep(NA,ncol(cluster.ind)))
for( k in 1:ncol(cluster.ind) ){
	logHCG.clust[[k]] <- logHCG[ ID %in% (1:n.pat)[cluster.ind[,k]] ]
	Times.clust[[k]] <- Times[ ID %in% (1:n.pat)[cluster.ind[,k]] ]
	fhat.clust[[k]] <- 0*Times.clust[[k]]
	ID.clust[[k]] <- ID[ ID %in% (1:n.pat)[cluster.ind[,k]] ]
	Corr.matrix[[k]] <- rho[k]^( abs(array(Times.clust[[k]], 
		c(length(Times.clust[[k]]),length(Times.clust[[k]])))-
		t(array(Times.clust[[k]], c(length(Times.clust[[k]]),
		length(Times.clust[[k]])))))/corr.unit) * (1*(array(ID.clust[[k]], 
		c(length(ID.clust[[k]]),length(ID.clust[[k]])))==
		t(array(ID.clust[[k]], c(length(ID.clust[[k]]),
		length(ID.clust[[k]]))))))
}

theta[,1] <- as.numeric( lapply(logHCG.clust,max) ) *1.01
theta[,2:3] <- t(sapply(1:ncol(cluster.ind), function(s)
	-as.numeric(lm(log(theta[s,1]/logHCG.clust[[s]]-1) ~ 
		Times.clust[[s]])$coef)[2:1]))
MH.mean <- t(theta)
theta <- apply(theta,1:2, function(s) rnorm(1,s,.5) )

for( k in 1:length(fhat.clust) ){
	fhat.clust[[k]] <- theta.function(Times.clust[[k]],theta[k,])
}

MH.counter <- rep(0,8)
	# MH.counter stores the number of times the MH steps accept the move
	# MH.counter[1:2] reports accepting theta for control and disease, respectively
	# MH.counter[3:4] reports accepting sig2 for control and disease, respectively
	# MH.counter[5:6] reports accepting rho for control and disease, respectively
	# MH.counter[7:8] reports accepting gam2 for control and disease, respectively
test.prob <- rep(0,n.pat)


keep.phi <- array(0, dim=c(dim(phi),floor(niter)))
keep.theta <- array(0, dim=c(dim(theta),floor(niter)))
keep.sigma2 <- array(0, dim=c(length(sigma2),floor(niter)))
keep.gamma2 <- array(0, dim=c(length(gamma2),floor(niter)))
keep.rho <- array(0, dim=c(length(rho),floor(niter)))
keep.MHcounter <- array(0, dim=c(length(MH.counter),floor(niter)))
keep.adpMCMC <- array(0, dim=c(4,floor(niter)))
keep.test <- array(0, dim=c(n.pat,niter))

set.seed(MCMC.seed)
time1 <- Sys.time()
for( iter in 1:(burn.in + niter*keep.prop) ){


##################################################
### Sample phi

phi <- rbeta(1, sum(Disease, na.rm=T)+phi.alpha, 
	sum(1-Disease, na.rm=T)+phi.beta)




##################################################
### Sample theta

for( k in 1:nrow(theta) ){
	theta.cand <- rmvnorm(1, theta[k,], MH.scale[k]*MH.cov[,,k])
	fhat.cand <- theta.function(Times.clust[[k]],theta.cand)

	MH.prob <- -.5 * tr( solve(sigma2[k]*Corr.matrix[[k]] + 
		gamma2[k]*(Corr.matrix[[k]]!=0)) %*% 
		( (logHCG.clust[[k]] - fhat.cand) %*%
		t(logHCG.clust[[k]] - fhat.cand) -
		(logHCG.clust[[k]] - fhat.clust[[k]]) %*%
		t(logHCG.clust[[k]] - fhat.clust[[k]])) ) +
		dmvnorm(theta.cand, theta.mean, theta.cov, log=T) -
		dmvnorm(theta[k,], theta.mean, theta.cov, log=T)
	if( runif(1) < exp(min(MH.prob,0)) ){
		theta[k,] <- theta.cand
		MH.counter[k] <- MH.counter[k]+1
		fhat.clust[[k]] <- fhat.cand
	}
	MH.mean[,k] <- MH.mean[,k] + MH.ss(iter)*(theta[k,] - MH.mean[,k])
	MH.cov[,,k] <- MH.cov[,,k] + MH.ss(iter)*( 
		(theta[k,] - MH.mean[,k])%*%t(theta[k,] - MH.mean[,k]) - 
		MH.cov[,,k])
	MH.scale[k] <- exp( log(MH.scale[k]) + MH.ss(iter)*(
		exp(min(MH.prob,0)) - alpha.star ) )
}




##################################################
### Sample sigma2

for( k in 1:nrow(theta) ){
	sig.cand <- rlnorm(1, log(sigma2[k]), MH.sig.sd)
	MH.prob <- dlnorm(sigma2[k], log(sig.cand), MH.sig.sd, log=T) - 
		dlnorm(sig.cand, log(sigma2[k]), MH.sig.sd, log=T) + 
		log(dinvgamma(sig.cand, sigma2.shape, sigma2.scale)) - 
		log(dinvgamma(sigma2[k], sigma2.shape, sigma2.scale))
	MH.prob <- MH.prob - 
		.5*tr( (solve(sig.cand*Corr.matrix[[k]]+
		gamma2[k]*(Corr.matrix[[k]]!=0)) - 
		solve(sigma2[k]*Corr.matrix[[k]]+
		gamma2[k]*(Corr.matrix[[k]]!=0))) %*%
		(logHCG.clust[[k]] - fhat.clust[[k]] ) %*% 
		t(logHCG.clust[[k]] - fhat.clust[[k]]) ) - 
		.5*as.numeric(determinant(sig.cand*Corr.matrix[[k]]+
		gamma2[k]*(Corr.matrix[[k]]!=0), log=T)$mod) + 
		.5*as.numeric(determinant(sigma2[k]*Corr.matrix[[k]]+
		gamma2[k]*(Corr.matrix[[k]]!=0),log=T)$mod)
	if( runif(1) < exp(min(MH.prob,0)) ){
		sigma2[k] <- sig.cand
		MH.counter[k+2] <- MH.counter[k+2]+1
	}
}




##################################################
### Sample rho

for( k in 1:nrow(theta) ){
	rho.cand <- rnorm(1, rho[k], MH.rho.sd)
	if( rho.cand>1 ){ rho.cand <- 2 - rho.cand }
	if( rho.cand<0 ){ rho.cand <- -rho.cand }
	Corr.temp <- rho.cand^( abs(array(Times.clust[[k]], 
		c(length(Times.clust[[k]]),length(Times.clust[[k]])))-
		t(array(Times.clust[[k]], c(length(Times.clust[[k]]),
		length(Times.clust[[k]])))))/corr.unit) * (1*(array(ID.clust[[k]], 
		c(length(ID.clust[[k]]),length(ID.clust[[k]])))==
		t(array(ID.clust[[k]], c(length(ID.clust[[k]]),
		length(ID.clust[[k]]))))))
	MH.prob <- - .5*tr( (solve(sigma2[k]*Corr.temp +
		gamma2[k]*(Corr.matrix[[k]]!=0)) - 
		solve(sigma2[k]*Corr.matrix[[k]]+
		gamma2[k]*(Corr.matrix[[k]]!=0))) %*%
		(logHCG.clust[[k]] - fhat.clust[[k]] ) %*% 
		t(logHCG.clust[[k]] - fhat.clust[[k]]) ) - 
		.5*as.numeric(determinant(sigma2[k]*Corr.temp +
		gamma2[k]*(Corr.matrix[[k]]!=0),log=T)$mod) + 
		.5*as.numeric(determinant(sigma2[k]*Corr.matrix[[k]]+
		gamma2[k]*(Corr.matrix[[k]]!=0),log=T)$mod)
	if( (runif(1) < exp(min(MH.prob,0))) & (rho.cand<1) & (rho.cand>0) ){
		rho[k] <- rho.cand
		Corr.matrix[[k]] <- Corr.temp
		MH.counter[4+k] <- MH.counter[4+k]+1
	}
}




##################################################
### Sample gamma2

for( k in 1:nrow(theta) ){
	gam.cand <- rlnorm(1, log(gamma2[k]), MH.gam.sd)
	MH.prob <- dlnorm(gamma2[k], log(gam.cand), MH.gam.sd, log=T) - 
		dlnorm(gam.cand, log(gamma2[k]), MH.gam.sd, log=T) + 
		log(dinvgamma(gam.cand, gamma2.shape, gamma2.scale)) - 
		log(dinvgamma(gamma2[k], gamma2.shape, gamma2.scale))
	MH.prob <- MH.prob - 
		.5*tr( (solve(sigma2[k]*Corr.matrix[[k]]+
		gam.cand*(Corr.matrix[[k]]!=0)) - 
		solve(sigma2[k]*Corr.matrix[[k]]+
		gamma2[k]*(Corr.matrix[[k]]!=0))) %*%
		(logHCG.clust[[k]] - fhat.clust[[k]] ) %*% 
		t(logHCG.clust[[k]] - fhat.clust[[k]]) ) - 
		.5*as.numeric(determinant(sigma2[k]*Corr.matrix[[k]]+
		gam.cand*(Corr.matrix[[k]]!=0), log=T)$mod) + 
		.5*as.numeric(determinant(sigma2[k]*Corr.matrix[[k]]+
		gamma2[k]*(Corr.matrix[[k]]!=0),log=T)$mod)
	if( runif(1) < exp(min(MH.prob,0)) ){
		gamma2[k] <- gam.cand
		MH.counter[k+6] <- MH.counter[k+6]+1
	}
}




##################################################
### Keep variables

if( (floor((iter-burn.in)/keep.prop)==(iter-burn.in)/keep.prop) * 
	(iter > burn.in) ){

# Computes the predictive probability for each patient
for( i in 1:n.pat ){
	set <- (i==ID)
	loglik <- rep(NA,nrow(theta))
	loglik <- log(c(1-phi,phi))
	for( k in 1:nrow(theta) ){
		Var.temp <- sigma2[k]*rho[k]^( abs(array(Times[set], 
			c(length(Times[set]),length(Times[set])))-
			t(array(Times[set], c(length(Times[set]),
			length(Times[set])))))/corr.unit) + gamma2[k]
		loglik[k] <- loglik[k] + dmvnorm(logHCG[set],
			theta.function(Times[set], theta[k,]),
			Var.temp, log=T)
	}
	loglik <- exp(loglik - max(loglik))
	loglik <- loglik/sum(loglik)
	test.prob[i] <- loglik[2]
}


	q <- (iter - burn.in) %/% keep.prop
	keep.phi[q] <- phi
	keep.theta[,,q] <- theta
	keep.sigma2[,q] <- sigma2
	keep.gamma2[,q] <- gamma2
	keep.rho[,q] <- rho
	keep.MHcounter[,q] <- MH.counter
	keep.adpMCMC[,q] <- c(log(apply(MH.cov,3,det)),MH.scale)
	keep.test[,q] <- test.prob
	MH.counter <- 0*MH.counter
}
if( iter %in% seq(10000, burn.in+niter*keep.prop, by=10000) ){
	print(paste('Iteration Number',iter,' out of',burn.in+niter*keep.prop,
		';    ', Sys.time() )) 
	flush.console()
}

}
finish.time <- Sys.time() - time1
print(finish.time)
warnings()




save(burn.in, niter, keep.prop, MCMC.seed, file.out, phi.alpha, 
	phi.beta, theta.mean, theta.cov, sigma2.shape, sigma2.scale, 
	gamma2.shape, gamma2.scale, corr.unit,
	theta.function, tr, MH.cov, MH.mean, MH.scale, MH.ss, alpha.star, 
	MH.sig.sd, MH.rho.sd, MH.gam.sd, data, Disease, n.pat, 
	keep.phi, keep.theta, keep.sigma2, keep.gamma2, keep.rho, 
	keep.MHcounter, keep.adpMCMC, keep.test, finish.time,  
	file=file.out)
