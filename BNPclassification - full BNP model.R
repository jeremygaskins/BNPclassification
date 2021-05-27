##################################################
### Read in output from 2-component model to initialize adaptive MH choices

load('BNPclassification - two component model.RData')

MH.cov <- mean(Disease==0, na.rm=T)*MH.cov[,,1] + 
	mean(Disease==1, na.rm=T)*MH.cov[,,2]
MH.scale <- exp( mean(Disease==0, na.rm=T)*log(MH.scale[1]) + 
	mean(Disease==1, na.rm=T)*log(MH.scale[2]) )
rm( list=setdiff(ls(), c('MH.cov','MH.scale','B')) )


burn.in <- 5000		#number of burn-in iterations
niter <- 2000		#number of iterations to store
keep.prop <- 50		#proportion of post-burn-in iterations to store
				#total number of iterations is (burn.in + niter*keep.prop)
H.max <- 30			#maximum number of clusters
MCMC.seed <- 5111984
file.out <- 'BNPclassification - full BNP model.RData'




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
# if not reading in 2 component model fit, need to initialize the following
# MH.cov <- diag(0.01,3)
# MH.scale <- 1
logMH.scale <- rep(log(MH.scale),H.max)
MH.ss <- function(s){ (max(burn.in,iter))^(-.5) }
alpha.star <- 0.234

# MH proposal parameters for the covariance terms
MH.rho.sd <- 0.3
MH.sig.sd <- 0.3
MH.gam.sd <- 0.5




##################################################
### Prior Choices

phi.alpha <- .5
phi.beta <- .5
theta.mean <- rep(1,3)
theta.prior.cov <- diag(100, 3)
sigma2.shape <- 0.1
sigma2.scale <- 0.1
alpha.shape <- 1
alpha.rate <- 1
theta.cov.df <- 5
theta.cov.scale <- diag(3)
gamma2.shape <- 0.1
gamma2.scale <- 0.1

alpha.sample <- TRUE
	#if FALSE, the DP concentration parameter will not be updated and will remain at alpha.DP
alpha.DP <- 1

corr.unit <- 7
# rho corresponds to the correlation between observations corr.unit apart




##################################################
### Load Data

load('BNCLD_sim1_1.RData')
attach(data)




##################################################
### Initial MCMC

phi <- rep(.5, H.max)
theta <- array(rnorm(3*H.max), c(H.max,3))
sigma2 <- 1
rho <- .5
gamma2 <- 0.5

V <- rep(.5, H.max)
V[H.max] <- 1
cluster.prob <- V*c(1,cumprod(1-V)[-H.max])

h.max <- min(10,H.max)
cluster.ind <- t(rmultinom(n.pat, 1, prob=rep(1/h.max,h.max)))
if( h.max != H.max ){
	cluster.ind <- cbind( cluster.ind, array(0, c(n.pat, H.max-h.max)) )
}
# cluster.ind[test.ID,] <- 0
cluster.ind <- apply(cluster.ind, 1:2, as.logical)


logHCG.clust <- Times.clust <- fhat.clust <- ID.clust <- Corr.matrix <-
	as.list(rep(NA,ncol(cluster.ind)))
for( k in 1:ncol(cluster.ind) ){
	logHCG.clust[[k]] <- logHCG[ ID %in% (1:n.pat)[cluster.ind[,k]] ]
	Times.clust[[k]] <- Times[ ID %in% (1:n.pat)[cluster.ind[,k]] ]
	fhat.clust[[k]] <- 0*Times.clust[[k]]
	ID.clust[[k]] <- ID[ ID %in% (1:n.pat)[cluster.ind[,k]] ]
	Corr.matrix[[k]] <- rho^( abs(array(Times.clust[[k]], 
		c(length(Times.clust[[k]]),length(Times.clust[[k]])))-
		t(array(Times.clust[[k]], c(length(Times.clust[[k]]),
		length(Times.clust[[k]])))))/corr.unit) * (1*(array(ID.clust[[k]], 
		c(length(ID.clust[[k]]),length(ID.clust[[k]])))==
		t(array(ID.clust[[k]], c(length(ID.clust[[k]]),
		length(ID.clust[[k]]))))))
}

theta[1:h.max,1] <- as.numeric( lapply(logHCG.clust[1:h.max],max) ) *1.01
theta[1:h.max,2:3] <- t(sapply(1:h.max, function(s)
	-as.numeric(lm(log(theta[s,1]/logHCG.clust[[s]]-1) ~ 
		Times.clust[[s]])$coef)[2:1]))

	
theta[is.na(theta)] <- 0
theta <- apply(theta,1:2, function(s) rnorm(1,s,.5) )
theta.star <- apply(theta,2,mean)
theta.cov <- cov(theta)

for( k in 1:length(fhat.clust) ){
	fhat.clust[[k]] <- theta.function(Times.clust[[k]],theta[k,])
}

MH.counter <- rep(0,4)
	# MH.counter stores the number of times the MH steps accept the move
	# MH.counter[1] reports accepting theta per number of non-empty clusters
	# MH.counter[2] reports accepting sig2
	# MH.counter[3] reports accepting gam2
	# MH.counter[4] reports accepting rho
test.prob <- rep(0,n.pat)

keep.clust.prob <- array(0, c(H.max,niter))
keep.alpha <- rep(0, niter)
keep.clust <- array(0, c(n.pat,niter))
keep.phi <- array(0, dim=c(length(phi),floor(niter)))
keep.theta <- array(0, dim=c(dim(theta),floor(niter)))
keep.sigma2 <- rep(0, niter)
keep.gamma2 <- rep(0, niter)
keep.th.star <- array(0, c(length(theta.star),niter))
keep.th.cov <- array(0, c(dim(theta.cov),niter))
keep.rho <- rep(0,niter)
keep.MHcounter <- array(0, dim=c(length(MH.counter),floor(niter)))
keep.adpMCMC <- array(0, dim=c(H.max,floor(niter)))
keep.test <- array(0, dim=c(n.pat,niter))

set.seed(MCMC.seed)
time1 <- Sys.time()
for( iter in 1:(burn.in + niter*keep.prop) ){


##################################################
### Sample cluster.prob

for( k in 1:(H.max-1) ){
	V[k] <- min(1-1e-8,rbeta( 1, 1+sum(cluster.ind[,k]), 
		alpha.DP + sum(cluster.ind[,(k+1):H.max])))
}
cluster.prob <- V*c(1,cumprod(1-V)[-H.max])


if( alpha.sample ){
	alpha.DP <- rgamma(1, H.max-1+alpha.shape, 
		rate=alpha.rate-sum(log(1-V[-H.max])))
}



##################################################
### Sample cluster membership

for( i in (1:n.pat) ){
	Corr.temp <- rho^( abs(array(Times[ID==i], c(sum(ID==i),sum(ID==i))) - 
		t(array(Times[ID==i], c(sum(ID==i),sum(ID==i))))) )
	cluster.mean <- array(apply(theta, 1, function(s) 
		theta.function( Times[ID==i], s) ), c(sum(ID==i),H.max))
	cluster.Gibbs <- log(cluster.prob) + dbinom(Disease[i], 1, phi, log=T) +
		sapply(1:H.max, function(s) dmvnorm( logHCG[ID==i], 
		cluster.mean[,s], sigma2*Corr.temp + gamma2, log=T))
	cluster.Gibbs <- exp( cluster.Gibbs - max(cluster.Gibbs) )
	cluster.ind[i,] <- rmultinom(1,1,cluster.Gibbs)
}
cluster.ind <- apply(cluster.ind, 1:2, as.logical)
for( k in 1:ncol(cluster.ind) ){
	logHCG.clust[[k]] <- logHCG[ ID %in% (1:n.pat)[cluster.ind[,k]] ]
	Times.clust[[k]] <- Times[ ID %in% (1:n.pat)[cluster.ind[,k]] ]
	fhat.clust[[k]] <- theta.function(Times.clust[[k]],theta[k,])
	ID.clust[[k]] <- ID[ ID %in% (1:n.pat)[cluster.ind[,k]] ]
	Corr.matrix[[k]] <- rho^( abs(array(Times.clust[[k]], 
		c(length(Times.clust[[k]]),length(Times.clust[[k]])))-
		t(array(Times.clust[[k]], c(length(Times.clust[[k]]),
		length(Times.clust[[k]])))))/corr.unit) * (1*(array(ID.clust[[k]], 
		c(length(ID.clust[[k]]),length(ID.clust[[k]])))==
		t(array(ID.clust[[k]], c(length(ID.clust[[k]]),
		length(ID.clust[[k]]))))))
}
clust.nonempty <- sum(apply(cluster.ind,2,sum)>0)
logMH.scale[apply(cluster.ind,2,sum)==0] <- mean(logMH.scale)



##################################################
### Sample phi

phi <- rbeta(H.max, 
	phi.alpha + sapply(1:H.max, function(s) sum(Disease[cluster.ind[,s]]) ),
	phi.beta + apply(cluster.ind,2,sum) - 
		sapply(1:H.max, function(s) sum(Disease[cluster.ind[,s]]) ))



##################################################
### Sample theta

for( k in (1:H.max)[apply(cluster.ind,2,sum)==0] ){
	theta[k,] <- rmvnorm(1, theta.star, theta.cov)
}
for( k in (1:H.max)[apply(cluster.ind,2,sum)>0] ){
	theta.cand <- rmvnorm(1, theta[k,], MH.cov * 
		exp(logMH.scale[k])/sqrt(length(ID.clust[[k]])) )
	fhat.cand <- theta.function(Times.clust[[k]],theta.cand)

	MH.prob <- -.5 * tr( solve(sigma2*Corr.matrix[[k]] + 
		gamma2*(Corr.matrix[[k]]!=0)) %*% 
		( (logHCG.clust[[k]] - fhat.cand) %*%
		t(logHCG.clust[[k]] - fhat.cand) -
		(logHCG.clust[[k]] - fhat.clust[[k]]) %*%
		t(logHCG.clust[[k]] - fhat.clust[[k]])) ) +
		dmvnorm(theta.cand, theta.star, theta.cov, log=T) -
		dmvnorm(theta[k,], theta.star, theta.cov, log=T)
	if( runif(1) < exp(min(MH.prob,0)) ){
		theta[k,] <- theta.cand
		MH.counter[1] <- MH.counter[1] + 1/clust.nonempty
		fhat.clust[[k]] <- fhat.cand
	}
	logMH.scale[k] <- logMH.scale[k] + MH.ss(iter)*(
		exp(min(MH.prob,0)) - alpha.star )
}



##################################################
### Sample theta.star, theta.cov

theta.star.cov <- solve( solve(theta.prior.cov)+H.max*solve(theta.cov) )
theta.star.mean <- theta.star.cov %*% (
	solve(theta.prior.cov)%*%theta.mean + 
	solve(theta.cov)%*%apply(theta,2,sum) )
theta.star <- rmvnorm(1, theta.star.mean, theta.star.cov)

theta.cov <- riwish( H.max + theta.cov.df,
	theta.cov.scale + t(theta-t(array(theta.star, c(3,H.max)))) %*% 
	(theta-t(array(theta.star, c(3,H.max)))) )



##################################################
### Sample sigma2

sig.cand <- rlnorm(1, log(sigma2), MH.sig.sd)
MH.prob <- dlnorm(sigma2, log(sig.cand), MH.sig.sd, log=T) - 
	dlnorm(sig.cand, log(sigma2), MH.sig.sd, log=T) + 
	log(dinvgamma(sig.cand, sigma2.shape, sigma2.scale)) - 
	log(dinvgamma(sigma2, sigma2.shape, sigma2.scale))
for( k in (1:ncol(cluster.ind))[apply(cluster.ind,2,sum)>0] ){
	MH.prob <- MH.prob - 
		.5*tr( (solve(sig.cand*Corr.matrix[[k]]+
		gamma2*(Corr.matrix[[k]]!=0)) - 
		solve(sigma2*Corr.matrix[[k]]+
		gamma2*(Corr.matrix[[k]]!=0))) %*%
		(logHCG.clust[[k]] - fhat.clust[[k]] ) %*% 
		t(logHCG.clust[[k]] - fhat.clust[[k]]) ) - 
		.5*as.numeric(determinant(sig.cand*Corr.matrix[[k]]+
		gamma2*(Corr.matrix[[k]]!=0), log=T)$mod) + 
		.5*as.numeric(determinant(sigma2*Corr.matrix[[k]]+
		gamma2*(Corr.matrix[[k]]!=0),log=T)$mod)
}
if( runif(1) < exp(min(MH.prob,0)) ){
		sigma2 <- sig.cand
		MH.counter[2] <- MH.counter[2]+1
}




##################################################
### Sample gamma2

gam.cand <- rlnorm(1, log(gamma2), MH.gam.sd)
MH.prob <- dlnorm(gamma2, log(gam.cand), MH.gam.sd, log=T) - 
	dlnorm(gam.cand, log(gamma2), MH.gam.sd, log=T) + 
	log(dinvgamma(gam.cand, gamma2.shape, gamma2.scale)) - 
	log(dinvgamma(gamma2, gamma2.shape, gamma2.scale))
for( k in (1:ncol(cluster.ind))[apply(cluster.ind,2,sum)>0] ){
	MH.prob <- MH.prob - 
		.5*tr( (solve(sigma2*Corr.matrix[[k]]+
		gam.cand*(Corr.matrix[[k]]!=0)) - 
		solve(sigma2*Corr.matrix[[k]]+
		gamma2*(Corr.matrix[[k]]!=0))) %*%
		(logHCG.clust[[k]] - fhat.clust[[k]] ) %*% 
		t(logHCG.clust[[k]] - fhat.clust[[k]]) ) - 
		.5*as.numeric(determinant(sigma2*Corr.matrix[[k]]+
		gam.cand*(Corr.matrix[[k]]!=0), log=T)$mod) + 
		.5*as.numeric(determinant(sigma2*Corr.matrix[[k]]+
		gamma2*(Corr.matrix[[k]]!=0),log=T)$mod)
}
if( runif(1) < exp(min(MH.prob,0)) ){
		gamma2 <- gam.cand
		MH.counter[3] <- MH.counter[3]+1
}




##################################################
### Sample rho

MH.prob <- 0
rho.cand <- rnorm(1, rho, MH.rho.sd)
if( rho.cand < 0 ){ rho.cand <- abs(rho.cand) }
if( rho.cand > 1 ){ rho.cand <- 1 - abs(rho.cand-1) }
if( rho.cand < 0 ){ rho.cand <- abs(rho.cand) }
if( rho.cand > 1 ){ rho.cand <- 1 - abs(rho.cand-1) }
Corr.temp <- Corr.matrix
for( k in (1:ncol(cluster.ind))[apply(cluster.ind,2,sum)>0] ){
	Corr.temp[[k]] <- rho.cand^( 
		abs(array(Times.clust[[k]], 
		c(length(Times.clust[[k]]),length(Times.clust[[k]])))-
		t(array(Times.clust[[k]], c(length(Times.clust[[k]]),
		length(Times.clust[[k]])))))/corr.unit) * (1*(array(ID.clust[[k]], 
		c(length(ID.clust[[k]]),length(ID.clust[[k]])))==
		t(array(ID.clust[[k]], c(length(ID.clust[[k]]),
		length(ID.clust[[k]]))))))
	MH.prob <- MH.prob - 
		.5*tr( (solve(sigma2*Corr.temp[[k]]+
		gamma2*(Corr.matrix[[k]]!=0)) - 
		solve(sigma2*Corr.matrix[[k]]+
		gamma2*(Corr.matrix[[k]]!=0))) %*%
		(logHCG.clust[[k]] - fhat.clust[[k]] ) %*% 
		t(logHCG.clust[[k]] - fhat.clust[[k]]) ) - 
		.5*as.numeric(determinant(sigma2*Corr.temp[[k]]+
		gamma2*(Corr.matrix[[k]]!=0),log=T)$mod) + 
		.5*as.numeric(determinant(sigma2*Corr.matrix[[k]]+
		gamma2*(Corr.matrix[[k]]!=0),log=T)$mod)
}
if( (runif(1) < exp(min(MH.prob,0))) & (rho.cand>0) & (rho.cand<1) ){
		rho <- rho.cand
		Corr.matrix <- Corr.temp
		MH.counter[4] <- MH.counter[4]+1
}



##################################################
### Keep variables

if( (floor((iter-burn.in)/keep.prop)==(iter-burn.in)/keep.prop) * 
	(iter > burn.in) ){
	q <- (iter - burn.in) %/% keep.prop

# Computes the predictive probability for each patient
for( i in 1:n.pat ){
	set <- (i==ID)
	loglik <- rep(NA,nrow(theta))
	loglik <- log(cluster.prob)
	Var.temp <- sigma2*rho^( abs(array(Times[set], 
		c(length(Times[set]),length(Times[set])))-
		t(array(Times[set], c(length(Times[set]),
		length(Times[set])))))/corr.unit) + gamma2
	for( k in 1:nrow(theta) ){
		loglik[k] <- loglik[k] + dmvnorm(logHCG[set],
			theta.function(Times[set], theta[k,]),
			Var.temp, log=T) 
	}
	loglik <- exp(loglik - max(loglik))
	loglik <- loglik/sum(loglik)
	test.prob[i] <- sum(phi*loglik)
}

	keep.clust.prob[,q] <- cluster.prob
	keep.alpha[q] <- alpha.DP
	keep.clust[,q] <- apply(cluster.ind,1,which.max)
	keep.phi[,q] <- phi
	keep.theta[,,q] <- theta
	keep.sigma2[q] <- sigma2
	keep.gamma2[q] <- gamma2
	keep.th.star[,q] <- theta.star
	keep.th.cov[,,q] <- theta.cov
	keep.rho[q] <- rho
	keep.MHcounter[,q] <- MH.counter
	keep.adpMCMC[,q] <- logMH.scale
	keep.test[,q] <- test.prob
	MH.counter <- 0*MH.counter
}
if( iter %in% seq(1000, burn.in+niter*keep.prop, by=1000) ){
	print(paste('Iteration Number',iter,' out of',burn.in+niter*keep.prop,
		';    ', Sys.time() )) 
	flush.console()
}

}
finish.time <- Sys.time() - time1
print(finish.time)
warnings()





save(MH.cov, burn.in, niter, keep.prop, H.max, MCMC.seed, file.out, 
	phi.alpha, phi.beta, theta.mean, theta.prior.cov, sigma2.shape, 
		sigma2.scale, alpha.shape, alpha.rate, theta.cov.df, 
		theta.cov.scale, gamma2.shape, gamma2.scale, alpha.sample, corr.unit,
	theta.function, tr, logMH.scale, MH.ss, alpha.star, MH.rho.sd, 
		MH.sig.sd, MH.gam.sd,
	data, Disease, n.pat, 
	keep.clust.prob, keep.alpha, keep.clust, keep.phi, keep.theta, 
		keep.sigma2, keep.gamma2, keep.th.star, keep.th.cov, keep.rho, 
		keep.MHcounter, keep.adpMCMC, keep.test, finish.time, 
	file=file.out)


