##################################################
### Read in output from full BNP model to estimate clustering

optimal.method <- 1; optimal.method
	#1=Dahl; 2=HCL Average, h=.75; 3=HCL Avg, k=median number of clusters
	#4=HCL Ward, k=median; 6=Avg Silhouette, 7=Ward Silhouette
	#8=Avg Gamma, 9=Avg Tau, 10=Ward Gamma, 11=Ward Tau
	#12=VI, 13=Binder 



load('BNPclassification - full BNP model.RData')




# Find optimal partition
cluster.ind <- array(FALSE, c(n.pat, H.max))
k <- median(apply(keep.clust,2,function(s) length(unique(s)) ))
MCMC.seed <- 5111984
clust.opt <- sample(1:niter,1)
if( optimal.method %in% c(1:4,6:7) ){
	clust.pairs <- array(NA, c(n.pat, n.pat, niter))
	for( i in 1:n.pat ){ for( j in i:n.pat ){
		clust.pairs[i,j,] <- clust.pairs[j,i,] <- 
			(keep.clust[i,]==keep.clust[j,])
	}}
	clust.hat <- apply(clust.pairs,1:2,mean)
	if( optimal.method == 1){
		clust.opt <- which.min(apply( clust.pairs, 3, 
			function(s) sum((s-clust.hat)^2) ))
		cluster.ind <- sapply(1:H.max, function(s) 
			s==keep.clust[,clust.opt])
	}
	if( optimal.method == 2){
		temp <- cutree(hclust( as.dist(1-clust.hat), 'average' ), h=.75)
		for( j in 1:max(temp) ){
			cluster.ind[,j] <- (temp==j)
		}
	}
	if( optimal.method == 3){
		temp <- cutree(hclust( as.dist(1-clust.hat), 'average' ), k=k)
		for( j in 1:max(temp) ){
			cluster.ind[,j] <- (temp==j)
		}
	}
	if( optimal.method == 4){
		library(cluster)
		temp <- cutree(agnes( as.dist(1-clust.hat), diss=T, 'ward' ), k=k)
		for( j in 1:max(temp) ){
			cluster.ind[,j] <- (temp==j)
		}
	}
	if( optimal.method == 6 ){
		library(cluster); library(NbClust)
		diss_matrix <- as.dist(1-clust.hat)
		temp <- NbClust(diss=diss_matrix, distance=NULL, min.nc=2,
			max.nc=12, method='average', index='silhouette')$Best.partition
		for( j in 1:max(temp) ){
			cluster.ind[,j] <- (temp==j)
		}
	}
	if( optimal.method == 7 ){
		library(cluster); library(NbClust)
		diss_matrix <- as.dist(1-clust.hat)
		temp <- NbClust(diss=diss_matrix, distance=NULL, min.nc=2,
			max.nc=12, method='ward.D2', index='silhouette')$Best.partition
		for( j in 1:max(temp) ){
			cluster.ind[,j] <- (temp==j)
		}
	}
}
if( optimal.method %in% 8:11 ){
	clust.pairs <- array(NA, c(n.pat, n.pat, niter))
	for( i in 1:n.pat ){ for( j in i:n.pat ){
		clust.pairs[i,j,] <- clust.pairs[j,i,] <- 
			(keep.clust[i,]==keep.clust[j,])
	}}
	clust.hat <- apply(clust.pairs,1:2,mean)
	if( optimal.method %in% 8:9 ){ fit <- hclust( as.dist(1-clust.hat), 'average' ) }
	if( optimal.method %in% 10:11 ){ fit <- hclust( as.dist(1-clust.hat), 'ward.D2' ) }
	
	K.max <- 12
	s.plus <- rep(0,K.max)
	s.minus <- rep(0,K.max)
	t.keep <- rep(0,K.max)
	time <- Sys.time()
	for( k in 2:K.max ){
		labels <- cutree(fit, k=k)
		for( i1 in 1:(n.pat-1) ){ for( i2 in i1:n.pat ){
		for( j1 in 1:(n.pat-1) ){ for( j2 in j1:n.pat ){
		if( (labels[i1]==labels[i2])&(labels[j1]!=labels[j2]) ){
			s.plus[k] <- s.plus[k] + ( clust.hat[i1,i2] > clust.hat[j1,j2] )
			s.minus[k] <- s.minus[k] + ( clust.hat[i1,i2] < clust.hat[j1,j2] )
		}
		if( (labels[i1]==labels[i2])&(labels[j1]==labels[j2]) ){
			t.keep[k] <- t.keep[k] + 1
		}
		if( (labels[i1]!=labels[i2])&(labels[j1]!=labels[j2]) ){
			t.keep[k] <- t.keep[k] + 1
		}
		}}}}
	}
	print(Sys.time() - time)

	if( optimal.method %in% c(8,10) ){
		Gamma <- (s.plus - s.minus)/(s.plus + s.minus)
		temp <- cutree(fit, k=which.max(Gamma))
	}
	if( optimal.method %in% c(9,11) ){
		Nt <- .5*n.pat*(n.pat+1)
		Tau <- (s.plus - s.minus)/sqrt( ( Nt*(Nt-1)/2 - t.keep/2)*Nt*(Nt-1)/2 )
		temp <- cutree(fit, k=which.max(Tau))
	}
	for( j in 1:max(temp) ){
		cluster.ind[,j] <- (temp==j)
	}
}
if( optimal.method %in% c(12,13) ){
	# devtools::install_github('sarawade/mcclust.ext')
	library(mcclust.ext)
	skip.clust <- keep.clust[,]
	psm <- comp.psm( t(skip.clust) )
	if( optimal.method==13 ){
		temp <- minbinder.ext( psm, t(skip.clust), method='draws')
		temp <- minbinder.ext( psm, t(skip.clust), method='greedy',
			start.cl.greedy=temp$cl)
	}
	if( optimal.method==12 ){
		set.seed(314159)
		temp <- minVI( psm, t(skip.clust), method='avg')
		temp <- minVI( psm, t(skip.clust), method='greedy', start.cl=temp$cl)
	}
	if( max(temp$cl) > ncol(cluster.ind) ){
		cluster.ind <- array(FALSE, c( nrow(cluster.ind), max(temp$cl) ))
	}
	for( j in 1:max(temp$cl) ){
		cluster.ind[,j] <- (temp$cl==j)
	}
}
if( optimal.method %in% c(5) ){ clust.opt <- niter+1 }

# Set initial values
phi <- keep.phi[,clust.opt]
theta <- keep.theta[,,clust.opt]
sigma2 <- keep.sigma2[clust.opt]
gamma2 <- keep.gamma2[clust.opt]
theta.star <- keep.th.star[,clust.opt]
theta.cov <- keep.th.cov[,,clust.opt]
rho <- keep.rho[clust.opt]

# Remove empty clusters
empty <- (1:H.max)[apply(cluster.ind,2,sum)==0]
if( length(empty) > 0 ){
	cluster.ind <- cluster.ind[,-empty]
	phi <- phi[-empty]
	theta <- theta[-empty,]
	logMH.scale <- logMH.scale[-empty]
	H.max <- length(phi)
}

if( length(phi) < ncol(cluster.ind) ){
	H.max <- ncol(cluster.ind)
	theta <- rbind(theta, 
		t(array(theta[1,], c(3,H.max-length(phi)))))
	logMH.scale <- c(logMH.scale, rep(logMH.scale[1], H.max-length(phi) ))
	phi <- c(phi, rep(phi[1], H.max-length(phi) ))
}


rm( list=setdiff(ls(), c('cluster.ind', 'phi', 'theta', 'sigma2', 'gamma2',
	'theta.star', 'theta.cov', 'rho', 'MH.cov', 'logMH.scale', 'H.max',
	'optimal.method', 'file.out')) )
	



##################################################
### print out of estimated cluster

optimal.method
H.max
cluster.ind
apply(cluster.ind,2,sum)




##################################################
### MCMC configuration for stage 2 fit

burn.in <- 2000		#number of burn-in iterations
niter <- 2000		#number of iterations to store
keep.prop <- 20		#proportion of post-burn-in iterations to store
				#total number of iterations is (burn.in + niter*keep.prop
MCMC.seed <- 5111984
file.out <- substr( file.out, 1, nchar(file.out)-6)
file.out <- paste(file.out, '_method', optimal.method, '.RData', sep='')




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
MH.ss <- function(s){ (max(burn.in,iter))^(-.5) }
alpha.star <- 0.234

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
dir.prior <- 1


corr.unit <- 7
# rho corresponds to the correlation between observations corr.unit apart





##################################################
### Load Data

load('BNCLD_sim1_1.RData')
attach(data)




##################################################
### Initial MCMC


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


MH.counter <- rep(0,4)
	# MH.counter stores the number of times the MH steps accept the move
	# MH.counter[1] reports accepting theta per number of clusters
	# MH.counter[2] reports accepting sig2
	# MH.counter[3] reports accepting gam2
	# MH.counter[4] reports accepting rho
test.prob <- rep(0,n.pat)

keep.clust.prob <- array(0, c(H.max,niter))
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

cluster.prob <- rdirichlet(1, dir.prior + apply(cluster.ind,2,sum) )



##################################################
### Sample cluster membership



##################################################
### Sample phi

phi <- rbeta(H.max, 
	phi.alpha + sapply(1:H.max, function(s) sum(Disease[cluster.ind[,s]]) ),
	phi.beta + apply(cluster.ind,2,sum) - 
		sapply(1:H.max, function(s) sum(Disease[cluster.ind[,s]]) ))



##################################################
### Sample theta

for( k in (1:H.max) ){
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
		MH.counter[1] <- MH.counter[1] + 1/H.max
		fhat.clust[[k]] <- fhat.cand
	}
	logMH.scale[k] <- logMH.scale[k] + MH.ss(iter)*(
		exp(min(MH.prob,0)) - alpha.star )
}




##################################################
### Sample sigma2, theta.star, theta.cov

sig.cand <- rlnorm(1, log(sigma2), MH.sig.sd)
MH.prob <- dlnorm(sigma2, log(sig.cand), MH.sig.sd, log=T) - 
	dlnorm(sig.cand, log(sigma2), MH.sig.sd, log=T) + 
	log(dinvgamma(sig.cand, sigma2.shape, sigma2.scale)) - 
	log(dinvgamma(sigma2, sigma2.shape, sigma2.scale))
for( k in 1:H.max ){
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

theta.star.cov <- solve( solve(theta.prior.cov)+H.max*solve(theta.cov) )
theta.star.mean <- theta.star.cov %*% (
	solve(theta.prior.cov)%*%theta.mean + 
	solve(theta.cov)%*%apply(theta,2,sum) )
theta.star <- rmvnorm(1, theta.star.mean, theta.star.cov)

theta.cov <- riwish( H.max + theta.cov.df,
	theta.cov.scale + t(theta-t(array(theta.star, c(3,H.max)))) %*% 
	(theta-t(array(theta.star, c(3,H.max)))) )




##################################################
### Sample gamma2

gam.cand <- rlnorm(1, log(gamma2), MH.gam.sd)
MH.prob <- dlnorm(gamma2, log(gam.cand), MH.gam.sd, log=T) - 
	dlnorm(gam.cand, log(gamma2), MH.gam.sd, log=T) + 
	log(dinvgamma(gam.cand, gamma2.shape, gamma2.scale)) - 
	log(dinvgamma(gamma2, gamma2.shape, gamma2.scale))
for( k in 1:H.max ){
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
for( k in 1:H.max ){
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
	optimal.method, cluster.ind, phi.alpha, phi.beta, theta.mean, 
		theta.prior.cov, sigma2.shape, sigma2.scale, theta.cov.df, 
		theta.cov.scale, gamma2.shape, gamma2.scale, dir.prior, corr.unit,
	theta.function, tr, logMH.scale, MH.ss, alpha.star, MH.rho.sd, 
		MH.sig.sd, MH.gam.sd, 
	data, Disease, n.pat, 
	keep.clust.prob, keep.phi, keep.theta, 
		keep.sigma2, keep.gamma2, keep.th.star, keep.th.cov, keep.rho, 
		keep.MHcounter, keep.adpMCMC, keep.test, finish.time,  
	file=paste(file.out, '.RData', sep=''))

