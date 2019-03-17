# Thompson without covariates
DtchoiceThompson=function(Y,D, #outcomes and treatments thus far
                          k, #number of treatments
                          C, #vector of treatment cost
                          Nt){ # number of observations for period t
  
  SS=tapply(Y,D,sum, default=0) #vector of successes
  NN=tapply(Y,D,length, default=0) #vector of trials
  A=1+SS
  B=1+NN-SS
  
  Dt=rep(0,Nt)
  previousD=-Inf # auxiliary variable to avoid repeat assignments of same D
  
  for (i  in 1:Nt) {
    thetadraw=sapply(1:k, function(j) rbeta(1, A[j], B[j]))
    Dt[i]=which.max(thetadraw-C)
    previousD = Dt[i]
  }
  
  factor(Dt, levels=1:k)
}



DtchoiceThompsonProbabilities=function(Y,D, #outcomes and treatments thus far
                                       k, #number of treatments
                                       C=rep(0,k), #vector of treatment cost
                                       RR=10000){ #number of replication draws
  
  # Repeat Thompson sampling RR times
  DtRR=DtchoiceThompson(Y,D,k,C, RR)
  P_Dt=table(DtRR) / RR #average count for each treatment value and covariate value, replicated sample
  
  P_Dt=as.tibble(matrix(P_Dt, 1,k))
  colnames(P_Dt)=paste(1:k)
  P_Dt
}

################################
# Thompson with covariates
# helper functions for MCMC

# prior for hyperparameters governing distribution of theta across strata within each treatment arm 
log.prior = function(alpha,beta) {
  -2.5*log(max(alpha + beta,0))
}

# sampling theta vector from posterior, given hyperparameters, for a given treatment arm
draw.thetas = function(alpha,beta, NNd, SSd, nx) {
  rbeta(nx,alpha+SSd,beta+NNd-SSd)
}

# sampling from posterior for hyperparameters, given theta vector, for a given treatment arm 
draw.alpha = function(alpha,beta,theta,prop.sd,nx) {
  alpha.star = rnorm(1,alpha,prop.sd)
  if (alpha.star > .5) { #restricting the support of hyperparameters
    num = nx*(lgamma(alpha.star+beta) - lgamma(alpha.star)) +
      alpha.star*sum(log(theta)) + log.prior(alpha.star,beta)
  } else num=-Inf
  den = nx*(lgamma(alpha+beta)      - lgamma(alpha)) +
    alpha     *sum(log(theta)) + log.prior(alpha,beta)
  ifelse((log(runif(1))<=num - den),alpha.star,alpha)
}

draw.beta = function(alpha,beta,theta,prop.sd,nx) {
  beta.star = rnorm(1,beta,prop.sd)
  if (beta.star > .5) { #restricting the support of hyperparameters
    num = nx*(lgamma(alpha+beta.star) - lgamma(beta.star)) +
      beta.star*sum(log(1-theta)) + log.prior(alpha,beta.star)
  } else num=-Inf 
  den = nx*(lgamma(alpha+beta)      - lgamma(beta)) +
    beta     *sum(log(1-theta)) + log.prior(alpha,beta)
  ifelse((log(runif(1))<=num - den),beta.star,beta)
}

sample.theta.d = function(NNd, SSd, nx, 
                          RR=10000) { #sampling period
  B = 1000 #burn in period
  MM = B + RR
  # Metropolis tuning parameters
  alpha.prop.sd =  1
  beta.prop.sd =   1
  
  alpha = rep(0,MM)
  beta = alpha
  theta = matrix(0,MM,nx)
  
  # Initial values for the chain
  alpha[1] = 2
  beta[1] = 2
  theta[1,] = draw.thetas(alpha[1],beta[1], NNd, SSd,nx)
 
  # MCMC simulation
  for (m in 2:MM) {
    alpha[m] = draw.alpha(alpha[m-1],beta[m-1],theta[m-1,],alpha.prop.sd,nx)
    beta[m] = draw.beta(alpha[m],beta[m-1],theta[m-1,],beta.prop.sd,nx)
    theta[m,] = draw.thetas(alpha[m],beta[m], NNd, SSd,nx)
  }

  theta[(B+1):MM,]
}

DtchoiceMCMCProbabilities=function(Y,D,X, #outcomes, treatments, and covariates thus far
                                   k,nx, #number of treatments and number of strata
                                   C=rep(0,k), #vector of treatment cost
                                   RR=10000){ #number of replication draws
  
  SS=tapply(Y,list(D,X),sum, default=0) #matrix of successes
  NN=tapply(Y,list(D,X),length, default=0) #matrix of trials

  P_Dt=matrix(0,nx,k)
  thetadraws=list()

  for (d in 1:k) {
    thetadraws[[d]]=sample.theta.d(NN[d,], SS[d,], nx, RR)
  }
  
  
  
  Dt_x=factor(rep(0,RR), levels=1:k)
  thetaxdraw=rep(0,k)
  for (x in 1:nx) {
    for (r in 1:RR) {
      for (d in 1:k) thetaxdraw[d]=thetadraws[[d]][r,x]
      Dt_x[r]=which.max(thetaxdraw-C)
    }
    P_Dt[x,]=table(Dt_x)/RR
  }
  
  P_Dt=as_tibble(P_Dt)
  colnames(P_Dt)=paste(1:k)
  P_Dt
}

