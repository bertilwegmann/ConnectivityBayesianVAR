rm(list=ls())
### This script estimates models using HMC for the case of subject-specific 
### Sigma_s and the analytically derived posterior for common Sigma.

### load packages
library(R.matlab)
library(LaplacesDemon)
library(MCMCpack)
library(MASS)
library(mvnfast)
library(Matrix)
library(matrixcalc)
library(CholWishart)
library(PEIP)
library(base)
library(invgamma)

### set working directory to the source file location ###
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

NLag <- 1

for (Control_Ind in 1:2){ # select 1 = data for Controls, 2 = data for ASDS
  print("Control_Ind = ")
  print(Control_Ind)
    for (nLag in 1:NLag){
      print("nLag = ")
      print(nLag)

### Read data
DataCat <- getwd()
if (Control_Ind==1){
  yMatData <- readMat(paste(DataCat,"/CONTROLS_NYU_cpac_nofilt_noglobal_rois_cc200.mat",sep=""))
  ToFileName <- "Control"
} else {
  yMatData <- readMat(paste(DataCat,"/ASDS_NYU_cpac_nofilt_noglobal_rois_cc200.mat",sep=""))
  ToFileName <- "ASDS"
}
yMatOrig <- yMatData$allData
T <- dim(yMatOrig)[3] # number of time points

### Initial settings

# Set subjects
S <- 20 # number of subjects is set to 20 in all models
#Subj_Ind <- sort( sample(1:dim(yMatOrig)[1], S) )
if (Control_Ind==1){ # Control
  Subj_Ind <- c(2,9,10,14,22,28,30,34,35,36,46,56,61,65,70,76,79,85,87,90)
} else { # ASDS
  Subj_Ind <- c(2,7,11,16,20,22,25,26,33,39,43,47,48,54,56,57,62,63,66,73)
}

vmPFC <- c(5,22,55,79);vlPFC <- c(141,144);post_cingulate <- c(6,76);precuneus <- c(136,163)
DMN <- c(vmPFC,vlPFC,post_cingulate,precuneus) # Default Mode Network (DMN)
parietal <- c(7,82,132,156);temporal <- c(11,49);insula <- c(121,137);SMA <- c(123,161)
SMN <- c(parietal,temporal,insula,SMA) # Sensory Motor Network (SMN)
p <- 20 # the number of brain regions
p_Ind <- sort(c(DMN,SMN)) # Select which brain regions to analyze

q <-  p*nLag # number of covariates, no intercept.
qp <- q*p # number of VAR coefficients, qp = q*p.
CurrT <- T-nLag # CurrT number of observations T minus nLag.
N <-  S*CurrT # number of subjects S times number of actual time points T

### Data arrangements according to nLag
Y <- array( 0 , dim = c(CurrT,p,S) )
X <- array( 0 , dim = c(CurrT,q,S) )
for (ii in 1:S){
  SubjData <- matrix(t(yMatOrig[Subj_Ind[ii],p_Ind,]),T,p) # Data for subject ii
  
  i_start <- nLag + 1
  i_end <- T
  Y[,,ii] <- SubjData[i_start:i_end,]
  for (jj in 1:nLag){
    i_start <- i_start - 1
    i_end <- i_end - 1
    Low <- 1 + (jj-1)*p
    High <- jj*p
    X[,Low:High,ii] <- SubjData[i_start:i_end,]
  }
}

# Optimizing over kappa_s and lambda

# Shrinkage for larger lags
l2_vec <- matrix(1,q,1)
for (iLag in 1:nLag){
  Lag_Ind <- (1+(iLag-1)*p):(iLag*p)
  l2_vec[Lag_Ind] <- iLag**2
}
l2_mat <- diag(q)
diag(l2_mat) <- l2_vec

# mean of the subjects' region-specific variances
var_mat <- matrix(0,S,p)
for (ii in 1:S){
  var_mat[ii,] <- apply(Y[,,ii],2,var) # estimated variances for the p variables
}
s2_mat <- diag(q)
diag(s2_mat) <- rep(apply(var_mat,2,mean),q/p)

# Prior settings
nu_0 <- p+2
Psi_0 <- diag(p)
diag(Psi_0) <- nu_0*apply(var_mat,2,max) # Note that Psi_0 is nu_0*Psi_0 in the paper.

fr <- function(kappa_lambda){
  
  Case1 <- 0
  Case2 <- 0
  Case3 <- 0
  if (length(kappa_lambda)==1){
    Case1 <- 1
    #print(kappa_lambda)
  }
  if (length(kappa_lambda)==2){
    Case2 <- 1
    #print(kappa_lambda)
  }
  if (length(kappa_lambda)>2){
    Case3 <- 1
    #print(kappa_lambda[21])
  }
  
  Log_c_k <- 0
  P_D <- matrix(0,q,q)
  Sum_Qs_Es <- matrix(0,q,p)
  Psi_D <- matrix(0,p,p)
  
  # Matrices for the model with Sigma_s
  R_s <- array(matrix(0,p,p),dim = c(S,p,p))
  E_s <- array(matrix(0,q,p),dim = c(S,q,p))
  Q_s_inv_mat <- array(matrix(0,q,q),dim = c(S,q,q))
  P_s_mat <- array(matrix(0,q,q),dim = c(S,q,q))
  #
  
  for (ii in 1:S){
    
    Ys <- Y[,,ii]
    Xs <- X[,,ii]
    if (Case1 == 1){
      P_s <- kappa_lambda**2 * s2_mat * l2_mat
      P_s_mat[ii,,] <- P_s
    }
    if (Case2 == 1){
      P_s <- kappa_lambda[1]**2 * s2_mat * l2_mat
      P_s_mat[ii,,] <- P_s
    }
    if (Case3 == 1){
      P_s <- kappa_lambda[ii]**2 * s2_mat * l2_mat
      P_s_mat[ii,,] <- P_s
    }
    
    K1_s <- solve(P_s + t(Xs) %*% Xs)
    K2_s <- Xs %*% K1_s %*% P_s
    K3_s <- diag(q) - K1_s %*% P_s
    K4_s <- Ys - Xs%*%K1_s%*%t(Xs)%*%Ys
    K5_s <- t(Ys)%*%Xs%*%K1_s
    
    Q_s_inv <- t(K2_s)%*%K2_s + t(K3_s)%*%P_s%*%K3_s
    P_D <- P_D + Q_s_inv
    
    Qsinv_Es <- t(K2_s)%*%K4_s + t(K3_s)%*%t(K2_s)%*%Ys
    Sum_Qs_Es <- Sum_Qs_Es + Qsinv_Es
    
    Psi_D_s <- t(K4_s)%*%K4_s + K5_s%*%P_s%*%t(K5_s)
    Psi_D <- Psi_D + Psi_D_s
    
    # Matrices for the model with Sigma_s
    Temp <- solve(Q_s_inv) %*% Qsinv_Es
    E_s[ii,,] <- Temp
    R_s[ii,,] <- Psi_D_s - t(Temp)%*%Q_s_inv%*%Temp
    Q_s_inv_mat[ii,,] <- Q_s_inv
    #
    
    Det_Ps_1 <- determinant(P_s)
    Det_Ps_2 <- determinant(P_s + t(Xs) %*% Xs)
    Log_c_k <- Log_c_k + 0.5*p*Det_Ps_1$modulus - 0.5*p*Det_Ps_2$modulus
    
  }
  
  if (Case1 == 1){
    P0 <- kappa_lambda**2 * s2_mat * l2_mat
  }
  if (Case2 == 1){
    P0 <- kappa_lambda[2]**2 * s2_mat * l2_mat
  }
  if (Case3 == 1){
    P0 <- kappa_lambda[21]**2 * s2_mat * l2_mat
  }
  P_n <- P0 + P_D
  B_n <- solve(P_n) %*% Sum_Qs_Es
  Psi_n <- Psi_0 + Psi_D - t(B_n)%*%Sum_Qs_Es
  nu_n <- nu_0 + N
  
  Det_P0 <- determinant(P0)
  Det_P_n <- determinant(P_n)
  Det_Psi_n <- determinant(Psi_n)
  
  if (MakeOpt==1){
    Log_c_k + 0.5*p*Det_P0$modulus - 0.5*p*Det_P_n$modulus - 0.5*(S*CurrT-q+nu_0)*Det_Psi_n$modulus
  } else {
    ResMat <- list("Log_c_k"=Log_c_k,"P0"=P0,"P_D"=P_D,"P_n"=P_n,"B_n"=B_n,"Psi_n"=Psi_n,"nu_n"=nu_n,
                   "R_s"=R_s,"E_s"=E_s,"Q_s_inv_mat"=Q_s_inv_mat,"P_s_mat"=P_s_mat)
    return(ResMat)
  }
  
}

MakeOpt <- 1
Opt_kappa_lambda <- optimize(f=fr,interval=c(0,100),maximum=TRUE)
print(Opt_kappa_lambda)

Opt_kappa_lambda <- optim(matrix(Opt_kappa_lambda$maximum,2,1),fr,control=list(fnscale=-1,maxit=10000))
print(Opt_kappa_lambda)

Opt_kappa_lambda <- optim(c(matrix(Opt_kappa_lambda$par[1],20,1),Opt_kappa_lambda$par[2]),fr,control=list(fnscale=-1,maxit=10000))
print(Opt_kappa_lambda)

### Inference settings ###
Niter <- 700
Nwarmup <- 200
Nchains <- 2
Ndraws <- Nchains*(Niter - Nwarmup)

DropBoxCat <- getwd()

### Make inference for the model with common dense Sigma ###
MakeOpt <- 0
ResMat <- fr(Opt_kappa_lambda$par)

P_0 <- ResMat$P0  # Prior precision for B
P_s_mat <- ResMat$P_s_mat # Prior precision for B_s. To be used to calculate WAIC.
P_D <- ResMat$P_D # Data precision
P_n <- ResMat$P_n # Posterior precision
P_n_inv <- as.symmetric.matrix(solve(P_n))
B_n <- ResMat$B_n # Posterior mean
Psi_n <- ResMat$Psi_n
nu_n <- ResMat$nu_n

ModelType <- 3
# Sample beta and Sigma
Sigma <- array(matrix(0,p,p),dim = c(Ndraws,p,p))
B <- array(matrix(0,q,p),dim = c(Ndraws,q,p))
for (ii in 1:Ndraws){
  SigmaCurr <- matrix(0,p,p)
  diag(SigmaCurr) <- rinvgamma(n=p,shape=0.5*(nu_n-2*p),rate=0.5*diag(Psi_n))
  Sigma[ii,,] <- as.symmetric.matrix(SigmaCurr)
  B[ii,,] <- rmatrixnorm(B_n, P_n_inv, SigmaCurr)
}
saveRDS(object=list(B,Sigma,nu_n),file = paste(DropBoxCat,"/",ModelType,"_",ToFileName,"_nLag_",nLag,".RDS",sep="")) 
###############################

ModelType <- 2
# Sample beta, Sigma and Rho
Sigma <- array(matrix(0,p,p),dim = c(Ndraws,p,p))
B <- array(matrix(0,q,p),dim = c(Ndraws,q,p))
Rho <- array(matrix(0,p,p),dim = c(Ndraws,p,p)) # sample correlations Rho
for (ii in 1:Ndraws){
  SigmaCurr <- as.symmetric.matrix(riwish(nu_n,Psi_n))
  Sigma[ii,,] <- SigmaCurr
  B[ii,,] <- rmatrixnorm(B_n, P_n_inv, SigmaCurr)
  DiagMat <- diag(sqrt(1/diag(SigmaCurr)),p,p)
  Rho[ii,,] <- DiagMat %*% SigmaCurr %*% DiagMat
}
saveRDS(object=list(B,Sigma,Rho,nu_n),file = paste(DropBoxCat,"/",ModelType,"_",ToFileName,"_nLag_",nLag,".RDS",sep=""))
##########################################

R_s <- ResMat$R_s
E_s <- ResMat$E_s
Q_s_inv <- ResMat$Q_s_inv_mat

ModelType <- 1
if (ModelType==1){
  ### Estimate the model with Sigma_s using HMC in Stan ###
  ### Estimation is also over unknown quantity nu ###
  I_Mat <- diag(qp) # Identity matrix
  B_0_spec <- as.vector(matrix(0,1,qp)) # Needs to be only zeros
  Chol_Cov_B <- chol(solve(P_0))# cholesky decomposition of prior covariance matrix
  
  data.list <- list(p=p,q=q,S=S,qp=qp,T=CurrT,R_s=R_s,E_s=E_s,Q_s_inv=Q_s_inv,
                    B_0_spec=B_0_spec,Chol_Cov_B=Chol_Cov_B,nu_0=nu_0,Psi_0=Psi_0,I_Mat=I_Mat)
  
  # Initial values to HMC
  initf1 <- function() {
    list( B = colMeans(B), Sigma = as.symmetric.matrix(colMeans(Sigma)) )
  }
  
  Est_Mod <- 
    stan(file=paste(DropBoxCat,"/","New_Est_Model.stan",sep=""),
                  data=data.list,init=initf1,chains=Nchains,cores=Nchains,iter=Niter,warmup=Nwarmup,refresh=10)
  
  PostDraws <- extract.samples(Est_Mod)
  B <- PostDraws$B
  Sigma <- PostDraws$Sigma
  nu <- PostDraws$nu
  Rho <- array(matrix(0,p,p),dim = c(Ndraws,p,p))
  for (jj in 1:Ndraws){
    CurrSigma <- Sigma[jj,,]
    DiagMat <- diag(sqrt(1/diag(CurrSigma)),p,p)
    Rho[jj,,] <- DiagMat %*% CurrSigma %*% DiagMat
  }
  saveRDS(object=list(Est_Mod,B,Sigma,Rho,nu),
          file = paste(DropBoxCat,"/",ModelType,"_",ToFileName,"_nLag_",nLag,".RDS",sep=""))
}

}
}




