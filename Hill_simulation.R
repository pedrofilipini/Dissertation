######################################
##Pedro Henrique Filipini dos Santos##
###University of São Paulo - Brazil###
######################################


#######################################
##Simulation study based on Real Data##
#######################################

#Libraries
library(BART)
library(bcf)
library(doSNOW)
library(foreach)

niters <- 1000
pthin <- 500
thin <- 100
burn <- 50000
npost <- 1000
ncores <- 20

#loading file
load(file="sim.data")

#excluding the controlled part of the study
obs <- imp1[,-c(imp1$treat)]

#only continuos variables will be used
covs.cont.n=c("bw","b.head","preterm","birth.o","nnhealth","momage")

#Number of observations
n <- nrow(obs)

#Creating training matrix
x <- t(obs[,c(covs.cont.n)])

#Standardizing the continuous variables for generation of Y
xmat <- obs[,c(covs.cont.n)]

xmat[,covs.cont.n] <- as.data.frame(t((t(xmat[,covs.cont.n])
                                       -unlist(lapply(xmat[,covs.cont.n],mean)))
                                      /sqrt(unlist(lapply(xmat[covs.cont.n],var)))))
xmat <- as.matrix(xmat)

#Rule to generate the propensity score
q <- -1*(xmat[,1]>(xmat[,2])) + 1*(xmat[,1]<(xmat[,2]))

#Generating the true p-score
pi <- pnorm(q)

#Generating treatment variables
set.seed(99)
z <- rbinom(n,1,pi)

#Alpha its the true treatment effect
alpha =  0.5*(xmat[,3] > -3/4) + 0.25*(xmat[,3] > 0) + 0.25*(xmat[,3]>3/4)

#RMSE ITE function
rmseITE <- function(ite, x){
  r <- numeric(0)
  n <- length(ite)
  for(i in 1:n){
    if(((x[3,(i)]-mean(x[3,]))/sd(x[3,]))>(-3/4) && ((x[3,(i)]-mean(x[3,]))/sd(x[3,]))<0){
      r[i] <- ((0.5-ite[i])^2)
    }else if(((x[3,(i)]-mean(x[3,]))/sd(x[3,]))<(3/4) && ((x[3,(i)]-mean(x[3,]))/sd(x[3,]))>0){
      r[i] <- ((0.75-ite[i])^2)
    }else if(((x[3,(i)]-mean(x[3,]))/sd(x[3,]))>(3/4)){
      r[i] <- ((1-ite[i])^2)
    }else{
      r[i] <- ((0-ite[i])^2)
    }
  }
  sqrt(mean(r))
}

#RMSE ATE function
rmseATE <- function(ate, alpha){
  aux <- rowMeans(ate)
  sqrt(mean((aux-mean(alpha))^2))
}

#starting RMSE vectors
ite_rmse_vanilla <- numeric(0)
ite_rmse_oracle <- numeric(0)
ite_rmse_ps <- numeric(0)
ite_rmse_psglm <- numeric(0)
ite_rmse_rand <- numeric(0)
ite_rmse_oraclebcf <- numeric(0)
ite_rmse_bartbcf <- numeric(0)
ite_rmse_glmbcf <- numeric(0)
ite_rmse_randbcf <- numeric(0)

ate_rmse_vanilla <- numeric(0)
ate_rmse_oracle <- numeric(0)
ate_rmse_ps <- numeric(0)
ate_rmse_psglm <- numeric(0)
ate_rmse_rand <- numeric(0)
ate_rmse_oraclebcf <- numeric(0)
ate_rmse_bartbcf <- numeric(0)
ate_rmse_glmbcf <- numeric(0)
ate_rmse_randbcf <- numeric(0)

rmse_pbart <- numeric(0)
pehe_tt_vanilla <- numeric(0)
pehe_tt_oracle <- numeric(0)
pehe_tt_ps <- numeric(0)
pehe_tt_psglm <- numeric(0)
pehe_tt_rand <- numeric(0)
pehe_tt_oraclebcf <- numeric(0)
pehe_tt_bartbcf <- numeric(0)
pehe_tt_glmbcf <- numeric(0)
pehe_tt_randbcf <- numeric(0)

pehe_tc_vanilla <- numeric(0)
pehe_tc_oracle <- numeric(0)
pehe_tc_ps <- numeric(0)
pehe_tc_psglm <- numeric(0)
pehe_tc_rand <- numeric(0)
pehe_tc_oraclebcf <- numeric(0)
pehe_tc_bartbcf <- numeric(0)
pehe_tc_glmbcf <- numeric(0)
pehe_tc_randbcf <- numeric(0)

att_rmse_vanilla <- numeric(0)
att_rmse_oracle <- numeric(0)
att_rmse_ps <- numeric(0)
att_rmse_psglm <- numeric(0)
att_rmse_rand <- numeric(0)
att_rmse_oraclebcf <- numeric(0)
att_rmse_bartbcf <- numeric(0)
att_rmse_glmbcf <- numeric(0)
att_rmse_randbcf <- numeric(0)

atc_rmse_vanilla <- numeric(0)
atc_rmse_oracle <- numeric(0)
atc_rmse_ps <- numeric(0)
atc_rmse_psglm <- numeric(0)
atc_rmse_rand <- numeric(0)
atc_rmse_oraclebcf <- numeric(0)
atc_rmse_bartbcf <- numeric(0)
atc_rmse_glmbcf <- numeric(0)
atc_rmse_randbcf <- numeric(0)

cluster = makeCluster(ncores, type = "SOCK")
registerDoSNOW(cluster)

time <- Sys.time()
results <- foreach(i=1:niters) %dopar%{
  library("BART")
  library("bcf")
  #setting different seeds for different iterations
  seedaux <- (7565 + i*5)
  
  set.seed(seedaux)
  
  #Generate data coefficients
  beta = sample(c(4, 3, 2, 1, 0),ncol(xmat),
                replace=TRUE,prob=c(.5,.2,.15,.1,.05))
  
  #Multiplying x matrix by the coefficients
  yhat = xmat %*% beta
  
  #Generating Y
  y <- yhat + q + z*alpha + 0.5*rnorm(n)
  
  #pihat is the true p-score
  pihat <- pi
  
  #Test Matrix
  X = cbind(rbind(cbind(t(x),pihat),cbind(t(x),pihat)),c(rep(1,n),rep(0,n)))
  colnames(X)[ncol(X)]="z"
  
  #VanillaBART
  set.seed(seedaux)
  vanilla = wbart(cbind(t(x),z),y,X[,-(ncol(X)-1)],nskip = burn,ndpost = npost, keepevery = thin)
  
  #OracleBART
  set.seed(seedaux)
  oracle = wbart(cbind(t(x),pihat,z),y,X,nskip = burn,ndpost = npost, keepevery = thin)
  
  #pihat with pbart
  set.seed(seedaux)
  fitz = pbart(t(x),z,nskip = burn,ndpost = npost, keepevery = pthin)
  pihat2 = fitz$prob.train.mean #Could use median instead
  X2 = cbind(rbind(cbind(t(x),pihat2),cbind(t(x),pihat2)),c(rep(1,n),rep(0,n)))
  colnames(X2)[ncol(X2)]="z"
  
  #PSBART
  set.seed(seedaux)
  ps = wbart(cbind(t(x),pihat2,z), y, X2, nskip = burn,ndpost = npost, keepevery = thin)
  
  #Calculating estimated ATE
  ate_est_vanilla = vanilla$yhat.test[,1:n] - vanilla$yhat.test[,(n+1):(2*n)]
  ate_est_oracle = oracle$yhat.test[,1:n] - oracle$yhat.test[,(n+1):(2*n)]
  ate_est_ps = ps$yhat.test[,1:n] - ps$yhat.test[,(n+1):(2*n)]
  
  #Calculating estimated ITE
  ite_est_vanilla = colMeans(vanilla$yhat.test[,1:n] - vanilla$yhat.test[,(n+1):(2*n)])
  ite_est_oracle = colMeans(oracle$yhat.test[,1:n] - oracle$yhat.test[,(n+1):(2*n)])
  ite_est_ps = colMeans(ps$yhat.test[,1:n] - ps$yhat.test[,(n+1):(2*n)])
  
  
  #pihat estimated by GLM
  fitz = glm(z ~ t(x), family = binomial())
  pihat4 = fitz$fitted.values
  X4 = cbind(rbind(cbind(t(x),pihat4),cbind(t(x),pihat4)),c(rep(1,n),rep(0,n)))
  colnames(X4)[ncol(X4)]="z"
  
  #PSGLM
  set.seed(seedaux)
  ps3 = wbart(cbind(t(x),pihat4,z),y,X4,nskip = burn,ndpost = npost, keepevery = thin)
  
  #Estimated Treatment Effects
  ate_est_psglm = ps3$yhat.test[,1:n] - ps3$yhat.test[,(n+1):(2*n)]
  ite_est_psglm = colMeans(ps3$yhat.test[,1:n] - ps3$yhat.test[,(n+1):(2*n)])
  
  
  #pihat randomly generated
  pihat5 = runif(n)
  X5 = cbind(rbind(cbind(t(x),pihat5),cbind(t(x),pihat5)),c(rep(1,n),rep(0,n)))
  colnames(X5)[ncol(X5)]="z"
  
  #PSRAND
  set.seed(seedaux)
  ps4 = wbart(cbind(t(x),pihat5,z),y,X5,nskip = burn,ndpost = npost, keepevery = thin)
  
  #Estimated Treatment Effects
  ate_est_psrand = ps4$yhat.test[,1:n] - ps4$yhat.test[,(n+1):(2*n)]
  ite_est_psrand = colMeans(ps4$yhat.test[,1:n] - ps4$yhat.test[,(n+1):(2*n)])
  
  
  #BCF with BART (pihat2)
  set.seed(seedaux)
  fitbcf2 = bcf(y, z, t(x), t(x), pihat2, nburn = burn, nsim = npost, nthin = thin, include_pi = "both")
  
  #Estimated Treatment Effects
  ate_est_bartbcf = fitbcf2$tau
  ite_est_bartbcf = colMeans(fitbcf2$tau)
  
  
  #Oracle BCF (pihat)
  set.seed(seedaux)
  fitbcf3 = bcf(y, z, t(x), t(x), pi, nburn = burn, nsim = npost, nthin = thin, include_pi = "both")
  
  #Estimated Treatment Effects
  ate_est_oraclebcf = fitbcf3$tau
  ite_est_oraclebcf = colMeans(fitbcf3$tau)
  
  
  #BCF with GLM (pihat4)
  set.seed(seedaux)
  fitbcf4 = bcf(y, z, t(x), t(x), pihat4, nburn = burn, nsim = npost, nthin = thin, include_pi = "both")
  
  #Estimated Treatment Effects
  ate_est_glmbcf = fitbcf4$tau
  ite_est_glmbcf = colMeans(fitbcf4$tau)
  
  
  #BCF with random (pihat5)
  set.seed(seedaux)
  fitbcf5 = bcf(y, z, t(x), t(x), pihat5, nburn = burn, nsim = npost, nthin = thin, include_pi = "both")
  
  #Estimated Treatment Effects
  ate_est_randbcf = fitbcf5$tau
  ite_est_randbcf = colMeans(fitbcf5$tau)
  
  
  #Calculating RMSE for this iteration
  ite_rmse_vanilla[i] <- rmseITE(ite_est_vanilla,x)
  pehe_tt_vanilla[i] <- rmseITE(ite_est_vanilla[(z==1)],x[,(z==1)])
  pehe_tc_vanilla[i] <- rmseITE(ite_est_vanilla[(z==0)],x[,(z==0)])
  ate_rmse_vanilla[i] <- rmseATE(ate_est_vanilla,alpha)
  att_rmse_vanilla[i] <- rmseATE(ate_est_vanilla[,(z==1)],alpha[(z==1)])
  atc_rmse_vanilla[i] <- rmseATE(ate_est_vanilla[,(z==0)],alpha[(z==0)])
  
  ite_rmse_oracle[i] <- rmseITE(ite_est_oracle,x)
  pehe_tt_oracle[i] <- rmseITE(ite_est_oracle[(z==1)],x[,(z==1)])
  pehe_tc_oracle[i] <- rmseITE(ite_est_oracle[(z==0)],x[,(z==0)])
  ate_rmse_oracle[i] <- rmseATE(ate_est_oracle,alpha)
  att_rmse_oracle[i] <- rmseATE(ate_est_oracle[,(z==1)],alpha[(z==1)])
  atc_rmse_oracle[i] <- rmseATE(ate_est_oracle[,(z==0)],alpha[(z==0)])
  
  rmse_pbart[i] <- sqrt(mean((pihat-pihat2)^2))
  
  ite_rmse_ps[i] <- rmseITE(ite_est_ps,x)
  pehe_tt_ps[i] <- rmseITE(ite_est_ps[(z==1)],x[,(z==1)])
  pehe_tc_ps[i] <- rmseITE(ite_est_ps[(z==0)],x[,(z==0)])
  ate_rmse_ps[i] <- rmseATE(ate_est_ps,alpha)
  att_rmse_ps[i] <- rmseATE(ate_est_ps[,(z==1)],alpha[(z==1)])
  atc_rmse_ps[i] <- rmseATE(ate_est_ps[,(z==0)],alpha[(z==0)])
  
  ite_rmse_psglm[i] <- rmseITE(ite_est_psglm,x)
  pehe_tt_psglm[i] <- rmseITE(ite_est_psglm[(z==1)],x[,(z==1)])
  pehe_tc_psglm[i] <- rmseITE(ite_est_psglm[(z==0)],x[,(z==0)])
  ate_rmse_psglm[i] <- rmseATE(ate_est_psglm,alpha)
  att_rmse_psglm[i] <- rmseATE(ate_est_psglm[,(z==1)],alpha[(z==1)])
  atc_rmse_psglm[i] <- rmseATE(ate_est_psglm[,(z==0)],alpha[(z==0)])
  
  ite_rmse_rand[i] <- rmseITE(ite_est_psrand,x)
  pehe_tt_rand[i] <- rmseITE(ite_est_psrand[(z==1)],x[,(z==1)])
  pehe_tc_rand[i] <- rmseITE(ite_est_psrand[(z==0)],x[,(z==0)])
  ate_rmse_rand[i] <- rmseATE(ate_est_psrand,alpha)
  att_rmse_rand[i] <- rmseATE(ate_est_psrand[,(z==1)],alpha[(z==1)])
  atc_rmse_rand[i] <- rmseATE(ate_est_psrand[,(z==0)],alpha[(z==0)])
  
  ite_rmse_oraclebcf[i] <- rmseITE(ite_est_oraclebcf,x)
  pehe_tt_oraclebcf[i] <- rmseITE(ite_est_oraclebcf[(z==1)],x[,(z==1)])
  pehe_tc_oraclebcf[i] <- rmseITE(ite_est_oraclebcf[(z==0)],x[,(z==0)])
  ate_rmse_oraclebcf[i] <- rmseATE(ate_est_oraclebcf,alpha)
  att_rmse_oraclebcf[i] <- rmseATE(ate_est_oraclebcf[,(z==1)],alpha[(z==1)])
  atc_rmse_oraclebcf[i] <- rmseATE(ate_est_oraclebcf[,(z==0)],alpha[(z==0)])
  
  ite_rmse_bartbcf[i] <- rmseITE(ite_est_bartbcf,x)
  pehe_tt_bartbcf[i] <- rmseITE(ite_est_bartbcf[(z==1)],x[,(z==1)])
  pehe_tc_bartbcf[i] <- rmseITE(ite_est_bartbcf[(z==0)],x[,(z==0)])
  ate_rmse_bartbcf[i] <- rmseATE(ate_est_bartbcf,alpha)
  att_rmse_bartbcf[i] <- rmseATE(ate_est_bartbcf[,(z==1)],alpha[(z==1)])
  atc_rmse_bartbcf[i] <- rmseATE(ate_est_bartbcf[,(z==0)],alpha[(z==0)])
  
  ite_rmse_glmbcf[i] <- rmseITE(ite_est_glmbcf,x)
  pehe_tt_glmbcf[i] <- rmseITE(ite_est_glmbcf[(z==1)],x[,(z==1)])
  pehe_tc_glmbcf[i] <- rmseITE(ite_est_glmbcf[(z==0)],x[,(z==0)])
  ate_rmse_glmbcf[i] <- rmseATE(ate_est_glmbcf,alpha)
  att_rmse_glmbcf[i] <- rmseATE(ate_est_glmbcf[,(z==1)],alpha[(z==1)])
  atc_rmse_glmbcf[i] <- rmseATE(ate_est_glmbcf[,(z==0)],alpha[(z==0)])
  
  ite_rmse_randbcf[i] <- rmseITE(ite_est_randbcf,x)
  pehe_tt_randbcf[i] <- rmseITE(ite_est_randbcf[(z==1)],x[,(z==1)])
  pehe_tc_randbcf[i] <- rmseITE(ite_est_randbcf[(z==0)],x[,(z==0)])
  ate_rmse_randbcf[i] <- rmseATE(ate_est_randbcf,alpha)
  att_rmse_randbcf[i] <- rmseATE(ate_est_randbcf[,(z==1)],alpha[(z==1)])
  atc_rmse_randbcf[i] <- rmseATE(ate_est_randbcf[,(z==0)],alpha[(z==0)])
  
  ########
  #Saving#
  ########
  
  aux <- list(ite_rmse_vanilla[i], #rmseITE(ite_est_vanilla,x)
              pehe_tt_vanilla[i], #rmseITE(ite_est_vanilla[(z==1)],x[,(z==1)])
              pehe_tc_vanilla[i], #rmseITE(ite_est_vanilla[(z==0)],x[,(z==0)])
              ate_rmse_vanilla[i], #rmseATE(ate_est_vanilla,alpha)
              att_rmse_vanilla[i], #rmseATE(ate_est_vanilla[(z==1)],alpha[(z==1)])
              atc_rmse_vanilla[i], #rmseATE(ate_est_vanilla[(z==0)],alpha[(z==0)])
              
              ite_rmse_oracle[i], #rmseITE(ite_est_oracle,x)
              pehe_tt_oracle[i], #rmseITE(ite_est_oracle[(z==1)],x[,(z==1)])
              pehe_tc_oracle[i], #rmseITE(ite_est_oracle[(z==0)],x[,(z==0)])
              ate_rmse_oracle[i], #rmseATE(ate_est_oracle,alpha)
              att_rmse_oracle[i], #rmseATE(ate_est_oracle[(z==1)],alpha[(z==1)])
              atc_rmse_oracle[i], #rmseATE(ate_est_oracle[(z==0)],alpha[(z==0)])
              
              rmse_pbart[i], #sqrt(mean((pihat-pihat2)^2))
              
              ite_rmse_ps[i], #rmseITE(ite_est_ps,x)
              pehe_tt_ps[i], #rmseITE(ite_est_ps[(z==1)],x[,(z==1)])
              pehe_tc_ps[i], #rmseITE(ite_est_ps[(z==0)],x[,(z==0)])
              ate_rmse_ps[i], #rmseATE(ate_est_ps,alpha)
              att_rmse_ps[i], #rmseATE(ate_est_ps[(z==1)],alpha[(z==1)])
              atc_rmse_ps[i], #rmseATE(ate_est_ps[(z==0)],alpha[(z==0)])
              
              ite_rmse_psglm[i], #rmseITE(ite_est_psglm,x)
              pehe_tt_psglm[i], #rmseITE(ite_est_psglm[(z==1)],x[,(z==1)])
              pehe_tc_psglm[i], #rmseITE(ite_est_psglm[(z==0)],x[,(z==0)])
              ate_rmse_psglm[i], #rmseATE(ate_est_psglm,alpha)
              att_rmse_psglm[i], #rmseATE(ate_est_psglm[(z==1)],alpha[(z==1)])
              atc_rmse_psglm[i], #rmseATE(ate_est_psglm[(z==0)],alpha[(z==0)])
              
              ite_rmse_rand[i], #rmseITE(ite_est_psrand,x)
              pehe_tt_rand[i], #rmseITE(ite_est_psrand[(z==1)],x[,(z==1)])
              pehe_tc_rand[i], #rmseITE(ite_est_psrand[(z==0)],x[,(z==0)])
              ate_rmse_rand[i], #rmseATE(ate_est_psrand,alpha)
              att_rmse_rand[i], #rmseATE(ate_est_psrand[(z==1)],alpha[(z==1)])
              atc_rmse_rand[i], #rmseATE(ate_est_psrand[(z==0)],alpha[(z==0)])
              
              ite_rmse_oraclebcf[i], #rmseITE(ite_est_oraclebcf,x)
              pehe_tt_oraclebcf[i], #rmseITE(ite_est_oraclebcf[(z==1)],x[,(z==1)])
              pehe_tc_oraclebcf[i], #rmseITE(ite_est_oraclebcf[(z==0)],x[,(z==0)])
              ate_rmse_oraclebcf[i], #rmseATE(ate_est_oraclebcf,alpha)
              att_rmse_oraclebcf[i], #rmseATE(ate_est_oraclebcf[(z==1)],alpha[(z==1)])
              atc_rmse_oraclebcf[i], #rmseATE(ate_est_oraclebcf[(z==0)],alpha[(z==0)])
              
              ite_rmse_bartbcf[i], #rmseITE(ite_est_bartbcf,x)
              pehe_tt_bartbcf[i], #rmseITE(ite_est_bartbcf[(z==1)],x[,(z==1)])
              pehe_tc_bartbcf[i], #rmseITE(ite_est_bartbcf[(z==0)],x[,(z==0)])
              ate_rmse_bartbcf[i], #rmseATE(ate_est_bartbcf,alpha)
              att_rmse_bartbcf[i], #rmseATE(ate_est_bartbcf[(z==1)],alpha[(z==1)])
              atc_rmse_bartbcf[i], #rmseATE(ate_est_bartbcf[(z==0)],alpha[(z==0)])
              
              ite_rmse_glmbcf[i], #rmseITE(ite_est_glmbcf,x)
              pehe_tt_glmbcf[i], #rmseITE(ite_est_glmbcf[(z==1)],x[,(z==1)])
              pehe_tc_glmbcf[i], #rmseITE(ite_est_glmbcf[(z==0)],x[,(z==0)])
              ate_rmse_glmbcf[i], #rmseATE(ate_est_glmbcf,alpha)
              att_rmse_glmbcf[i], #rmseATE(ate_est_glmbcf[(z==1)],alpha[(z==1)])
              atc_rmse_glmbcf[i], #rmseATE(ate_est_glmbcf[(z==0)],alpha[(z==0)])
              
              ite_rmse_randbcf[i], #rmseITE(ite_est_randbcf,x)
              pehe_tt_randbcf[i], #rmseITE(ite_est_randbcf[(z==1)],x[,(z==1)])
              pehe_tc_randbcf[i], #rmseITE(ite_est_randbcf[(z==0)],x[,(z==0)])
              ate_rmse_randbcf[i], #rmseATE(ate_est_randbcf,alpha)
              att_rmse_randbcf[i], #rmseATE(ate_est_randbcf[(z==1)],alpha[(z==1)])
              atc_rmse_randbcf[i]) #rmseATE(ate_est_randbcf[(z==0)],alpha[(z==0)]))
  
  save(aux, file = paste0("results_sim_new_hill",i,".RData"))
  
  return(aux)
  
}

save(results, file = paste0("results_sim_new_hill.RData"))

stopCluster(cluster)
