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
library(ICEbox)
#library(doSNOW)
#library(foreach)

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

#cluster = makeCluster(ncores, type = "SOCK")
#registerDoSNOW(cluster)

#time <- Sys.time()
#results <- foreach(i=1:niters) %dopar%{
#  library("BART")
#  library("bcf")
for(i in c(1)){
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
  fitzglm = glm(z ~ t(x), family = binomial())
  pihat4 = fitzglm$fitted.values
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
  
  #save(aux, file = paste0("results_sim_new_hill",i,".RData"))
  
  return(aux)
  
}

#save(results, file = paste0("results_sim_new_hill.RData"))

#stopCluster(cluster)



##########
##ICEbox##
##########

#Just plot the loaded file to see the ICEbox plot

pred.BCF.ATE <- function(object, newdata){
  pred.aux <- predict(object, newdata)-object$mu
  pred.aux <- as.matrix(pred.aux)/as.vector(object$vars)
  ite <- colMeans(pred.aux)
  return(ite)
}

pred.BCF.ATE.0025 <- function(object, newdata){
  pred.aux <- predict(object, newdata)-object$mu
  pred.aux <- as.matrix(pred.aux)/as.vector(object$vars)
  ite <- colMeans(pred.aux)
  ite <- apply(pred.aux, 2, quantile, 0.025)
  return(ite)
}

pred.BCF.ATE.0975 <- function(object, newdata){
  pred.aux <- predict(object, newdata)-object$mu
  pred.aux <- as.matrix(pred.aux)/as.vector(object$vars)
  ite <- colMeans(pred.aux)
  ite <- apply(pred.aux, 2, quantile, 0.975)
  return(ite)
}


#BARTBCF
fitbcf2aux <- vanilla
fitbcf2aux$treedraws <- fitbcf2$treedraws
fitbcf2aux$mu <- mean(y)
fitbcf2aux$vars <- (predict(fitbcf2aux, newdata = cbind(t(x),pihat2))[,1]-mean(y))/(fitbcf2$tau[,1])

Xice = cbind(t(x),pihat2)

ice.bcf.bart.x1 <- ice(fitbcf2aux, Xice, predictor = "bw", 
                       predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.x1, file="ice_bcf.bart.x1.RData")

ice.bcf.bart.x1.0025 <- ice(fitbcf2aux, Xice, predictor = "bw", 
                       predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.x1.0025, file="ice_bcf.bart.x1.0025.RData")

ice.bcf.bart.x1.0975 <- ice(fitbcf2aux, Xice, predictor = "bw", 
                            predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.x1.0975, file="ice_bcf.bart.x1.0975.RData")

###


ice.bcf.bart.x2 <- ice(fitbcf2aux, Xice, predictor = "b.head", 
                       predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.x2, file="ice_bcf.bart.x2.RData")

ice.bcf.bart.x2.0025 <- ice(fitbcf2aux, Xice, predictor = "b.head", 
                            predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.x2.0025, file="ice_bcf.bart.x2.0025.RData")

ice.bcf.bart.x2.0975 <- ice(fitbcf2aux, Xice, predictor = "b.head", 
                            predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.x2.0975, file="ice_bcf.bart.x2.0975.RData")

###

ice.bcf.bart.x3 <- ice(fitbcf2aux, Xice, predictor = "preterm", 
                       predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.x3, file="ice_bcf.bart.x3.RData")

ice.bcf.bart.x3.0025 <- ice(fitbcf2aux, Xice, predictor = "preterm", 
                            predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.x3.0025, file="ice_bcf.bart.x3.0025.RData")

ice.bcf.bart.x3.0975 <- ice(fitbcf2aux, Xice, predictor = "preterm", 
                            predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.x3.0975, file="ice_bcf.bart.x3.0975.RData")

###

ice.bcf.bart.x4 <- ice(fitbcf2aux, Xice, predictor = "birth.o", 
                       predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.x4, file="ice_bcf.bart.x4.RData")

ice.bcf.bart.x4.0025 <- ice(fitbcf2aux, Xice, predictor = "birth.o", 
                            predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.x4.0025, file="ice_bcf.bart.x4.0025.RData")

ice.bcf.bart.x4.0975 <- ice(fitbcf2aux, Xice, predictor = "birth.o", 
                            predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.x4.0975, file="ice_bcf.bart.x4.0975.RData")

###

ice.bcf.bart.x5 <- ice(fitbcf2aux, Xice, predictor = "nnhealth", 
                       predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.x5, file="ice_bcf.bart.x5.RData")

ice.bcf.bart.x5.0025 <- ice(fitbcf2aux, Xice, predictor = "nnhealth", 
                            predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.x5.0025, file="ice_bcf.bart.x5.0025.RData")

ice.bcf.bart.x5.0975 <- ice(fitbcf2aux, Xice, predictor = "nnhealth", 
                            predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.x5.0975, file="ice_bcf.bart.x5.0975.RData")

###

ice.bcf.bart.x6 <- ice(fitbcf2aux, Xice, predictor = "momage", 
                       predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.x6, file="ice_bcf.bart.x6.RData")

ice.bcf.bart.x6.0025 <- ice(fitbcf2aux, Xice, predictor = "momage", 
                            predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.x6.0025, file="ice_bcf.bart.x6.0025.RData")

ice.bcf.bart.x6.0975 <- ice(fitbcf2aux, Xice, predictor = "momage", 
                            predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.x6.0975, file="ice_bcf.bart.x6.0975.RData")

###

ice.bcf.bart.pihat <- ice(fitbcf2aux, Xice, predictor = "pihat2", 
                          predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.pihat, file="ice_bcf.bart.pihat.RData")

ice.bcf.bart.pihat.0025 <- ice(fitbcf2aux, Xice, predictor = "pihat2", 
                            predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.pihat.0025, file="ice_bcf.bart.pihat.0025.RData")

ice.bcf.bart.pihat.0975 <- ice(fitbcf2aux, Xice, predictor = "pihat2", 
                            predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.pihat.0975, file="ice_bcf.bart.pihat.0975.RData")

###

#OracleBCF
fitbcf3aux <- vanilla
fitbcf3aux$treedraws <- fitbcf3$treedraws
fitbcf3aux$mu <- mean(y)
fitbcf3aux$vars <- (predict(fitbcf3aux, newdata = cbind(t(x),pihat))[,1]-mean(y))/(fitbcf3$tau[,1])

Xice = cbind(t(x),pihat)

ice.bcf.oracle.x1 <- ice(fitbcf3aux, Xice, predictor = "bw", 
                         predictfcn = pred.BCF.ATE)
save(ice.bcf.oracle.x1, file="ice_bcf.oracle.x1.RData")

ice.bcf.oracle.x1.0025 <- ice(fitbcf3aux, Xice, predictor = "bw", 
                         predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.oracle.x1.0025, file="ice_bcf.oracle.x1.0025.RData")

ice.bcf.oracle.x1.0975 <- ice(fitbcf3aux, Xice, predictor = "bw", 
                              predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.oracle.x1.0975, file="ice_bcf.oracle.x1.0975.RData")

###

ice.bcf.oracle.x2 <- ice(fitbcf3aux, Xice, predictor = "b.head", 
                         predictfcn = pred.BCF.ATE)
save(ice.bcf.oracle.x2, file="ice_bcf.oracle.x2.RData")

ice.bcf.oracle.x2.0025 <- ice(fitbcf3aux, Xice, predictor = "b.head", 
                              predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.oracle.x2.0025, file="ice_bcf.oracle.x2.0025.RData")

ice.bcf.oracle.x2.0975 <- ice(fitbcf3aux, Xice, predictor = "b.head", 
                              predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.oracle.x2.0975, file="ice_bcf.oracle.x2.0975.RData")

###


ice.bcf.oracle.x3 <- ice(fitbcf3aux, Xice, predictor = "preterm", 
                         predictfcn = pred.BCF.ATE)
save(ice.bcf.oracle.x3, file="ice_bcf.oracle.x3.RData")

ice.bcf.oracle.x3.0025 <- ice(fitbcf3aux, Xice, predictor = "preterm", 
                              predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.oracle.x3.0025, file="ice_bcf.oracle.x3.0025.RData")

ice.bcf.oracle.x3.0975 <- ice(fitbcf3aux, Xice, predictor = "preterm", 
                              predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.oracle.x3.0975, file="ice_bcf.oracle.x3.0975.RData")

###

ice.bcf.oracle.x4 <- ice(fitbcf3aux, Xice, predictor = "birth.o", 
                         predictfcn = pred.BCF.ATE)
save(ice.bcf.oracle.x4, file="ice_bcf.oracle.x4.RData")

ice.bcf.oracle.x4.0025 <- ice(fitbcf3aux, Xice, predictor = "birth.o", 
                              predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.oracle.x4.0025, file="ice_bcf.oracle.x4.0025.RData")

ice.bcf.oracle.x4.0975 <- ice(fitbcf3aux, Xice, predictor = "birth.o", 
                              predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.oracle.x4.0975, file="ice_bcf.oracle.x4.0975.RData")

###

ice.bcf.oracle.x5 <- ice(fitbcf3aux, Xice, predictor = "nnhealth", 
                         predictfcn = pred.BCF.ATE)
save(ice.bcf.oracle.x5, file="ice_bcf.oracle.x5.RData")

ice.bcf.oracle.x5.0025 <- ice(fitbcf3aux, Xice, predictor = "nnhealth", 
                              predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.oracle.x5.0025, file="ice_bcf.oracle.x5.0025.RData")

ice.bcf.oracle.x5.0975 <- ice(fitbcf3aux, Xice, predictor = "nnhealth", 
                              predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.oracle.x5.0975, file="ice_bcf.oracle.x5.0975.RData")

###

ice.bcf.oracle.x6 <- ice(fitbcf3aux, Xice, predictor = "momage", 
                         predictfcn = pred.BCF.ATE)
save(ice.bcf.oracle.x6, file="ice_bcf.oracle.x6.RData")

ice.bcf.oracle.x6.0025 <- ice(fitbcf3aux, Xice, predictor = "momage", 
                              predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.oracle.x6.0025, file="ice_bcf.oracle.x6.0025.RData")

ice.bcf.oracle.x6.0975 <- ice(fitbcf3aux, Xice, predictor = "momage", 
                              predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.oracle.x6.0975, file="ice_bcf.oracle.x6.0975.RData")

###

ice.bcf.oracle.pihat <- ice(fitbcf3aux, Xice, predictor = "pihat", 
                            predictfcn = pred.BCF.ATE)
save(ice.bcf.oracle.pihat, file="ice_bcf.oracle.pihat.RData")

ice.bcf.oracle.pihat.0025 <- ice(fitbcf3aux, Xice, predictor = "pihat", 
                              predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.oracle.pihat.0025, file="ice_bcf.oracle.pihat.0025.RData")

ice.bcf.oracle.pihat.0975 <- ice(fitbcf3aux, Xice, predictor = "pihat", 
                              predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.oracle.pihat.0975, file="ice_bcf.oracle.pihat.0975.RData")

###

#GLMBCF
fitbcf4aux <- vanilla
fitbcf4aux$treedraws <- fitbcf4$treedraws
fitbcf4aux$mu <- mean(y)
fitbcf4aux$vars <- (predict(fitbcf4aux, newdata = cbind(t(x),pihat4))[,1]-mean(y))/(fitbcf4$tau[,1])

Xice = cbind(t(x),pihat4)

ice.bcf.glm.x1 <- ice(fitbcf4aux, Xice, predictor = "bw", 
                      predictfcn = pred.BCF.ATE)
save(ice.bcf.glm.x1, file="ice_bcf.glm.x1.RData")

ice.bcf.glm.x1.0025 <- ice(fitbcf4aux, Xice, predictor = "bw", 
                      predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.glm.x1.0025, file="ice_bcf.glm.x1.0025.RData")

ice.bcf.glm.x1.0975 <- ice(fitbcf4aux, Xice, predictor = "bw", 
                           predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.glm.x1.0975, file="ice_bcf.glm.x1.0975.RData")

###

ice.bcf.glm.x2 <- ice(fitbcf4aux, Xice, predictor = "b.head", 
                      predictfcn = pred.BCF.ATE)
save(ice.bcf.glm.x2, file="ice_bcf.glm.x2.RData")

ice.bcf.glm.x2.0025 <- ice(fitbcf4aux, Xice, predictor = "b.head", 
                           predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.glm.x2.0025, file="ice_bcf.glm.x2.0025.RData")

ice.bcf.glm.x2.0975 <- ice(fitbcf4aux, Xice, predictor = "b.head", 
                           predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.glm.x2.0975, file="ice_bcf.glm.x2.0975.RData")

###

ice.bcf.glm.x3 <- ice(fitbcf4aux, Xice, predictor = "preterm", 
                      predictfcn = pred.BCF.ATE)
save(ice.bcf.glm.x3, file="ice_bcf.glm.x3.RData")

ice.bcf.glm.x3.0025 <- ice(fitbcf4aux, Xice, predictor = "preterm", 
                           predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.glm.x3.0025, file="ice_bcf.glm.x3.0025.RData")

ice.bcf.glm.x3.0975 <- ice(fitbcf4aux, Xice, predictor = "preterm", 
                           predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.glm.x3.0975, file="ice_bcf.glm.x3.0975.RData")

###

ice.bcf.glm.x4 <- ice(fitbcf4aux, Xice, predictor = "birth.o", 
                      predictfcn = pred.BCF.ATE)
save(ice.bcf.glm.x4, file="ice_bcf.glm.x4.RData")

ice.bcf.glm.x4.0025 <- ice(fitbcf4aux, Xice, predictor = "birth.o", 
                           predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.glm.x4.0025, file="ice_bcf.glm.x4.0025.RData")

ice.bcf.glm.x4.0975 <- ice(fitbcf4aux, Xice, predictor = "birth.o", 
                           predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.glm.x4.0975, file="ice_bcf.glm.x4.0975.RData")

###

ice.bcf.glm.x5 <- ice(fitbcf4aux, Xice, predictor = "nnhealth", 
                      predictfcn = pred.BCF.ATE)
save(ice.bcf.glm.x5, file="ice_bcf.glm.x5.RData")

ice.bcf.glm.x5.0025 <- ice(fitbcf4aux, Xice, predictor = "nnhealth", 
                           predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.glm.x5.0025, file="ice_bcf.glm.x5.0025.RData")

ice.bcf.glm.x5.0975 <- ice(fitbcf4aux, Xice, predictor = "nnhealth", 
                           predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.glm.x5.0975, file="ice_bcf.glm.x5.0975.RData")

###

ice.bcf.glm.x6 <- ice(fitbcf4aux, Xice, predictor = "momage", 
                      predictfcn = pred.BCF.ATE)
save(ice.bcf.glm.x6, file="ice_bcf.glm.x6.RData")

ice.bcf.glm.x6.0025 <- ice(fitbcf4aux, Xice, predictor = "momage", 
                           predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.glm.x6.0025, file="ice_bcf.glm.x6.0025.RData")

ice.bcf.glm.x6.0975 <- ice(fitbcf4aux, Xice, predictor = "momage", 
                           predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.glm.x6.0975, file="ice_bcf.glm.x6.0975.RData")

###

ice.bcf.glm.pihat <- ice(fitbcf4aux, Xice, predictor = "pihat4", 
                         predictfcn = pred.BCF.ATE)
save(ice.bcf.glm.pihat, file="ice_bcf.glm.pihat.RData")

ice.bcf.glm.pihat.0025 <- ice(fitbcf4aux, Xice, predictor = "pihat4", 
                           predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.glm.pihat.0025, file="ice_bcf.glm.pihat.0025.RData")

ice.bcf.glm.pihat.0975 <- ice(fitbcf4aux, Xice, predictor = "pihat4", 
                           predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.glm.pihat.0975, file="ice_bcf.glm.pihat.0975.RData")

###

#RandBCF
fitbcf5aux <- vanilla
fitbcf5aux$treedraws <- fitbcf5$treedraws
fitbcf5aux$mu <- mean(y)
fitbcf5aux$vars <- (predict(fitbcf5aux, newdata = cbind(t(x),pihat5))[,1]-mean(y))/(fitbcf5$tau[,1])

Xice = cbind(t(x),pihat5)

ice.bcf.rand.x1 <- ice(fitbcf5aux, Xice, predictor = "bw", 
                       predictfcn = pred.BCF.ATE)
save(ice.bcf.rand.x1, file="ice_bcf.rand.x1.RData")

ice.bcf.rand.x1.0025 <- ice(fitbcf5aux, Xice, predictor = "bw", 
                       predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.rand.x1.0025, file="ice_bcf.rand.x1.0025.RData")

ice.bcf.rand.x1.0975 <- ice(fitbcf5aux, Xice, predictor = "bw", 
                            predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.rand.x1.0975, file="ice_bcf.rand.x1.0975.RData")

###

ice.bcf.rand.x2 <- ice(fitbcf5aux, Xice, predictor = "b.head", 
                       predictfcn = pred.BCF.ATE)
save(ice.bcf.rand.x2, file="ice_bcf.rand.x2.RData")

ice.bcf.rand.x2.0025 <- ice(fitbcf5aux, Xice, predictor = "b.head", 
                            predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.rand.x2.0025, file="ice_bcf.rand.x2.0025.RData")

ice.bcf.rand.x2.0975 <- ice(fitbcf5aux, Xice, predictor = "b.head", 
                            predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.rand.x2.0975, file="ice_bcf.rand.x2.0975.RData")

###


ice.bcf.rand.x3 <- ice(fitbcf5aux, Xice, predictor = "preterm", 
                       predictfcn = pred.BCF.ATE)
save(ice.bcf.rand.x3, file="ice_bcf.rand.x3.RData")

ice.bcf.rand.x3.0025 <- ice(fitbcf5aux, Xice, predictor = "preterm", 
                            predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.rand.x3.0025, file="ice_bcf.rand.x3.0025.RData")

ice.bcf.rand.x3.0975 <- ice(fitbcf5aux, Xice, predictor = "preterm", 
                            predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.rand.x3.0975, file="ice_bcf.rand.x3.0975.RData")

###

ice.bcf.rand.x4 <- ice(fitbcf5aux, Xice, predictor = "birth.o", 
                       predictfcn = pred.BCF.ATE)
save(ice.bcf.rand.x4, file="ice_bcf.rand.x4.RData")

ice.bcf.rand.x4.0025 <- ice(fitbcf5aux, Xice, predictor = "birth.o", 
                            predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.rand.x4.0025, file="ice_bcf.rand.x4.0025.RData")

ice.bcf.rand.x4.0975 <- ice(fitbcf5aux, Xice, predictor = "birth.o", 
                            predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.rand.x4.0975, file="ice_bcf.rand.x4.0975.RData")

###

ice.bcf.rand.x5 <- ice(fitbcf5aux, Xice, predictor = "nnhealth", 
                       predictfcn = pred.BCF.ATE)
save(ice.bcf.rand.x5, file="ice_bcf.rand.x5.RData")

ice.bcf.rand.x5.0025 <- ice(fitbcf5aux, Xice, predictor = "nnhealth", 
                            predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.rand.x5.0025, file="ice_bcf.rand.x5.0025.RData")

ice.bcf.rand.x5.0975 <- ice(fitbcf5aux, Xice, predictor = "nnhealth", 
                            predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.rand.x5.0975, file="ice_bcf.rand.x5.0975.RData")

###

ice.bcf.rand.x6 <- ice(fitbcf5aux, Xice, predictor = "momage", 
                       predictfcn = pred.BCF.ATE)
save(ice.bcf.rand.x6, file="ice_bcf.rand.x6.RData")

ice.bcf.rand.x6.0025 <- ice(fitbcf5aux, Xice, predictor = "momage", 
                            predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.rand.x6.0025, file="ice_bcf.rand.x6.0025.RData")

ice.bcf.rand.x6.0975 <- ice(fitbcf5aux, Xice, predictor = "momage", 
                            predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.rand.x6.0975, file="ice_bcf.rand.x6.0975.RData")

###


ice.bcf.rand.pihat <- ice(fitbcf5aux, Xice, predictor = "pihat5", 
                          predictfcn = pred.BCF.ATE)
save(ice.bcf.rand.pihat, file="ice_bcf.rand.pihat.RData")

ice.bcf.rand.pihat.0025 <- ice(fitbcf5aux, Xice, predictor = "pihat5", 
                            predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.rand.pihat.0025, file="ice_bcf.rand.pihat.0025.RData")

ice.bcf.rand.pihat.0975 <- ice(fitbcf5aux, Xice, predictor = "pihat5", 
                            predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.rand.pihat.0975, file="ice_bcf.rand.pihat.0975.RData")

###




#ATE prediction function
pred.BART.ATE <- function(object, newdata){
  X.aux <- cbind(rbind(newdata[,-(ncol(newdata))],newdata[,-(ncol(newdata))]),
                 c(rep(1,n),rep(0,n)))
  colnames(X.aux)[(ncol(X.aux))]="z"
  pred.aux <- predict(object, X.aux)
  size <- length(object$yhat.train.mean)
  ite <- colMeans(pred.aux[,1:size] - pred.aux[,(size+1):(2*size)])
}

#Credible Interval 2,5%
pred.BART.ATE.0025 <- function(object, newdata){
  X.aux <- cbind(rbind(newdata[,-(ncol(newdata))],newdata[,-(ncol(newdata))]),
                 c(rep(1,n),rep(0,n)))
  colnames(X.aux)[(ncol(X.aux))]="z"
  pred.aux <- predict(object, X.aux)
  size <- length(object$yhat.train.mean)
  ite <- apply((pred.aux[,1:size] - pred.aux[,(size+1):(2*size)]), 2, quantile, 0.025)
}

#Credible Interval 97,5%
pred.BART.ATE.0975 <- function(object, newdata){
  X.aux <- cbind(rbind(newdata[,-(ncol(newdata))],newdata[,-(ncol(newdata))]),
                 c(rep(1,n),rep(0,n)))
  colnames(X.aux)[(ncol(X.aux))]="z"
  pred.aux <- predict(object, X.aux)
  size <- length(object$yhat.train.mean)
  ite <- apply((pred.aux[,1:size] - pred.aux[,(size+1):(2*size)]), 2, quantile, 0.975)
}


#Vanilla
Xice = cbind(t(x),z)

ice.bart.vanilla.x1 <- ice(vanilla, Xice, predictor = "bw", 
                           predictfcn = pred.BART.ATE)
save(ice.bart.vanilla.x1, file="ice_vanilla_x1.RData")

ice.bart.vanilla.x1.0025 <- ice(vanilla, Xice, predictor = "bw", 
                                predictfcn = pred.BART.ATE.0025)
save(ice.bart.vanilla.x1.0025, file="ice_vanilla_x1.0025.RData")

ice.bart.vanilla.x1.0975 <- ice(vanilla, Xice, predictor = "bw", 
                                predictfcn = pred.BART.ATE.0975)
save(ice.bart.vanilla.x1.0975, file="ice_vanilla_x1.0975.RData")

####

ice.bart.vanilla.x2 <- ice(vanilla, Xice, predictor = "b.head", 
                           predictfcn = pred.BART.ATE)
save(ice.bart.vanilla.x2, file="ice_vanilla_x2.RData")

ice.bart.vanilla.x2.0025 <- ice(vanilla, Xice, predictor = "b.head", 
                                predictfcn = pred.BART.ATE.0025)
save(ice.bart.vanilla.x2.0025, file="ice_vanilla_x2.0025.RData")

ice.bart.vanilla.x2.0975 <- ice(vanilla, Xice, predictor = "b.head", 
                                predictfcn = pred.BART.ATE.0975)
save(ice.bart.vanilla.x2.0975, file="ice_vanilla_x2.0975.RData")

###

ice.bart.vanilla.x3 <- ice(vanilla, Xice, predictor = "preterm", 
                           predictfcn = pred.BART.ATE)
save(ice.bart.vanilla.x3, file="ice_vanilla_x3.RData")

ice.bart.vanilla.x3.0025 <- ice(vanilla, Xice, predictor = "preterm", 
                                predictfcn = pred.BART.ATE.0025)
save(ice.bart.vanilla.x3.0025, file="ice_vanilla_x3.0025.RData")

ice.bart.vanilla.x3.0975 <- ice(vanilla, Xice, predictor = "preterm", 
                                predictfcn = pred.BART.ATE.0975)
save(ice.bart.vanilla.x3.0975, file="ice_vanilla_x3.0975.RData")

###


ice.bart.vanilla.x4 <- ice(vanilla, Xice, predictor = "birth.o", 
                           predictfcn = pred.BART.ATE)
save(ice.bart.vanilla.x4, file="ice_vanilla_x4.RData")

ice.bart.vanilla.x4.0025 <- ice(vanilla, Xice, predictor = "birth.o", 
                                predictfcn = pred.BART.ATE.0025)
save(ice.bart.vanilla.x4.0025, file="ice_vanilla_x4.0025.RData")

ice.bart.vanilla.x4.0975 <- ice(vanilla, Xice, predictor = "birth.o", 
                                predictfcn = pred.BART.ATE.0975)
save(ice.bart.vanilla.x4.0975, file="ice_vanilla_x4.0975.RData")

###

ice.bart.vanilla.x5 <- ice(vanilla, Xice, predictor = "nnhealth", 
                           predictfcn = pred.BART.ATE)
save(ice.bart.vanilla.x5, file="ice_vanilla_x5.RData")

ice.bart.vanilla.x5.0025 <- ice(vanilla, Xice, predictor = "nnhealth", 
                                predictfcn = pred.BART.ATE.0025)
save(ice.bart.vanilla.x5.0025, file="ice_vanilla_x5.0025.RData")

ice.bart.vanilla.x5.0975 <- ice(vanilla, Xice, predictor = "nnhealth", 
                                predictfcn = pred.BART.ATE.0975)
save(ice.bart.vanilla.x5.0975, file="ice_vanilla_x5.0975.RData")

###

ice.bart.vanilla.x6 <- ice(vanilla, Xice, predictor = "momage", 
                           predictfcn = pred.BART.ATE)
save(ice.bart.vanilla.x6, file="ice_vanilla_x6.RData")

ice.bart.vanilla.x6.0025 <- ice(vanilla, Xice, predictor = "momage", 
                                predictfcn = pred.BART.ATE.0025)
save(ice.bart.vanilla.x6.0025, file="ice_vanilla_x6.0025.RData")

ice.bart.vanilla.x6.0975 <- ice(vanilla, Xice, predictor = "momage", 
                                predictfcn = pred.BART.ATE.0975)
save(ice.bart.vanilla.x6.0975, file="ice_vanilla_x6.0975.RData")

###

#Oracle
Xice = cbind(t(x),pihat,z)

ice.bart.oracle.x1 <- ice(oracle, Xice, predictor = "bw", 
                          predictfcn = pred.BART.ATE)
save(ice.bart.oracle.x1, file="ice_oracle_x1.RData")

ice.bart.oracle.x1.0025 <- ice(oracle, Xice, predictor = "bw", 
                               predictfcn = pred.BART.ATE.0025)
save(ice.bart.oracle.x1.0025, file="ice_oracle_x1.0025.RData")

ice.bart.oracle.x1.0975 <- ice(oracle, Xice, predictor = "bw", 
                               predictfcn = pred.BART.ATE.0975)
save(ice.bart.oracle.x1.0975, file="ice_oracle_x1.0975.RData")

####

ice.bart.oracle.x2 <- ice(oracle, Xice, predictor = "b.head", 
                          predictfcn = pred.BART.ATE)
save(ice.bart.oracle.x2, file="ice_oracle_x2.RData")

ice.bart.oracle.x2.0025 <- ice(oracle, Xice, predictor = "b.head", 
                               predictfcn = pred.BART.ATE.0025)
save(ice.bart.oracle.x2.0025, file="ice_oracle_x2.0025.RData")

ice.bart.oracle.x2.0975 <- ice(oracle, Xice, predictor = "b.head", 
                               predictfcn = pred.BART.ATE.0975)
save(ice.bart.oracle.x2.0975, file="ice_oracle_x2.0975.RData")

###

ice.bart.oracle.x3 <- ice(oracle, Xice, predictor = "preterm", 
                          predictfcn = pred.BART.ATE)
save(ice.bart.oracle.x3, file="ice_oracle_x3.RData")

ice.bart.oracle.x3.0025 <- ice(oracle, Xice, predictor = "preterm", 
                               predictfcn = pred.BART.ATE.0025)
save(ice.bart.oracle.x3.0025, file="ice_oracle_x3.0025.RData")

ice.bart.oracle.x3.0975 <- ice(oracle, Xice, predictor = "preterm", 
                               predictfcn = pred.BART.ATE.0975)
save(ice.bart.oracle.x3.0975, file="ice_oracle_x3.0975.RData")

###

ice.bart.oracle.x4 <- ice(oracle, Xice, predictor = "birth.o", 
                          predictfcn = pred.BART.ATE)
save(ice.bart.oracle.x4, file="ice_oracle_x4.RData")

ice.bart.oracle.x4.0025 <- ice(oracle, Xice, predictor = "birth.o", 
                               predictfcn = pred.BART.ATE.0025)
save(ice.bart.oracle.x4.0025, file="ice_oracle_x4.0025.RData")

ice.bart.oracle.x4.0975 <- ice(oracle, Xice, predictor = "birth.o", 
                               predictfcn = pred.BART.ATE.0975)
save(ice.bart.oracle.x4.0975, file="ice_oracle_x4.0975.RData")

###

ice.bart.oracle.x5 <- ice(oracle, Xice, predictor = "nnhealth", 
                          predictfcn = pred.BART.ATE)
save(ice.bart.oracle.x5, file="ice_oracle_x5.RData")

ice.bart.oracle.x5.0025 <- ice(oracle, Xice, predictor = "nnhealth", 
                               predictfcn = pred.BART.ATE.0025)
save(ice.bart.oracle.x5.0025, file="ice_oracle_x5.0025.RData")

ice.bart.oracle.x5.0975 <- ice(oracle, Xice, predictor = "nnhealth", 
                               predictfcn = pred.BART.ATE.0975)
save(ice.bart.oracle.x5.0975, file="ice_oracle_x5.0975.RData")

###

ice.bart.oracle.x6 <- ice(oracle, Xice, predictor = "momage", 
                          predictfcn = pred.BART.ATE)
save(ice.bart.oracle.x6, file="ice_oracle_x6.RData")

ice.bart.oracle.x6.0025 <- ice(oracle, Xice, predictor = "momage", 
                               predictfcn = pred.BART.ATE.0025)
save(ice.bart.oracle.x6.0025, file="ice_oracle_x6.0025.RData")

ice.bart.oracle.x6.0975 <- ice(oracle, Xice, predictor = "momage", 
                               predictfcn = pred.BART.ATE.0975)
save(ice.bart.oracle.x6.0975, file="ice_oracle_x6.0975.RData")

###

ice.bart.oracle.pihat <- ice(oracle, Xice, predictor = "pihat", 
                             predictfcn = pred.BART.ATE)
save(ice.bart.oracle.pihat, file="ice_oracle_pihat.RData")

ice.bart.oracle.pihat.0025 <- ice(oracle, Xice, predictor = "pihat", 
                                  predictfcn = pred.BART.ATE.0025)
save(ice.bart.oracle.pihat.0025, file="ice_oracle_pihat.0025.RData")

ice.bart.oracle.pihat.0975 <- ice(oracle, Xice, predictor = "pihat", 
                                  predictfcn = pred.BART.ATE.0975)
save(ice.bart.oracle.pihat.0975, file="ice_oracle_pihat.0975.RData")

###


#PSBART
Xice = cbind(t(x),pihat2,z)

ice.bart.ps.x1 <- ice(ps, Xice, predictor = "bw", 
                      predictfcn = pred.BART.ATE)
save(ice.bart.ps.x1, file="ice_ps_x1.RData")

ice.bart.ps.x1.0025 <- ice(ps, Xice, predictor = "bw", 
                           predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps.x1.0025, file="ice_ps_x1.0025.RData")

ice.bart.ps.x1.0975 <- ice(ps, Xice, predictor = "bw", 
                           predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps.x1.0975, file="ice_ps_x1.0975.RData")

####

ice.bart.ps.x2 <- ice(ps, Xice, predictor = "b.head", 
                      predictfcn = pred.BART.ATE)
save(ice.bart.ps.x2, file="ice_ps_x2.RData")

ice.bart.ps.x2.0025 <- ice(ps, Xice, predictor = "b.head", 
                           predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps.x2.0025, file="ice_ps_x2.0025.RData")

ice.bart.ps.x2.0975 <- ice(ps, Xice, predictor = "b.head", 
                           predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps.x2.0975, file="ice_ps_x2.0975.RData")

###

ice.bart.ps.x3 <- ice(ps, Xice, predictor = "preterm", 
                      predictfcn = pred.BART.ATE)
save(ice.bart.ps.x3, file="ice_ps_x3.RData")

ice.bart.ps.x3.0025 <- ice(ps, Xice, predictor = "preterm", 
                           predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps.x3.0025, file="ice_ps_x3.0025.RData")

ice.bart.ps.x3.0975 <- ice(ps, Xice, predictor = "preterm", 
                           predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps.x3.0975, file="ice_ps_x3.0975.RData")

###

ice.bart.ps.x4 <- ice(ps, Xice, predictor = "birth.o", 
                      predictfcn = pred.BART.ATE)
save(ice.bart.ps.x4, file="ice_ps_x4.RData")

ice.bart.ps.x4.0025 <- ice(ps, Xice, predictor = "birth.o", 
                           predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps.x4.0025, file="ice_ps_x4.0025.RData")

ice.bart.ps.x4.0975 <- ice(ps, Xice, predictor = "birth.o", 
                           predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps.x4.0975, file="ice_ps_x4.0975.RData")

###

ice.bart.ps.x5 <- ice(ps, Xice, predictor = "nnhealth", 
                      predictfcn = pred.BART.ATE)
save(ice.bart.ps.x5, file="ice_ps_x5.RData")

ice.bart.ps.x5.0025 <- ice(ps, Xice, predictor = "nnhealth", 
                           predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps.x5.0025, file="ice_ps_x5.0025.RData")

ice.bart.ps.x5.0975 <- ice(ps, Xice, predictor = "nnhealth", 
                           predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps.x5.0975, file="ice_ps_x5.0975.RData")

###

ice.bart.ps.x6 <- ice(ps, Xice, predictor = "momage", 
                      predictfcn = pred.BART.ATE)
save(ice.bart.ps.x6, file="ice_ps_x6.RData")

ice.bart.ps.x6.0025 <- ice(ps, Xice, predictor = "momage", 
                           predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps.x6.0025, file="ice_ps_x6.0025.RData")

ice.bart.ps.x6.0975 <- ice(ps, Xice, predictor = "momage", 
                           predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps.x6.0975, file="ice_ps_x6.0975.RData")

###

ice.bart.ps.pihat <- ice(ps, Xice, predictor = "pihat2", 
                         predictfcn = pred.BART.ATE)
save(ice.bart.ps.pihat, file="ice_ps_pihat.RData")

ice.bart.ps.pihat.0025 <- ice(ps, Xice, predictor = "pihat2", 
                              predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps.pihat.0025, file="ice_ps_pihat.0025.RData")

ice.bart.ps.pihat.0975 <- ice(ps, Xice, predictor = "pihat2", 
                              predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps.pihat.0975, file="ice_ps_pihat.0975.RData")

###



#GLM-BART
Xice = cbind(t(x),pihat4,z)

ice.bart.ps3.x1 <- ice(ps3, Xice, predictor = "bw", 
                       predictfcn = pred.BART.ATE)
save(ice.bart.ps3.x1, file="ice_ps3_x1.RData")

ice.bart.ps3.x1.0025 <- ice(ps3, Xice, predictor = "bw", 
                            predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps3.x1.0025, file="ice_ps3_x1.0025.RData")

ice.bart.ps3.x1.0975 <- ice(ps3, Xice, predictor = "bw", 
                            predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps3.x1.0975, file="ice_ps3_x1.0975.RData")

####

ice.bart.ps3.x2 <- ice(ps3, Xice, predictor = "b.head", 
                       predictfcn = pred.BART.ATE)
save(ice.bart.ps3.x2, file="ice_ps3_x2.RData")

ice.bart.ps3.x2.0025 <- ice(ps3, Xice, predictor = "b.head", 
                            predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps3.x2.0025, file="ice_ps3_x2.0025.RData")

ice.bart.ps3.x2.0975 <- ice(ps3, Xice, predictor = "b.head", 
                            predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps3.x2.0975, file="ice_ps3_x2.0975.RData")

###

ice.bart.ps3.x3 <- ice(ps3, Xice, predictor = "preterm", 
                       predictfcn = pred.BART.ATE)
save(ice.bart.ps3.x3, file="ice_ps3_x3.RData")

ice.bart.ps3.x3.0025 <- ice(ps3, Xice, predictor = "preterm", 
                            predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps3.x3.0025, file="ice_ps3_x3.0025.RData")

ice.bart.ps3.x3.0975 <- ice(ps3, Xice, predictor = "preterm", 
                            predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps3.x3.0975, file="ice_ps3_x3.0975.RData")

###

ice.bart.ps3.x4 <- ice(ps3, Xice, predictor = "birth.o", 
                       predictfcn = pred.BART.ATE)
save(ice.bart.ps3.x4, file="ice_ps3_x4.RData")

ice.bart.ps3.x4.0025 <- ice(ps3, Xice, predictor = "birth.o", 
                            predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps3.x4.0025, file="ice_ps3_x4.0025.RData")

ice.bart.ps3.x4.0975 <- ice(ps3, Xice, predictor = "birth.o", 
                            predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps3.x4.0975, file="ice_ps3_x4.0975.RData")

###

ice.bart.ps3.x5 <- ice(ps3, Xice, predictor = "nnhealth", 
                       predictfcn = pred.BART.ATE)
save(ice.bart.ps3.x5, file="ice_ps3_x5.RData")

ice.bart.ps3.x5.0025 <- ice(ps3, Xice, predictor = "nnhealth", 
                            predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps3.x5.0025, file="ice_ps3_x5.0025.RData")

ice.bart.ps3.x5.0975 <- ice(ps3, Xice, predictor = "nnhealth", 
                            predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps3.x5.0975, file="ice_ps3_x5.0975.RData")

###

ice.bart.ps3.x6 <- ice(ps3, Xice, predictor = "momage", 
                       predictfcn = pred.BART.ATE)
save(ice.bart.ps3.x6, file="ice_ps3_x6.RData")

ice.bart.ps3.x6.0025 <- ice(ps3, Xice, predictor = "momage", 
                            predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps3.x6.0025, file="ice_ps3_x6.0025.RData")

ice.bart.ps3.x6.0975 <- ice(ps3, Xice, predictor = "momage", 
                            predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps3.x6.0975, file="ice_ps3_x6.0975.RData")

###

ice.bart.ps3.pihat <- ice(ps3, Xice, predictor = "pihat4", 
                          predictfcn = pred.BART.ATE)
save(ice.bart.ps3.pihat, file="ice_ps3_pihat.RData")

ice.bart.ps3.pihat.0025 <- ice(ps3, Xice, predictor = "pihat4", 
                               predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps3.pihat.0025, file="ice_ps3_pihat.0025.RData")

ice.bart.ps3.pihat.0975 <- ice(ps3, Xice, predictor = "pihat4", 
                               predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps3.pihat.0975, file="ice_ps3_pihat.0975.RData")

###




#Random-BART
Xice = cbind(t(x),pihat5,z)

ice.bart.ps4.x1 <- ice(ps4, Xice, predictor = "bw", 
                       predictfcn = pred.BART.ATE)
save(ice.bart.ps4.x1, file="ice_ps4_x1.RData")

ice.bart.ps4.x1.0025 <- ice(ps4, Xice, predictor = "bw", 
                            predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps4.x1.0025, file="ice_ps4_x1.0025.RData")

ice.bart.ps4.x1.0975 <- ice(ps4, Xice, predictor = "bw", 
                            predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps4.x1.0975, file="ice_ps4_x1.0975.RData")

####

ice.bart.ps4.x2 <- ice(ps4, Xice, predictor = "b.head", 
                       predictfcn = pred.BART.ATE)
save(ice.bart.ps4.x2, file="ice_ps4_x2.RData")

ice.bart.ps4.x2.0025 <- ice(ps4, Xice, predictor = "b.head", 
                            predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps4.x2.0025, file="ice_ps4_x2.0025.RData")

ice.bart.ps4.x2.0975 <- ice(ps4, Xice, predictor = "b.head", 
                            predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps4.x2.0975, file="ice_ps4_x2.0975.RData")

###

ice.bart.ps4.x3 <- ice(ps4, Xice, predictor = "preterm", 
                       predictfcn = pred.BART.ATE)
save(ice.bart.ps4.x3, file="ice_ps4_x3.RData")

ice.bart.ps4.x3.0025 <- ice(ps4, Xice, predictor = "preterm", 
                            predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps4.x3.0025, file="ice_ps4_x3.0025.RData")

ice.bart.ps4.x3.0975 <- ice(ps4, Xice, predictor = "preterm", 
                            predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps4.x3.0975, file="ice_ps4_x3.0975.RData")

###

ice.bart.ps4.x4 <- ice(ps4, Xice, predictor = "birth.o", 
                       predictfcn = pred.BART.ATE)
save(ice.bart.ps4.x4, file="ice_ps4_x4.RData")

ice.bart.ps4.x4.0025 <- ice(ps4, Xice, predictor = "birth.o", 
                            predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps4.x4.0025, file="ice_ps4_x4.0025.RData")

ice.bart.ps4.x4.0975 <- ice(ps4, Xice, predictor = "birth.o", 
                            predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps4.x4.0975, file="ice_ps4_x4.0975.RData")

###

ice.bart.ps4.x5 <- ice(ps4, Xice, predictor = "nnhealth", 
                       predictfcn = pred.BART.ATE)
save(ice.bart.ps4.x5, file="ice_ps4_x5.RData")

ice.bart.ps4.x5.0025 <- ice(ps4, Xice, predictor = "nnhealth", 
                            predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps4.x5.0025, file="ice_ps4_x5.0025.RData")

ice.bart.ps4.x5.0975 <- ice(ps4, Xice, predictor = "nnhealth", 
                            predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps4.x5.0975, file="ice_ps4_x5.0975.RData")

###

ice.bart.ps4.x6 <- ice(ps4, Xice, predictor = "momage", 
                       predictfcn = pred.BART.ATE)
save(ice.bart.ps4.x6, file="ice_ps4_x6.RData")

ice.bart.ps4.x6.0025 <- ice(ps4, Xice, predictor = "momage", 
                            predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps4.x6.0025, file="ice_ps4_x6.0025.RData")

ice.bart.ps4.x6.0975 <- ice(ps4, Xice, predictor = "momage", 
                            predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps4.x6.0975, file="ice_ps4_x6.0975.RData")

###

ice.bart.ps4.pihat <- ice(ps4, Xice, predictor = "pihat5", 
                          predictfcn = pred.BART.ATE)
save(ice.bart.ps4.pihat, file="ice_ps4_pihat.RData")

ice.bart.ps4.pihat.0025 <- ice(ps4, Xice, predictor = "pihat5", 
                               predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps4.pihat.0025, file="ice_ps4_pihat.0025.RData")

ice.bart.ps4.pihat.0975 <- ice(ps4, Xice, predictor = "pihat5", 
                               predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps4.pihat.0975, file="ice_ps4_pihat.0975.RData")

###

#############################################################



#######################
#######Graphics########
#######################

ATE <- data.frame(rowMeans(ate_est_oracle),  rowMeans(ate_est_oraclebcf),
                  rowMeans(ate_est_vanilla), 
                  rowMeans(ate_est_ps), rowMeans(ate_est_bartbcf),
                  rowMeans(ate_est_psglm), rowMeans(ate_est_glmbcf), 
                  rowMeans(ate_est_psrand), rowMeans(ate_est_randbcf))
png("hill01.png", width = 7, height = 7, units = 'in', res = 100)
boxplot(ATE, las = 2,
        ylab = "CATE",
        names = c("Oracle", "Oracle-BCF",
                  "Vanilla", 
                  "PS-BART", "BART-BCF", 
                  "GLM-BART", "GLM-BCF",
                  "Rand-BART", "Rand-BCF"), 
        cex = 0.5,
        cex.lab = 0.8,
        cex.axis = 0.7)
#True ATE
abline(h = mean(alpha), col = "red", lwd = 2)
abline(v = 2.5, col = "gray10", lwd = 2, lty = 2)
dev.off()

#################################
##Assuring Convergence of pbart##
#################################
p <- 6
i <- floor(seq(1, n, length.out=20))

auto.corr <- acf(fitz$yhat.train[ , i], plot=FALSE)
max.lag <- max(auto.corr$lag[ , 1, 1])

j <- seq(-0.5, 0.4, length.out=10)
png("hill02.png", width = 7, height = 7, units = 'in', res = 100)
for(h in 1:10) {
  if(h==1)
    plot(1:max.lag+j[h], auto.corr$acf[1+(1:max.lag), h, h],
         type='h', xlim=c(0, max.lag+1), ylim=c(-1, 1),
         ylab='ACF', xlab='Lag')
  else
    lines(1:max.lag+j[h], auto.corr$acf[1+(1:max.lag), h, h],
          type='h', col=h)
}
dev.off()

png("hill03.png", width = 7, height = 7, units = 'in', res = 100)
for(j in 1:10) {
  if(j==1)
    plot(pnorm(fitz$yhat.train[ , i[j]]),
         type='l', ylim=c(0, 1),
         sub=paste0('N:', n, ', p:', p, ', thin:', pthin),
         ylab=expression(Phi(f(x))), xlab='m')
  else
    lines(pnorm(fitz$yhat.train[ , i[j]]),
          type='l', col=j)
}
dev.off()

geweke <- gewekediag(fitz$yhat.train)

j <- -80
png("hill04.png", width = 7, height = 7, units = 'in', res = 100)
plot(geweke$z, pch='.', cex=2, ylab='z', xlab='i',
     sub=paste0('N:', n, ', p:', p, ', thin:', pthin),
     xlim=c(j, n), ylim=c(-5, 5))
lines(1:n, rep(-1.96, n), type='l', col=6)
lines(1:n, rep(+1.96, n), type='l', col=6)
lines(1:n, rep(-2.576, n), type='l', col=5)
lines(1:n, rep(+2.576, n), type='l', col=5)
lines(1:n, rep(-3.291, n), type='l', col=4)
lines(1:n, rep(+3.291, n), type='l', col=4)
lines(1:n, rep(-3.891, n), type='l', col=3)
lines(1:n, rep(+3.891, n), type='l', col=3)
lines(1:n, rep(-4.417, n), type='l', col=2)
lines(1:n, rep(+4.417, n), type='l', col=2)
text(c(1, 1), c(-1.96, 1.96), pos=2, cex=0.6, labels='0.95')
text(c(1, 1), c(-2.576, 2.576), pos=2, cex=0.6, labels='0.99')
text(c(1, 1), c(-3.291, 3.291), pos=2, cex=0.6, labels='0.999')
text(c(1, 1), c(-3.891, 3.891), pos=2, cex=0.6, labels='0.9999')
text(c(1, 1), c(-4.417, 4.417), pos=2, cex=0.6, labels='0.99999')
dev.off()

length(geweke$z[geweke$z>-1.96 & geweke$z<1.96])/length(geweke$z)
length(geweke$z[geweke$z>-2.576 & geweke$z<2.576])/length(geweke$z)
length(geweke$z[geweke$z>-3.291 & geweke$z<3.291])/length(geweke$z)
length(geweke$z[geweke$z>-3.891 & geweke$z<3.891])/length(geweke$z)
length(geweke$z[geweke$z>-4.417 & geweke$z<4.417])/length(geweke$z)

########################
########################

png("hillvanilla_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(vanilla$sigma[seq(50000,150000, by = 100)], type = "l",
     ylab = expression(sigma),
     xlab = "Draw",
     ylim= c(0.45, 0.65))
dev.off()
png("hillvanilla_2.png", width = 7, height = 7, units = 'in', res = 100)
acf(vanilla$sigma[seq(50000,150000, by = 100)], main = "")
dev.off()

png("hilloracle_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(oracle$sigma[seq(50000,150000, by = 100)], type = "l",
     ylab = expression(sigma),
     xlab = "Draw",
     ylim= c(0.45, 0.65))
dev.off()
png("hilloracle_2.png", width = 7, height = 7, units = 'in', res = 100)
acf(oracle$sigma[seq(50000,150000, by = 100)], main = "")
dev.off()

png("hillps_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ps$sigma[seq(50000,150000, by = 100)], type = "l",
     ylab = expression(sigma),
     xlab = "Draw",
     ylim= c(0.45, 0.65))
dev.off()
png("hillps_2.png", width = 7, height = 7, units = 'in', res = 100)
acf(ps$sigma[seq(50000,150000, by = 100)], main = "")
dev.off()

png("hillglm_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ps3$sigma[seq(50000,150000, by = 100)], type = "l",
     ylab = expression(sigma),
     xlab = "Draw",
     ylim= c(0.45, 0.65))
dev.off()
png("hillglm_2.png", width = 7, height = 7, units = 'in', res = 100)
acf(ps3$sigma[seq(50000,150000, by = 100)], main = "")
dev.off()

png("hillrand_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ps4$sigma[seq(50000,150000, by = 100)], type = "l",
     ylab = expression(sigma),
     xlab = "Draw",
     ylim= c(0.45, 0.65))
dev.off()
png("hillrand_2.png", width = 7, height = 7, units = 'in', res = 100)
acf(ps4$sigma[seq(50000,150000, by = 100)], main = "")
dev.off()

png("hillbartbcf_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(fitbcf2$sigma, type = "l",
     ylab = expression(sigma),
     xlab = "Draw",
     ylim= c(0.45, 0.65))
dev.off()
png("hillbartbcf_2.png", width = 7, height = 7, units = 'in', res = 100)
acf(fitbcf2$sigma, main = "")
dev.off()

png("hilloraclebcf_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(fitbcf3$sigma, type = "l",
     ylab = expression(sigma),
     xlab = "Draw",
     ylim= c(0.45, 0.65))
dev.off()
png("hilloraclebcf_2.png", width = 7, height = 7, units = 'in', res = 100)
acf(fitbcf3$sigma, main = "")
dev.off()

png("hillglmbcf_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(fitbcf4$sigma, type = "l",
     ylab = expression(sigma),
     xlab = "Draw",
     ylim= c(0.45, 0.65))
dev.off()
png("hillglmbcf_2.png", width = 7, height = 7, units = 'in', res = 100)
acf(fitbcf4$sigma, main = "")
dev.off()

png("hillrandbcf_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(fitbcf5$sigma, type = "l",
     ylab = expression(sigma),
     xlab = "Draw",
     ylim= c(0.45, 0.65))
dev.off()
png("hillrandbcf_2.png", width = 7, height = 7, units = 'in', res = 100)
acf(fitbcf5$sigma, main = "")
dev.off()



par(cex=0.5)

#################
#####Vanilla#####
#################

png("hill.vanilla.icex1_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x1,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[1]))
lines(as.numeric(names(apply(ice.bart.vanilla.x1.0975$ice_curves,2,mean))),
      apply(ice.bart.vanilla.x1.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.vanilla.x1.0025$ice_curves,2,mean))),
      apply(ice.bart.vanilla.x1.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.vanilla.icex1_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x1.0025,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("hill.vanilla.icex1_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x1.0975,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("hill.vanilla.icex1_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x1,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[1]))
dev.off()

png("hill.vanilla.icex1_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.vanilla.x1),
     ylab = "Derivative Y",
     xlab = expression(x[1]))
dev.off()


png("hill.vanilla.icex2_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x2,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[2]))
lines(as.numeric(names(apply(ice.bart.vanilla.x2.0975$ice_curves,2,mean))),
      apply(ice.bart.vanilla.x2.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.vanilla.x2.0025$ice_curves,2,mean))),
      apply(ice.bart.vanilla.x2.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.vanilla.icex2_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x2.0025,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("hill.vanilla.icex2_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x2.0975,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("hill.vanilla.icex2_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x2,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[2]))
dev.off()

png("hill.vanilla.icex2_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.vanilla.x2),
     ylab = "Derivative Y",
     xlab = expression(x[2]))
dev.off()


png("hill.vanilla.icex3_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x3,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[3]))
lines(as.numeric(names(apply(ice.bart.vanilla.x3.0975$ice_curves,2,mean))),
      apply(ice.bart.vanilla.x3.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.vanilla.x3.0025$ice_curves,2,mean))),
      apply(ice.bart.vanilla.x3.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
{segments(x0 = -3, y0 = 0, x1 = -3/4*sd(x[3,])+mean(x[3,]), y1 = 0, col = "red", lwd = 3)
  segments(x0 = -3/4*sd(x[3,])+mean(x[3,]), y0 = 0, x1 = -3/4*sd(x[3,])+mean(x[3,]), y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = -3/4*sd(x[3,])+mean(x[3,]), y0 = 0.5, x1 = 0*sd(x[3,])+mean(x[3,]), y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 0*sd(x[3,])+mean(x[3,]), y0 = 0.5, x1 = 0*sd(x[3,])+mean(x[3,]), y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 0*sd(x[3,])+mean(x[3,]), y0 = 0.75, x1 = 3/4*sd(x[3,])+mean(x[3,]), y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 3/4*sd(x[3,])+mean(x[3,]), y0 = 0.75, x1 = 3/4*sd(x[3,])+mean(x[3,]), y1 = 1, col = "red", lwd = 3)
  segments(x0 = 3/4*sd(x[3,])+mean(x[3,]), y0 = 1, x1 = 15*sd(x[3,])+mean(x[3,]), y1 = 1, col = "red", lwd = 3)}

dev.off()


png("hill.vanilla.icex3_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x3.0025,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("hill.vanilla.icex3_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x3.0975,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("hill.vanilla.icex3_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x3,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[3]))
dev.off()

png("hill.vanilla.icex3_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.vanilla.x3),
     ylab = "Derivative Y",
     xlab = expression(x[3]))
dev.off()


png("hill.vanilla.icex4_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x4,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[4]))
lines(as.numeric(names(apply(ice.bart.vanilla.x4.0975$ice_curves,2,mean))),
      apply(ice.bart.vanilla.x4.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.vanilla.x4.0025$ice_curves,2,mean))),
      apply(ice.bart.vanilla.x4.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.vanilla.icex4_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x4.0025,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("hill.vanilla.icex4_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x4.0975,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("hill.vanilla.icex4_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x4,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[4]))
dev.off()

png("hill.vanilla.icex4_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.vanilla.x4),
     ylab = "Derivative Y",
     xlab = expression(x[4]))
dev.off()

png("hill.vanilla.icex5_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x5,
     ylim = c(-5,3),
     ylab = "Partial Y",
     xlab = expression(x[5]))
lines(as.numeric(names(apply(ice.bart.vanilla.x5.0975$ice_curves,2,mean))),
      apply(ice.bart.vanilla.x5.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.vanilla.x5.0025$ice_curves,2,mean))),
      apply(ice.bart.vanilla.x5.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.vanilla.icex5_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x5.0025,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("hill.vanilla.icex5_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x5.0975,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("hill.vanilla.icex5_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x5,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[5]))
dev.off()

png("hill.vanilla.icex5_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.vanilla.x5),
     ylab = "Derivative Y",
     xlab = expression(x[5]))
dev.off()

png("hill.vanilla.icex6_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x6,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[6]))
lines(as.numeric(names(apply(ice.bart.vanilla.x6.0975$ice_curves,2,mean))),
      apply(ice.bart.vanilla.x6.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.vanilla.x6.0025$ice_curves,2,mean))),
      apply(ice.bart.vanilla.x6.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.vanilla.icex6_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x6.0025,
     ylab = "Partial Y",
     xlab = expression(x[6]))
dev.off()

png("hill.vanilla.icex6_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x6.0975,
     ylab = "Partial Y",
     xlab = expression(x[6]))
dev.off()

png("hill.vanilla.icex6_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x6,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[6]))
dev.off()

png("hill.vanilla.icex6_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.vanilla.x6),
     ylab = "Derivative Y",
     xlab = expression(x[6]))
dev.off()

###################
######Oracle#######
###################

png("hill.oracle.icex1_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x1,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[1]))
lines(as.numeric(names(apply(ice.bart.oracle.x1.0975$ice_curves,2,mean))),
      apply(ice.bart.oracle.x1.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.oracle.x1.0025$ice_curves,2,mean))),
      apply(ice.bart.oracle.x1.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.oracle.icex1_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x1.0025,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("hill.oracle.icex1_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x1.0975,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("hill.oracle.icex1_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x1,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[1]))
dev.off()

png("hill.oracle.icex1_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.oracle.x1),
     ylab = "Derivative Y",
     xlab = expression(x[1]))
dev.off()


png("hill.oracle.icex2_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x2,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[2]))
lines(as.numeric(names(apply(ice.bart.oracle.x2.0975$ice_curves,2,mean))),
      apply(ice.bart.oracle.x2.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.oracle.x2.0025$ice_curves,2,mean))),
      apply(ice.bart.oracle.x2.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.oracle.icex2_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x2.0025,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("hill.oracle.icex2_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x2.0975,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("hill.oracle.icex2_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x2,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[2]))
dev.off()

png("hill.oracle.icex2_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.oracle.x2),
     ylab = "Derivative Y",
     xlab = expression(x[2]))
dev.off()


png("hill.oracle.icex3_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x3,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[3]))
lines(as.numeric(names(apply(ice.bart.oracle.x3.0975$ice_curves,2,mean))),
      apply(ice.bart.oracle.x3.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.oracle.x3.0025$ice_curves,2,mean))),
      apply(ice.bart.oracle.x3.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
{segments(x0 = -3, y0 = 0, x1 = -3/4*sd(x[3,])+mean(x[3,]), y1 = 0, col = "red", lwd = 3)
  segments(x0 = -3/4*sd(x[3,])+mean(x[3,]), y0 = 0, x1 = -3/4*sd(x[3,])+mean(x[3,]), y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = -3/4*sd(x[3,])+mean(x[3,]), y0 = 0.5, x1 = 0*sd(x[3,])+mean(x[3,]), y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 0*sd(x[3,])+mean(x[3,]), y0 = 0.5, x1 = 0*sd(x[3,])+mean(x[3,]), y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 0*sd(x[3,])+mean(x[3,]), y0 = 0.75, x1 = 3/4*sd(x[3,])+mean(x[3,]), y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 3/4*sd(x[3,])+mean(x[3,]), y0 = 0.75, x1 = 3/4*sd(x[3,])+mean(x[3,]), y1 = 1, col = "red", lwd = 3)
  segments(x0 = 3/4*sd(x[3,])+mean(x[3,]), y0 = 1, x1 = 15*sd(x[3,])+mean(x[3,]), y1 = 1, col = "red", lwd = 3)}

dev.off()


png("hill.oracle.icex3_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x3.0025,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("hill.oracle.icex3_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x3.0975,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("hill.oracle.icex3_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x3,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[3]))
dev.off()

png("hill.oracle.icex3_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.oracle.x3),
     ylab = "Derivative Y",
     xlab = expression(x[3]))
dev.off()

png("hill.oracle.icepihat_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.pihat,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
lines(as.numeric(names(apply(ice.bart.oracle.pihat.0975$ice_curves,2,mean))),
      apply(ice.bart.oracle.pihat.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.oracle.pihat.0025$ice_curves,2,mean))),
      apply(ice.bart.oracle.pihat.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()

png("hill.oracle.icepihat_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.pihat.0025,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("hill.oracle.icepihat_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.pihat.0975,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("hill.oracle.icepihat_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.pihat,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("hill.oracle.icepihat_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.oracle.pihat),
     ylab = "Derivative Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("hill.oracle.icex4_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x4,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[4]))
lines(as.numeric(names(apply(ice.bart.oracle.x4.0975$ice_curves,2,mean))),
      apply(ice.bart.oracle.x4.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.oracle.x4.0025$ice_curves,2,mean))),
      apply(ice.bart.oracle.x4.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.oracle.icex4_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x4.0025,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("hill.oracle.icex4_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x4.0975,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("hill.oracle.icex4_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x4,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[4]))
dev.off()

png("hill.oracle.icex4_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.oracle.x4),
     ylab = "Derivative Y",
     xlab = expression(x[4]))
dev.off()

png("hill.oracle.icex5_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x5,
     ylim = c(-5,3),
     ylab = "Partial Y",
     xlab = expression(x[5]))
lines(as.numeric(names(apply(ice.bart.oracle.x5.0975$ice_curves,2,mean))),
      apply(ice.bart.oracle.x5.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.oracle.x5.0025$ice_curves,2,mean))),
      apply(ice.bart.oracle.x5.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.oracle.icex5_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x5.0025,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("hill.oracle.icex5_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x5.0975,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("hill.oracle.icex5_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x5,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[5]))
dev.off()

png("hill.oracle.icex5_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.oracle.x5),
     ylab = "Derivative Y",
     xlab = expression(x[5]))
dev.off()

png("hill.oracle.icex6_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x6,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[6]))
lines(as.numeric(names(apply(ice.bart.oracle.x6.0975$ice_curves,2,mean))),
      apply(ice.bart.oracle.x6.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.oracle.x6.0025$ice_curves,2,mean))),
      apply(ice.bart.oracle.x6.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.oracle.icex6_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x6.0025,
     ylab = "Partial Y",
     xlab = expression(x[6]))
dev.off()

png("hill.oracle.icex6_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x6.0975,
     ylab = "Partial Y",
     xlab = expression(x[6]))
dev.off()

png("hill.oracle.icex6_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x6,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[6]))
dev.off()

png("hill.oracle.icex6_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.oracle.x6),
     ylab = "Derivative Y",
     xlab = expression(x[6]))
dev.off()

###########
#####PS####
###########

png("hill.ps.icex1_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x1,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[1]))
lines(as.numeric(names(apply(ice.bart.ps.x1.0975$ice_curves,2,mean))),
      apply(ice.bart.ps.x1.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.ps.x1.0025$ice_curves,2,mean))),
      apply(ice.bart.ps.x1.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.ps.icex1_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x1.0025,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("hill.ps.icex1_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x1.0975,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("hill.ps.icex1_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x1,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[1]))
dev.off()

png("hill.ps.icex1_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps.x1),
     ylab = "Derivative Y",
     xlab = expression(x[1]))
dev.off()


png("hill.ps.icex2_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x2,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[2]))
lines(as.numeric(names(apply(ice.bart.ps.x2.0975$ice_curves,2,mean))),
      apply(ice.bart.ps.x2.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.ps.x2.0025$ice_curves,2,mean))),
      apply(ice.bart.ps.x2.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.ps.icex2_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x2.0025,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("hill.ps.icex2_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x2.0975,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("hill.ps.icex2_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x2,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[2]))
dev.off()

png("hill.ps.icex2_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps.x2),
     ylab = "Derivative Y",
     xlab = expression(x[2]))
dev.off()


png("hill.ps.icex3_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x3,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[3]))
lines(as.numeric(names(apply(ice.bart.ps.x3.0975$ice_curves,2,mean))),
      apply(ice.bart.ps.x3.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.ps.x3.0025$ice_curves,2,mean))),
      apply(ice.bart.ps.x3.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
{segments(x0 = -3, y0 = 0, x1 = -3/4*sd(x[3,])+mean(x[3,]), y1 = 0, col = "red", lwd = 3)
  segments(x0 = -3/4*sd(x[3,])+mean(x[3,]), y0 = 0, x1 = -3/4*sd(x[3,])+mean(x[3,]), y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = -3/4*sd(x[3,])+mean(x[3,]), y0 = 0.5, x1 = 0*sd(x[3,])+mean(x[3,]), y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 0*sd(x[3,])+mean(x[3,]), y0 = 0.5, x1 = 0*sd(x[3,])+mean(x[3,]), y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 0*sd(x[3,])+mean(x[3,]), y0 = 0.75, x1 = 3/4*sd(x[3,])+mean(x[3,]), y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 3/4*sd(x[3,])+mean(x[3,]), y0 = 0.75, x1 = 3/4*sd(x[3,])+mean(x[3,]), y1 = 1, col = "red", lwd = 3)
  segments(x0 = 3/4*sd(x[3,])+mean(x[3,]), y0 = 1, x1 = 15*sd(x[3,])+mean(x[3,]), y1 = 1, col = "red", lwd = 3)}

dev.off()


png("hill.ps.icex3_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x3.0025,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("hill.ps.icex3_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x3.0975,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("hill.ps.icex3_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x3,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[3]))
dev.off()

png("hill.ps.icex3_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps.x3),
     ylab = "Derivative Y",
     xlab = expression(x[3]))
dev.off()

png("hill.ps.icepihat_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.pihat,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
lines(as.numeric(names(apply(ice.bart.ps.pihat.0975$ice_curves,2,mean))),
      apply(ice.bart.ps.pihat.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.ps.pihat.0025$ice_curves,2,mean))),
      apply(ice.bart.ps.pihat.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.ps.icepihat_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.pihat.0025,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("hill.ps.icepihat_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.pihat.0975,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("hill.ps.icepihat_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.pihat,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("hill.ps.icepihat_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps.pihat),
     ylab = "Derivative Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("hill.ps.icex4_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x4,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[4]))
lines(as.numeric(names(apply(ice.bart.ps.x4.0975$ice_curves,2,mean))),
      apply(ice.bart.ps.x4.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.ps.x4.0025$ice_curves,2,mean))),
      apply(ice.bart.ps.x4.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.ps.icex4_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x4.0025,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("hill.ps.icex4_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x4.0975,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("hill.ps.icex4_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x4,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[4]))
dev.off()

png("hill.ps.icex4_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps.x4),
     ylab = "Derivative Y",
     xlab = expression(x[4]))
dev.off()

png("hill.ps.icex5_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x5,
     ylim = c(-5,3),
     ylab = "Partial Y",
     xlab = expression(x[5]))
lines(as.numeric(names(apply(ice.bart.ps.x5.0975$ice_curves,2,mean))),
      apply(ice.bart.ps.x5.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.ps.x5.0025$ice_curves,2,mean))),
      apply(ice.bart.ps.x5.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.ps.icex5_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x5.0025,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("hill.ps.icex5_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x5.0975,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("hill.ps.icex5_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x5,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[5]))
dev.off()

png("hill.ps.icex5_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps.x5),
     ylab = "Derivative Y",
     xlab = expression(x[5]))
dev.off()

png("hill.ps.icex6_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x6,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[6]))
lines(as.numeric(names(apply(ice.bart.ps.x6.0975$ice_curves,2,mean))),
      apply(ice.bart.ps.x6.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.ps.x6.0025$ice_curves,2,mean))),
      apply(ice.bart.ps.x6.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.ps.icex6_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x6.0025,
     ylab = "Partial Y",
     xlab = expression(x[6]))
dev.off()

png("hill.ps.icex6_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x6.0975,
     ylab = "Partial Y",
     xlab = expression(x[6]))
dev.off()

png("hill.ps.icex6_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x6,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[6]))
dev.off()

png("hill.ps.icex6_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps.x6),
     ylab = "Derivative Y",
     xlab = expression(x[6]))
dev.off()

#########
###PS3###
#########

png("hill.ps3.icex1_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x1,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[1]))
lines(as.numeric(names(apply(ice.bart.ps3.x1.0975$ice_curves,2,mean))),
      apply(ice.bart.ps3.x1.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.ps3.x1.0025$ice_curves,2,mean))),
      apply(ice.bart.ps3.x1.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.ps3.icex1_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x1.0025,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("hill.ps3.icex1_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x1.0975,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("hill.ps3.icex1_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x1,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[1]))
dev.off()

png("hill.ps3.icex1_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps3.x1),
     ylab = "Derivative Y",
     xlab = expression(x[1]))
dev.off()


png("hill.ps3.icex2_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x2,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[2]))
lines(as.numeric(names(apply(ice.bart.ps3.x2.0975$ice_curves,2,mean))),
      apply(ice.bart.ps3.x2.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.ps3.x2.0025$ice_curves,2,mean))),
      apply(ice.bart.ps3.x2.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.ps3.icex2_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x2.0025,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("hill.ps3.icex2_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x2.0975,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("hill.ps3.icex2_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x2,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[2]))
dev.off()

png("hill.ps3.icex2_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps3.x2),
     ylab = "Derivative Y",
     xlab = expression(x[2]))
dev.off()


png("hill.ps3.icex3_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x3,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[3]))
lines(as.numeric(names(apply(ice.bart.ps3.x3.0975$ice_curves,2,mean))),
      apply(ice.bart.ps3.x3.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.ps3.x3.0025$ice_curves,2,mean))),
      apply(ice.bart.ps3.x3.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
{segments(x0 = -3, y0 = 0, x1 = -3/4*sd(x[3,])+mean(x[3,]), y1 = 0, col = "red", lwd = 3)
  segments(x0 = -3/4*sd(x[3,])+mean(x[3,]), y0 = 0, x1 = -3/4*sd(x[3,])+mean(x[3,]), y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = -3/4*sd(x[3,])+mean(x[3,]), y0 = 0.5, x1 = 0*sd(x[3,])+mean(x[3,]), y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 0*sd(x[3,])+mean(x[3,]), y0 = 0.5, x1 = 0*sd(x[3,])+mean(x[3,]), y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 0*sd(x[3,])+mean(x[3,]), y0 = 0.75, x1 = 3/4*sd(x[3,])+mean(x[3,]), y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 3/4*sd(x[3,])+mean(x[3,]), y0 = 0.75, x1 = 3/4*sd(x[3,])+mean(x[3,]), y1 = 1, col = "red", lwd = 3)
  segments(x0 = 3/4*sd(x[3,])+mean(x[3,]), y0 = 1, x1 = 15*sd(x[3,])+mean(x[3,]), y1 = 1, col = "red", lwd = 3)}

dev.off()


png("hill.ps3.icex3_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x3.0025,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("hill.ps3.icex3_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x3.0975,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("hill.ps3.icex3_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x3,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[3]))
dev.off()

png("hill.ps3.icex3_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps3.x3),
     ylab = "Derivative Y",
     xlab = expression(x[3]))
dev.off()

png("hill.ps3.icepihat_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.pihat,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
lines(as.numeric(names(apply(ice.bart.ps3.pihat.0975$ice_curves,2,mean))),
      apply(ice.bart.ps3.pihat.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.ps3.pihat.0025$ice_curves,2,mean))),
      apply(ice.bart.ps3.pihat.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.ps3.icepihat_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.pihat.0025,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("hill.ps3.icepihat_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.pihat.0975,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("hill.ps3.icepihat_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.pihat,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("hill.ps3.icepihat_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps3.pihat),
     ylab = "Derivative Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("hill.ps3.icex4_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x4,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[4]))
lines(as.numeric(names(apply(ice.bart.ps3.x4.0975$ice_curves,2,mean))),
      apply(ice.bart.ps3.x4.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.ps3.x4.0025$ice_curves,2,mean))),
      apply(ice.bart.ps3.x4.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.ps3.icex4_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x4.0025,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("hill.ps3.icex4_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x4.0975,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("hill.ps3.icex4_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x4,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[4]))
dev.off()

png("hill.ps3.icex4_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps3.x4),
     ylab = "Derivative Y",
     xlab = expression(x[4]))
dev.off()

png("hill.ps3.icex5_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x5,
     ylim = c(-5,3),
     ylab = "Partial Y",
     xlab = expression(x[5]))
lines(as.numeric(names(apply(ice.bart.ps3.x5.0975$ice_curves,2,mean))),
      apply(ice.bart.ps3.x5.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.ps3.x5.0025$ice_curves,2,mean))),
      apply(ice.bart.ps3.x5.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.ps3.icex5_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x5.0025,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("hill.ps3.icex5_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x5.0975,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("hill.ps3.icex5_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x5,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[5]))
dev.off()

png("hill.ps3.icex5_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps3.x5),
     ylab = "Derivative Y",
     xlab = expression(x[5]))
dev.off()

png("hill.ps3.icex6_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x6,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[6]))
lines(as.numeric(names(apply(ice.bart.ps3.x6.0975$ice_curves,2,mean))),
      apply(ice.bart.ps3.x6.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.ps3.x6.0025$ice_curves,2,mean))),
      apply(ice.bart.ps3.x6.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.ps3.icex6_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x6.0025,
     ylab = "Partial Y",
     xlab = expression(x[6]))
dev.off()

png("hill.ps3.icex6_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x6.0975,
     ylab = "Partial Y",
     xlab = expression(x[6]))
dev.off()

png("hill.ps3.icex6_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x6,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[6]))
dev.off()

png("hill.ps3.icex6_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps3.x6),
     ylab = "Derivative Y",
     xlab = expression(x[6]))
dev.off()

#######
##PS4##
#######

png("hill.ps4.icex1_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x1,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[1]))
lines(as.numeric(names(apply(ice.bart.ps4.x1.0975$ice_curves,2,mean))),
      apply(ice.bart.ps4.x1.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.ps4.x1.0025$ice_curves,2,mean))),
      apply(ice.bart.ps4.x1.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.ps4.icex1_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x1.0025,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("hill.ps4.icex1_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x1.0975,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("hill.ps4.icex1_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x1,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[1]))
dev.off()

png("hill.ps4.icex1_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps4.x1),
     ylab = "Derivative Y",
     xlab = expression(x[1]))
dev.off()


png("hill.ps4.icex2_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x2,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[2]))
lines(as.numeric(names(apply(ice.bart.ps4.x2.0975$ice_curves,2,mean))),
      apply(ice.bart.ps4.x2.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.ps4.x2.0025$ice_curves,2,mean))),
      apply(ice.bart.ps4.x2.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.ps4.icex2_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x2.0025,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("hill.ps4.icex2_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x2.0975,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("hill.ps4.icex2_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x2,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[2]))
dev.off()

png("hill.ps4.icex2_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps4.x2),
     ylab = "Derivative Y",
     xlab = expression(x[2]))
dev.off()


png("hill.ps4.icex3_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x3,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[3]))
lines(as.numeric(names(apply(ice.bart.ps4.x3.0975$ice_curves,2,mean))),
      apply(ice.bart.ps4.x3.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.ps4.x3.0025$ice_curves,2,mean))),
      apply(ice.bart.ps4.x3.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
{segments(x0 = -3, y0 = 0, x1 = -3/4*sd(x[3,])+mean(x[3,]), y1 = 0, col = "red", lwd = 3)
  segments(x0 = -3/4*sd(x[3,])+mean(x[3,]), y0 = 0, x1 = -3/4*sd(x[3,])+mean(x[3,]), y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = -3/4*sd(x[3,])+mean(x[3,]), y0 = 0.5, x1 = 0*sd(x[3,])+mean(x[3,]), y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 0*sd(x[3,])+mean(x[3,]), y0 = 0.5, x1 = 0*sd(x[3,])+mean(x[3,]), y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 0*sd(x[3,])+mean(x[3,]), y0 = 0.75, x1 = 3/4*sd(x[3,])+mean(x[3,]), y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 3/4*sd(x[3,])+mean(x[3,]), y0 = 0.75, x1 = 3/4*sd(x[3,])+mean(x[3,]), y1 = 1, col = "red", lwd = 3)
  segments(x0 = 3/4*sd(x[3,])+mean(x[3,]), y0 = 1, x1 = 15*sd(x[3,])+mean(x[3,]), y1 = 1, col = "red", lwd = 3)}

dev.off()


png("hill.ps4.icex3_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x3.0025,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("hill.ps4.icex3_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x3.0975,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("hill.ps4.icex3_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x3,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[3]))
dev.off()

png("hill.ps4.icex3_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps4.x3),
     ylab = "Derivative Y",
     xlab = expression(x[3]))
dev.off()

png("hill.ps4.icepihat_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.pihat,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
lines(as.numeric(names(apply(ice.bart.ps4.pihat.0975$ice_curves,2,mean))),
      apply(ice.bart.ps4.pihat.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.ps4.pihat.0025$ice_curves,2,mean))),
      apply(ice.bart.ps4.pihat.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.ps4.icepihat_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.pihat.0025,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("hill.ps4.icepihat_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.pihat.0975,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("hill.ps4.icepihat_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.pihat,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("hill.ps4.icepihat_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps4.pihat),
     ylab = "Derivative Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("hill.ps4.icex4_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x4,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[4]))
lines(as.numeric(names(apply(ice.bart.ps4.x4.0975$ice_curves,2,mean))),
      apply(ice.bart.ps4.x4.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.ps4.x4.0025$ice_curves,2,mean))),
      apply(ice.bart.ps4.x4.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.ps4.icex4_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x4.0025,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("hill.ps4.icex4_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x4.0975,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("hill.ps4.icex4_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x4,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[4]))
dev.off()

png("hill.ps4.icex4_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps4.x4),
     ylab = "Derivative Y",
     xlab = expression(x[4]))
dev.off()

png("hill.ps4.icex5_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x5,
     ylim = c(-5,3),
     ylab = "Partial Y",
     xlab = expression(x[5]))
lines(as.numeric(names(apply(ice.bart.ps4.x5.0975$ice_curves,2,mean))),
      apply(ice.bart.ps4.x5.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.ps4.x5.0025$ice_curves,2,mean))),
      apply(ice.bart.ps4.x5.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.ps4.icex5_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x5.0025,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("hill.ps4.icex5_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x5.0975,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("hill.ps4.icex5_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x5,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[5]))
dev.off()

png("hill.ps4.icex5_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps4.x5),
     ylab = "Derivative Y",
     xlab = expression(x[5]))
dev.off()

png("hill.ps4.icex6_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x6,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[6]))
lines(as.numeric(names(apply(ice.bart.ps4.x6.0975$ice_curves,2,mean))),
      apply(ice.bart.ps4.x6.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.ps4.x6.0025$ice_curves,2,mean))),
      apply(ice.bart.ps4.x6.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.ps4.icex6_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x6.0025,
     ylab = "Partial Y",
     xlab = expression(x[6]))
dev.off()

png("hill.ps4.icex6_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x6.0975,
     ylab = "Partial Y",
     xlab = expression(x[6]))
dev.off()

png("hill.ps4.icex6_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x6,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[6]))
dev.off()

png("hill.ps4.icex6_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps4.x6),
     ylab = "Derivative Y",
     xlab = expression(x[6]))
dev.off()

#####################################
############## BCF ##################
#####################################

###############
####Oracle#####
###############


png("hill.bcf.oracle.icex1_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x1,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[1]))
lines(as.numeric(names(apply(ice.bcf.oracle.x1.0975$ice_curves,2,mean))),
      apply(ice.bcf.oracle.x1.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.oracle.x1.0025$ice_curves,2,mean))),
      apply(ice.bcf.oracle.x1.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.bcf.oracle.icex1_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x1.0025,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("hill.bcf.oracle.icex1_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x1.0975,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("hill.bcf.oracle.icex1_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x1,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[1]))
dev.off()

png("hill.bcf.oracle.icex1_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.oracle.x1),
     ylab = "Derivative Y",
     xlab = expression(x[1]))
dev.off()


png("hill.bcf.oracle.icex2_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x2,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[2]))
lines(as.numeric(names(apply(ice.bcf.oracle.x2.0975$ice_curves,2,mean))),
      apply(ice.bcf.oracle.x2.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.oracle.x2.0025$ice_curves,2,mean))),
      apply(ice.bcf.oracle.x2.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.bcf.oracle.icex2_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x2.0025,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("hill.bcf.oracle.icex2_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x2.0975,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("hill.bcf.oracle.icex2_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x2,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[2]))
dev.off()

png("hill.bcf.oracle.icex2_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.oracle.x2),
     ylab = "Derivative Y",
     xlab = expression(x[2]))
dev.off()


png("hill.bcf.oracle.icex3_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x3,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[3]))
lines(as.numeric(names(apply(ice.bcf.oracle.x3.0975$ice_curves,2,mean))),
      apply(ice.bcf.oracle.x3.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.oracle.x3.0025$ice_curves,2,mean))),
      apply(ice.bcf.oracle.x3.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
{segments(x0 = -3, y0 = 0, x1 = -3/4*sd(x[3,])+mean(x[3,]), y1 = 0, col = "red", lwd = 3)
  segments(x0 = -3/4*sd(x[3,])+mean(x[3,]), y0 = 0, x1 = -3/4*sd(x[3,])+mean(x[3,]), y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = -3/4*sd(x[3,])+mean(x[3,]), y0 = 0.5, x1 = 0*sd(x[3,])+mean(x[3,]), y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 0*sd(x[3,])+mean(x[3,]), y0 = 0.5, x1 = 0*sd(x[3,])+mean(x[3,]), y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 0*sd(x[3,])+mean(x[3,]), y0 = 0.75, x1 = 3/4*sd(x[3,])+mean(x[3,]), y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 3/4*sd(x[3,])+mean(x[3,]), y0 = 0.75, x1 = 3/4*sd(x[3,])+mean(x[3,]), y1 = 1, col = "red", lwd = 3)
  segments(x0 = 3/4*sd(x[3,])+mean(x[3,]), y0 = 1, x1 = 15*sd(x[3,])+mean(x[3,]), y1 = 1, col = "red", lwd = 3)}

dev.off()


png("hill.bcf.oracle.icex3_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x3.0025,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("hill.bcf.oracle.icex3_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x3.0975,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("hill.bcf.oracle.icex3_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x3,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[3]))
dev.off()

png("hill.bcf.oracle.icex3_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.oracle.x3),
     ylab = "Derivative Y",
     xlab = expression(x[3]))
dev.off()

png("hill.bcf.oracle.icepihat_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.pihat,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
lines(as.numeric(names(apply(ice.bcf.oracle.pihat.0975$ice_curves,2,mean))),
      apply(ice.bcf.oracle.pihat.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.oracle.pihat.0025$ice_curves,2,mean))),
      apply(ice.bcf.oracle.pihat.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.bcf.oracle.icepihat_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.pihat.0025,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("hill.bcf.oracle.icepihat_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.pihat.0975,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("hill.bcf.oracle.icepihat_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.pihat,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("hill.bcf.oracle.icepihat_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.oracle.pihat),
     ylab = "Derivative Y",
     xlab = expression(hat(Pi)(x)))
dev.off()


png("hill.bcf.oracle.icex4_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x4,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[4]))
lines(as.numeric(names(apply(ice.bcf.oracle.x4.0975$ice_curves,2,mean))),
      apply(ice.bcf.oracle.x4.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.oracle.x4.0025$ice_curves,2,mean))),
      apply(ice.bcf.oracle.x4.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.bcf.oracle.icex4_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x4.0025,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("hill.bcf.oracle.icex4_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x4.0975,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("hill.bcf.oracle.icex4_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x4,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[4]))
dev.off()

png("hill.bcf.oracle.icex4_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.oracle.x4),
     ylab = "Derivative Y",
     xlab = expression(x[4]))
dev.off()

png("hill.bcf.oracle.icex5_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x5,
     ylim = c(-5,3),
     ylab = "Partial Y",
     xlab = expression(x[5]))
lines(as.numeric(names(apply(ice.bcf.oracle.x5.0975$ice_curves,2,mean))),
      apply(ice.bcf.oracle.x5.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.oracle.x5.0025$ice_curves,2,mean))),
      apply(ice.bcf.oracle.x5.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.bcf.oracle.icex5_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x5.0025,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("hill.bcf.oracle.icex5_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x5.0975,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("hill.bcf.oracle.icex5_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x5,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[5]))
dev.off()

png("hill.bcf.oracle.icex5_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.oracle.x5),
     ylab = "Derivative Y",
     xlab = expression(x[5]))
dev.off()

png("hill.bcf.oracle.icex6_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x6,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[6]))
lines(as.numeric(names(apply(ice.bcf.oracle.x6.0975$ice_curves,2,mean))),
      apply(ice.bcf.oracle.x6.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.oracle.x6.0025$ice_curves,2,mean))),
      apply(ice.bcf.oracle.x6.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.bcf.oracle.icex6_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x6.0025,
     ylab = "Partial Y",
     xlab = expression(x[6]))
dev.off()

png("hill.bcf.oracle.icex6_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x6.0975,
     ylab = "Partial Y",
     xlab = expression(x[6]))
dev.off()

png("hill.bcf.oracle.icex6_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x6,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[6]))
dev.off()

png("hill.bcf.oracle.icex6_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.oracle.x6),
     ylab = "Derivative Y",
     xlab = expression(x[6]))
dev.off()

#############
#####Bart####
#############

png("hill.bcf.bart.icex1_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x1,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[1]))
lines(as.numeric(names(apply(ice.bcf.bart.x1.0975$ice_curves,2,mean))),
      apply(ice.bcf.bart.x1.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.bart.x1.0025$ice_curves,2,mean))),
      apply(ice.bcf.bart.x1.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.bcf.bart.icex1_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x1.0025,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("hill.bcf.bart.icex1_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x1.0975,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("hill.bcf.bart.icex1_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x1,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[1]))
dev.off()

png("hill.bcf.bart.icex1_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.bart.x1),
     ylab = "Derivative Y",
     xlab = expression(x[1]))
dev.off()


png("hill.bcf.bart.icex2_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x2,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[2]))
lines(as.numeric(names(apply(ice.bcf.bart.x2.0975$ice_curves,2,mean))),
      apply(ice.bcf.bart.x2.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.bart.x2.0025$ice_curves,2,mean))),
      apply(ice.bcf.bart.x2.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.bcf.bart.icex2_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x2.0025,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("hill.bcf.bart.icex2_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x2.0975,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("hill.bcf.bart.icex2_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x2,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[2]))
dev.off()

png("hill.bcf.bart.icex2_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.bart.x2),
     ylab = "Derivative Y",
     xlab = expression(x[2]))
dev.off()


png("hill.bcf.bart.icex3_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x3,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[3]))
lines(as.numeric(names(apply(ice.bcf.bart.x3.0975$ice_curves,2,mean))),
      apply(ice.bcf.bart.x3.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.bart.x3.0025$ice_curves,2,mean))),
      apply(ice.bcf.bart.x3.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
{segments(x0 = -3, y0 = 0, x1 = -3/4*sd(x[3,])+mean(x[3,]), y1 = 0, col = "red", lwd = 3)
  segments(x0 = -3/4*sd(x[3,])+mean(x[3,]), y0 = 0, x1 = -3/4*sd(x[3,])+mean(x[3,]), y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = -3/4*sd(x[3,])+mean(x[3,]), y0 = 0.5, x1 = 0*sd(x[3,])+mean(x[3,]), y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 0*sd(x[3,])+mean(x[3,]), y0 = 0.5, x1 = 0*sd(x[3,])+mean(x[3,]), y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 0*sd(x[3,])+mean(x[3,]), y0 = 0.75, x1 = 3/4*sd(x[3,])+mean(x[3,]), y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 3/4*sd(x[3,])+mean(x[3,]), y0 = 0.75, x1 = 3/4*sd(x[3,])+mean(x[3,]), y1 = 1, col = "red", lwd = 3)
  segments(x0 = 3/4*sd(x[3,])+mean(x[3,]), y0 = 1, x1 = 15*sd(x[3,])+mean(x[3,]), y1 = 1, col = "red", lwd = 3)}

dev.off()


png("hill.bcf.bart.icex3_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x3.0025,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("hill.bcf.bart.icex3_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x3.0975,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("hill.bcf.bart.icex3_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x3,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[3]))
dev.off()

png("hill.bcf.bart.icex3_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.bart.x3),
     ylab = "Derivative Y",
     xlab = expression(x[3]))
dev.off()

png("hill.bcf.bart.icepihat_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.pihat,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
lines(as.numeric(names(apply(ice.bcf.bart.pihat.0975$ice_curves,2,mean))),
      apply(ice.bcf.bart.pihat.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.bart.pihat.0025$ice_curves,2,mean))),
      apply(ice.bcf.bart.pihat.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.bcf.bart.icepihat_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.pihat.0025,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("hill.bcf.bart.icepihat_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.pihat.0975,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("hill.bcf.bart.icepihat_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.pihat,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("hill.bcf.bart.icepihat_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.bart.pihat),
     ylab = "Derivative Y",
     xlab = expression(hat(Pi)(x)))
dev.off()


png("hill.bcf.bart.icex4_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x4,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[4]))
lines(as.numeric(names(apply(ice.bcf.bart.x4.0975$ice_curves,2,mean))),
      apply(ice.bcf.bart.x4.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.bart.x4.0025$ice_curves,2,mean))),
      apply(ice.bcf.bart.x4.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.bcf.bart.icex4_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x4.0025,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("hill.bcf.bart.icex4_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x4.0975,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("hill.bcf.bart.icex4_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x4,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[4]))
dev.off()

png("hill.bcf.bart.icex4_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.bart.x4),
     ylab = "Derivative Y",
     xlab = expression(x[4]))
dev.off()

png("hill.bcf.bart.icex5_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x5,
     ylim = c(-5,3),
     ylab = "Partial Y",
     xlab = expression(x[5]))
lines(as.numeric(names(apply(ice.bcf.bart.x5.0975$ice_curves,2,mean))),
      apply(ice.bcf.bart.x5.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.bart.x5.0025$ice_curves,2,mean))),
      apply(ice.bcf.bart.x5.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.bcf.bart.icex5_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x5.0025,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("hill.bcf.bart.icex5_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x5.0975,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("hill.bcf.bart.icex5_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x5,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[5]))
dev.off()

png("hill.bcf.bart.icex5_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.bart.x5),
     ylab = "Derivative Y",
     xlab = expression(x[5]))
dev.off()

png("hill.bcf.bart.icex6_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x6,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[6]))
lines(as.numeric(names(apply(ice.bcf.bart.x6.0975$ice_curves,2,mean))),
      apply(ice.bcf.bart.x6.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.bart.x6.0025$ice_curves,2,mean))),
      apply(ice.bcf.bart.x6.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.bcf.bart.icex6_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x6.0025,
     ylab = "Partial Y",
     xlab = expression(x[6]))
dev.off()

png("hill.bcf.bart.icex6_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x6.0975,
     ylab = "Partial Y",
     xlab = expression(x[6]))
dev.off()

png("hill.bcf.bart.icex6_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x6,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[6]))
dev.off()

png("hill.bcf.bart.icex6_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.bart.x6),
     ylab = "Derivative Y",
     xlab = expression(x[6]))
dev.off()


#########
###glm###
#########

png("hill.bcf.glm.icex1_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x1,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[1]))
lines(as.numeric(names(apply(ice.bcf.glm.x1.0975$ice_curves,2,mean))),
      apply(ice.bcf.glm.x1.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.glm.x1.0025$ice_curves,2,mean))),
      apply(ice.bcf.glm.x1.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.bcf.glm.icex1_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x1.0025,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("hill.bcf.glm.icex1_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x1.0975,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("hill.bcf.glm.icex1_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x1,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[1]))
dev.off()

png("hill.bcf.glm.icex1_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.glm.x1),
     ylab = "Derivative Y",
     xlab = expression(x[1]))
dev.off()


png("hill.bcf.glm.icex2_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x2,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[2]))
lines(as.numeric(names(apply(ice.bcf.glm.x2.0975$ice_curves,2,mean))),
      apply(ice.bcf.glm.x2.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.glm.x2.0025$ice_curves,2,mean))),
      apply(ice.bcf.glm.x2.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.bcf.glm.icex2_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x2.0025,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("hill.bcf.glm.icex2_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x2.0975,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("hill.bcf.glm.icex2_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x2,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[2]))
dev.off()

png("hill.bcf.glm.icex2_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.glm.x2),
     ylab = "Derivative Y",
     xlab = expression(x[2]))
dev.off()


png("hill.bcf.glm.icex3_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x3,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[3]))
lines(as.numeric(names(apply(ice.bcf.glm.x3.0975$ice_curves,2,mean))),
      apply(ice.bcf.glm.x3.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.glm.x3.0025$ice_curves,2,mean))),
      apply(ice.bcf.glm.x3.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
{segments(x0 = -3, y0 = 0, x1 = -3/4*sd(x[3,])+mean(x[3,]), y1 = 0, col = "red", lwd = 3)
  segments(x0 = -3/4*sd(x[3,])+mean(x[3,]), y0 = 0, x1 = -3/4*sd(x[3,])+mean(x[3,]), y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = -3/4*sd(x[3,])+mean(x[3,]), y0 = 0.5, x1 = 0*sd(x[3,])+mean(x[3,]), y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 0*sd(x[3,])+mean(x[3,]), y0 = 0.5, x1 = 0*sd(x[3,])+mean(x[3,]), y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 0*sd(x[3,])+mean(x[3,]), y0 = 0.75, x1 = 3/4*sd(x[3,])+mean(x[3,]), y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 3/4*sd(x[3,])+mean(x[3,]), y0 = 0.75, x1 = 3/4*sd(x[3,])+mean(x[3,]), y1 = 1, col = "red", lwd = 3)
  segments(x0 = 3/4*sd(x[3,])+mean(x[3,]), y0 = 1, x1 = 15*sd(x[3,])+mean(x[3,]), y1 = 1, col = "red", lwd = 3)}

dev.off()


png("hill.bcf.glm.icex3_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x3.0025,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("hill.bcf.glm.icex3_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x3.0975,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("hill.bcf.glm.icex3_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x3,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[3]))
dev.off()

png("hill.bcf.glm.icex3_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.glm.x3),
     ylab = "Derivative Y",
     xlab = expression(x[3]))
dev.off()

png("hill.bcf.glm.icepihat_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.pihat,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
lines(as.numeric(names(apply(ice.bcf.glm.pihat.0975$ice_curves,2,mean))),
      apply(ice.bcf.glm.pihat.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.glm.pihat.0025$ice_curves,2,mean))),
      apply(ice.bcf.glm.pihat.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.bcf.glm.icepihat_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.pihat.0025,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("hill.bcf.glm.icepihat_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.pihat.0975,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("hill.bcf.glm.icepihat_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.pihat,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("hill.bcf.glm.icepihat_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.glm.pihat),
     ylab = "Derivative Y",
     xlab = expression(hat(Pi)(x)))
dev.off()


png("hill.bcf.glm.icex4_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x4,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[4]))
lines(as.numeric(names(apply(ice.bcf.glm.x4.0975$ice_curves,2,mean))),
      apply(ice.bcf.glm.x4.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.glm.x4.0025$ice_curves,2,mean))),
      apply(ice.bcf.glm.x4.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.bcf.glm.icex4_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x4.0025,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("hill.bcf.glm.icex4_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x4.0975,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("hill.bcf.glm.icex4_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x4,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[4]))
dev.off()

png("hill.bcf.glm.icex4_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.glm.x4),
     ylab = "Derivative Y",
     xlab = expression(x[4]))
dev.off()

png("hill.bcf.glm.icex5_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x5,
     ylim = c(-5,3),
     ylab = "Partial Y",
     xlab = expression(x[5]))
lines(as.numeric(names(apply(ice.bcf.glm.x5.0975$ice_curves,2,mean))),
      apply(ice.bcf.glm.x5.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.glm.x5.0025$ice_curves,2,mean))),
      apply(ice.bcf.glm.x5.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.bcf.glm.icex5_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x5.0025,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("hill.bcf.glm.icex5_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x5.0975,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("hill.bcf.glm.icex5_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x5,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[5]))
dev.off()

png("hill.bcf.glm.icex5_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.glm.x5),
     ylab = "Derivative Y",
     xlab = expression(x[5]))
dev.off()

png("hill.bcf.glm.icex6_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x6,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[6]))
lines(as.numeric(names(apply(ice.bcf.glm.x6.0975$ice_curves,2,mean))),
      apply(ice.bcf.glm.x6.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.glm.x6.0025$ice_curves,2,mean))),
      apply(ice.bcf.glm.x6.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.bcf.glm.icex6_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x6.0025,
     ylab = "Partial Y",
     xlab = expression(x[6]))
dev.off()

png("hill.bcf.glm.icex6_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x6.0975,
     ylab = "Partial Y",
     xlab = expression(x[6]))
dev.off()

png("hill.bcf.glm.icex6_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x6,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[6]))
dev.off()

png("hill.bcf.glm.icex6_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.glm.x6),
     ylab = "Derivative Y",
     xlab = expression(x[6]))
dev.off()


#######
##rand##
#######

png("hill.bcf.rand.icex1_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x1,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[1]))
lines(as.numeric(names(apply(ice.bcf.rand.x1.0975$ice_curves,2,mean))),
      apply(ice.bcf.rand.x1.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.rand.x1.0025$ice_curves,2,mean))),
      apply(ice.bcf.rand.x1.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.bcf.rand.icex1_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x1.0025,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("hill.bcf.rand.icex1_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x1.0975,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("hill.bcf.rand.icex1_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x1,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[1]))
dev.off()

png("hill.bcf.rand.icex1_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.rand.x1),
     ylab = "Derivative Y",
     xlab = expression(x[1]))
dev.off()


png("hill.bcf.rand.icex2_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x2,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[2]))
lines(as.numeric(names(apply(ice.bcf.rand.x2.0975$ice_curves,2,mean))),
      apply(ice.bcf.rand.x2.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.rand.x2.0025$ice_curves,2,mean))),
      apply(ice.bcf.rand.x2.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.bcf.rand.icex2_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x2.0025,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("hill.bcf.rand.icex2_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x2.0975,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("hill.bcf.rand.icex2_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x2,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[2]))
dev.off()

png("hill.bcf.rand.icex2_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.rand.x2),
     ylab = "Derivative Y",
     xlab = expression(x[2]))
dev.off()


png("hill.bcf.rand.icex3_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x3,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[3]))
lines(as.numeric(names(apply(ice.bcf.rand.x3.0975$ice_curves,2,mean))),
      apply(ice.bcf.rand.x3.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.rand.x3.0025$ice_curves,2,mean))),
      apply(ice.bcf.rand.x3.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
{segments(x0 = -3, y0 = 0, x1 = -3/4*sd(x[3,])+mean(x[3,]), y1 = 0, col = "red", lwd = 3)
  segments(x0 = -3/4*sd(x[3,])+mean(x[3,]), y0 = 0, x1 = -3/4*sd(x[3,])+mean(x[3,]), y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = -3/4*sd(x[3,])+mean(x[3,]), y0 = 0.5, x1 = 0*sd(x[3,])+mean(x[3,]), y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 0*sd(x[3,])+mean(x[3,]), y0 = 0.5, x1 = 0*sd(x[3,])+mean(x[3,]), y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 0*sd(x[3,])+mean(x[3,]), y0 = 0.75, x1 = 3/4*sd(x[3,])+mean(x[3,]), y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 3/4*sd(x[3,])+mean(x[3,]), y0 = 0.75, x1 = 3/4*sd(x[3,])+mean(x[3,]), y1 = 1, col = "red", lwd = 3)
  segments(x0 = 3/4*sd(x[3,])+mean(x[3,]), y0 = 1, x1 = 15*sd(x[3,])+mean(x[3,]), y1 = 1, col = "red", lwd = 3)}
dev.off()


png("hill.bcf.rand.icex3_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x3.0025,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("hill.bcf.rand.icex3_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x3.0975,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("hill.bcf.rand.icex3_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x3,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[3]))
dev.off()

png("hill.bcf.rand.icex3_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.rand.x3),
     ylab = "Derivative Y",
     xlab = expression(x[3]))
dev.off()

png("hill.bcf.rand.icepihat_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.pihat,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
lines(as.numeric(names(apply(ice.bcf.rand.pihat.0975$ice_curves,2,mean))),
      apply(ice.bcf.rand.pihat.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.rand.pihat.0025$ice_curves,2,mean))),
      apply(ice.bcf.rand.pihat.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.bcf.rand.icepihat_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.pihat.0025,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("hill.bcf.rand.icepihat_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.pihat.0975,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("hill.bcf.rand.icepihat_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.pihat,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("hill.bcf.rand.icepihat_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.rand.pihat),
     ylab = "Derivative Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("hill.bcf.rand.icex4_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x4,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[4]))
lines(as.numeric(names(apply(ice.bcf.rand.x4.0975$ice_curves,2,mean))),
      apply(ice.bcf.rand.x4.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.rand.x4.0025$ice_curves,2,mean))),
      apply(ice.bcf.rand.x4.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.bcf.rand.icex4_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x4.0025,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("hill.bcf.rand.icex4_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x4.0975,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("hill.bcf.rand.icex4_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x4,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[4]))
dev.off()

png("hill.bcf.rand.icex4_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.rand.x4),
     ylab = "Derivative Y",
     xlab = expression(x[4]))
dev.off()

png("hill.bcf.rand.icex5_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x5,
     ylim = c(-5,3),
     ylab = "Partial Y",
     xlab = expression(x[5]))
lines(as.numeric(names(apply(ice.bcf.rand.x5.0975$ice_curves,2,mean))),
      apply(ice.bcf.rand.x5.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.rand.x5.0025$ice_curves,2,mean))),
      apply(ice.bcf.rand.x5.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.bcf.rand.icex5_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x5.0025,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("hill.bcf.rand.icex5_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x5.0975,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("hill.bcf.rand.icex5_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x5,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[5]))
dev.off()

png("hill.bcf.rand.icex5_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.rand.x5),
     ylab = "Derivative Y",
     xlab = expression(x[5]))
dev.off()

png("hill.bcf.rand.icex6_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x6,
     ylim = c(-4,3),
     ylab = "Partial Y",
     xlab = expression(x[6]))
lines(as.numeric(names(apply(ice.bcf.rand.x6.0975$ice_curves,2,mean))),
      apply(ice.bcf.rand.x6.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.rand.x6.0025$ice_curves,2,mean))),
      apply(ice.bcf.rand.x6.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("hill.bcf.rand.icex6_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x6.0025,
     ylab = "Partial Y",
     xlab = expression(x[6]))
dev.off()

png("hill.bcf.rand.icex6_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x6.0975,
     ylab = "Partial Y",
     xlab = expression(x[6]))
dev.off()

png("hill.bcf.rand.icex6_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x6,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[6]))
dev.off()

png("hill.bcf.rand.icex6_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.rand.x6),
     ylab = "Derivative Y",
     xlab = expression(x[6]))
dev.off()


