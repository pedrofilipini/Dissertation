######################################
##Pedro Henrique Filipini dos Santos##
###University of São Paulo - Brazil###
######################################


##############################################
##Simulation study based on Freidman Example##
##############################################

#Libraries
library(BART)
library(ICEbox)
#library(doSNOW)
#library(foreach)

niters <- 1000
pthin <- 500
thin <- 100
burn <- 50000
npost <- 1000
ncores <- 20

#RMSE ITE function
rmseITE <- function(ite,x){
  r <- numeric(0)
  n <- length(ite)
  for(i in 1:n){
    if(((x[3,(i)])>(1/4)) && ((x[3,(i)])<(2/4))){
      r[i] <- (0.5-ite[i])^2
    }
    if(((x[3,(i)])<(2/4)) && ((x[3,(i)])>(3/4))){
      r[i] <- (0.75-ite[i])^2
    }
    if((x[3,(i)])>(3/4)){
      r[i] <- (1-ite[i])^2
    }else{
      r[i] <- (0-ite[i])^2
    }
  }
  sqrt(mean(r))
}

#RMSE ATE function
rmseATE <- function(ate, alpha){
  aux <- rowMeans(ate)
  sqrt(mean((aux-mean(alpha))^2))
}

#Function to create DFs of different sizes
joinDF <- function(Var, x, i){
  if(is.matrix(vanillaVar)==0){
    Var <- as.matrix(x)
  }else{
    Var <- cbind(Var, as.matrix(x))
  }
}

#Posterior Inclusion Probability
PIP <- function(fit){
  colMeans(fit$varcount > 0)
}

#starting RMSE vectors
ite_rmse_vanilla <- numeric(0)
ite_rmse_oracle <- numeric(0)
ite_rmse_ps <- numeric(0)
ite_rmse_psglm <- numeric(0)
ite_rmse_rand <- numeric(0)
dite_rmse_vanilla <- numeric(0)
dite_rmse_oracle <- numeric(0)
dite_rmse_ps <- numeric(0)
dite_rmse_psglm <- numeric(0)
dite_rmse_rand <- numeric(0)

ate_rmse_vanilla <- numeric(0)
ate_rmse_oracle <- numeric(0)
ate_rmse_ps <- numeric(0)
ate_rmse_psglm <- numeric(0)
ate_rmse_rand <- numeric(0)
date_rmse_vanilla <- numeric(0)
date_rmse_oracle <- numeric(0)
date_rmse_ps <- numeric(0)
date_rmse_psglm <- numeric(0)
date_rmse_rand <- numeric(0)

vanillaVar <- vector("list", length= niters)
oracleVar <- vector("list", length= niters)
fitzVar <- vector("list", length= niters)
psVar <- vector("list", length= niters)
ps3Var <- vector("list", length= niters)
ps4Var <- vector("list", length= niters)

vanillaCount <- vector("list", length= niters)
oracleCount <- vector("list", length= niters)
fitzCount <- vector("list", length= niters)
psCount <- vector("list", length= niters)
ps3Count <- vector("list", length= niters)
ps4Count <- vector("list", length= niters)
dvanillaCount <- vector("list", length= niters)
doracleCount <- vector("list", length= niters)
dfitzCount <- vector("list", length= niters)
dpsCount <- vector("list", length= niters)
dps3Count <- vector("list", length= niters)
dps4Count <- vector("list", length= niters)

rmse_pbart <- numeric(0)
rmse_dpbart <- numeric(0)

pehe_tt_vanilla <- numeric(0)
pehe_tt_oracle <- numeric(0)
pehe_tt_ps <- numeric(0)
pehe_tt_psglm <- numeric(0)
pehe_tt_rand <- numeric(0)
pehe_tt_dvanilla <- numeric(0)
pehe_tt_doracle <- numeric(0)
pehe_tt_dps <- numeric(0)
pehe_tt_dpsglm <- numeric(0)
pehe_tt_drand <- numeric(0)

pehe_tc_vanilla <- numeric(0)
pehe_tc_oracle <- numeric(0)
pehe_tc_ps <- numeric(0)
pehe_tc_psglm <- numeric(0)
pehe_tc_rand <- numeric(0)
pehe_tc_dvanilla <- numeric(0)
pehe_tc_doracle <- numeric(0)
pehe_tc_dps <- numeric(0)
pehe_tc_dpsglm <- numeric(0)
pehe_tc_drand <- numeric(0)

att_rmse_vanilla <- numeric(0)
att_rmse_oracle <- numeric(0)
att_rmse_ps <- numeric(0)
att_rmse_psglm <- numeric(0)
att_rmse_rand <- numeric(0)
datt_rmse_vanilla <- numeric(0)
datt_rmse_oracle <- numeric(0)
datt_rmse_ps <- numeric(0)
datt_rmse_psglm <- numeric(0)
datt_rmse_rand <- numeric(0)

atc_rmse_vanilla <- numeric(0)
atc_rmse_oracle <- numeric(0)
atc_rmse_ps <- numeric(0)
atc_rmse_psglm <- numeric(0)
atc_rmse_rand <- numeric(0)
datc_rmse_vanilla <- numeric(0)
datc_rmse_oracle <- numeric(0)
datc_rmse_ps <- numeric(0)
datc_rmse_psglm <- numeric(0)
datc_rmse_rand <- numeric(0)


if(.Platform$OS.type=="unix"){
  mc.cores.detected <- 1
}


#cluster = makeCluster(ncores, type = "SOCK")
#registerDoSNOW(cluster)

#time <- Sys.time()
#results <- foreach(i=1:niters) %dopar%{
#  library("BART")

for(i in c(1)){
  #setting different seeds for different iterations
  seedaux <- (7565 + i*5)
  
  set.seed(seedaux)
  #Gerando os dados
  p = 98
  n = 1000
  x = t(data.frame(matrix(runif(p*n),n,p)))
  
  #Gerando u[i]
  q = -1*(x[1,]>(x[2,])) + 1*(x[1,]<(x[2,]))
  
  #Gerando tratamento
  pi = pnorm(q)
  z = rbinom(n,1,pi)
  
  #Alpha ? o verdadeiro treatment effect
  alpha = 0.5*(x[3,] > 1/4) + 0.25*(x[3,] > 2/4) + 0.25*(x[3,]>3/4)
  
  #Cria um sigma para o erro
  sigma = diff(range(q + alpha*pi))/8
  
  ###############################
  #####Heterogeneous Effects#####
  ###############################
  
  #Y[i] = q[i] + Z[i]*a[i] + e[i]
  
  #Gerando Y
  mu = q + alpha*z
  
  #Hahn Example
  #y = mu + 0.1*x[1,] + 0.1*x[2,] + sigma*rnorm(n)
  
  #Friedman Example
  y <- 10*sin(pi*x[1,]*x[2,]) + 20*((x[3,] - 0.5)^2) + 10*x[4,] + 5*x[5,] + mu + sigma*rnorm(n)
  
  pihat = pi
  
  
  #Test Matrix
  X = cbind(rbind(cbind(t(x),pihat),cbind(t(x),pihat)),c(rep(1,n),rep(0,n)))
  colnames(X)[ncol(X)]="z"
  
  #VanillaBART
  set.seed(seedaux)
  
  if(.Platform$OS.type=="unix"){
    vanilla = mc.wbart(cbind(t(x),z),y,X[,-(ncol(X)-1)],nskip = burn,ndpost = npost,keepevery = thin, 
                       seed = seedaux, mc.cores = mc.cores.detected)
  }else{
    vanilla = wbart(cbind(t(x),z),y,X[,-(ncol(X)-1)],nskip = burn,ndpost = npost,keepevery = thin)
  }  
  
  #Calculating estimated ATE
  ate_est_vanilla = vanilla$yhat.test[,1:n] - vanilla$yhat.test[,(n+1):(2*n)]
  ite_est_vanilla = colMeans(vanilla$yhat.test[,1:n] - vanilla$yhat.test[,(n+1):(2*n)])
  
  #Calculating RMSE ITE for this iteration
  ite_rmse_vanilla[i] <- rmseITE(ite_est_vanilla,x)
  pehe_tt_vanilla[i] <- rmseITE(ite_est_vanilla[(z==1)],x[,(z==1)])
  pehe_tc_vanilla[i] <- rmseITE(ite_est_vanilla[(z==0)],x[,(z==0)])
  ate_rmse_vanilla[i] <- rmseATE(ate_est_vanilla,alpha)
  att_rmse_vanilla[i] <- rmseATE(ate_est_vanilla[,(z==1)],alpha[(z==1)])
  atc_rmse_vanilla[i] <- rmseATE(ate_est_vanilla[,(z==0)],alpha[(z==0)])
  vanillaCount[[i]] <- joinDF(vanillaCount,PIP(vanilla),i)
  
  #OracleBART
  set.seed(seedaux)
  
  if(.Platform$OS.type=="unix"){
    oracle = mc.wbart(cbind(t(x),pihat,z),y,X,nskip = burn,ndpost = npost,keepevery = thin, 
                      seed = seedaux, mc.cores = mc.cores.detected)
  }else{
    oracle = wbart(cbind(t(x),pihat,z),y,X,nskip = burn,ndpost = npost,keepevery = thin)
  }  
  
  #Calculating estimated ATE
  ate_est_oracle = oracle$yhat.test[,1:n] - oracle$yhat.test[,(n+1):(2*n)]
  ite_est_oracle = colMeans(oracle$yhat.test[,1:n] - oracle$yhat.test[,(n+1):(2*n)])
  
  
  #Calculating RMSE ITE for this iteration
  ite_rmse_oracle[i] <- rmseITE(ite_est_oracle,x)
  pehe_tt_oracle[i] <- rmseITE(ite_est_oracle[(z==1)],x[,(z==1)])
  pehe_tc_oracle[i] <- rmseITE(ite_est_oracle[(z==0)],x[,(z==0)])
  ate_rmse_oracle[i] <- rmseATE(ate_est_oracle,alpha)
  att_rmse_oracle[i] <- rmseATE(ate_est_oracle[,(z==1)],alpha[(z==1)])
  atc_rmse_oracle[i] <- rmseATE(ate_est_oracle[,(z==0)],alpha[(z==0)])
  oracleCount[[i]] <- joinDF(oracleCount,PIP(oracle),i)
  
  #pihat with pbart
  set.seed(seedaux)
  
  if(.Platform$OS.type=="unix"){
    fitz = mc.pbart(cbind(t(x)),z,nskip = burn,ndpost = npost,keepevery = pthin, 
                    seed = seedaux, mc.cores = mc.cores.detected)
  }else{
    fitz = pbart(cbind(t(x)),z,nskip = burn,ndpost = npost,keepevery = pthin)
  }
  
  pihat2 = fitz$prob.train.mean #Could use median instead
  X2 = cbind(rbind(cbind(t(x),pihat2),cbind(t(x),pihat2)),c(rep(1,n),rep(0,n)))
  colnames(X2)[ncol(X2)]="z"
  
  rmse_pbart[i] <- sqrt(mean((pihat-pihat2)^2))
  fitzCount[[i]] <- joinDF(fitzCount,PIP(fitz),i)
  
  
  #PSBART
  set.seed(seedaux)
  
  if(.Platform$OS.type=="unix"){
    ps = mc.wbart(cbind(t(x),pihat2,z), y, X2, nskip = burn,ndpost = npost,keepevery = thin, 
                  seed = seedaux, mc.cores = mc.cores.detected)
  }else{
    ps = wbart(cbind(t(x),pihat2,z), y, X2, nskip = burn,ndpost = npost,keepevery = thin)
  }
  
  #Calculating estimated ATE
  ate_est_ps = ps$yhat.test[,1:n] - ps$yhat.test[,(n+1):(2*n)]
  
  #Calculating estimated ITE
  ite_est_ps = colMeans(ps$yhat.test[,1:n] - ps$yhat.test[,(n+1):(2*n)])
  
  
  #Calculating RMSE ITE for this iteration
  ite_rmse_ps[i] <- rmseITE(ite_est_ps,x)
  pehe_tt_ps[i] <- rmseITE(ite_est_ps[(z==1)],x[,(z==1)])
  pehe_tc_ps[i] <- rmseITE(ite_est_ps[(z==0)],x[,(z==0)])
  ate_rmse_ps[i] <- rmseATE(ate_est_ps,alpha)
  att_rmse_ps[i] <- rmseATE(ate_est_ps[,(z==1)],alpha[(z==1)])
  atc_rmse_ps[i] <- rmseATE(ate_est_ps[,(z==0)],alpha[(z==0)])
  psCount[[i]] <- joinDF(psCount,PIP(ps),i)
  
  
  
  #pihat estimated by GLM
  fitzglm = glm(z ~ cbind(t(x)), family = binomial())
  pihat4 = fitzglm$fitted.values
  X4 = cbind(rbind(cbind(t(x),pihat4),cbind(t(x),pihat4)),c(rep(1,n),rep(0,n)))
  colnames(X4)[ncol(X4)]="z"
  
  #PSGLM
  set.seed(seedaux)
  
  if(.Platform$OS.type=="unix"){
    ps3 = mc.wbart(cbind(t(x),pihat4,z),y,X4,nskip = burn,ndpost = npost,keepevery = thin, 
                   seed = seedaux, mc.cores = mc.cores.detected)
  }else{
    ps3 = wbart(cbind(t(x),pihat4,z),y,X4,nskip = burn,ndpost = npost,keepevery = thin)
  }
  #Estimated Treatment Effects
  ate_est_psglm = ps3$yhat.test[,1:n] - ps3$yhat.test[,(n+1):(2*n)]
  ite_est_psglm = colMeans(ps3$yhat.test[,1:n] - ps3$yhat.test[,(n+1):(2*n)])
  
  
  #Calculating RMSE ITE for this iteration
  ite_rmse_psglm[i] <- rmseITE(ite_est_psglm,x)
  pehe_tt_psglm[i] <- rmseITE(ite_est_psglm[(z==1)],x[,(z==1)])
  pehe_tc_psglm[i] <- rmseITE(ite_est_psglm[(z==0)],x[,(z==0)])
  ate_rmse_psglm[i] <- rmseATE(ate_est_psglm,alpha)
  att_rmse_psglm[i] <- rmseATE(ate_est_psglm[,(z==1)],alpha[(z==1)])
  atc_rmse_psglm[i] <- rmseATE(ate_est_psglm[,(z==0)],alpha[(z==0)])
  ps3Count[[i]] <- joinDF(ps3Count,PIP(ps3),i)
  
  
  #pihat randomly generated
  pihat5 = runif(n)
  X5 = cbind(rbind(cbind(t(x),pihat5),cbind(t(x),pihat5)),c(rep(1,n),rep(0,n)))
  colnames(X5)[ncol(X5)]="z"
  
  #PSRAND
  set.seed(seedaux)
  
  if(.Platform$OS.type=="unix"){
    ps4 = mc.wbart(cbind(t(x),pihat5,z),y,X5,nskip = burn,ndpost = npost,keepevery = thin, 
                   seed = seedaux, mc.cores = mc.cores.detected)
  }else{
    ps4 = wbart(cbind(t(x),pihat5,z),y,X5,nskip = burn,ndpost = npost,keepevery = thin)
  }
  
  #Estimated Treatment Effects
  ate_est_psrand = ps4$yhat.test[,1:n] - ps4$yhat.test[,(n+1):(2*n)]
  ite_est_psrand = colMeans(ps4$yhat.test[,1:n] - ps4$yhat.test[,(n+1):(2*n)])
  
  
  #Calculating RMSE ITE for this iteration
  ite_rmse_rand[i] <- rmseITE(ite_est_psrand,x)
  pehe_tt_rand[i] <- rmseITE(ite_est_psrand[(z==1)],x[,(z==1)])
  pehe_tc_rand[i] <- rmseITE(ite_est_psrand[(z==0)],x[,(z==0)])
  ate_rmse_rand[i] <- rmseATE(ate_est_psrand,alpha)
  att_rmse_rand[i] <- rmseATE(ate_est_psrand[,(z==1)],alpha[(z==1)])
  atc_rmse_rand[i] <- rmseATE(ate_est_psrand[,(z==0)],alpha[(z==0)])
  ps4Count[[i]] <- joinDF(ps4Count,PIP(ps4),i)
  
  
  ####################################
  #BART Models with sparse setting on#
  ####################################
  
  #VanillaBART
  set.seed(seedaux)
  
  if(.Platform$OS.type=="unix"){
    dvanilla = mc.wbart(cbind(t(x),z),y,X[,-(ncol(X)-1)], sparse = T, nskip = burn,ndpost = npost,keepevery = thin, 
                        seed = seedaux, mc.cores = mc.cores.detected)
  }else{
    dvanilla = wbart(cbind(t(x),z),y,X[,-(ncol(X)-1)], sparse = T, nskip = burn,ndpost = npost,keepevery = thin)
  }
  
  #Calculating estimated ATE
  date_est_vanilla = dvanilla$yhat.test[,1:n] - dvanilla$yhat.test[,(n+1):(2*n)]
  dite_est_vanilla = colMeans(dvanilla$yhat.test[,1:n] - dvanilla$yhat.test[,(n+1):(2*n)])
  
  #Calculating RMSE ITE for this iteration
  dite_rmse_vanilla[i] <- rmseITE(dite_est_vanilla,x)
  pehe_tt_dvanilla[i] <- rmseITE(dite_est_vanilla[(z==1)],x[,(z==1)])
  pehe_tc_dvanilla[i] <- rmseITE(dite_est_vanilla[(z==0)],x[,(z==0)])
  date_rmse_vanilla[i] <- rmseATE(date_est_vanilla,alpha)
  datt_rmse_vanilla[i] <- rmseATE(date_est_vanilla[,(z==1)],alpha[(z==1)])
  datc_rmse_vanilla[i] <- rmseATE(date_est_vanilla[,(z==0)],alpha[(z==0)])
  vanillaVar[[i]] <- joinDF(vanillaVar,dvanilla$varprob.mean,i)
  dvanillaCount[[i]] <- joinDF(dvanillaCount,PIP(dvanilla),i)
  
  #OracleBART
  set.seed(seedaux)
  
  if(.Platform$OS.type=="unix"){
    doracle = mc.wbart(cbind(t(x),pihat,z),y,X,sparse = T, nskip = burn,ndpost = npost,keepevery = thin, 
                       seed = seedaux, mc.cores = mc.cores.detected)
  }else{
    doracle = wbart(cbind(t(x),pihat,z),y,X,sparse = T, nskip = burn,ndpost = npost,keepevery = thin)
  }
  
  
  #Calculating estimated ATE
  date_est_oracle = doracle$yhat.test[,1:n] - doracle$yhat.test[,(n+1):(2*n)]
  dite_est_oracle = colMeans(doracle$yhat.test[,1:n] - doracle$yhat.test[,(n+1):(2*n)])
  
  
  #Calculating RMSE ITE for this iteration
  dite_rmse_oracle[i] <- rmseITE(dite_est_oracle,x)
  pehe_tt_doracle[i] <- rmseITE(dite_est_oracle[(z==1)],x[,(z==1)])
  pehe_tc_doracle[i] <- rmseITE(dite_est_oracle[(z==0)],x[,(z==0)])
  date_rmse_oracle[i] <- rmseATE(date_est_oracle,alpha)
  datt_rmse_oracle[i] <- rmseATE(date_est_oracle[,(z==1)],alpha[(z==1)])
  datc_rmse_oracle[i] <- rmseATE(date_est_oracle[,(z==0)],alpha[(z==0)])
  oracleVar[[i]] <- joinDF(oracleVar,doracle$varprob.mean,i)
  doracleCount[[i]] <- joinDF(doracleCount,PIP(doracle),i)
  
  #pihat with pbart
  set.seed(seedaux)
  if(.Platform$OS.type=="unix"){
    dfitz = mc.pbart(cbind(t(x)),z,sparse = T, nskip = burn,ndpost = npost,keepevery = pthin, 
                     seed = seedaux, mc.cores = mc.cores.detected)
  }else{
    dfitz = pbart(cbind(t(x)),z,sparse = T, nskip = burn,ndpost = npost,keepevery = pthin)
  }
  
  dpihat2 = dfitz$prob.train.mean #Could use median instead
  
  #Calculating Number of times Variable was used
  rmse_dpbart[i] <- sqrt(mean((pihat-dpihat2)^2))
  fitzVar <- joinDF(fitzVar,dfitz$varprob.mean,i)
  dfitzCount <- joinDF(dfitzCount,PIP(dfitz),i)
  
  dX2 = cbind(rbind(cbind(t(x),dpihat2),cbind(t(x),dpihat2)),c(rep(1,n),rep(0,n)))
  colnames(dX2)[ncol(dX2)]="z"
  
  #PSBART
  set.seed(seedaux)
  
  if(.Platform$OS.type=="unix"){
    dps = mc.wbart(cbind(t(x),dpihat2,z), y, dX2, sparse = T, nskip = burn,ndpost = npost,keepevery = thin, 
                   seed = seedaux, mc.cores = mc.cores.detected)
  }else{
    dps = wbart(cbind(t(x),dpihat2,z), y, dX2, sparse = T, nskip = burn,ndpost = npost,keepevery = thin)
  }
  
  #Calculating estimated ATE
  date_est_ps = dps$yhat.test[,1:n] - dps$yhat.test[,(n+1):(2*n)]
  
  #Calculating estimated ITE
  dite_est_ps = colMeans(dps$yhat.test[,1:n] - dps$yhat.test[,(n+1):(2*n)])
  
  
  
  #Calculating RMSE ITE for this iteration
  dite_rmse_ps[i] <- rmseITE(dite_est_ps,x)
  pehe_tt_dps[i] <- rmseITE(dite_est_ps[(z==1)],x[,(z==1)])
  pehe_tc_dps[i] <- rmseITE(dite_est_ps[(z==0)],x[,(z==0)])
  date_rmse_ps[i] <- rmseATE(date_est_ps,alpha)
  datt_rmse_ps[i] <- rmseATE(date_est_ps[,(z==1)],alpha[(z==1)])
  datc_rmse_ps[i] <- rmseATE(date_est_ps[,(z==0)],alpha[(z==0)])
  psVar[[i]] <- joinDF(psVar,dps$varprob.mean,i)
  dpsCount[[i]] <- joinDF(dpsCount,PIP(dps),i)
  
  #PSGLM
  set.seed(seedaux)
  
  if(.Platform$OS.type=="unix"){
    dps3 = mc.wbart(cbind(t(x),pihat4,z),y,sparse = T, X4,nskip = burn,ndpost = npost,keepevery = thin, 
                    seed = seedaux, mc.cores = mc.cores.detected)
  }else{
    dps3 = wbart(cbind(t(x),pihat4,z),y,sparse = T, X4,nskip = burn,ndpost = npost,keepevery = thin)
  }
  #Estimated Treatment Effects
  date_est_psglm = dps3$yhat.test[,1:n] - dps3$yhat.test[,(n+1):(2*n)]
  dite_est_psglm = colMeans(dps3$yhat.test[,1:n] - dps3$yhat.test[,(n+1):(2*n)])
  
  
  #Calculating RMSE ITE for this iteration
  dite_rmse_psglm[i] <- rmseITE(dite_est_psglm,x)
  pehe_tt_dpsglm[i] <- rmseITE(dite_est_psglm[(z==1)],x[,(z==1)])
  pehe_tc_dpsglm[i] <- rmseITE(dite_est_psglm[(z==0)],x[,(z==0)])
  date_rmse_psglm[i] <- rmseATE(date_est_psglm,alpha)
  datt_rmse_psglm[i] <- rmseATE(date_est_psglm[,(z==1)],alpha[(z==1)])
  datc_rmse_psglm[i] <- rmseATE(date_est_psglm[,(z==0)],alpha[(z==0)])
  ps3Var[[i]] <- joinDF(ps3Var,dps3$varprob.mean,i)
  dps3Count[[i]] <- joinDF(dps3Count,PIP(dps3),i)
  
  #PSRAND
  set.seed(seedaux)
  
  if(.Platform$OS.type=="unix"){
    dps4 = mc.wbart(cbind(t(x),pihat5,z),y,X5,sparse = T, nskip = burn,ndpost = npost,keepevery = thin, 
                    seed = seedaux, mc.cores = mc.cores.detected)
  }else{
    dps4 = wbart(cbind(t(x),pihat5,z),y,X5,sparse = T, nskip = burn,ndpost = npost,keepevery = thin)
  }
  
  #Estimated Treatment Effects
  date_est_psrand = dps4$yhat.test[,1:n] - dps4$yhat.test[,(n+1):(2*n)]
  dite_est_psrand = colMeans(dps4$yhat.test[,1:n] - dps4$yhat.test[,(n+1):(2*n)])
  
  #Calculating RMSE ITE for this iteration
  dite_rmse_rand[i] <- rmseITE(dite_est_psrand,x)
  pehe_tt_drand[i] <- rmseITE(dite_est_psrand[(z==1)],x[,(z==1)])
  pehe_tc_drand[i] <- rmseITE(dite_est_psrand[(z==0)],x[,(z==0)])
  date_rmse_rand[i] <- rmseATE(date_est_psrand,alpha)
  datt_rmse_rand[i] <- rmseATE(date_est_psrand[,(z==1)],alpha[(z==1)])
  datc_rmse_rand[i] <- rmseATE(date_est_psrand[,(z==0)],alpha[(z==0)])
  ps4Var[[i]] <- joinDF(ps4Var,dps4$varprob.mean,i)
  dps4Count[[i]] <- joinDF(dps4Count,PIP(dps4),i)
  
  ########
  #Saving#
  ########
  
  aux <- list(ite_rmse_vanilla[i],
              ite_rmse_oracle[i],
              ite_rmse_ps[i],
              ite_rmse_psglm[i],
              ite_rmse_rand[i],
              dite_rmse_vanilla[i],
              dite_rmse_oracle[i],
              dite_rmse_ps[i],
              dite_rmse_psglm[i],
              dite_rmse_rand[i],
              ate_rmse_vanilla[i],
              ate_rmse_oracle[i],
              ate_rmse_ps[i],
              ate_rmse_psglm[i],
              ate_rmse_rand[i],
              date_rmse_vanilla[i],
              date_rmse_oracle[i],
              date_rmse_ps[i],
              date_rmse_psglm[i],
              date_rmse_rand[i],
              vanillaVar[[i]],
              oracleVar[[i]],
              fitzVar[[i]],
              psVar[[i]],
              ps3Var[[i]],
              ps4Var[[i]],
              vanillaCount[[i]],
              oracleCount[[i]],
              fitzCount[[i]],
              psCount[[i]],
              ps3Count[[i]],
              ps4Count[[i]],
              dvanillaCount[[i]],
              doracleCount[[i]],
              dfitzCount[[i]],
              dpsCount[[i]],
              dps3Count[[i]],
              dps4Count[[i]],
              rmse_pbart[i],
              pehe_tt_vanilla[i],
              pehe_tt_oracle[i],
              pehe_tt_ps[i],
              pehe_tt_psglm[i],
              pehe_tt_rand[i],
              pehe_tt_dvanilla[i],
              pehe_tt_doracle[i],
              pehe_tt_dps[i],
              pehe_tt_dpsglm[i],
              pehe_tt_drand[i],
              pehe_tc_vanilla[i],
              pehe_tc_oracle[i],
              pehe_tc_ps[i],
              pehe_tc_psglm[i],
              pehe_tc_rand[i],
              pehe_tc_dvanilla[i],
              pehe_tc_doracle[i],
              pehe_tc_dps[i],
              pehe_tc_dpsglm[i],
              pehe_tc_drand[i],
              att_rmse_vanilla[i],
              att_rmse_oracle[i],
              att_rmse_ps[i],
              att_rmse_psglm[i],
              att_rmse_rand[i],
              datt_rmse_vanilla[i],
              datt_rmse_oracle[i],
              datt_rmse_ps[i],
              datt_rmse_psglm[i],
              datt_rmse_rand[i],
              atc_rmse_vanilla[i],
              atc_rmse_oracle[i],
              atc_rmse_ps[i],
              atc_rmse_psglm[i],
              atc_rmse_rand[i],
              datc_rmse_vanilla[i],
              datc_rmse_oracle[i],
              datc_rmse_ps[i],
              datc_rmse_psglm[i],
              datc_rmse_rand[i],
              rmse_dpbart[i])
  
#  save(aux, file = paste0("results_sim_new_fried",i,".RData"))
  
  gc()
  
  return(aux)
  
}

#save(results, file = paste0("results_sim_new_fried.RData"))

#stopCluster(cluster)


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

ice.bart.vanilla.x1 <- ice(vanilla, Xice, predictor = "X1", 
                           predictfcn = pred.BART.ATE)
save(ice.bart.vanilla.x1, file="ice_vanilla_x1.RData")

ice.bart.vanilla.x1.0025 <- ice(vanilla, Xice, predictor = "X1", 
                           predictfcn = pred.BART.ATE.0025)
save(ice.bart.vanilla.x1.0025, file="ice_vanilla_x1.0025.RData")

ice.bart.vanilla.x1.0975 <- ice(vanilla, Xice, predictor = "X1", 
                                predictfcn = pred.BART.ATE.0975)
save(ice.bart.vanilla.x1.0975, file="ice_vanilla_x1.0975.RData")

####

ice.bart.vanilla.x2 <- ice(vanilla, Xice, predictor = "X2", 
                           predictfcn = pred.BART.ATE)
save(ice.bart.vanilla.x2, file="ice_vanilla_x2.RData")

ice.bart.vanilla.x2.0025 <- ice(vanilla, Xice, predictor = "X2", 
                           predictfcn = pred.BART.ATE.0025)
save(ice.bart.vanilla.x2.0025, file="ice_vanilla_x2.0025.RData")

ice.bart.vanilla.x2.0975 <- ice(vanilla, Xice, predictor = "X2", 
                                predictfcn = pred.BART.ATE.0975)
save(ice.bart.vanilla.x2.0975, file="ice_vanilla_x2.0975.RData")

###

ice.bart.vanilla.x3 <- ice(vanilla, Xice, predictor = "X3", 
                           predictfcn = pred.BART.ATE)
save(ice.bart.vanilla.x3, file="ice_vanilla_x3.RData")

ice.bart.vanilla.x3.0025 <- ice(vanilla, Xice, predictor = "X3", 
                           predictfcn = pred.BART.ATE.0025)
save(ice.bart.vanilla.x3.0025, file="ice_vanilla_x3.0025.RData")

ice.bart.vanilla.x3.0975 <- ice(vanilla, Xice, predictor = "X3", 
                                predictfcn = pred.BART.ATE.0975)
save(ice.bart.vanilla.x3.0975, file="ice_vanilla_x3.0975.RData")

###

ice.bart.vanilla.x4 <- ice(vanilla, Xice, predictor = "X4", 
                           predictfcn = pred.BART.ATE)
save(ice.bart.vanilla.x4, file="ice_vanilla_x4.RData")

ice.bart.vanilla.x4.0025 <- ice(vanilla, Xice, predictor = "X4", 
                                predictfcn = pred.BART.ATE.0025)
save(ice.bart.vanilla.x4.0025, file="ice_vanilla_x4.0025.RData")

ice.bart.vanilla.x4.0975 <- ice(vanilla, Xice, predictor = "X4", 
                                predictfcn = pred.BART.ATE.0975)
save(ice.bart.vanilla.x4.0975, file="ice_vanilla_x4.0975.RData")

###

ice.bart.vanilla.x5 <- ice(vanilla, Xice, predictor = "X5", 
                           predictfcn = pred.BART.ATE)
save(ice.bart.vanilla.x5, file="ice_vanilla_x5.RData")

ice.bart.vanilla.x5.0025 <- ice(vanilla, Xice, predictor = "X5", 
                                predictfcn = pred.BART.ATE.0025)
save(ice.bart.vanilla.x5.0025, file="ice_vanilla_x5.0025.RData")

ice.bart.vanilla.x5.0975 <- ice(vanilla, Xice, predictor = "X5", 
                                predictfcn = pred.BART.ATE.0975)
save(ice.bart.vanilla.x5.0975, file="ice_vanilla_x5.0975.RData")

###

#Oracle
Xice = cbind(t(x),pihat,z)

ice.bart.oracle.x1 <- ice(oracle, Xice, predictor = "X1", 
                           predictfcn = pred.BART.ATE)
save(ice.bart.oracle.x1, file="ice_oracle_x1.RData")

ice.bart.oracle.x1.0025 <- ice(oracle, Xice, predictor = "X1", 
                                predictfcn = pred.BART.ATE.0025)
save(ice.bart.oracle.x1.0025, file="ice_oracle_x1.0025.RData")

ice.bart.oracle.x1.0975 <- ice(oracle, Xice, predictor = "X1", 
                                predictfcn = pred.BART.ATE.0975)
save(ice.bart.oracle.x1.0975, file="ice_oracle_x1.0975.RData")

####

ice.bart.oracle.x2 <- ice(oracle, Xice, predictor = "X2", 
                           predictfcn = pred.BART.ATE)
save(ice.bart.oracle.x2, file="ice_oracle_x2.RData")

ice.bart.oracle.x2.0025 <- ice(oracle, Xice, predictor = "X2", 
                                predictfcn = pred.BART.ATE.0025)
save(ice.bart.oracle.x2.0025, file="ice_oracle_x2.0025.RData")

ice.bart.oracle.x2.0975 <- ice(oracle, Xice, predictor = "X2", 
                                predictfcn = pred.BART.ATE.0975)
save(ice.bart.oracle.x2.0975, file="ice_oracle_x2.0975.RData")

###

ice.bart.oracle.x3 <- ice(oracle, Xice, predictor = "X3", 
                           predictfcn = pred.BART.ATE)
save(ice.bart.oracle.x3, file="ice_oracle_x3.RData")

ice.bart.oracle.x3.0025 <- ice(oracle, Xice, predictor = "X3", 
                                predictfcn = pred.BART.ATE.0025)
save(ice.bart.oracle.x3.0025, file="ice_oracle_x3.0025.RData")

ice.bart.oracle.x3.0975 <- ice(oracle, Xice, predictor = "X3", 
                                predictfcn = pred.BART.ATE.0975)
save(ice.bart.oracle.x3.0975, file="ice_oracle_x3.0975.RData")

###

ice.bart.oracle.x4 <- ice(oracle, Xice, predictor = "X4", 
                           predictfcn = pred.BART.ATE)
save(ice.bart.oracle.x4, file="ice_oracle_x4.RData")

ice.bart.oracle.x4.0025 <- ice(oracle, Xice, predictor = "X4", 
                                predictfcn = pred.BART.ATE.0025)
save(ice.bart.oracle.x4.0025, file="ice_oracle_x4.0025.RData")

ice.bart.oracle.x4.0975 <- ice(oracle, Xice, predictor = "X4", 
                                predictfcn = pred.BART.ATE.0975)
save(ice.bart.oracle.x4.0975, file="ice_oracle_x4.0975.RData")

###

ice.bart.oracle.x5 <- ice(oracle, Xice, predictor = "X5", 
                           predictfcn = pred.BART.ATE)
save(ice.bart.oracle.x5, file="ice_oracle_x5.RData")

ice.bart.oracle.x5.0025 <- ice(oracle, Xice, predictor = "X5", 
                                predictfcn = pred.BART.ATE.0025)
save(ice.bart.oracle.x5.0025, file="ice_oracle_x5.0025.RData")

ice.bart.oracle.x5.0975 <- ice(oracle, Xice, predictor = "X5", 
                                predictfcn = pred.BART.ATE.0975)
save(ice.bart.oracle.x5.0975, file="ice_oracle_x5.0975.RData")

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

ice.bart.ps.x1 <- ice(ps, Xice, predictor = "X1", 
                          predictfcn = pred.BART.ATE)
save(ice.bart.ps.x1, file="ice_ps_x1.RData")

ice.bart.ps.x1.0025 <- ice(ps, Xice, predictor = "X1", 
                               predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps.x1.0025, file="ice_ps_x1.0025.RData")

ice.bart.ps.x1.0975 <- ice(ps, Xice, predictor = "X1", 
                               predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps.x1.0975, file="ice_ps_x1.0975.RData")

####

ice.bart.ps.x2 <- ice(ps, Xice, predictor = "X2", 
                          predictfcn = pred.BART.ATE)
save(ice.bart.ps.x2, file="ice_ps_x2.RData")

ice.bart.ps.x2.0025 <- ice(ps, Xice, predictor = "X2", 
                               predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps.x2.0025, file="ice_ps_x2.0025.RData")

ice.bart.ps.x2.0975 <- ice(ps, Xice, predictor = "X2", 
                               predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps.x2.0975, file="ice_ps_x2.0975.RData")

###

ice.bart.ps.x3 <- ice(ps, Xice, predictor = "X3", 
                          predictfcn = pred.BART.ATE)
save(ice.bart.ps.x3, file="ice_ps_x3.RData")

ice.bart.ps.x3.0025 <- ice(ps, Xice, predictor = "X3", 
                               predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps.x3.0025, file="ice_ps_x3.0025.RData")

ice.bart.ps.x3.0975 <- ice(ps, Xice, predictor = "X3", 
                               predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps.x3.0975, file="ice_ps_x3.0975.RData")

###

ice.bart.ps.x4 <- ice(ps, Xice, predictor = "X4", 
                          predictfcn = pred.BART.ATE)
save(ice.bart.ps.x4, file="ice_ps_x4.RData")

ice.bart.ps.x4.0025 <- ice(ps, Xice, predictor = "X4", 
                               predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps.x4.0025, file="ice_ps_x4.0025.RData")

ice.bart.ps.x4.0975 <- ice(ps, Xice, predictor = "X4", 
                               predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps.x4.0975, file="ice_ps_x4.0975.RData")

###

ice.bart.ps.x5 <- ice(ps, Xice, predictor = "X5", 
                          predictfcn = pred.BART.ATE)
save(ice.bart.ps.x5, file="ice_ps_x5.RData")

ice.bart.ps.x5.0025 <- ice(ps, Xice, predictor = "X5", 
                               predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps.x5.0025, file="ice_ps_x5.0025.RData")

ice.bart.ps.x5.0975 <- ice(ps, Xice, predictor = "X5", 
                               predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps.x5.0975, file="ice_ps_x5.0975.RData")

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

ice.bart.ps3.x1 <- ice(ps3, Xice, predictor = "X1", 
                      predictfcn = pred.BART.ATE)
save(ice.bart.ps3.x1, file="ice_ps3_x1.RData")

ice.bart.ps3.x1.0025 <- ice(ps3, Xice, predictor = "X1", 
                           predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps3.x1.0025, file="ice_ps3_x1.0025.RData")

ice.bart.ps3.x1.0975 <- ice(ps3, Xice, predictor = "X1", 
                           predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps3.x1.0975, file="ice_ps3_x1.0975.RData")

####

ice.bart.ps3.x2 <- ice(ps3, Xice, predictor = "X2", 
                      predictfcn = pred.BART.ATE)
save(ice.bart.ps3.x2, file="ice_ps3_x2.RData")

ice.bart.ps3.x2.0025 <- ice(ps3, Xice, predictor = "X2", 
                           predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps3.x2.0025, file="ice_ps3_x2.0025.RData")

ice.bart.ps3.x2.0975 <- ice(ps3, Xice, predictor = "X2", 
                           predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps3.x2.0975, file="ice_ps3_x2.0975.RData")

###

ice.bart.ps3.x3 <- ice(ps3, Xice, predictor = "X3", 
                      predictfcn = pred.BART.ATE)
save(ice.bart.ps3.x3, file="ice_ps3_x3.RData")

ice.bart.ps3.x3.0025 <- ice(ps3, Xice, predictor = "X3", 
                           predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps3.x3.0025, file="ice_ps3_x3.0025.RData")

ice.bart.ps3.x3.0975 <- ice(ps3, Xice, predictor = "X3", 
                           predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps3.x3.0975, file="ice_ps3_x3.0975.RData")

###

ice.bart.ps3.x4 <- ice(ps3, Xice, predictor = "X4", 
                      predictfcn = pred.BART.ATE)
save(ice.bart.ps3.x4, file="ice_ps3_x4.RData")

ice.bart.ps3.x4.0025 <- ice(ps3, Xice, predictor = "X4", 
                           predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps3.x4.0025, file="ice_ps3_x4.0025.RData")

ice.bart.ps3.x4.0975 <- ice(ps3, Xice, predictor = "X4", 
                           predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps3.x4.0975, file="ice_ps3_x4.0975.RData")

###

ice.bart.ps3.x5 <- ice(ps3, Xice, predictor = "X5", 
                      predictfcn = pred.BART.ATE)
save(ice.bart.ps3.x5, file="ice_ps3_x5.RData")

ice.bart.ps3.x5.0025 <- ice(ps3, Xice, predictor = "X5", 
                           predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps3.x5.0025, file="ice_ps3_x5.0025.RData")

ice.bart.ps3.x5.0975 <- ice(ps3, Xice, predictor = "X5", 
                           predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps3.x5.0975, file="ice_ps3_x5.0975.RData")

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

ice.bart.ps4.x1 <- ice(ps4, Xice, predictor = "X1", 
                       predictfcn = pred.BART.ATE)
save(ice.bart.ps4.x1, file="ice_ps4_x1.RData")

ice.bart.ps4.x1.0025 <- ice(ps4, Xice, predictor = "X1", 
                            predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps4.x1.0025, file="ice_ps4_x1.0025.RData")

ice.bart.ps4.x1.0975 <- ice(ps4, Xice, predictor = "X1", 
                            predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps4.x1.0975, file="ice_ps4_x1.0975.RData")

####

ice.bart.ps4.x2 <- ice(ps4, Xice, predictor = "X2", 
                       predictfcn = pred.BART.ATE)
save(ice.bart.ps4.x2, file="ice_ps4_x2.RData")

ice.bart.ps4.x2.0025 <- ice(ps4, Xice, predictor = "X2", 
                            predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps4.x2.0025, file="ice_ps4_x2.0025.RData")

ice.bart.ps4.x2.0975 <- ice(ps4, Xice, predictor = "X2", 
                            predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps4.x2.0975, file="ice_ps4_x2.0975.RData")

###

ice.bart.ps4.x3 <- ice(ps4, Xice, predictor = "X3", 
                       predictfcn = pred.BART.ATE)
save(ice.bart.ps4.x3, file="ice_ps4_x3.RData")

ice.bart.ps4.x3.0025 <- ice(ps4, Xice, predictor = "X3", 
                            predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps4.x3.0025, file="ice_ps4_x3.0025.RData")

ice.bart.ps4.x3.0975 <- ice(ps4, Xice, predictor = "X3", 
                            predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps4.x3.0975, file="ice_ps4_x3.0975.RData")

###

ice.bart.ps4.x4 <- ice(ps4, Xice, predictor = "X4", 
                       predictfcn = pred.BART.ATE)
save(ice.bart.ps4.x4, file="ice_ps4_x4.RData")

ice.bart.ps4.x4.0025 <- ice(ps4, Xice, predictor = "X4", 
                            predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps4.x4.0025, file="ice_ps4_x4.0025.RData")

ice.bart.ps4.x4.0975 <- ice(ps4, Xice, predictor = "X4", 
                            predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps4.x4.0975, file="ice_ps4_x4.0975.RData")

###

ice.bart.ps4.x5 <- ice(ps4, Xice, predictor = "X5", 
                       predictfcn = pred.BART.ATE)
save(ice.bart.ps4.x5, file="ice_ps4_x5.RData")

ice.bart.ps4.x5.0025 <- ice(ps4, Xice, predictor = "X5", 
                            predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps4.x5.0025, file="ice_ps4_x5.0025.RData")

ice.bart.ps4.x5.0975 <- ice(ps4, Xice, predictor = "X5", 
                            predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps4.x5.0975, file="ice_ps4_x5.0975.RData")

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

#DART

#Vanilla
Xice = cbind(t(x),z)

ice.bart.dvanilla.x1 <- ice(dvanilla, Xice, predictor = "X1", 
                           predictfcn = pred.BART.ATE)
save(ice.bart.dvanilla.x1, file="ice_dvanilla_x1.RData")

ice.bart.dvanilla.x1.0025 <- ice(dvanilla, Xice, predictor = "X1", 
                                predictfcn = pred.BART.ATE.0025)
save(ice.bart.dvanilla.x1.0025, file="ice_dvanilla_x1.0025.RData")

ice.bart.dvanilla.x1.0975 <- ice(dvanilla, Xice, predictor = "X1", 
                                predictfcn = pred.BART.ATE.0975)
save(ice.bart.dvanilla.x1.0975, file="ice_dvanilla_x1.0975.RData")

####

ice.bart.dvanilla.x2 <- ice(dvanilla, Xice, predictor = "X2", 
                           predictfcn = pred.BART.ATE)
save(ice.bart.dvanilla.x2, file="ice_dvanilla_x2.RData")

ice.bart.dvanilla.x2.0025 <- ice(dvanilla, Xice, predictor = "X2", 
                                predictfcn = pred.BART.ATE.0025)
save(ice.bart.dvanilla.x2.0025, file="ice_dvanilla_x2.0025.RData")

ice.bart.dvanilla.x2.0975 <- ice(dvanilla, Xice, predictor = "X2", 
                                predictfcn = pred.BART.ATE.0975)
save(ice.bart.dvanilla.x2.0975, file="ice_dvanilla_x2.0975.RData")

###

ice.bart.dvanilla.x3 <- ice(dvanilla, Xice, predictor = "X3", 
                           predictfcn = pred.BART.ATE)
save(ice.bart.dvanilla.x3, file="ice_dvanilla_x3.RData")

ice.bart.dvanilla.x3.0025 <- ice(dvanilla, Xice, predictor = "X3", 
                                predictfcn = pred.BART.ATE.0025)
save(ice.bart.dvanilla.x3.0025, file="ice_dvanilla_x3.0025.RData")

ice.bart.dvanilla.x3.0975 <- ice(dvanilla, Xice, predictor = "X3", 
                                predictfcn = pred.BART.ATE.0975)
save(ice.bart.dvanilla.x3.0975, file="ice_dvanilla_x3.0975.RData")

###

ice.bart.dvanilla.x4 <- ice(dvanilla, Xice, predictor = "X4", 
                           predictfcn = pred.BART.ATE)
save(ice.bart.dvanilla.x4, file="ice_dvanilla_x4.RData")

ice.bart.dvanilla.x4.0025 <- ice(dvanilla, Xice, predictor = "X4", 
                                predictfcn = pred.BART.ATE.0025)
save(ice.bart.dvanilla.x4.0025, file="ice_dvanilla_x4.0025.RData")

ice.bart.dvanilla.x4.0975 <- ice(dvanilla, Xice, predictor = "X4", 
                                predictfcn = pred.BART.ATE.0975)
save(ice.bart.dvanilla.x4.0975, file="ice_dvanilla_x4.0975.RData")

###

ice.bart.dvanilla.x5 <- ice(dvanilla, Xice, predictor = "X5", 
                           predictfcn = pred.BART.ATE)
save(ice.bart.dvanilla.x5, file="ice_dvanilla_x5.RData")

ice.bart.dvanilla.x5.0025 <- ice(dvanilla, Xice, predictor = "X5", 
                                predictfcn = pred.BART.ATE.0025)
save(ice.bart.dvanilla.x5.0025, file="ice_dvanilla_x5.0025.RData")

ice.bart.dvanilla.x5.0975 <- ice(dvanilla, Xice, predictor = "X5", 
                                predictfcn = pred.BART.ATE.0975)
save(ice.bart.dvanilla.x5.0975, file="ice_dvanilla_x5.0975.RData")

###

#Oracle
Xice = cbind(t(x),pihat,z)

ice.bart.doracle.x1 <- ice(doracle, Xice, predictor = "X1", 
                          predictfcn = pred.BART.ATE)
save(ice.bart.doracle.x1, file="ice_doracle_x1.RData")

ice.bart.doracle.x1.0025 <- ice(doracle, Xice, predictor = "X1", 
                               predictfcn = pred.BART.ATE.0025)
save(ice.bart.doracle.x1.0025, file="ice_doracle_x1.0025.RData")

ice.bart.doracle.x1.0975 <- ice(doracle, Xice, predictor = "X1", 
                               predictfcn = pred.BART.ATE.0975)
save(ice.bart.doracle.x1.0975, file="ice_doracle_x1.0975.RData")

####

ice.bart.doracle.x2 <- ice(doracle, Xice, predictor = "X2", 
                          predictfcn = pred.BART.ATE)
save(ice.bart.doracle.x2, file="ice_doracle_x2.RData")

ice.bart.doracle.x2.0025 <- ice(doracle, Xice, predictor = "X2", 
                               predictfcn = pred.BART.ATE.0025)
save(ice.bart.doracle.x2.0025, file="ice_doracle_x2.0025.RData")

ice.bart.doracle.x2.0975 <- ice(doracle, Xice, predictor = "X2", 
                               predictfcn = pred.BART.ATE.0975)
save(ice.bart.doracle.x2.0975, file="ice_doracle_x2.0975.RData")

###

ice.bart.doracle.x3 <- ice(doracle, Xice, predictor = "X3", 
                          predictfcn = pred.BART.ATE)
save(ice.bart.doracle.x3, file="ice_doracle_x3.RData")

ice.bart.doracle.x3.0025 <- ice(doracle, Xice, predictor = "X3", 
                               predictfcn = pred.BART.ATE.0025)
save(ice.bart.doracle.x3.0025, file="ice_doracle_x3.0025.RData")

ice.bart.doracle.x3.0975 <- ice(doracle, Xice, predictor = "X3", 
                               predictfcn = pred.BART.ATE.0975)
save(ice.bart.doracle.x3.0975, file="ice_doracle_x3.0975.RData")

###

ice.bart.doracle.x4 <- ice(doracle, Xice, predictor = "X4", 
                          predictfcn = pred.BART.ATE)
save(ice.bart.doracle.x4, file="ice_doracle_x4.RData")

ice.bart.doracle.x4.0025 <- ice(doracle, Xice, predictor = "X4", 
                               predictfcn = pred.BART.ATE.0025)
save(ice.bart.doracle.x4.0025, file="ice_doracle_x4.0025.RData")

ice.bart.doracle.x4.0975 <- ice(doracle, Xice, predictor = "X4", 
                               predictfcn = pred.BART.ATE.0975)
save(ice.bart.doracle.x4.0975, file="ice_doracle_x4.0975.RData")

###

ice.bart.doracle.x5 <- ice(doracle, Xice, predictor = "X5", 
                          predictfcn = pred.BART.ATE)
save(ice.bart.doracle.x5, file="ice_doracle_x5.RData")

ice.bart.doracle.x5.0025 <- ice(doracle, Xice, predictor = "X5", 
                               predictfcn = pred.BART.ATE.0025)
save(ice.bart.doracle.x5.0025, file="ice_doracle_x5.0025.RData")

ice.bart.doracle.x5.0975 <- ice(doracle, Xice, predictor = "X5", 
                               predictfcn = pred.BART.ATE.0975)
save(ice.bart.doracle.x5.0975, file="ice_doracle_x5.0975.RData")

###

ice.bart.doracle.pihat <- ice(doracle, Xice, predictor = "pihat", 
                             predictfcn = pred.BART.ATE)
save(ice.bart.doracle.pihat, file="ice_doracle_pihat.RData")

ice.bart.doracle.pihat.0025 <- ice(doracle, Xice, predictor = "pihat", 
                                  predictfcn = pred.BART.ATE.0025)
save(ice.bart.doracle.pihat.0025, file="ice_doracle_pihat.0025.RData")

ice.bart.doracle.pihat.0975 <- ice(doracle, Xice, predictor = "pihat", 
                                  predictfcn = pred.BART.ATE.0975)
save(ice.bart.doracle.pihat.0975, file="ice_doracle_pihat.0975.RData")

###


#PSBART
Xice = cbind(t(x),dpihat2,z)

ice.bart.dps.x1 <- ice(dps, Xice, predictor = "X1", 
                      predictfcn = pred.BART.ATE)
save(ice.bart.dps.x1, file="ice_dps_x1.RData")

ice.bart.dps.x1.0025 <- ice(dps, Xice, predictor = "X1", 
                           predictfcn = pred.BART.ATE.0025)
save(ice.bart.dps.x1.0025, file="ice_dps_x1.0025.RData")

ice.bart.dps.x1.0975 <- ice(dps, Xice, predictor = "X1", 
                           predictfcn = pred.BART.ATE.0975)
save(ice.bart.dps.x1.0975, file="ice_dps_x1.0975.RData")

####

ice.bart.dps.x2 <- ice(dps, Xice, predictor = "X2", 
                      predictfcn = pred.BART.ATE)
save(ice.bart.dps.x2, file="ice_dps_x2.RData")

ice.bart.dps.x2.0025 <- ice(dps, Xice, predictor = "X2", 
                           predictfcn = pred.BART.ATE.0025)
save(ice.bart.dps.x2.0025, file="ice_dps_x2.0025.RData")

ice.bart.dps.x2.0975 <- ice(dps, Xice, predictor = "X2", 
                           predictfcn = pred.BART.ATE.0975)
save(ice.bart.dps.x2.0975, file="ice_dps_x2.0975.RData")

###

ice.bart.dps.x3 <- ice(dps, Xice, predictor = "X3", 
                      predictfcn = pred.BART.ATE)
save(ice.bart.dps.x3, file="ice_dps_x3.RData")

ice.bart.dps.x3.0025 <- ice(dps, Xice, predictor = "X3", 
                           predictfcn = pred.BART.ATE.0025)
save(ice.bart.dps.x3.0025, file="ice_dps_x3.0025.RData")

ice.bart.dps.x3.0975 <- ice(dps, Xice, predictor = "X3", 
                           predictfcn = pred.BART.ATE.0975)
save(ice.bart.dps.x3.0975, file="ice_dps_x3.0975.RData")

###

ice.bart.dps.x4 <- ice(dps, Xice, predictor = "X4", 
                      predictfcn = pred.BART.ATE)
save(ice.bart.dps.x4, file="ice_dps_x4.RData")

ice.bart.dps.x4.0025 <- ice(dps, Xice, predictor = "X4", 
                           predictfcn = pred.BART.ATE.0025)
save(ice.bart.dps.x4.0025, file="ice_dps_x4.0025.RData")

ice.bart.dps.x4.0975 <- ice(dps, Xice, predictor = "X4", 
                           predictfcn = pred.BART.ATE.0975)
save(ice.bart.dps.x4.0975, file="ice_dps_x4.0975.RData")

###

ice.bart.dps.x5 <- ice(dps, Xice, predictor = "X5", 
                      predictfcn = pred.BART.ATE)
save(ice.bart.dps.x5, file="ice_dps_x5.RData")

ice.bart.dps.x5.0025 <- ice(dps, Xice, predictor = "X5", 
                           predictfcn = pred.BART.ATE.0025)
save(ice.bart.dps.x5.0025, file="ice_dps_x5.0025.RData")

ice.bart.dps.x5.0975 <- ice(dps, Xice, predictor = "X5", 
                           predictfcn = pred.BART.ATE.0975)
save(ice.bart.dps.x5.0975, file="ice_dps_x5.0975.RData")

###

ice.bart.dps.pihat <- ice(dps, Xice, predictor = "dpihat2", 
                         predictfcn = pred.BART.ATE)
save(ice.bart.dps.pihat, file="ice_dps_pihat.RData")

ice.bart.dps.pihat.0025 <- ice(dps, Xice, predictor = "dpihat2", 
                              predictfcn = pred.BART.ATE.0025)
save(ice.bart.dps.pihat.0025, file="ice_dps_pihat.0025.RData")

ice.bart.dps.pihat.0975 <- ice(dps, Xice, predictor = "dpihat2", 
                              predictfcn = pred.BART.ATE.0975)
save(ice.bart.dps.pihat.0975, file="ice_dps_pihat.0975.RData")

###



#GLM-BART
Xice = cbind(t(x),pihat4,z)

ice.bart.dps3.x1 <- ice(dps3, Xice, predictor = "X1", 
                       predictfcn = pred.BART.ATE)
save(ice.bart.dps3.x1, file="ice_dps3_x1.RData")

ice.bart.dps3.x1.0025 <- ice(dps3, Xice, predictor = "X1", 
                            predictfcn = pred.BART.ATE.0025)
save(ice.bart.dps3.x1.0025, file="ice_dps3_x1.0025.RData")

ice.bart.dps3.x1.0975 <- ice(dps3, Xice, predictor = "X1", 
                            predictfcn = pred.BART.ATE.0975)
save(ice.bart.dps3.x1.0975, file="ice_dps3_x1.0975.RData")

####

ice.bart.dps3.x2 <- ice(dps3, Xice, predictor = "X2", 
                       predictfcn = pred.BART.ATE)
save(ice.bart.dps3.x2, file="ice_dps3_x2.RData")

ice.bart.dps3.x2.0025 <- ice(dps3, Xice, predictor = "X2", 
                            predictfcn = pred.BART.ATE.0025)
save(ice.bart.dps3.x2.0025, file="ice_dps3_x2.0025.RData")

ice.bart.dps3.x2.0975 <- ice(dps3, Xice, predictor = "X2", 
                            predictfcn = pred.BART.ATE.0975)
save(ice.bart.dps3.x2.0975, file="ice_dps3_x2.0975.RData")

###

ice.bart.dps3.x3 <- ice(dps3, Xice, predictor = "X3", 
                       predictfcn = pred.BART.ATE)
save(ice.bart.dps3.x3, file="ice_dps3_x3.RData")

ice.bart.dps3.x3.0025 <- ice(dps3, Xice, predictor = "X3", 
                            predictfcn = pred.BART.ATE.0025)
save(ice.bart.dps3.x3.0025, file="ice_dps3_x3.0025.RData")

ice.bart.dps3.x3.0975 <- ice(dps3, Xice, predictor = "X3", 
                            predictfcn = pred.BART.ATE.0975)
save(ice.bart.dps3.x3.0975, file="ice_dps3_x3.0975.RData")

###

ice.bart.dps3.x4 <- ice(dps3, Xice, predictor = "X4", 
                       predictfcn = pred.BART.ATE)
save(ice.bart.dps3.x4, file="ice_dps3_x4.RData")

ice.bart.dps3.x4.0025 <- ice(dps3, Xice, predictor = "X4", 
                            predictfcn = pred.BART.ATE.0025)
save(ice.bart.dps3.x4.0025, file="ice_dps3_x4.0025.RData")

ice.bart.dps3.x4.0975 <- ice(dps3, Xice, predictor = "X4", 
                            predictfcn = pred.BART.ATE.0975)
save(ice.bart.dps3.x4.0975, file="ice_dps3_x4.0975.RData")

###

ice.bart.dps3.x5 <- ice(dps3, Xice, predictor = "X5", 
                       predictfcn = pred.BART.ATE)
save(ice.bart.dps3.x5, file="ice_dps3_x5.RData")

ice.bart.dps3.x5.0025 <- ice(dps3, Xice, predictor = "X5", 
                            predictfcn = pred.BART.ATE.0025)
save(ice.bart.dps3.x5.0025, file="ice_dps3_x5.0025.RData")

ice.bart.dps3.x5.0975 <- ice(dps3, Xice, predictor = "X5", 
                            predictfcn = pred.BART.ATE.0975)
save(ice.bart.dps3.x5.0975, file="ice_dps3_x5.0975.RData")

###

ice.bart.dps3.pihat <- ice(dps3, Xice, predictor = "pihat4", 
                          predictfcn = pred.BART.ATE)
save(ice.bart.dps3.pihat, file="ice_dps3_pihat.RData")

ice.bart.dps3.pihat.0025 <- ice(dps3, Xice, predictor = "pihat4", 
                               predictfcn = pred.BART.ATE.0025)
save(ice.bart.dps3.pihat.0025, file="ice_dps3_pihat.0025.RData")

ice.bart.dps3.pihat.0975 <- ice(dps3, Xice, predictor = "pihat4", 
                               predictfcn = pred.BART.ATE.0975)
save(ice.bart.dps3.pihat.0975, file="ice_dps3_pihat.0975.RData")

###




#Random-BART
Xice = cbind(t(x),pihat5,z)

ice.bart.dps4.x1 <- ice(dps4, Xice, predictor = "X1", 
                       predictfcn = pred.BART.ATE)
save(ice.bart.dps4.x1, file="ice_dps4_x1.RData")

ice.bart.dps4.x1.0025 <- ice(dps4, Xice, predictor = "X1", 
                            predictfcn = pred.BART.ATE.0025)
save(ice.bart.dps4.x1.0025, file="ice_dps4_x1.0025.RData")

ice.bart.dps4.x1.0975 <- ice(dps4, Xice, predictor = "X1", 
                            predictfcn = pred.BART.ATE.0975)
save(ice.bart.dps4.x1.0975, file="ice_dps4_x1.0975.RData")

####

ice.bart.dps4.x2 <- ice(dps4, Xice, predictor = "X2", 
                       predictfcn = pred.BART.ATE)
save(ice.bart.dps4.x2, file="ice_dps4_x2.RData")

ice.bart.dps4.x2.0025 <- ice(dps4, Xice, predictor = "X2", 
                            predictfcn = pred.BART.ATE.0025)
save(ice.bart.dps4.x2.0025, file="ice_dps4_x2.0025.RData")

ice.bart.dps4.x2.0975 <- ice(dps4, Xice, predictor = "X2", 
                            predictfcn = pred.BART.ATE.0975)
save(ice.bart.dps4.x2.0975, file="ice_dps4_x2.0975.RData")

###

ice.bart.dps4.x3 <- ice(dps4, Xice, predictor = "X3", 
                       predictfcn = pred.BART.ATE)
save(ice.bart.dps4.x3, file="ice_dps4_x3.RData")

ice.bart.dps4.x3.0025 <- ice(dps4, Xice, predictor = "X3", 
                            predictfcn = pred.BART.ATE.0025)
save(ice.bart.dps4.x3.0025, file="ice_dps4_x3.0025.RData")

ice.bart.dps4.x3.0975 <- ice(dps4, Xice, predictor = "X3", 
                            predictfcn = pred.BART.ATE.0975)
save(ice.bart.dps4.x3.0975, file="ice_dps4_x3.0975.RData")

###

ice.bart.dps4.x4 <- ice(dps4, Xice, predictor = "X4", 
                       predictfcn = pred.BART.ATE)
save(ice.bart.dps4.x4, file="ice_dps4_x4.RData")

ice.bart.dps4.x4.0025 <- ice(dps4, Xice, predictor = "X4", 
                            predictfcn = pred.BART.ATE.0025)
save(ice.bart.dps4.x4.0025, file="ice_dps4_x4.0025.RData")

ice.bart.dps4.x4.0975 <- ice(dps4, Xice, predictor = "X4", 
                            predictfcn = pred.BART.ATE.0975)
save(ice.bart.dps4.x4.0975, file="ice_dps4_x4.0975.RData")

###

ice.bart.dps4.x5 <- ice(dps4, Xice, predictor = "X5", 
                       predictfcn = pred.BART.ATE)
save(ice.bart.dps4.x5, file="ice_dps4_x5.RData")

ice.bart.dps4.x5.0025 <- ice(dps4, Xice, predictor = "X5", 
                            predictfcn = pred.BART.ATE.0025)
save(ice.bart.dps4.x5.0025, file="ice_dps4_x5.0025.RData")

ice.bart.dps4.x5.0975 <- ice(dps4, Xice, predictor = "X5", 
                            predictfcn = pred.BART.ATE.0975)
save(ice.bart.dps4.x5.0975, file="ice_dps4_x5.0975.RData")

###

ice.bart.dps4.pihat <- ice(dps4, Xice, predictor = "pihat5", 
                          predictfcn = pred.BART.ATE)
save(ice.bart.dps4.pihat, file="ice_dps4_pihat.RData")

ice.bart.dps4.pihat.0025 <- ice(dps4, Xice, predictor = "pihat5", 
                               predictfcn = pred.BART.ATE.0025)
save(ice.bart.dps4.pihat.0025, file="ice_dps4_pihat.0025.RData")

ice.bart.dps4.pihat.0975 <- ice(dps4, Xice, predictor = "pihat5", 
                               predictfcn = pred.BART.ATE.0975)
save(ice.bart.dps4.pihat.0975, file="ice_dps4_pihat.0975.RData")

###




######
#Plot#
######
ATE <- data.frame(rowMeans(ate_est_oracle),  rowMeans(date_est_oracle),
                  rowMeans(ate_est_vanilla), rowMeans(date_est_vanilla), 
                  rowMeans(ate_est_ps), rowMeans(date_est_ps),
                  rowMeans(ate_est_psglm), rowMeans(date_est_psglm),
                  rowMeans(ate_est_psrand), rowMeans(date_est_psrand))
png("boxplot_fried.png", width = 7, height = 7, units = 'in', res = 100)
boxplot(ATE, las = 2,
        ylab = "CATE",
        names = c("Oracle", "Oracle-DART",
                  "Vanilla", "Vanilla-DART",
                  "PS-BART", "PS-DART", 
                  "GLM-BART", "GLM-DART",
                  "Rand-BART", "Rand-DART"), 
        cex = 0.5,
        cex.lab = 0.8,
        cex.axis = 0.7)
#True ATE
abline(h = mean(alpha), col = "red", lwd = 2)
abline(v = 2.5, col = "gray10", lwd = 2, lty = 2)
dev.off()

##########################

png("fried_dirich_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(dvanilla$varprob.mean, ylim = c(0,1),
     ylab = "Vanilla-DART - Dirichlet Hyperprior Draws", 
     xlab = "Variables",
     col = c("red", "red", "red", "red", "red", rep("black",93), "red"))
for(i in 1:length(dvanilla$varprob.mean)){
  segments(x0 = i, y0 = as.numeric(dvanilla$varprob.mean)[i],
           x1 = i, y1 = as.numeric(apply(dvanilla$varprob,2,quantile, probs = 0.025))[i],
           col = "black", lwd = 2)
}
for(i in 1:length(dvanilla$varprob.mean)){
  segments(x0 = i, y0 = as.numeric(dvanilla$varprob.mean)[i],
           x1 = i, y1 = as.numeric(apply(dvanilla$varprob,2,quantile, probs = 0.975))[i],
           col = "black", lwd = 2)
}
for(i in c(1,2,3,4,5,99)){
  segments(x0 = i, y0 = as.numeric(dvanilla$varprob.mean)[i],
           x1 = i, y1 = as.numeric(apply(dvanilla$varprob,2,quantile, probs = 0.025))[i],
           col = "red", lwd = 2)
  segments(x0 = i, y0 = as.numeric(dvanilla$varprob.mean)[i],
           x1 = i, y1 = as.numeric(apply(dvanilla$varprob,2,quantile, probs = 0.975))[i],
           col = "red", lwd = 2)
}
dev.off()


png("fried_dirich_2.png", width = 7, height = 7, units = 'in', res = 100)
plot(doracle$varprob.mean, ylim = c(0,1), 
     ylab = "Oracle-DART - Dirichlet Hyperprior Draws", 
     xlab = "Variables",
     col = c("red", "red", "red", "red", "red", rep("black",93), "red", "red"))
for(i in 1:length(doracle$varprob.mean)){
  segments(x0 = i, y0 = as.numeric(doracle$varprob.mean)[i],
           x1 = i, y1 = as.numeric(apply(doracle$varprob,2,quantile, probs = 0.025))[i],
           col = "black", lwd = 2)
}
for(i in 1:length(doracle$varprob.mean)){
  segments(x0 = i, y0 = as.numeric(doracle$varprob.mean)[i],
           x1 = i, y1 = as.numeric(apply(doracle$varprob,2,quantile, probs = 0.975))[i],
           col = "black", lwd = 2)
}
for(i in c(1,2,3,4,5,99,100)){
  segments(x0 = i, y0 = as.numeric(doracle$varprob.mean)[i],
           x1 = i, y1 = as.numeric(apply(doracle$varprob,2,quantile, probs = 0.025))[i],
           col = "red", lwd = 2)
  segments(x0 = i, y0 = as.numeric(doracle$varprob.mean)[i],
           x1 = i, y1 = as.numeric(apply(doracle$varprob,2,quantile, probs = 0.975))[i],
           col = "red", lwd = 2)
}
dev.off()


png("fried_dirich_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(dps$varprob.mean, ylim = c(0,1), 
     ylab = "PS-DART - Dirichlet Hyperprior Draws", 
     xlab = "Variables",
     col = c("red", "red", "red", "red", "red", rep("black",93), "red", "red"))
for(i in 1:length(dps$varprob.mean)){
  segments(x0 = i, y0 = as.numeric(dps$varprob.mean)[i],
           x1 = i, y1 = as.numeric(apply(dps$varprob,2,quantile, probs = 0.025))[i],
           col = "black", lwd = 2)
}
for(i in 1:length(dps$varprob.mean)){
  segments(x0 = i, y0 = as.numeric(dps$varprob.mean)[i],
           x1 = i, y1 = as.numeric(apply(dps$varprob,2,quantile, probs = 0.975))[i],
           col = "black", lwd = 2)
}
for(i in c(1,2,3,4,5,99,100)){
  segments(x0 = i, y0 = as.numeric(dps$varprob.mean)[i],
           x1 = i, y1 = as.numeric(apply(dps$varprob,2,quantile, probs = 0.025))[i],
           col = "red", lwd = 2)
  segments(x0 = i, y0 = as.numeric(dps$varprob.mean)[i],
           x1 = i, y1 = as.numeric(apply(dps$varprob,2,quantile, probs = 0.975))[i],
           col = "red", lwd = 2)
}
dev.off()


png("fried_dirich_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dps3$varprob.mean, ylim = c(0,1), 
     ylab = "GLM-DART - Dirichlet Hyperprior Draws", 
     xlab = "Variables",
     col = c("red", "red", "red", "red", "red", rep("black",93), "red", "red"))
for(i in 1:length(dps3$varprob.mean)){
  segments(x0 = i, y0 = as.numeric(dps3$varprob.mean)[i],
           x1 = i, y1 = as.numeric(apply(dps3$varprob,2,quantile, probs = 0.025))[i],
           col = "black", lwd = 2)
}
for(i in 1:length(dps3$varprob.mean)){
  segments(x0 = i, y0 = as.numeric(dps3$varprob.mean)[i],
           x1 = i, y1 = as.numeric(apply(dps3$varprob,2,quantile, probs = 0.975))[i],
           col = "black", lwd = 2)
}
for(i in c(1,2,3,4,5,99,100)){
  segments(x0 = i, y0 = as.numeric(dps3$varprob.mean)[i],
           x1 = i, y1 = as.numeric(apply(dps3$varprob,2,quantile, probs = 0.025))[i],
           col = "red", lwd = 2)
  segments(x0 = i, y0 = as.numeric(dps3$varprob.mean)[i],
           x1 = i, y1 = as.numeric(apply(dps3$varprob,2,quantile, probs = 0.975))[i],
           col = "red", lwd = 2)
}
dev.off()


png("fried_dirich_5.png", width = 7, height = 7, units = 'in', res = 100)
plot(dps4$varprob.mean, ylim = c(0,1), 
     ylab = "Rand-DART - Dirichlet Hyperprior Draws", 
     xlab = "Variables",
     col = c("red", "red", "red", "red", "red", rep("black",93), "red", "red"))
for(i in 1:length(dps4$varprob.mean)){
  segments(x0 = i, y0 = as.numeric(dps4$varprob.mean)[i],
           x1 = i, y1 = as.numeric(apply(dps4$varprob,2,quantile, probs = 0.025))[i],
           col = "black", lwd = 2)
}
for(i in 1:length(dps4$varprob.mean)){
  segments(x0 = i, y0 = as.numeric(dps4$varprob.mean)[i],
           x1 = i, y1 = as.numeric(apply(dps4$varprob,2,quantile, probs = 0.975))[i],
           col = "black", lwd = 2)
}
for(i in c(1,2,3,4,5,99,100)){
  segments(x0 = i, y0 = as.numeric(dps4$varprob.mean)[i],
           x1 = i, y1 = as.numeric(apply(dps4$varprob,2,quantile, probs = 0.025))[i],
           col = "red", lwd = 2)
  segments(x0 = i, y0 = as.numeric(dps4$varprob.mean)[i],
           x1 = i, y1 = as.numeric(apply(dps4$varprob,2,quantile, probs = 0.975))[i],
           col = "red", lwd = 2)
}
dev.off()


###########################
png("fried_dirich_pdart.png", width = 7, height = 7, units = 'in', res = 100)
plot(dfitz$varprob.mean, ylim = c(0,1), 
     ylab = "Probit-DART - Dirichlet Hyperprior Draws", 
     xlab = "Variables",
     col = c("red", "red", rep("black",96)))
for(i in 1:length(dfitz$varprob.mean)){
  segments(x0 = i, y0 = as.numeric(dfitz$varprob.mean)[i],
           x1 = i, y1 = as.numeric(apply(dfitz$varprob,2,quantile, probs = 0.025))[i],
           col = "black", lwd = 2)
}
for(i in 1:length(dfitz$varprob.mean)){
  segments(x0 = i, y0 = as.numeric(dfitz$varprob.mean)[i],
           x1 = i, y1 = as.numeric(apply(dfitz$varprob,2,quantile, probs = 0.975))[i],
           col = "black", lwd = 2)
}
for(i in c(1,2)){
  segments(x0 = i, y0 = as.numeric(dfitz$varprob.mean)[i],
           x1 = i, y1 = as.numeric(apply(dfitz$varprob,2,quantile, probs = 0.025))[i],
           col = "red", lwd = 2)
  segments(x0 = i, y0 = as.numeric(dfitz$varprob.mean)[i],
           x1 = i, y1 = as.numeric(apply(dfitz$varprob,2,quantile, probs = 0.975))[i],
           col = "red", lwd = 2)
}
dev.off()

png("fried_varcount_pbart.png", width = 7, height = 7, units = 'in', res = 100)
plot(colMeans(fitz$varcount > 0), ylim=c(0,1),
     ylab = "Probit-BART - Posterior Inclusion Probabilities", 
     xlab = "Variables",
     col = c("red", "red",rep("black",96)))
abline(h=0.5, lwd = 2)
dev.off()

png("fried_varcount_pdart.png", width = 7, height = 7, units = 'in', res = 100)
plot(colMeans(dfitz$varcount > 0), ylim=c(0,1),
     ylab = "Probit-DART - Posterior Inclusion Probabilities", 
     xlab = "Variables",
     col = c("red", "red",rep("black",96)))
abline(h=0.5, lwd = 2)
dev.off()



png("fried_varcount_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(colMeans(vanilla$varcount > 0), ylim=c(0,1),
     ylab = "Vanilla - Posterior Inclusion Probabilities", 
     xlab = "Variables",
     col = c("red", "red", "red", "red", "red",rep("black",93), "red"))
abline(h=0.5, lwd = 2)
dev.off()

png("fried_varcount_2.png", width = 7, height = 7, units = 'in', res = 100)
plot(colMeans(dvanilla$varcount > 0), ylim=c(0,1),
     ylab = "Vanilla-DART - Posterior Inclusion Probabilities", 
     xlab = "Variables",
     col = c("red", "red", "red", "red", "red",rep("black",93), "red"))
abline(h=0.5, lwd = 2)
dev.off()

png("fried_varcount_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(colMeans(oracle$varcount > 0), ylim=c(0,1),
     ylab = "Oracle - Posterior Inclusion Probabilities", 
     xlab = "Variables",
     col = c("red", "red", "red", "red", "red",rep("black",93), "red", "red"))
abline(h=0.5, lwd = 2)
dev.off()

png("fried_varcount_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(colMeans(doracle$varcount > 0), ylim=c(0,1),
     ylab = "Oracle-DART - Posterior Inclusion Probabilities", 
     xlab = "Variables",
     col = c("red", "red", "red","red", "red", rep("black",93), "red", "red"))
abline(h=0.5, lwd = 2)
dev.off()

png("fried_varcount_5.png", width = 7, height = 7, units = 'in', res = 100)
plot(colMeans(ps$varcount > 0), ylim=c(0,1),
     ylab = "PS-BART - Posterior Inclusion Probabilities", 
     xlab = "Variables",
     col = c("red", "red", "red","red", "red", rep("black",93), "red", "red"))
abline(h=0.5, lwd = 2)
dev.off()

png("fried_varcount_6.png", width = 7, height = 7, units = 'in', res = 100)
plot(colMeans(dps$varcount > 0), ylim=c(0,1),
     ylab = "PS-DART - Posterior Inclusion Probabilities", 
     xlab = "Variables",
     col = c("red", "red", "red","red", "red", rep("black",93), "red", "red"))
abline(h=0.5, lwd = 2)
dev.off()

png("fried_varcount_7.png", width = 7, height = 7, units = 'in', res = 100)
plot(colMeans(ps3$varcount > 0), ylim=c(0,1),
     ylab = "GLM-BART - Posterior Inclusion Probabilities", 
     xlab = "Variables",
     col = c("red", "red", "red","red", "red", rep("black",93), "red", "red"))
abline(h=0.5, lwd = 2)
dev.off()

png("fried_varcount_8.png", width = 7, height = 7, units = 'in', res = 100)
plot(colMeans(dps3$varcount > 0), ylim=c(0,1),
     ylab = "GLM-DART - Posterior Inclusion Probabilities", 
     xlab = "Variables",
     col = c("red", "red", "red","red", "red", rep("black",93), "red", "red"))
abline(h=0.5, lwd = 2)
dev.off()

png("fried_varcount_9.png", width = 7, height = 7, units = 'in', res = 100)
plot(colMeans(ps4$varcount > 0), ylim=c(0,1),
     ylab = "Rand-BART - Posterior Inclusion Probabilities", 
     xlab = "Variables",
     col = c("red", "red", "red","red", "red", rep("black",93), "red", "red"))
abline(h=0.5, lwd = 2)
dev.off()

png("fried_varcount_10.png", width = 7, height = 7, units = 'in', res = 100)
plot(colMeans(dps4$varcount > 0), ylim=c(0,1),
     ylab = "Rand-DART - Posterior Inclusion Probabilities", 
     xlab = "Variables",
     col = c("red", "red", "red","red", "red", rep("black",93), "red", "red"))
abline(h=0.5, lwd = 2)
dev.off()


#################################
##Assuring Convergence of pbart##
#################################

#Fitz

i <- floor(seq(1, n, length.out=20))

auto.corr <- acf(fitz$yhat.train[ , i], plot=FALSE)
max.lag <- max(auto.corr$lag[ , 1, 1])

j <- seq(-0.5, 0.4, length.out=10)
png("fried_pbart_1.png", width = 7, height = 7, units = 'in', res = 100)
for(h in 1:10) {
  if(h==1)
    plot(1:max.lag+j[h], auto.corr$acf[1+(1:max.lag), h, h],
         type='h', xlim=c(0, max.lag+1), ylim=c(-1, 1),
         ylab='acf', xlab='lag')
  else
    lines(1:max.lag+j[h], auto.corr$acf[1+(1:max.lag), h, h],
          type='h', col=h)
}
dev.off()

png("fried_pbart_2.png", width = 7, height = 7, units = 'in', res = 100)
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

j <- -10^(log10(n)-1)
png("fried_pbart_3.png", width = 7, height = 7, units = 'in', res = 100)
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



#DFitz

i <- floor(seq(1, n, length.out=20))

auto.corr <- acf(dfitz$yhat.train[ , i], plot=FALSE)
max.lag <- max(auto.corr$lag[ , 1, 1])

j <- seq(-0.5, 0.4, length.out=10)
png("fried_pdart_1.png", width = 7, height = 7, units = 'in', res = 100)
for(h in 1:10) {
  if(h==1)
    plot(1:max.lag+j[h], auto.corr$acf[1+(1:max.lag), h, h],
         type='h', xlim=c(0, max.lag+1), ylim=c(-1, 1),
         ylab='acf', xlab='lag')
  else
    lines(1:max.lag+j[h], auto.corr$acf[1+(1:max.lag), h, h],
          type='h', col=h)
}
dev.off()

png("fried_pdart_2.png", width = 7, height = 7, units = 'in', res = 100)
for(j in 1:10) {
  if(j==1)
    plot(pnorm(dfitz$yhat.train[ , i[j]]),
         type='l', ylim=c(0, 1),
         sub=paste0('N:', n, ', p:', p, ', thin:', pthin),
         ylab=expression(Phi(f(x))), xlab='m')
  else
    lines(pnorm(dfitz$yhat.train[ , i[j]]),
          type='l', col=j)
}
dev.off()

geweke <- gewekediag(dfitz$yhat.train)

j <- -10^(log10(n)-1)
png("fried_pdart_3.png", width = 7, height = 7, units = 'in', res = 100)
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


########################
########################

png("friedvanilla_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(vanilla$sigma[seq(50000,150000, by = 100)], type = "l",
     ylab = expression(sigma),
     xlab = "Draw",
     ylim= c(0.2, 0.65))
dev.off()
png("friedvanilla_2.png", width = 7, height = 7, units = 'in', res = 100)
acf(vanilla$sigma[seq(50000,150000, by = 100)], main = "")
dev.off()

png("friedoracle_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(oracle$sigma[seq(50000,150000, by = 100)], type = "l",
     ylab = expression(sigma),
     xlab = "Draw",
     ylim= c(0.2, 0.65))
dev.off()
png("friedoracle_2.png", width = 7, height = 7, units = 'in', res = 100)
acf(oracle$sigma[seq(50000,150000, by = 100)], main = "")
dev.off()

png("friedps_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ps$sigma[seq(50000,150000, by = 100)], type = "l",
     ylab = expression(sigma),
     xlab = "Draw",
     ylim= c(0.2, 0.65))
dev.off()
png("friedps_2.png", width = 7, height = 7, units = 'in', res = 100)
acf(ps$sigma[seq(50000,150000, by = 100)], main = "")
dev.off()

png("friedglm_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ps3$sigma[seq(50000,150000, by = 100)], type = "l",
     ylab = expression(sigma),
     xlab = "Draw",
     ylim= c(0.2, 0.65))
dev.off()
png("friedglm_2.png", width = 7, height = 7, units = 'in', res = 100)
acf(ps3$sigma[seq(50000,150000, by = 100)], main = "")
dev.off()

png("friedrand_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ps4$sigma[seq(50000,150000, by = 100)], type = "l",
     ylab = expression(sigma),
     xlab = "Draw",
     ylim= c(0.2, 0.65))
dev.off()
png("friedrand_2.png", width = 7, height = 7, units = 'in', res = 100)
acf(ps4$sigma[seq(50000,150000, by = 100)], main = "")
dev.off()

png("frieddvanilla_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(dvanilla$sigma[seq(50000,150000, by = 100)], type = "l",
     ylab = expression(sigma),
     xlab = "Draw",
     ylim= c(0.2, 0.65))
dev.off()
png("frieddvanilla_2.png", width = 7, height = 7, units = 'in', res = 100)
acf(dvanilla$sigma[seq(50000,150000, by = 100)], main = "")
dev.off()

png("frieddoracle_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(doracle$sigma[seq(50000,150000, by = 100)], type = "l",
     ylab = expression(sigma),
     xlab = "Draw",
     ylim= c(0.2, 0.65))
dev.off()
png("frieddoracle_2.png", width = 7, height = 7, units = 'in', res = 100)
acf(doracle$sigma[seq(50000,150000, by = 100)], main = "")
dev.off()

png("frieddps_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(dps$sigma[seq(50000,150000, by = 100)], type = "l",
     ylab = expression(sigma),
     xlab = "Draw",
     ylim= c(0.2, 0.65))
dev.off()
png("frieddps_2.png", width = 7, height = 7, units = 'in', res = 100)
acf(dps$sigma[seq(50000,150000, by = 100)], main = "")
dev.off()

png("frieddps3_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(dps3$sigma[seq(50000,150000, by = 100)], type = "l",
     ylab = expression(sigma),
     xlab = "Draw",
     ylim= c(0.2, 0.65))
dev.off()
png("frieddps3_2.png", width = 7, height = 7, units = 'in', res = 100)
acf(dps3$sigma[seq(50000,150000, by = 100)], main = "")
dev.off()

png("frieddps4_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(dps4$sigma[seq(50000,150000, by = 100)], type = "l",
     ylab = expression(sigma),
     xlab = "Draw",
     ylim= c(0.2, 0.65))
dev.off()
png("frieddps4_2.png", width = 7, height = 7, units = 'in', res = 100)
acf(dps4$sigma[seq(50000,150000, by = 100)], main = "")
dev.off()



par(cex=0.5)

#################
#####Vanilla#####
#################

png("fried.vanilla.icex1_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x1,
     ylim = c(-4,6),
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


png("fried.vanilla.icex1_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x1.0025,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("fried.vanilla.icex1_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x1.0975,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("fried.vanilla.icex1_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x1,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[1]))
dev.off()

png("fried.vanilla.icex1_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.vanilla.x1),
     ylab = "Derivative Y",
     xlab = expression(x[1]))
dev.off()


png("fried.vanilla.icex2_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x2,
     ylim = c(-4,6),
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


png("fried.vanilla.icex2_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x2.0025,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("fried.vanilla.icex2_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x2.0975,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("fried.vanilla.icex2_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x2,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[2]))
dev.off()

png("fried.vanilla.icex2_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.vanilla.x2),
     ylab = "Derivative Y",
     xlab = expression(x[2]))
dev.off()


png("fried.vanilla.icex3_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x3,
     ylim = c(-4,6),
     ylab = "Partial Y",
     xlab = expression(x[3]))
lines(as.numeric(names(apply(ice.bart.vanilla.x3.0975$ice_curves,2,mean))),
      apply(ice.bart.vanilla.x3.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.vanilla.x3.0025$ice_curves,2,mean))),
      apply(ice.bart.vanilla.x3.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
{segments(x0 = -3, y0 = 0, x1 = 1/4, y1 = 0, col = "red", lwd = 3)
  segments(x0 = 1/4, y0 = 0, x1 = 1/4, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 1/4, y0 = 0.5, x1 = 2/4, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 2/4, y0 = 0.5, x1 = 2/4, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 2/4, y0 = 0.75, x1 = 3/4, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 0.75, x1 = 3/4, y1 = 1, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 1, x1 = 15*sd(x[3,])+mean(x[3,]), y1 = 1, col = "red", lwd = 3)}

dev.off()


png("fried.vanilla.icex3_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x3.0025,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("fried.vanilla.icex3_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x3.0975,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("fried.vanilla.icex3_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x3,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[3]))
dev.off()

png("fried.vanilla.icex3_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.vanilla.x3),
     ylab = "Derivative Y",
     xlab = expression(x[3]))
dev.off()


png("fried.vanilla.icex4_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x4,
     ylim = c(-4,6),
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


png("fried.vanilla.icex4_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x4.0025,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("fried.vanilla.icex4_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x4.0975,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("fried.vanilla.icex4_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x4,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[4]))
dev.off()

png("fried.vanilla.icex4_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.vanilla.x4),
     ylab = "Derivative Y",
     xlab = expression(x[4]))
dev.off()

png("fried.vanilla.icex5_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x5,
     ylim = c(-4,6),
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


png("fried.vanilla.icex5_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x5.0025,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("fried.vanilla.icex5_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x5.0975,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("fried.vanilla.icex5_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x5,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[5]))
dev.off()

png("fried.vanilla.icex5_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.vanilla.x5),
     ylab = "Derivative Y",
     xlab = expression(x[5]))
dev.off()


###################
######Oracle#######
###################

png("fried.oracle.icex1_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x1,
     ylim = c(-4,6),
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


png("fried.oracle.icex1_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x1.0025,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("fried.oracle.icex1_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x1.0975,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("fried.oracle.icex1_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x1,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[1]))
dev.off()

png("fried.oracle.icex1_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.oracle.x1),
     ylab = "Derivative Y",
     xlab = expression(x[1]))
dev.off()


png("fried.oracle.icex2_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x2,
     ylim = c(-4,6),
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


png("fried.oracle.icex2_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x2.0025,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("fried.oracle.icex2_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x2.0975,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("fried.oracle.icex2_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x2,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[2]))
dev.off()

png("fried.oracle.icex2_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.oracle.x2),
     ylab = "Derivative Y",
     xlab = expression(x[2]))
dev.off()


png("fried.oracle.icex3_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x3,
     ylim = c(-4,6),
     ylab = "Partial Y",
     xlab = expression(x[3]))
lines(as.numeric(names(apply(ice.bart.oracle.x3.0975$ice_curves,2,mean))),
      apply(ice.bart.oracle.x3.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.oracle.x3.0025$ice_curves,2,mean))),
      apply(ice.bart.oracle.x3.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
{segments(x0 = -3, y0 = 0, x1 = 1/4, y1 = 0, col = "red", lwd = 3)
  segments(x0 = 1/4, y0 = 0, x1 = 1/4, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 1/4, y0 = 0.5, x1 = 2/4, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 2/4, y0 = 0.5, x1 = 2/4, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 2/4, y0 = 0.75, x1 = 3/4, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 0.75, x1 = 3/4, y1 = 1, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 1, x1 = 15*sd(x[3,])+mean(x[3,]), y1 = 1, col = "red", lwd = 3)}

dev.off()


png("fried.oracle.icex3_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x3.0025,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("fried.oracle.icex3_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x3.0975,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("fried.oracle.icex3_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x3,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[3]))
dev.off()

png("fried.oracle.icex3_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.oracle.x3),
     ylab = "Derivative Y",
     xlab = expression(x[3]))
dev.off()

png("fried.oracle.icepihat_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.pihat,
     ylim = c(-4,6),
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

png("fried.oracle.icepihat_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.pihat.0025,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("fried.oracle.icepihat_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.pihat.0975,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("fried.oracle.icepihat_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.pihat,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("fried.oracle.icepihat_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.oracle.pihat),
     ylab = "Derivative Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("fried.oracle.icex4_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x4,
     ylim = c(-4,6),
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


png("fried.oracle.icex4_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x4.0025,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("fried.oracle.icex4_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x4.0975,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("fried.oracle.icex4_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x4,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[4]))
dev.off()

png("fried.oracle.icex4_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.oracle.x4),
     ylab = "Derivative Y",
     xlab = expression(x[4]))
dev.off()

png("fried.oracle.icex5_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x5,
     ylim = c(-4,6),
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


png("fried.oracle.icex5_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x5.0025,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("fried.oracle.icex5_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x5.0975,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("fried.oracle.icex5_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x5,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[5]))
dev.off()

png("fried.oracle.icex5_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.oracle.x5),
     ylab = "Derivative Y",
     xlab = expression(x[5]))
dev.off()


###########
#####PS####
###########

png("fried.ps.icex1_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x1,
     ylim = c(-4,6),
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


png("fried.ps.icex1_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x1.0025,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("fried.ps.icex1_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x1.0975,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("fried.ps.icex1_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x1,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[1]))
dev.off()

png("fried.ps.icex1_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps.x1),
     ylab = "Derivative Y",
     xlab = expression(x[1]))
dev.off()


png("fried.ps.icex2_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x2,
     ylim = c(-4,6),
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


png("fried.ps.icex2_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x2.0025,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("fried.ps.icex2_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x2.0975,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("fried.ps.icex2_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x2,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[2]))
dev.off()

png("fried.ps.icex2_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps.x2),
     ylab = "Derivative Y",
     xlab = expression(x[2]))
dev.off()


png("fried.ps.icex3_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x3,
     ylim = c(-4,6),
     ylab = "Partial Y",
     xlab = expression(x[3]))
lines(as.numeric(names(apply(ice.bart.ps.x3.0975$ice_curves,2,mean))),
      apply(ice.bart.ps.x3.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.ps.x3.0025$ice_curves,2,mean))),
      apply(ice.bart.ps.x3.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
{segments(x0 = -3, y0 = 0, x1 = 1/4, y1 = 0, col = "red", lwd = 3)
  segments(x0 = 1/4, y0 = 0, x1 = 1/4, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 1/4, y0 = 0.5, x1 = 2/4, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 2/4, y0 = 0.5, x1 = 2/4, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 2/4, y0 = 0.75, x1 = 3/4, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 0.75, x1 = 3/4, y1 = 1, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 1, x1 = 15*sd(x[3,])+mean(x[3,]), y1 = 1, col = "red", lwd = 3)}

dev.off()


png("fried.ps.icex3_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x3.0025,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("fried.ps.icex3_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x3.0975,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("fried.ps.icex3_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x3,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[3]))
dev.off()

png("fried.ps.icex3_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps.x3),
     ylab = "Derivative Y",
     xlab = expression(x[3]))
dev.off()

png("fried.ps.icepihat_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.pihat,
     ylim = c(-4,6),
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


png("fried.ps.icepihat_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.pihat.0025,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("fried.ps.icepihat_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.pihat.0975,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("fried.ps.icepihat_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.pihat,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("fried.ps.icepihat_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps.pihat),
     ylab = "Derivative Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("fried.ps.icex4_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x4,
     ylim = c(-4,6),
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


png("fried.ps.icex4_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x4.0025,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("fried.ps.icex4_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x4.0975,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("fried.ps.icex4_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x4,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[4]))
dev.off()

png("fried.ps.icex4_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps.x4),
     ylab = "Derivative Y",
     xlab = expression(x[4]))
dev.off()

png("fried.ps.icex5_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x5,
     ylim = c(-4,6),
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


png("fried.ps.icex5_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x5.0025,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("fried.ps.icex5_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x5.0975,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("fried.ps.icex5_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x5,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[5]))
dev.off()

png("fried.ps.icex5_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps.x5),
     ylab = "Derivative Y",
     xlab = expression(x[5]))
dev.off()


#########
###PS3###
#########

png("fried.ps3.icex1_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x1,
     ylim = c(-4,6),
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


png("fried.ps3.icex1_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x1.0025,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("fried.ps3.icex1_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x1.0975,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("fried.ps3.icex1_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x1,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[1]))
dev.off()

png("fried.ps3.icex1_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps3.x1),
     ylab = "Derivative Y",
     xlab = expression(x[1]))
dev.off()


png("fried.ps3.icex2_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x2,
     ylim = c(-4,6),
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


png("fried.ps3.icex2_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x2.0025,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("fried.ps3.icex2_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x2.0975,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("fried.ps3.icex2_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x2,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[2]))
dev.off()

png("fried.ps3.icex2_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps3.x2),
     ylab = "Derivative Y",
     xlab = expression(x[2]))
dev.off()


png("fried.ps3.icex3_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x3,
     ylim = c(-4,6),
     ylab = "Partial Y",
     xlab = expression(x[3]))
lines(as.numeric(names(apply(ice.bart.ps3.x3.0975$ice_curves,2,mean))),
      apply(ice.bart.ps3.x3.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.ps3.x3.0025$ice_curves,2,mean))),
      apply(ice.bart.ps3.x3.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
{segments(x0 = -3, y0 = 0, x1 = 1/4, y1 = 0, col = "red", lwd = 3)
  segments(x0 = 1/4, y0 = 0, x1 = 1/4, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 1/4, y0 = 0.5, x1 = 2/4, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 2/4, y0 = 0.5, x1 = 2/4, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 2/4, y0 = 0.75, x1 = 3/4, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 0.75, x1 = 3/4, y1 = 1, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 1, x1 = 15*sd(x[3,])+mean(x[3,]), y1 = 1, col = "red", lwd = 3)}

dev.off()


png("fried.ps3.icex3_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x3.0025,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("fried.ps3.icex3_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x3.0975,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("fried.ps3.icex3_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x3,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[3]))
dev.off()

png("fried.ps3.icex3_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps3.x3),
     ylab = "Derivative Y",
     xlab = expression(x[3]))
dev.off()

png("fried.ps3.icepihat_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.pihat,
     ylim = c(-4,6),
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


png("fried.ps3.icepihat_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.pihat.0025,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("fried.ps3.icepihat_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.pihat.0975,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("fried.ps3.icepihat_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.pihat,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("fried.ps3.icepihat_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps3.pihat),
     ylab = "Derivative Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("fried.ps3.icex4_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x4,
     ylim = c(-4,6),
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


png("fried.ps3.icex4_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x4.0025,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("fried.ps3.icex4_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x4.0975,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("fried.ps3.icex4_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x4,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[4]))
dev.off()

png("fried.ps3.icex4_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps3.x4),
     ylab = "Derivative Y",
     xlab = expression(x[4]))
dev.off()

png("fried.ps3.icex5_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x5,
     ylim = c(-4,6),
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


png("fried.ps3.icex5_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x5.0025,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("fried.ps3.icex5_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x5.0975,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("fried.ps3.icex5_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x5,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[5]))
dev.off()

png("fried.ps3.icex5_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps3.x5),
     ylab = "Derivative Y",
     xlab = expression(x[5]))
dev.off()


#######
##PS4##
#######

png("fried.ps4.icex1_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x1,
     ylim = c(-4,6),
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


png("fried.ps4.icex1_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x1.0025,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("fried.ps4.icex1_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x1.0975,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("fried.ps4.icex1_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x1,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[1]))
dev.off()

png("fried.ps4.icex1_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps4.x1),
     ylab = "Derivative Y",
     xlab = expression(x[1]))
dev.off()


png("fried.ps4.icex2_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x2,
     ylim = c(-4,6),
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


png("fried.ps4.icex2_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x2.0025,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("fried.ps4.icex2_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x2.0975,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("fried.ps4.icex2_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x2,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[2]))
dev.off()

png("fried.ps4.icex2_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps4.x2),
     ylab = "Derivative Y",
     xlab = expression(x[2]))
dev.off()


png("fried.ps4.icex3_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x3,
     ylim = c(-4,6),
     ylab = "Partial Y",
     xlab = expression(x[3]))
lines(as.numeric(names(apply(ice.bart.ps4.x3.0975$ice_curves,2,mean))),
      apply(ice.bart.ps4.x3.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.ps4.x3.0025$ice_curves,2,mean))),
      apply(ice.bart.ps4.x3.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
{segments(x0 = -3, y0 = 0, x1 = 1/4, y1 = 0, col = "red", lwd = 3)
  segments(x0 = 1/4, y0 = 0, x1 = 1/4, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 1/4, y0 = 0.5, x1 = 2/4, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 2/4, y0 = 0.5, x1 = 2/4, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 2/4, y0 = 0.75, x1 = 3/4, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 0.75, x1 = 3/4, y1 = 1, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 1, x1 = 15*sd(x[3,])+mean(x[3,]), y1 = 1, col = "red", lwd = 3)}

dev.off()


png("fried.ps4.icex3_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x3.0025,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("fried.ps4.icex3_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x3.0975,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("fried.ps4.icex3_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x3,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[3]))
dev.off()

png("fried.ps4.icex3_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps4.x3),
     ylab = "Derivative Y",
     xlab = expression(x[3]))
dev.off()

png("fried.ps4.icepihat_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.pihat,
     ylim = c(-4,6),
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


png("fried.ps4.icepihat_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.pihat.0025,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("fried.ps4.icepihat_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.pihat.0975,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("fried.ps4.icepihat_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.pihat,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("fried.ps4.icepihat_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps4.pihat),
     ylab = "Derivative Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("fried.ps4.icex4_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x4,
     ylim = c(-4,6),
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


png("fried.ps4.icex4_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x4.0025,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("fried.ps4.icex4_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x4.0975,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("fried.ps4.icex4_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x4,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[4]))
dev.off()

png("fried.ps4.icex4_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps4.x4),
     ylab = "Derivative Y",
     xlab = expression(x[4]))
dev.off()

png("fried.ps4.icex5_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x5,
     ylim = c(-4,6),
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


png("fried.ps4.icex5_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x5.0025,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("fried.ps4.icex5_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x5.0975,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("fried.ps4.icex5_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x5,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[5]))
dev.off()

png("fried.ps4.icex5_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps4.x5),
     ylab = "Derivative Y",
     xlab = expression(x[5]))
dev.off()


#####################################
############## DART #################
#####################################
###############
####dvanilla#####
###############


png("fried.bart.dvanilla.icex1_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dvanilla.x1,
     ylim = c(-4,6),
     ylab = "Partial Y",
     xlab = expression(x[1]))
lines(as.numeric(names(apply(ice.bart.dvanilla.x1.0975$ice_curves,2,mean))),
      apply(ice.bart.dvanilla.x1.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.dvanilla.x1.0025$ice_curves,2,mean))),
      apply(ice.bart.dvanilla.x1.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("fried.bart.dvanilla.icex1_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dvanilla.x1.0025,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("fried.bart.dvanilla.icex1_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dvanilla.x1.0975,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("fried.bart.dvanilla.icex1_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dvanilla.x1,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[1]))
dev.off()

png("fried.bart.dvanilla.icex1_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.dvanilla.x1),
     ylab = "Derivative Y",
     xlab = expression(x[1]))
dev.off()


png("fried.bart.dvanilla.icex2_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dvanilla.x2,
     ylim = c(-4,6),
     ylab = "Partial Y",
     xlab = expression(x[2]))
lines(as.numeric(names(apply(ice.bart.dvanilla.x2.0975$ice_curves,2,mean))),
      apply(ice.bart.dvanilla.x2.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.dvanilla.x2.0025$ice_curves,2,mean))),
      apply(ice.bart.dvanilla.x2.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("fried.bart.dvanilla.icex2_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dvanilla.x2.0025,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("fried.bart.dvanilla.icex2_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dvanilla.x2.0975,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("fried.bart.dvanilla.icex2_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dvanilla.x2,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[2]))
dev.off()

png("fried.bart.dvanilla.icex2_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.dvanilla.x2),
     ylab = "Derivative Y",
     xlab = expression(x[2]))
dev.off()


png("fried.bart.dvanilla.icex3_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dvanilla.x3,
     ylim = c(-4,6),
     ylab = "Partial Y",
     xlab = expression(x[3]))
lines(as.numeric(names(apply(ice.bart.dvanilla.x3.0975$ice_curves,2,mean))),
      apply(ice.bart.dvanilla.x3.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.dvanilla.x3.0025$ice_curves,2,mean))),
      apply(ice.bart.dvanilla.x3.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
{segments(x0 = -3, y0 = 0, x1 = 1/4, y1 = 0, col = "red", lwd = 3)
  segments(x0 = 1/4, y0 = 0, x1 = 1/4, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 1/4, y0 = 0.5, x1 = 2/4, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 2/4, y0 = 0.5, x1 = 2/4, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 2/4, y0 = 0.75, x1 = 3/4, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 0.75, x1 = 3/4, y1 = 1, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 1, x1 = 15*sd(x[3,])+mean(x[3,]), y1 = 1, col = "red", lwd = 3)}

dev.off()


png("fried.bart.dvanilla.icex3_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dvanilla.x3.0025,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("fried.bart.dvanilla.icex3_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dvanilla.x3.0975,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("fried.bart.dvanilla.icex3_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dvanilla.x3,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[3]))
dev.off()

png("fried.bart.dvanilla.icex3_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.dvanilla.x3),
     ylab = "Derivative Y",
     xlab = expression(x[3]))
dev.off()


png("fried.bart.dvanilla.icex4_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dvanilla.x4,
     ylim = c(-4,6),
     ylab = "Partial Y",
     xlab = expression(x[4]))
lines(as.numeric(names(apply(ice.bart.dvanilla.x4.0975$ice_curves,2,mean))),
      apply(ice.bart.dvanilla.x4.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.dvanilla.x4.0025$ice_curves,2,mean))),
      apply(ice.bart.dvanilla.x4.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("fried.bart.dvanilla.icex4_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dvanilla.x4.0025,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("fried.bart.dvanilla.icex4_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dvanilla.x4.0975,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("fried.bart.dvanilla.icex4_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dvanilla.x4,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[4]))
dev.off()

png("fried.bart.dvanilla.icex4_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.dvanilla.x4),
     ylab = "Derivative Y",
     xlab = expression(x[4]))
dev.off()

png("fried.bart.dvanilla.icex5_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dvanilla.x5,
     ylim = c(-4,6),
     ylab = "Partial Y",
     xlab = expression(x[5]))
lines(as.numeric(names(apply(ice.bart.dvanilla.x5.0975$ice_curves,2,mean))),
      apply(ice.bart.dvanilla.x5.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.dvanilla.x5.0025$ice_curves,2,mean))),
      apply(ice.bart.dvanilla.x5.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("fried.bart.dvanilla.icex5_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dvanilla.x5.0025,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("fried.bart.dvanilla.icex5_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dvanilla.x5.0975,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("fried.bart.dvanilla.icex5_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dvanilla.x5,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[5]))
dev.off()

png("fried.bart.dvanilla.icex5_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.dvanilla.x5),
     ylab = "Derivative Y",
     xlab = expression(x[5]))
dev.off()



###############
####doracle#####
###############


png("fried.bart.doracle.icex1_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.doracle.x1,
     ylim = c(-4,6),
     ylab = "Partial Y",
     xlab = expression(x[1]))
lines(as.numeric(names(apply(ice.bart.doracle.x1.0975$ice_curves,2,mean))),
      apply(ice.bart.doracle.x1.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.doracle.x1.0025$ice_curves,2,mean))),
      apply(ice.bart.doracle.x1.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("fried.bart.doracle.icex1_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.doracle.x1.0025,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("fried.bart.doracle.icex1_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.doracle.x1.0975,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("fried.bart.doracle.icex1_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.doracle.x1,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[1]))
dev.off()

png("fried.bart.doracle.icex1_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.doracle.x1),
     ylab = "Derivative Y",
     xlab = expression(x[1]))
dev.off()


png("fried.bart.doracle.icex2_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.doracle.x2,
     ylim = c(-4,6),
     ylab = "Partial Y",
     xlab = expression(x[2]))
lines(as.numeric(names(apply(ice.bart.doracle.x2.0975$ice_curves,2,mean))),
      apply(ice.bart.doracle.x2.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.doracle.x2.0025$ice_curves,2,mean))),
      apply(ice.bart.doracle.x2.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("fried.bart.doracle.icex2_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.doracle.x2.0025,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("fried.bart.doracle.icex2_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.doracle.x2.0975,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("fried.bart.doracle.icex2_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.doracle.x2,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[2]))
dev.off()

png("fried.bart.doracle.icex2_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.doracle.x2),
     ylab = "Derivative Y",
     xlab = expression(x[2]))
dev.off()


png("fried.bart.doracle.icex3_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.doracle.x3,
     ylim = c(-4,6),
     ylab = "Partial Y",
     xlab = expression(x[3]))
lines(as.numeric(names(apply(ice.bart.doracle.x3.0975$ice_curves,2,mean))),
      apply(ice.bart.doracle.x3.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.doracle.x3.0025$ice_curves,2,mean))),
      apply(ice.bart.doracle.x3.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
{segments(x0 = -3, y0 = 0, x1 = 1/4, y1 = 0, col = "red", lwd = 3)
  segments(x0 = 1/4, y0 = 0, x1 = 1/4, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 1/4, y0 = 0.5, x1 = 2/4, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 2/4, y0 = 0.5, x1 = 2/4, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 2/4, y0 = 0.75, x1 = 3/4, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 0.75, x1 = 3/4, y1 = 1, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 1, x1 = 15*sd(x[3,])+mean(x[3,]), y1 = 1, col = "red", lwd = 3)}

dev.off()


png("fried.bart.doracle.icex3_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.doracle.x3.0025,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("fried.bart.doracle.icex3_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.doracle.x3.0975,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("fried.bart.doracle.icex3_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.doracle.x3,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[3]))
dev.off()

png("fried.bart.doracle.icex3_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.doracle.x3),
     ylab = "Derivative Y",
     xlab = expression(x[3]))
dev.off()

png("fried.bart.doracle.icepihat_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.doracle.pihat,
     ylim = c(-4,6),
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
lines(as.numeric(names(apply(ice.bart.doracle.pihat.0975$ice_curves,2,mean))),
      apply(ice.bart.doracle.pihat.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.doracle.pihat.0025$ice_curves,2,mean))),
      apply(ice.bart.doracle.pihat.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("fried.bart.doracle.icepihat_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.doracle.pihat.0025,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("fried.bart.doracle.icepihat_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.doracle.pihat.0975,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("fried.bart.doracle.icepihat_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.doracle.pihat,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("fried.bart.doracle.icepihat_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.doracle.pihat),
     ylab = "Derivative Y",
     xlab = expression(hat(Pi)(x)))
dev.off()


png("fried.bart.doracle.icex4_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.doracle.x4,
     ylim = c(-4,6),
     ylab = "Partial Y",
     xlab = expression(x[4]))
lines(as.numeric(names(apply(ice.bart.doracle.x4.0975$ice_curves,2,mean))),
      apply(ice.bart.doracle.x4.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.doracle.x4.0025$ice_curves,2,mean))),
      apply(ice.bart.doracle.x4.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("fried.bart.doracle.icex4_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.doracle.x4.0025,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("fried.bart.doracle.icex4_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.doracle.x4.0975,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("fried.bart.doracle.icex4_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.doracle.x4,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[4]))
dev.off()

png("fried.bart.doracle.icex4_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.doracle.x4),
     ylab = "Derivative Y",
     xlab = expression(x[4]))
dev.off()

png("fried.bart.doracle.icex5_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.doracle.x5,
     ylim = c(-4,6),
     ylab = "Partial Y",
     xlab = expression(x[5]))
lines(as.numeric(names(apply(ice.bart.doracle.x5.0975$ice_curves,2,mean))),
      apply(ice.bart.doracle.x5.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.doracle.x5.0025$ice_curves,2,mean))),
      apply(ice.bart.doracle.x5.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("fried.bart.doracle.icex5_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.doracle.x5.0025,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("fried.bart.doracle.icex5_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.doracle.x5.0975,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("fried.bart.doracle.icex5_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.doracle.x5,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[5]))
dev.off()

png("fried.bart.doracle.icex5_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.doracle.x5),
     ylab = "Derivative Y",
     xlab = expression(x[5]))
dev.off()


#############
#####Bart####
#############

png("fried.bart.dps.icex1_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps.x1,
     ylim = c(-4,6),
     ylab = "Partial Y",
     xlab = expression(x[1]))
lines(as.numeric(names(apply(ice.bart.dps.x1.0975$ice_curves,2,mean))),
      apply(ice.bart.dps.x1.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.dps.x1.0025$ice_curves,2,mean))),
      apply(ice.bart.dps.x1.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("fried.bart.dps.icex1_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps.x1.0025,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("fried.bart.dps.icex1_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps.x1.0975,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("fried.bart.dps.icex1_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps.x1,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[1]))
dev.off()

png("fried.bart.dps.icex1_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.dps.x1),
     ylab = "Derivative Y",
     xlab = expression(x[1]))
dev.off()


png("fried.bart.dps.icex2_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps.x2,
     ylim = c(-4,6),
     ylab = "Partial Y",
     xlab = expression(x[2]))
lines(as.numeric(names(apply(ice.bart.dps.x2.0975$ice_curves,2,mean))),
      apply(ice.bart.dps.x2.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.dps.x2.0025$ice_curves,2,mean))),
      apply(ice.bart.dps.x2.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("fried.bart.dps.icex2_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps.x2.0025,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("fried.bart.dps.icex2_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps.x2.0975,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("fried.bart.dps.icex2_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps.x2,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[2]))
dev.off()

png("fried.bart.dps.icex2_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.dps.x2),
     ylab = "Derivative Y",
     xlab = expression(x[2]))
dev.off()


png("fried.bart.dps.icex3_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps.x3,
     ylim = c(-4,6),
     ylab = "Partial Y",
     xlab = expression(x[3]))
lines(as.numeric(names(apply(ice.bart.dps.x3.0975$ice_curves,2,mean))),
      apply(ice.bart.dps.x3.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.dps.x3.0025$ice_curves,2,mean))),
      apply(ice.bart.dps.x3.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
{segments(x0 = -3, y0 = 0, x1 = 1/4, y1 = 0, col = "red", lwd = 3)
  segments(x0 = 1/4, y0 = 0, x1 = 1/4, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 1/4, y0 = 0.5, x1 = 2/4, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 2/4, y0 = 0.5, x1 = 2/4, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 2/4, y0 = 0.75, x1 = 3/4, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 0.75, x1 = 3/4, y1 = 1, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 1, x1 = 15*sd(x[3,])+mean(x[3,]), y1 = 1, col = "red", lwd = 3)}

dev.off()


png("fried.bart.dps.icex3_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps.x3.0025,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("fried.bart.dps.icex3_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps.x3.0975,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("fried.bart.dps.icex3_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps.x3,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[3]))
dev.off()

png("fried.bart.dps.icex3_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.dps.x3),
     ylab = "Derivative Y",
     xlab = expression(x[3]))
dev.off()

png("fried.bart.dps.icepihat_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps.pihat,
     ylim = c(-4,6),
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
lines(as.numeric(names(apply(ice.bart.dps.pihat.0975$ice_curves,2,mean))),
      apply(ice.bart.dps.pihat.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.dps.pihat.0025$ice_curves,2,mean))),
      apply(ice.bart.dps.pihat.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("fried.bart.dps.icepihat_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps.pihat.0025,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("fried.bart.dps.icepihat_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps.pihat.0975,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("fried.bart.dps.icepihat_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps.pihat,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("fried.bart.dps.icepihat_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.dps.pihat),
     ylab = "Derivative Y",
     xlab = expression(hat(Pi)(x)))
dev.off()


png("fried.bart.dps.icex4_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps.x4,
     ylim = c(-4,6),
     ylab = "Partial Y",
     xlab = expression(x[4]))
lines(as.numeric(names(apply(ice.bart.dps.x4.0975$ice_curves,2,mean))),
      apply(ice.bart.dps.x4.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.dps.x4.0025$ice_curves,2,mean))),
      apply(ice.bart.dps.x4.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("fried.bart.dps.icex4_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps.x4.0025,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("fried.bart.dps.icex4_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps.x4.0975,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("fried.bart.dps.icex4_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps.x4,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[4]))
dev.off()

png("fried.bart.dps.icex4_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.dps.x4),
     ylab = "Derivative Y",
     xlab = expression(x[4]))
dev.off()

png("fried.bart.dps.icex5_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps.x5,
     ylim = c(-4,6),
     ylab = "Partial Y",
     xlab = expression(x[5]))
lines(as.numeric(names(apply(ice.bart.dps.x5.0975$ice_curves,2,mean))),
      apply(ice.bart.dps.x5.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.dps.x5.0025$ice_curves,2,mean))),
      apply(ice.bart.dps.x5.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("fried.bart.dps.icex5_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps.x5.0025,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("fried.bart.dps.icex5_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps.x5.0975,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("fried.bart.dps.icex5_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps.x5,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[5]))
dev.off()

png("fried.bart.dps.icex5_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.dps.x5),
     ylab = "Derivative Y",
     xlab = expression(x[5]))
dev.off()


#########
###dps3###
#########

png("fried.bart.dps3.icex1_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps3.x1,
     ylim = c(-4,6),
     ylab = "Partial Y",
     xlab = expression(x[1]))
lines(as.numeric(names(apply(ice.bart.dps3.x1.0975$ice_curves,2,mean))),
      apply(ice.bart.dps3.x1.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.dps3.x1.0025$ice_curves,2,mean))),
      apply(ice.bart.dps3.x1.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("fried.bart.dps3.icex1_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps3.x1.0025,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("fried.bart.dps3.icex1_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps3.x1.0975,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("fried.bart.dps3.icex1_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps3.x1,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[1]))
dev.off()

png("fried.bart.dps3.icex1_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.dps3.x1),
     ylab = "Derivative Y",
     xlab = expression(x[1]))
dev.off()


png("fried.bart.dps3.icex2_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps3.x2,
     ylim = c(-4,6),
     ylab = "Partial Y",
     xlab = expression(x[2]))
lines(as.numeric(names(apply(ice.bart.dps3.x2.0975$ice_curves,2,mean))),
      apply(ice.bart.dps3.x2.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.dps3.x2.0025$ice_curves,2,mean))),
      apply(ice.bart.dps3.x2.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("fried.bart.dps3.icex2_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps3.x2.0025,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("fried.bart.dps3.icex2_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps3.x2.0975,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("fried.bart.dps3.icex2_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps3.x2,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[2]))
dev.off()

png("fried.bart.dps3.icex2_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.dps3.x2),
     ylab = "Derivative Y",
     xlab = expression(x[2]))
dev.off()


png("fried.bart.dps3.icex3_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps3.x3,
     ylim = c(-4,6),
     ylab = "Partial Y",
     xlab = expression(x[3]))
lines(as.numeric(names(apply(ice.bart.dps3.x3.0975$ice_curves,2,mean))),
      apply(ice.bart.dps3.x3.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.dps3.x3.0025$ice_curves,2,mean))),
      apply(ice.bart.dps3.x3.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
{segments(x0 = -3, y0 = 0, x1 = 1/4, y1 = 0, col = "red", lwd = 3)
  segments(x0 = 1/4, y0 = 0, x1 = 1/4, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 1/4, y0 = 0.5, x1 = 2/4, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 2/4, y0 = 0.5, x1 = 2/4, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 2/4, y0 = 0.75, x1 = 3/4, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 0.75, x1 = 3/4, y1 = 1, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 1, x1 = 15*sd(x[3,])+mean(x[3,]), y1 = 1, col = "red", lwd = 3)}

dev.off()


png("fried.bart.dps3.icex3_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps3.x3.0025,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("fried.bart.dps3.icex3_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps3.x3.0975,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("fried.bart.dps3.icex3_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps3.x3,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[3]))
dev.off()

png("fried.bart.dps3.icex3_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.dps3.x3),
     ylab = "Derivative Y",
     xlab = expression(x[3]))
dev.off()

png("fried.bart.dps3.icepihat_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps3.pihat,
     ylim = c(-4,6),
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
lines(as.numeric(names(apply(ice.bart.dps3.pihat.0975$ice_curves,2,mean))),
      apply(ice.bart.dps3.pihat.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.dps3.pihat.0025$ice_curves,2,mean))),
      apply(ice.bart.dps3.pihat.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("fried.bart.dps3.icepihat_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps3.pihat.0025,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("fried.bart.dps3.icepihat_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps3.pihat.0975,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("fried.bart.dps3.icepihat_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps3.pihat,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("fried.bart.dps3.icepihat_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.dps3.pihat),
     ylab = "Derivative Y",
     xlab = expression(hat(Pi)(x)))
dev.off()


png("fried.bart.dps3.icex4_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps3.x4,
     ylim = c(-4,6),
     ylab = "Partial Y",
     xlab = expression(x[4]))
lines(as.numeric(names(apply(ice.bart.dps3.x4.0975$ice_curves,2,mean))),
      apply(ice.bart.dps3.x4.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.dps3.x4.0025$ice_curves,2,mean))),
      apply(ice.bart.dps3.x4.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("fried.bart.dps3.icex4_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps3.x4.0025,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("fried.bart.dps3.icex4_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps3.x4.0975,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("fried.bart.dps3.icex4_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps3.x4,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[4]))
dev.off()

png("fried.bart.dps3.icex4_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.dps3.x4),
     ylab = "Derivative Y",
     xlab = expression(x[4]))
dev.off()

png("fried.bart.dps3.icex5_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps3.x5,
     ylim = c(-4,6),
     ylab = "Partial Y",
     xlab = expression(x[5]))
lines(as.numeric(names(apply(ice.bart.dps3.x5.0975$ice_curves,2,mean))),
      apply(ice.bart.dps3.x5.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.dps3.x5.0025$ice_curves,2,mean))),
      apply(ice.bart.dps3.x5.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("fried.bart.dps3.icex5_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps3.x5.0025,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("fried.bart.dps3.icex5_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps3.x5.0975,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("fried.bart.dps3.icex5_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps3.x5,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[5]))
dev.off()

png("fried.bart.dps3.icex5_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.dps3.x5),
     ylab = "Derivative Y",
     xlab = expression(x[5]))
dev.off()


#######
##dps4##
#######

png("fried.bart.dps4.icex1_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps4.x1,
     ylim = c(-4,6),
     ylab = "Partial Y",
     xlab = expression(x[1]))
lines(as.numeric(names(apply(ice.bart.dps4.x1.0975$ice_curves,2,mean))),
      apply(ice.bart.dps4.x1.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.dps4.x1.0025$ice_curves,2,mean))),
      apply(ice.bart.dps4.x1.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("fried.bart.dps4.icex1_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps4.x1.0025,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("fried.bart.dps4.icex1_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps4.x1.0975,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("fried.bart.dps4.icex1_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps4.x1,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[1]))
dev.off()

png("fried.bart.dps4.icex1_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.dps4.x1),
     ylab = "Derivative Y",
     xlab = expression(x[1]))
dev.off()


png("fried.bart.dps4.icex2_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps4.x2,
     ylim = c(-4,6),
     ylab = "Partial Y",
     xlab = expression(x[2]))
lines(as.numeric(names(apply(ice.bart.dps4.x2.0975$ice_curves,2,mean))),
      apply(ice.bart.dps4.x2.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.dps4.x2.0025$ice_curves,2,mean))),
      apply(ice.bart.dps4.x2.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("fried.bart.dps4.icex2_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps4.x2.0025,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("fried.bart.dps4.icex2_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps4.x2.0975,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("fried.bart.dps4.icex2_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps4.x2,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[2]))
dev.off()

png("fried.bart.dps4.icex2_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.dps4.x2),
     ylab = "Derivative Y",
     xlab = expression(x[2]))
dev.off()


png("fried.bart.dps4.icex3_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps4.x3,
     ylim = c(-4,6),
     ylab = "Partial Y",
     xlab = expression(x[3]))
lines(as.numeric(names(apply(ice.bart.dps4.x3.0975$ice_curves,2,mean))),
      apply(ice.bart.dps4.x3.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.dps4.x3.0025$ice_curves,2,mean))),
      apply(ice.bart.dps4.x3.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
{segments(x0 = -3, y0 = 0, x1 = 1/4, y1 = 0, col = "red", lwd = 3)
  segments(x0 = 1/4, y0 = 0, x1 = 1/4, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 1/4, y0 = 0.5, x1 = 2/4, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 2/4, y0 = 0.5, x1 = 2/4, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 2/4, y0 = 0.75, x1 = 3/4, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 0.75, x1 = 3/4, y1 = 1, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 1, x1 = 15*sd(x[3,])+mean(x[3,]), y1 = 1, col = "red", lwd = 3)}
dev.off()


png("fried.bart.dps4.icex3_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps4.x3.0025,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("fried.bart.dps4.icex3_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps4.x3.0975,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("fried.bart.dps4.icex3_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps4.x3,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[3]))
dev.off()

png("fried.bart.dps4.icex3_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.dps4.x3),
     ylab = "Derivative Y",
     xlab = expression(x[3]))
dev.off()

png("fried.bart.dps4.icepihat_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps4.pihat,
     ylim = c(-4,6),
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
lines(as.numeric(names(apply(ice.bart.dps4.pihat.0975$ice_curves,2,mean))),
      apply(ice.bart.dps4.pihat.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.dps4.pihat.0025$ice_curves,2,mean))),
      apply(ice.bart.dps4.pihat.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("fried.bart.dps4.icepihat_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps4.pihat.0025,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("fried.bart.dps4.icepihat_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps4.pihat.0975,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("fried.bart.dps4.icepihat_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps4.pihat,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("fried.bart.dps4.icepihat_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.dps4.pihat),
     ylab = "Derivative Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("fried.bart.dps4.icex4_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps4.x4,
     ylim = c(-4,6),
     ylab = "Partial Y",
     xlab = expression(x[4]))
lines(as.numeric(names(apply(ice.bart.dps4.x4.0975$ice_curves,2,mean))),
      apply(ice.bart.dps4.x4.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.dps4.x4.0025$ice_curves,2,mean))),
      apply(ice.bart.dps4.x4.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("fried.bart.dps4.icex4_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps4.x4.0025,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("fried.bart.dps4.icex4_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps4.x4.0975,
     ylab = "Partial Y",
     xlab = expression(x[4]))
dev.off()

png("fried.bart.dps4.icex4_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps4.x4,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[4]))
dev.off()

png("fried.bart.dps4.icex4_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.dps4.x4),
     ylab = "Derivative Y",
     xlab = expression(x[4]))
dev.off()

png("fried.bart.dps4.icex5_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps4.x5,
     ylim = c(-4,6),
     ylab = "Partial Y",
     xlab = expression(x[5]))
lines(as.numeric(names(apply(ice.bart.dps4.x5.0975$ice_curves,2,mean))),
      apply(ice.bart.dps4.x5.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.dps4.x5.0025$ice_curves,2,mean))),
      apply(ice.bart.dps4.x5.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
abline(h = mean(alpha), col = "red", lwd = 3)
dev.off()


png("fried.bart.dps4.icex5_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps4.x5.0025,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("fried.bart.dps4.icex5_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps4.x5.0975,
     ylab = "Partial Y",
     xlab = expression(x[5]))
dev.off()

png("fried.bart.dps4.icex5_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.dps4.x5,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[5]))
dev.off()

png("fried.bart.dps4.icex5_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.dps4.x5),
     ylab = "Derivative Y",
     xlab = expression(x[5]))
dev.off()

