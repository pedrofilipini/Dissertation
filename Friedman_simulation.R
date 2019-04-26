######################################
##Pedro Henrique Filipini dos Santos##
###University of São Paulo - Brazil###
######################################


##############################################
##Simulation study based on Freidman Example##
##############################################

#Libraries
library(BART)
library(doSNOW)
library(foreach)

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
joinDF <- function(Var, x){
  if(is.matrix(Var)==0){
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


cluster = makeCluster(ncores, type = "SOCK")
registerDoSNOW(cluster)

time <- Sys.time()
results <- foreach(i=1:niters) %dopar%{
  library("BART")
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
  
  #  if(.Platform$OS.type=="unix"){
  #    vanilla = mc.wbart(cbind(t(x),z),y,X[,-(ncol(X)-1)],nskip = burn,ndpost = npost,keepevery = thin, 
  #                       seed = seedaux, mc.cores = mc.cores.detected)
  #  }else{
  vanilla = wbart(cbind(t(x),z),y,X[,-(ncol(X)-1)],nskip = burn,ndpost = npost,keepevery = thin)
  #  }  
  
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
  vanillaCount[[i]] <- joinDF(vanillaCount,PIP(vanilla))
  
  #OracleBART
  set.seed(seedaux)
  
  #  if(.Platform$OS.type=="unix"){
  #    oracle = mc.wbart(cbind(t(x),pihat,z),y,X,nskip = burn,ndpost = npost,keepevery = thin, 
  #                      seed = seedaux, mc.cores = mc.cores.detected)
  #  }else{
  oracle = wbart(cbind(t(x),pihat,z),y,X,nskip = burn,ndpost = npost,keepevery = thin)
  #  }  
  
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
  oracleCount[[i]] <- joinDF(oracleCount,PIP(oracle))
  
  #pihat with pbart
  set.seed(seedaux)
  
  #  if(.Platform$OS.type=="unix"){
  #    fitz = mc.pbart(cbind(t(x)),z,nskip = burn,ndpost = npost,keepevery = pthin, 
  #                    seed = seedaux, mc.cores = mc.cores.detected)
  #  }else{
  fitz = pbart(cbind(t(x)),z,nskip = burn,ndpost = npost,keepevery = pthin)
  #  }
  
  pihat2 = fitz$prob.train.mean #Could use median instead
  X2 = cbind(rbind(cbind(t(x),pihat2),cbind(t(x),pihat2)),c(rep(1,n),rep(0,n)))
  colnames(X2)[ncol(X2)]="z"
  
  rmse_pbart[i] <- sqrt(mean((pihat-pihat2)^2))
  fitzCount[[i]] <- joinDF(fitzCount,PIP(fitz))
  
  
  #PSBART
  set.seed(seedaux)
  
  #  if(.Platform$OS.type=="unix"){
  #    ps = mc.wbart(cbind(t(x),pihat2,z), y, X2, nskip = burn,ndpost = npost,keepevery = thin, 
  #                  seed = seedaux, mc.cores = mc.cores.detected)
  #  }else{
  ps = wbart(cbind(t(x),pihat2,z), y, X2, nskip = burn,ndpost = npost,keepevery = thin)
  #  }
  
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
  psCount[[i]] <- joinDF(psCount,PIP(ps))
  
  
  
  #pihat estimated by GLM
  fitz = glm(z ~ cbind(t(x)), family = binomial())
  pihat4 = fitz$fitted.values
  X4 = cbind(rbind(cbind(t(x),pihat4),cbind(t(x),pihat4)),c(rep(1,n),rep(0,n)))
  colnames(X4)[ncol(X4)]="z"
  
  #PSGLM
  set.seed(seedaux)
  
  #  if(.Platform$OS.type=="unix"){
  #    ps3 = mc.wbart(cbind(t(x),pihat4,z),y,X4,nskip = burn,ndpost = npost,keepevery = thin, 
  #                   seed = seedaux, mc.cores = mc.cores.detected)
  #  }else{
  ps3 = wbart(cbind(t(x),pihat4,z),y,X4,nskip = burn,ndpost = npost,keepevery = thin)
  #  }
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
  ps3Count[[i]] <- joinDF(ps3Count,PIP(ps3))
  
  
  #pihat randomly generated
  pihat5 = runif(n)
  X5 = cbind(rbind(cbind(t(x),pihat5),cbind(t(x),pihat5)),c(rep(1,n),rep(0,n)))
  colnames(X5)[ncol(X5)]="z"
  
  #PSRAND
  set.seed(seedaux)
  
  #  if(.Platform$OS.type=="unix"){
  #    ps4 = mc.wbart(cbind(t(x),pihat5,z),y,X5,nskip = burn,ndpost = npost,keepevery = thin, 
  #                   seed = seedaux, mc.cores = mc.cores.detected)
  #  }else{
  ps4 = wbart(cbind(t(x),pihat5,z),y,X5,nskip = burn,ndpost = npost,keepevery = thin)
  #  }
  
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
  ps4Count[[i]] <- joinDF(ps4Count,PIP(ps4))
  
  
  ####################################
  #BART Models with sparse setting on#
  ####################################
  
  #VanillaBART
  set.seed(seedaux)
  
  #  if(.Platform$OS.type=="unix"){
  #    dvanilla = mc.wbart(cbind(t(x),z),y,X[,-(ncol(X)-1)], sparse = T, nskip = burn,ndpost = npost,keepevery = thin, 
  #                        seed = seedaux, mc.cores = mc.cores.detected)
  #  }else{
  dvanilla = wbart(cbind(t(x),z),y,X[,-(ncol(X)-1)], sparse = T, nskip = burn,ndpost = npost,keepevery = thin)
  #  }
  
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
  vanillaVar[[i]] <- joinDF(vanillaVar,dvanilla$varprob.mean)
  dvanillaCount[[i]] <- joinDF(dvanillaCount,PIP(dvanilla))
  
  #OracleBART
  set.seed(seedaux)
  
  #  if(.Platform$OS.type=="unix"){
  #    doracle = mc.wbart(cbind(t(x),pihat,z),y,X,sparse = T, nskip = burn,ndpost = npost,keepevery = thin, 
  #                       seed = seedaux, mc.cores = mc.cores.detected)
  #  }else{
  doracle = wbart(cbind(t(x),pihat,z),y,X,sparse = T, nskip = burn,ndpost = npost,keepevery = thin)
  #  }
  
  
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
  oracleVar[[i]] <- joinDF(oracleVar,doracle$varprob.mean)
  doracleCount[[i]] <- joinDF(doracleCount,PIP(doracle))
  
  #pihat with pbart
  set.seed(seedaux)
  #  if(.Platform$OS.type=="unix"){
  #    dfitz = mc.pbart(cbind(t(x)),z,sparse = T, nskip = burn,ndpost = npost,keepevery = pthin, 
  #                     seed = seedaux, mc.cores = mc.cores.detected)
  #  }else{
  dfitz = pbart(cbind(t(x)),z,sparse = T, nskip = burn,ndpost = npost,keepevery = pthin)
  #  }
  
  dpihat2 = dfitz$prob.train.mean #Could use median instead
  
  #Calculating Number of times Variable was used
  rmse_dpbart[i] <- sqrt(mean((pihat-dpihat2)^2))
  fitzVar[[i]] <- joinDF(fitzVar,dfitz$varprob.mean)
  dfitzCount[[i]] <- joinDF(dfitzCount,PIP(dfitz))
  
  dX2 = cbind(rbind(cbind(t(x),dpihat2),cbind(t(x),dpihat2)),c(rep(1,n),rep(0,n)))
  colnames(dX2)[ncol(dX2)]="z"
  
  #PSBART
  set.seed(seedaux)
  
  #  if(.Platform$OS.type=="unix"){
  #    dps = mc.wbart(cbind(t(x),dpihat2,z), y, dX2, sparse = T, nskip = burn,ndpost = npost,keepevery = thin, 
  #                   seed = seedaux, mc.cores = mc.cores.detected)
  #  }else{
  dps = wbart(cbind(t(x),dpihat2,z), y, dX2, sparse = T, nskip = burn,ndpost = npost,keepevery = thin)
  #  }
  
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
  psVar[[i]] <- joinDF(psVar,dps$varprob.mean)
  dpsCount[[i]] <- joinDF(dpsCount,PIP(dps))
  
  #PSGLM
  set.seed(seedaux)
  
  #  if(.Platform$OS.type=="unix"){
  #    dps3 = mc.wbart(cbind(t(x),pihat4,z),y,sparse = T, X4,nskip = burn,ndpost = npost,keepevery = thin, 
  #                    seed = seedaux, mc.cores = mc.cores.detected)
  #  }else{
  dps3 = wbart(cbind(t(x),pihat4,z),y,sparse = T, X4,nskip = burn,ndpost = npost,keepevery = thin)
  #  }
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
  ps3Var[[i]] <- joinDF(ps3Var,dps3$varprob.mean)
  dps3Count[[i]] <- joinDF(dps3Count,PIP(dps3))
  
  #PSRAND
  set.seed(seedaux)
  
  #  if(.Platform$OS.type=="unix"){
  #    dps4 = mc.wbart(cbind(t(x),pihat5,z),y,X5,sparse = T, nskip = burn,ndpost = npost,keepevery = thin, 
  #                    seed = seedaux, mc.cores = mc.cores.detected)
  #  }else{
  dps4 = wbart(cbind(t(x),pihat5,z),y,X5,sparse = T, nskip = burn,ndpost = npost,keepevery = thin)
  #  }
  
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
  ps4Var[[i]] <- joinDF(ps4Var,dps4$varprob.mean)
  dps4Count[[i]] <- joinDF(dps4Count,PIP(dps4))
  
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
  
  save(aux, file = paste0("results_sim_new_fried",i,".RData"))
  
  gc()
  
  return(aux)
  
}

save(results, file = paste0("results_sim_new_fried.RData"))

stopCluster(cluster)
