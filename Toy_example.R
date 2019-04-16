######################################
##Pedro Henrique Filipini dos Santos##
###University of S?o Paulo - Brazil###
######################################


##########################
##Boxplot and ICE Images##
##########################

#Basically I used the first iteration of the loop as a base for my
#boxplot images and my ICEs plots.

#Libraries
library(BART)
library(bcf)
library(ICEbox)


#number of iterations (which I have set as 1 instead of 1000)
niters <- 1
burn <- 5000
post<- 500
thin <- 150

set.seed(99)

n = 200

q <- numeric(0)
z <- numeric(0)
ps <- numeric(0)
x1 <- numeric(0)
x2 <- numeric(0)
x3 <- numeric(0)
y <- numeric(0)

#Gerando vari?veis
x1 <- rnorm(n, 0, 1)
x2 <- rnorm(n, 0, 1)
x3 <- rnorm(n, 0, 1)

x <- t(data.frame(x1=x1,x2=x2,
                  x3=x3))


for(i in 1:n){
  #Rule to generate the propensity score
  if(x1[i]<x2[i]){ #indicator function (1 or -1)
    q[i] <- 1
  }else{
    q[i] <- -1
  }
  
  ps[i] <- pnorm(q[i]) #P(z[i]=1|x1[i],x2[i]) the real propensity score
  
  p <- runif(1, 0, 1) #runs an uniform in order to sample from probability
  
  if(p <= ps[i]){ #uses the sample to assign a value to z[i] (1 or 0)
    z[i] <- 1
  }else{
    z[i] <- 0
  }
  
}

#Alpha its the true treatment effect
alpha =  0.5*(t(x)[,3] > -3/4) + 0.25*(t(x)[,3] > 0) + 0.25*(t(x)[,3]>3/4)


#simulation starts
for(i in 1:niters){
  #setting different seeds for different iterations
  seedaux <- (7565 + i*5)
  
  set.seed(seedaux)
  
  #Generating the response Y
  y <- q + z*alpha + 0.5*rnorm(n)
  
  #pihat is the true p-score
  pihat <- ps
  
  #Test Matrix
  X = cbind(rbind(cbind(t(x),pihat),cbind(t(x),pihat)),c(rep(1,n),rep(0,n)))
  colnames(X)[ncol(X)]="z"
  
  #VanillaBART
  set.seed(seedaux)
  vanilla = wbart(cbind(t(x),z),y,X[,-(ncol(X)-1)],
                  nskip = burn,ndpost = post, keepevery = thin)
  
  #OracleBART
  set.seed(seedaux)
  oracle = wbart(cbind(t(x),pihat,z),y,X,nskip = burn,ndpost = post, keepevery = thin)
  
  
  #pihat with pbart
  fitz = pbart(t(x),z,nskip = burn,ndpost = post, keepevery = thin)
  pihat2 = fitz$prob.train.mean #Could use median instead
  X2 = cbind(rbind(cbind(t(x),pihat2),cbind(t(x),pihat2)),c(rep(1,n),rep(0,n)))
  colnames(X2)[ncol(X2)]="z"
  
  #PSBART
  set.seed(seedaux)
  ps = wbart(cbind(t(x),pihat2,z), y, X2, nskip = burn,ndpost = post, keepevery = thin)
  
  #Calculating estimated ATE
  ate_est_vanilla = rowMeans(vanilla$yhat.test[,1:n] - vanilla$yhat.test[,(n+1):(2*n)])
  ate_est_oracle = rowMeans(oracle$yhat.test[,1:n] - oracle$yhat.test[,(n+1):(2*n)])
  ate_est_ps = rowMeans(ps$yhat.test[,1:n] - ps$yhat.test[,(n+1):(2*n)])
  
  #Calculating estimated ITE
  ite_est_vanilla = vanilla$yhat.test[,1:n] - vanilla$yhat.test[,(n+1):(2*n)]
  ite_est_oracle = oracle$yhat.test[,1:n] - oracle$yhat.test[,(n+1):(2*n)]
  ite_est_ps = ps$yhat.test[,1:n] - ps$yhat.test[,(n+1):(2*n)]
  
  
  #pihat estimated by GLM
  fitzglm = glm(z ~ t(x), family = binomial())
  pihat4 = fitzglm$fitted.values
  X4 = cbind(rbind(cbind(t(x),pihat4),cbind(t(x),pihat4)),c(rep(1,n),rep(0,n)))
  colnames(X4)[ncol(X4)]="z"
  
  #PSGLM
  set.seed(seedaux)
  ps3 = wbart(cbind(t(x),pihat4,z),y,X4,nskip = burn,ndpost = post, keepevery = thin)
  
  #Estimated Treatment Effects
  ate_est_psglm = rowMeans(ps3$yhat.test[,1:n] - ps3$yhat.test[,(n+1):(2*n)])
  ite_est_psglm = ps3$yhat.test[,1:n] - ps3$yhat.test[,(n+1):(2*n)]
  
  
  #pihat randomly generated
  pihat5 = runif(n)
  X5 = cbind(rbind(cbind(t(x),pihat5),cbind(t(x),pihat5)),c(rep(1,n),rep(0,n)))
  colnames(X5)[ncol(X5)]="z"
  
  #PSRAND
  set.seed(seedaux)
  ps4 = wbart(cbind(t(x),pihat5,z),y,X5,nskip = burn,ndpost = post, keepevery = thin)
  
  #Estimated Treatment Effects
  ate_est_psrand = rowMeans(ps4$yhat.test[,1:n] - ps4$yhat.test[,(n+1):(2*n)])
  ite_est_psrand = ps4$yhat.test[,1:n] - ps4$yhat.test[,(n+1):(2*n)]
  
  
  
  
  
  #BCF with BART (pihat2)
  set.seed(seedaux)
  fitbcf2 = bcf(y, z, t(x), t(x), pihat2, burn, post, nthin = thin, include_pi="both")
  
  #Estimated Treatment Effects
  ate_est_bartbcf = rowMeans(fitbcf2$tau)
  ite_est_bartbcf = colMeans(fitbcf2$tau)
  
  
  #Oracle BCF (pihat)
  set.seed(seedaux)
  fitbcf3 = bcf(y, z, t(x), t(x), pihat, burn, post, nthin = thin, include_pi="both")
  
  #Estimated Treatment Effects
  ate_est_oraclebcf = rowMeans(fitbcf3$tau)
  ite_est_oraclebcf = colMeans(fitbcf3$tau)
  
  
  #BCF with GLM (pihat4)
  set.seed(seedaux)
  fitbcf4 = bcf(y, z, t(x), t(x), pihat4, burn, post, nthin = thin, include_pi="both")
  
  #Estimated Treatment Effects
  ate_est_glmbcf = rowMeans(fitbcf4$tau)
  ite_est_glmbcf = colMeans(fitbcf4$tau)
  
  
  #BCF with random (pihat5)
  set.seed(seedaux)
  fitbcf5 = bcf(y, z, t(x), t(x), pihat5, burn, post, nthin = thin, include_pi="both")
  
  #Estimated Treatment Effects
  ate_est_randbcf = rowMeans(fitbcf5$tau)
  ite_est_randbcf = colMeans(fitbcf5$tau)
  
  
  ######
  #Plot#
  ######
  ATE <- data.frame(ate_est_vanilla, ate_est_oracle, ate_est_ps,
                    ate_est_psglm, ate_est_psrand, ate_est_oraclebcf, 
                    ate_est_bartbcf, ate_est_glmbcf, ate_est_randbcf)
  #png(filename="boxplot_all.png")
  boxplot(ATE, las = 2,
          ylab = "Average Treatment Effect",
          names = c("Vanilla", "Oracle", "PS-BART", 
                    "GLM-BART","Rand-BART", "Oracle-BCF", 
                    "BART-BCF", "GLM-BCF", "Rand-BCF"))
  #True ATE
  abline(h = mean(alpha), col = "black", lwd = 2)
  #dev.off()
}


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

ice.bcf.bart.x1 <- ice(fitbcf2aux, Xice, predictor = "x1", 
                       predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.x1, file="ice_bcf.bart.x1.RData")

ice.bcf.bart.x1.0025 <- ice(fitbcf2aux, Xice, predictor = "x1", 
                            predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.x1.0025, file="ice_bcf.bart.x1.0025.RData")

ice.bcf.bart.x1.0975 <- ice(fitbcf2aux, Xice, predictor = "x1", 
                            predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.x1.0975, file="ice_bcf.bart.x1.0975.RData")

###


ice.bcf.bart.x2 <- ice(fitbcf2aux, Xice, predictor = "x2", 
                       predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.x2, file="ice_bcf.bart.x2.RData")

ice.bcf.bart.x2.0025 <- ice(fitbcf2aux, Xice, predictor = "x2", 
                            predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.x2.0025, file="ice_bcf.bart.x2.0025.RData")

ice.bcf.bart.x2.0975 <- ice(fitbcf2aux, Xice, predictor = "x2", 
                            predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.x2.0975, file="ice_bcf.bart.x2.0975.RData")

###

ice.bcf.bart.x3 <- ice(fitbcf2aux, Xice, predictor = "x3", 
                       predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.x3, file="ice_bcf.bart.x3.RData")

ice.bcf.bart.x3.0025 <- ice(fitbcf2aux, Xice, predictor = "x3", 
                            predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.x3.0025, file="ice_bcf.bart.x3.0025.RData")

ice.bcf.bart.x3.0975 <- ice(fitbcf2aux, Xice, predictor = "x3", 
                            predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.x3.0975, file="ice_bcf.bart.x3.0975.RData")

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

ice.bcf.oracle.x1 <- ice(fitbcf3aux, Xice, predictor = "x1", 
                         predictfcn = pred.BCF.ATE)
save(ice.bcf.oracle.x1, file="ice_bcf.oracle.x1.RData")

ice.bcf.oracle.x1.0025 <- ice(fitbcf3aux, Xice, predictor = "x1", 
                              predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.oracle.x1.0025, file="ice_bcf.oracle.x1.0025.RData")

ice.bcf.oracle.x1.0975 <- ice(fitbcf3aux, Xice, predictor = "x1", 
                              predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.oracle.x1.0975, file="ice_bcf.oracle.x1.0975.RData")

###

ice.bcf.oracle.x2 <- ice(fitbcf3aux, Xice, predictor = "x2", 
                         predictfcn = pred.BCF.ATE)
save(ice.bcf.oracle.x2, file="ice_bcf.oracle.x2.RData")

ice.bcf.oracle.x2.0025 <- ice(fitbcf3aux, Xice, predictor = "x2", 
                              predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.oracle.x2.0025, file="ice_bcf.oracle.x2.0025.RData")

ice.bcf.oracle.x2.0975 <- ice(fitbcf3aux, Xice, predictor = "x2", 
                              predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.oracle.x2.0975, file="ice_bcf.oracle.x2.0975.RData")

###


ice.bcf.oracle.x3 <- ice(fitbcf3aux, Xice, predictor = "x3", 
                         predictfcn = pred.BCF.ATE)
save(ice.bcf.oracle.x3, file="ice_bcf.oracle.x3.RData")

ice.bcf.oracle.x3.0025 <- ice(fitbcf3aux, Xice, predictor = "x3", 
                              predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.oracle.x3.0025, file="ice_bcf.oracle.x3.0025.RData")

ice.bcf.oracle.x3.0975 <- ice(fitbcf3aux, Xice, predictor = "x3", 
                              predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.oracle.x3.0975, file="ice_bcf.oracle.x3.0975.RData")

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

ice.bcf.glm.x1 <- ice(fitbcf4aux, Xice, predictor = "x1", 
                      predictfcn = pred.BCF.ATE)
save(ice.bcf.glm.x1, file="ice_bcf.glm.x1.RData")

ice.bcf.glm.x1.0025 <- ice(fitbcf4aux, Xice, predictor = "x1", 
                           predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.glm.x1.0025, file="ice_bcf.glm.x1.0025.RData")

ice.bcf.glm.x1.0975 <- ice(fitbcf4aux, Xice, predictor = "x1", 
                           predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.glm.x1.0975, file="ice_bcf.glm.x1.0975.RData")

###

ice.bcf.glm.x2 <- ice(fitbcf4aux, Xice, predictor = "x2", 
                      predictfcn = pred.BCF.ATE)
save(ice.bcf.glm.x2, file="ice_bcf.glm.x2.RData")

ice.bcf.glm.x2.0025 <- ice(fitbcf4aux, Xice, predictor = "x2", 
                           predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.glm.x2.0025, file="ice_bcf.glm.x2.0025.RData")

ice.bcf.glm.x2.0975 <- ice(fitbcf4aux, Xice, predictor = "x2", 
                           predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.glm.x2.0975, file="ice_bcf.glm.x2.0975.RData")

###

ice.bcf.glm.x3 <- ice(fitbcf4aux, Xice, predictor = "x3", 
                      predictfcn = pred.BCF.ATE)
save(ice.bcf.glm.x3, file="ice_bcf.glm.x3.RData")

ice.bcf.glm.x3.0025 <- ice(fitbcf4aux, Xice, predictor = "x3", 
                           predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.glm.x3.0025, file="ice_bcf.glm.x3.0025.RData")

ice.bcf.glm.x3.0975 <- ice(fitbcf4aux, Xice, predictor = "x3", 
                           predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.glm.x3.0975, file="ice_bcf.glm.x3.0975.RData")

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

ice.bcf.rand.x1 <- ice(fitbcf5aux, Xice, predictor = "x1", 
                       predictfcn = pred.BCF.ATE)
save(ice.bcf.rand.x1, file="ice_bcf.rand.x1.RData")

ice.bcf.rand.x1.0025 <- ice(fitbcf5aux, Xice, predictor = "x1", 
                            predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.rand.x1.0025, file="ice_bcf.rand.x1.0025.RData")

ice.bcf.rand.x1.0975 <- ice(fitbcf5aux, Xice, predictor = "x1", 
                            predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.rand.x1.0975, file="ice_bcf.rand.x1.0975.RData")

###

ice.bcf.rand.x2 <- ice(fitbcf5aux, Xice, predictor = "x2", 
                       predictfcn = pred.BCF.ATE)
save(ice.bcf.rand.x2, file="ice_bcf.rand.x2.RData")

ice.bcf.rand.x2.0025 <- ice(fitbcf5aux, Xice, predictor = "x2", 
                            predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.rand.x2.0025, file="ice_bcf.rand.x2.0025.RData")

ice.bcf.rand.x2.0975 <- ice(fitbcf5aux, Xice, predictor = "x2", 
                            predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.rand.x2.0975, file="ice_bcf.rand.x2.0975.RData")

###


ice.bcf.rand.x3 <- ice(fitbcf5aux, Xice, predictor = "x3", 
                       predictfcn = pred.BCF.ATE)
save(ice.bcf.rand.x3, file="ice_bcf.rand.x3.RData")

ice.bcf.rand.x3.0025 <- ice(fitbcf5aux, Xice, predictor = "x3", 
                            predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.rand.x3.0025, file="ice_bcf.rand.x3.0025.RData")

ice.bcf.rand.x3.0975 <- ice(fitbcf5aux, Xice, predictor = "x3", 
                            predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.rand.x3.0975, file="ice_bcf.rand.x3.0975.RData")

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

ice.bart.vanilla.x1 <- ice(vanilla, Xice, predictor = "x1", 
                           predictfcn = pred.BART.ATE)
save(ice.bart.vanilla.x1, file="ice_vanilla_x1.RData")

ice.bart.vanilla.x1.0025 <- ice(vanilla, Xice, predictor = "x1", 
                                predictfcn = pred.BART.ATE.0025)
save(ice.bart.vanilla.x1.0025, file="ice_vanilla_x1.0025.RData")

ice.bart.vanilla.x1.0975 <- ice(vanilla, Xice, predictor = "x1", 
                                predictfcn = pred.BART.ATE.0975)
save(ice.bart.vanilla.x1.0975, file="ice_vanilla_x1.0975.RData")

####

ice.bart.vanilla.x2 <- ice(vanilla, Xice, predictor = "x2", 
                           predictfcn = pred.BART.ATE)
save(ice.bart.vanilla.x2, file="ice_vanilla_x2.RData")

ice.bart.vanilla.x2.0025 <- ice(vanilla, Xice, predictor = "x2", 
                                predictfcn = pred.BART.ATE.0025)
save(ice.bart.vanilla.x2.0025, file="ice_vanilla_x2.0025.RData")

ice.bart.vanilla.x2.0975 <- ice(vanilla, Xice, predictor = "x2", 
                                predictfcn = pred.BART.ATE.0975)
save(ice.bart.vanilla.x2.0975, file="ice_vanilla_x2.0975.RData")

###

ice.bart.vanilla.x3 <- ice(vanilla, Xice, predictor = "x3", 
                           predictfcn = pred.BART.ATE)
save(ice.bart.vanilla.x3, file="ice_vanilla_x3.RData")

ice.bart.vanilla.x3.0025 <- ice(vanilla, Xice, predictor = "x3", 
                                predictfcn = pred.BART.ATE.0025)
save(ice.bart.vanilla.x3.0025, file="ice_vanilla_x3.0025.RData")

ice.bart.vanilla.x3.0975 <- ice(vanilla, Xice, predictor = "x3", 
                                predictfcn = pred.BART.ATE.0975)
save(ice.bart.vanilla.x3.0975, file="ice_vanilla_x3.0975.RData")

###

#Oracle
Xice = cbind(t(x),pihat,z)

ice.bart.oracle.x1 <- ice(oracle, Xice, predictor = "x1", 
                          predictfcn = pred.BART.ATE)
save(ice.bart.oracle.x1, file="ice_oracle_x1.RData")

ice.bart.oracle.x1.0025 <- ice(oracle, Xice, predictor = "x1", 
                               predictfcn = pred.BART.ATE.0025)
save(ice.bart.oracle.x1.0025, file="ice_oracle_x1.0025.RData")

ice.bart.oracle.x1.0975 <- ice(oracle, Xice, predictor = "x1", 
                               predictfcn = pred.BART.ATE.0975)
save(ice.bart.oracle.x1.0975, file="ice_oracle_x1.0975.RData")

####

ice.bart.oracle.x2 <- ice(oracle, Xice, predictor = "x2", 
                          predictfcn = pred.BART.ATE)
save(ice.bart.oracle.x2, file="ice_oracle_x2.RData")

ice.bart.oracle.x2.0025 <- ice(oracle, Xice, predictor = "x2", 
                               predictfcn = pred.BART.ATE.0025)
save(ice.bart.oracle.x2.0025, file="ice_oracle_x2.0025.RData")

ice.bart.oracle.x2.0975 <- ice(oracle, Xice, predictor = "x2", 
                               predictfcn = pred.BART.ATE.0975)
save(ice.bart.oracle.x2.0975, file="ice_oracle_x2.0975.RData")

###

ice.bart.oracle.x3 <- ice(oracle, Xice, predictor = "x3", 
                          predictfcn = pred.BART.ATE)
save(ice.bart.oracle.x3, file="ice_oracle_x3.RData")

ice.bart.oracle.x3.0025 <- ice(oracle, Xice, predictor = "x3", 
                               predictfcn = pred.BART.ATE.0025)
save(ice.bart.oracle.x3.0025, file="ice_oracle_x3.0025.RData")

ice.bart.oracle.x3.0975 <- ice(oracle, Xice, predictor = "x3", 
                               predictfcn = pred.BART.ATE.0975)
save(ice.bart.oracle.x3.0975, file="ice_oracle_x3.0975.RData")

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

ice.bart.ps.x1 <- ice(ps, Xice, predictor = "x1", 
                      predictfcn = pred.BART.ATE)
save(ice.bart.ps.x1, file="ice_ps_x1.RData")

ice.bart.ps.x1.0025 <- ice(ps, Xice, predictor = "x1", 
                           predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps.x1.0025, file="ice_ps_x1.0025.RData")

ice.bart.ps.x1.0975 <- ice(ps, Xice, predictor = "x1", 
                           predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps.x1.0975, file="ice_ps_x1.0975.RData")

####

ice.bart.ps.x2 <- ice(ps, Xice, predictor = "x2", 
                      predictfcn = pred.BART.ATE)
save(ice.bart.ps.x2, file="ice_ps_x2.RData")

ice.bart.ps.x2.0025 <- ice(ps, Xice, predictor = "x2", 
                           predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps.x2.0025, file="ice_ps_x2.0025.RData")

ice.bart.ps.x2.0975 <- ice(ps, Xice, predictor = "x2", 
                           predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps.x2.0975, file="ice_ps_x2.0975.RData")

###

ice.bart.ps.x3 <- ice(ps, Xice, predictor = "x3", 
                      predictfcn = pred.BART.ATE)
save(ice.bart.ps.x3, file="ice_ps_x3.RData")

ice.bart.ps.x3.0025 <- ice(ps, Xice, predictor = "x3", 
                           predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps.x3.0025, file="ice_ps_x3.0025.RData")

ice.bart.ps.x3.0975 <- ice(ps, Xice, predictor = "x3", 
                           predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps.x3.0975, file="ice_ps_x3.0975.RData")

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

ice.bart.ps3.x1 <- ice(ps3, Xice, predictor = "x1", 
                       predictfcn = pred.BART.ATE)
save(ice.bart.ps3.x1, file="ice_ps3_x1.RData")

ice.bart.ps3.x1.0025 <- ice(ps3, Xice, predictor = "x1", 
                            predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps3.x1.0025, file="ice_ps3_x1.0025.RData")

ice.bart.ps3.x1.0975 <- ice(ps3, Xice, predictor = "x1", 
                            predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps3.x1.0975, file="ice_ps3_x1.0975.RData")

####

ice.bart.ps3.x2 <- ice(ps3, Xice, predictor = "x2", 
                       predictfcn = pred.BART.ATE)
save(ice.bart.ps3.x2, file="ice_ps3_x2.RData")

ice.bart.ps3.x2.0025 <- ice(ps3, Xice, predictor = "x2", 
                            predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps3.x2.0025, file="ice_ps3_x2.0025.RData")

ice.bart.ps3.x2.0975 <- ice(ps3, Xice, predictor = "x2", 
                            predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps3.x2.0975, file="ice_ps3_x2.0975.RData")

###

ice.bart.ps3.x3 <- ice(ps3, Xice, predictor = "x3", 
                       predictfcn = pred.BART.ATE)
save(ice.bart.ps3.x3, file="ice_ps3_x3.RData")

ice.bart.ps3.x3.0025 <- ice(ps3, Xice, predictor = "x3", 
                            predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps3.x3.0025, file="ice_ps3_x3.0025.RData")

ice.bart.ps3.x3.0975 <- ice(ps3, Xice, predictor = "x3", 
                            predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps3.x3.0975, file="ice_ps3_x3.0975.RData")

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

ice.bart.ps4.x1 <- ice(ps4, Xice, predictor = "x1", 
                       predictfcn = pred.BART.ATE)
save(ice.bart.ps4.x1, file="ice_ps4_x1.RData")

ice.bart.ps4.x1.0025 <- ice(ps4, Xice, predictor = "x1", 
                            predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps4.x1.0025, file="ice_ps4_x1.0025.RData")

ice.bart.ps4.x1.0975 <- ice(ps4, Xice, predictor = "x1", 
                            predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps4.x1.0975, file="ice_ps4_x1.0975.RData")

####

ice.bart.ps4.x2 <- ice(ps4, Xice, predictor = "x2", 
                       predictfcn = pred.BART.ATE)
save(ice.bart.ps4.x2, file="ice_ps4_x2.RData")

ice.bart.ps4.x2.0025 <- ice(ps4, Xice, predictor = "x2", 
                            predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps4.x2.0025, file="ice_ps4_x2.0025.RData")

ice.bart.ps4.x2.0975 <- ice(ps4, Xice, predictor = "x2", 
                            predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps4.x2.0975, file="ice_ps4_x2.0975.RData")

###

ice.bart.ps4.x3 <- ice(ps4, Xice, predictor = "x3", 
                       predictfcn = pred.BART.ATE)
save(ice.bart.ps4.x3, file="ice_ps4_x3.RData")

ice.bart.ps4.x3.0025 <- ice(ps4, Xice, predictor = "x3", 
                            predictfcn = pred.BART.ATE.0025)
save(ice.bart.ps4.x3.0025, file="ice_ps4_x3.0025.RData")

ice.bart.ps4.x3.0975 <- ice(ps4, Xice, predictor = "x3", 
                            predictfcn = pred.BART.ATE.0975)
save(ice.bart.ps4.x3.0975, file="ice_ps4_x3.0975.RData")

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


#Calculating kappa function

ke <- function(model){
  #salva os cutpoints
  cuts <- model$treedraws$cutpoints 
  
  #Putting trees in a list
  trees <- model$treedraws$trees
  trees <- strsplit(trees, "\n")
  trees <- strsplit(trees[[1]], "\\s+")
  trees <- lapply(trees, function(y) unlist(lapply(y, as.numeric)))
  
  #rm(model)
  gc()
  
  #Getting the main information
  post <- trees[[1]][1] #posterior size
  m <- trees[[1]][2] #number of trees
  p <- trees[[1]][3] #number of variables
  
  ###################################
  ###################################
  #This part is temporary
  
  #Verify if the rule 0 0 has been used
  r <- 0 #tree counter
  s <- numeric(0) #final node counter
  n <- numeric(0) #inner node counter
  t <- numeric(0) #node counter
  
  
  leafs <- vector("list", m*post) #list of final nodes
  nodes <- vector("list", m*post) #list of inner nodes
  predaux <- vector("list", post) #list of positions os inner nodes
  predaux <- lapply(predaux, function(x){vector("list", p)})
  countpost <- 1 #inicializing variable countpost
  
  for(i in 1:length(trees)){ #Runs through every posterior tree
    
    if(length(trees[[i]]) == 1){ #If it is the beginning of a tree
      
      r <- r+1 #Tree counter
      s[r] <- 0 #Inicializing final node counter
      n[r] <- 0 #Inicializing inner node counter
      t[r] <- trees[[i]] #Total node counter
      
      for(j in (i+1):(i+(trees[[i]]))){ #Visiting nodes
        
        if((trees[[j]][2]==0) & (trees[[j]][3]==0)){ #Condition to final node
          trees[[j]][5] <- 1
          s[r] <- s[r]+1 #Counter of final node by tree
          leafs[[r]][(s[r])] <- trees[[j]][1] #Save the tree IDs in the format of a list
        }else{
          trees[[j]][5] <- 0
          n[r] <- n[r]+1 #Counter of inner node by tree
          nodes[[r]][(n[r])] <- trees[[j]][1] #Save the inner IDs in a list
          predaux[[countpost]][[((trees[[j]][2])+1)]] <- c(predaux[[countpost]][[((trees[[j]][2])+1)]],cuts[[((trees[[j]][2])+1)]][((trees[[j]][3])+1)]) #Save the local of the IDs
        }
      }
      i <- i+trees[[i]] #Make the variable i go to the begin of the next tree
      if(r%%m==0){
        countpost <- countpost + 1
      }
    }
  }
  
  predaux <- lapply(predaux,function(x){lapply(x,function(y){sort(unique(y))})}) #Arruma lista de regras
  
  for(a in 1:length(predaux)){
    for(b in 1:(lengths(predaux)[1])){
      if(is.null(predaux[[a]][[b]])){
        predaux[[a]][[b]] <- 0 #If the variable was not used, include it here
      }else{
        predaux[[a]][[b]] <- c((predaux[[a]][[b]][1]-1),predaux[[a]][[b]]) #Adiciona elemento final na lista
      }
    }
  }
  
  
  use <- if(sum((t+1)/2-s)==0){
    use <- 1 #Can use the rule
  }else{
    use <- 0 #Cannot use the rule (fix temporary part of the code)
  }
  
  
  #Function to test if it is final node
  
  isleaf <- function(id,i,trees){
    nodes <- trees[[i]] #number of nodes on tree i
    local <- 0 #id location
    
    for(q in 1:nodes){
      if(trees[[i+q]][1]==id){ #walking through the tree
        local <- i+q #local where the id is
      }
    }
    
    if(local==0){
      return("Error: Could not find node") #Wrong argument: Node does not exist
    }
    
    if(((trees[[local]][2]==0) & (trees[[local]][3]==0)) & (use==1)){ #if it is final node
      return(1)
    }else if(((trees[[local]][2]==0) & (trees[[local]][3]==0)) & (isid(2*id,i,trees)==0)){
      return(1)
    }else{
      return(0)
    }
    
  }
  
  
  
  ###################################
  ###################################
  
  #Verify if id exists
  isid <- function(id,i,trees){
    nodes <- trees[[i]] #number of nodes on tree
    local <- 0 #local of the id
    
    for(q in 1:nodes){
      if(trees[[i+q]][1]==id){ #walks through the tree
        local <- i+q #local where id is
      }
    }
    
    if(local==0){ #If id does not exist
      return(0) 
    }else{
      return(1)
    }
  }
  
  #Verifica a regra
  isleft <- function(id,i,trees,pred,cuts){
    nodes <- trees[[i]] #number of nodes on tree
    local <- 0 #location of the id
    
    for(q in 1:nodes){
      if(trees[[i+q]][1]==id){ #walks through the tree
        local <- i+q #place where the id is
      }
    }
    
    if(local==0){
      return("Error: Could not find node") #Some error in the arguments
    }
    
    if(pred[((trees[[local]][2])+1)]<cuts[[((trees[[local]][2])+1)]][(trees[[local]][3]+1)]){ #verify rule
      return(1)
    }else{
      return(0)
    }
    
  }
  
  treeprev <- function(id,i,trees){
    nodes <- trees[[i]] #n?mero de n?s na ?rvore
    local <- 0 #local que o id est?
    
    for(q in 1:nodes){
      if(trees[[i+q]][1]==id){ #caminha pela ?rvore
        local <- i+q #lugar onde est? o id em quest?o
      }
    }
    
    if(local==0){
      return("Error: Could not find node") #caso tenha colocado argumento errado
    }
    
    return(trees[[local]][4])
  }
  
  ke <- vector("numeric",post) #vetor com a posteriori que quero
  
  compareNodes <- function(x,k,trees,cuts,leafs,h){
    
    id <- 1 #id do in?cio da ?rvore
    r <- k+1
    
    while(trees[[r]][5]==0){
      l <- isleft(id,k,trees,x,cuts) #verifica a regra
      if(l==1){
        id = 2*id #manda pro n? da esquerda
      }else if(l==0){
        id = 2*id +1 #manda pro n? da direita
      }
      
      nodes <- trees[[k]] #n?mero de n?s na ?rvore
      local <- 0 #location of the id
      
      for(q in 1:nodes){
        if(trees[[k+q]][1]==id){ #walks through the tree
          local <- k+q #place where the id is
        }
      }
      
      if(local==0){
        return("Error: Could not find node") #Some error in the arguments
      }
      r <- local
    }
    
    return(list(as.numeric(leafs[[(m*(h-1)+tr)]]==id))) #Vejo em qual n? estou
    #return(treeprev(id,k,trees))
  }
  
  
  for(h in 1:post){
    
    if(h==1){
      print(paste0("MCMC"))
      
      k <- 2 #lugar onde est? a primeira ?rvore
    }
    
    pred <- expand.grid(predaux[[h]])
    print(nrow(pred))
    
    tr <- 1 #conta as ?rvores e reseta a cada posteriori
    
    Raux <- list(0) #Cria matriz auxiliar para parti??es
    
    while(tr<=m){
      
      if(tr%%10==0){
        print(paste0("Post: ",h," out of ",post," ; Tree: ",tr," out of ",m))
      }
      
      Raux[[tr]] <- apply(pred, 1, compareNodes,k=k,trees=trees,cuts=cuts,leafs=leafs,h=h)
      
      k <- trees[[k]]+k+1 #vai pro come?o da pr?xima ?rvore
      
      tr <- tr+1 #vai pra pr?xima ?rvore
      
      if(tr%%10==0){
        gc()
      }
      
    }
    
    #Reduce(`+`,Raux) #Somando quando estou estimando yhat
    #Reduce(`+`,Raux)+mean(y)-predict(model, newdata = pred)[3,]
    
    Raux <- lapply(Raux, function(x){data.frame(matrix(unlist(x), ncol=lengths(x[[1]]), byrow=T))})
    R <- data.frame(matrix(unlist(Raux), nrow=nrow(pred), byrow=F))
    R <- unique(R)
    R <- svd(as.matrix(R))$d
    
    ke[h] <- max(R)/min(R[abs(R)>0.0000001])
    
    rm(R)
    
    #if(h%%100==0){
    print(paste0("done ", h," (out of ",post,")")) #Printar de 100 em 100
    #}
    rm(pred)
    
  }
  
  #boxplot(ke) #boxplot do que interessa
  return(ke)
}
kvanilla <- ke(vanilla)
koracle <- ke(oracle)
kfitbcf2 <- ke(kfitbcf2)
kfitbcf3 <- ke(kfitbcf3)
kfitbcf4 <- ke(kfitbcf4)
kfitbcf5 <- ke(kfitbcf5)





#######################
#######Graphics########
#######################

ATE <- data.frame(ate_est_oracle,  ate_est_oraclebcf,
                  ate_est_vanilla, ate_est_ps, ate_est_bartbcf,
                  ate_est_psglm, ate_est_glmbcf, ate_est_psrand,
                  ate_est_randbcf)
png("toy01.png", width = 7, height = 7, units = 'in', res = 100)
boxplot(ATE, las = 2,
        ylab = "CATE",
        names = c("Oracle", "Oracle-BCF",
                  "Vanilla", "PS-BART", "BART-BCF",
                  "GLM-BART", "GLM-BCF", "Rand-BART",
                  "Rand-BCF"), 
        cex = 0.5,
        cex.lab = 0.8,
        cex.axis = 0.7)
#True ATE
abline(h = mean(alpha), col = "black", lwd = 2)
abline(v = 2.5, col = "gray10", lwd = 2, lty = 2)
dev.off()

#################################
##Assuring Convergence of pbart##
#################################
p <- 3
i <- floor(seq(1, n, length.out=20))

auto.corr <- acf(fitz$yhat.train[ , i], plot=FALSE)
max.lag <- max(auto.corr$lag[ , 1, 1])

j <- seq(-0.5, 0.4, length.out=10)
png("toy02.png", width = 7, height = 7, units = 'in', res = 100)
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

png("toy03.png", width = 7, height = 7, units = 'in', res = 100)
for(j in 1:10) {
  if(j==1)
    plot(pnorm(fitz$yhat.train[ , i[j]]),
         type='l', ylim=c(0, 1),
         sub=paste0('N:', n, ', p:', p, ', thin:', thin),
         ylab=expression(Phi(f(x))), xlab='m')
  else
    lines(pnorm(fitz$yhat.train[ , i[j]]),
          type='l', col=j)
}
dev.off()

geweke <- gewekediag(fitz$yhat.train)

j <- -40
png("toy04.png", width = 7, height = 7, units = 'in', res = 100)
plot(geweke$z, pch='.', cex=2, ylab='z', xlab='i',
     sub=paste0('N:', n, ', p:', p, ', thin:', thin),
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

png("toyvanilla_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(vanilla$sigma[seq(5000,80000, by = 150)], type = "l",
     ylab = expression(sigma),
     xlab = "Draw",
     ylim= c(0.40, 0.80))
dev.off()
png("toyvanilla_2.png", width = 7, height = 7, units = 'in', res = 100)
acf(vanilla$sigma[seq(5000,80000, by = 150)], main = "")
dev.off()

png("toyoracle_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(oracle$sigma[seq(5000,80000, by = 150)], type = "l",
     ylab = expression(sigma),
     xlab = "Draw",
     ylim= c(0.40, 0.80))
dev.off()
png("toyoracle_2.png", width = 7, height = 7, units = 'in', res = 100)
acf(oracle$sigma[seq(5000,80000, by = 150)], main = "")
dev.off()

png("toyps_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ps$sigma[seq(5000,80000, by = 150)], type = "l",
     ylab = expression(sigma),
     xlab = "Draw",
     ylim= c(0.40, 0.80))
dev.off()
png("toyps_2.png", width = 7, height = 7, units = 'in', res = 100)
acf(ps$sigma[seq(5000,80000, by = 150)], main = "")
dev.off()

png("toyglm_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ps3$sigma[seq(5000,80000, by = 150)], type = "l",
     ylab = expression(sigma),
     xlab = "Draw",
     ylim= c(0.40, 0.80))
dev.off()
png("toyglm_2.png", width = 7, height = 7, units = 'in', res = 100)
acf(ps3$sigma[seq(5000,80000, by = 150)], main = "")
dev.off()

png("toyrand_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ps4$sigma[seq(5000,80000, by = 150)], type = "l",
     ylab = expression(sigma),
     xlab = "Draw",
     ylim= c(0.40, 0.80))
dev.off()
png("toyrand_2.png", width = 7, height = 7, units = 'in', res = 100)
acf(ps4$sigma[seq(5000,80000, by = 150)], main = "")
dev.off()

png("toybartbcf_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(fitbcf2$sigma, type = "l",
     ylab = expression(sigma),
     xlab = "Draw",
     ylim= c(0.40, 0.80))
dev.off()
png("toybartbcf_2.png", width = 7, height = 7, units = 'in', res = 100)
acf(fitbcf2$sigma, main = "")
dev.off()

png("toyoraclebcf_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(fitbcf3$sigma, type = "l",
     ylab = expression(sigma),
     xlab = "Draw",
     ylim= c(0.40, 0.80))
dev.off()
png("toyoraclebcf_2.png", width = 7, height = 7, units = 'in', res = 100)
acf(fitbcf3$sigma, main = "")
dev.off()

png("toyglmbcf_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(fitbcf4$sigma, type = "l",
     ylab = expression(sigma),
     xlab = "Draw",
     ylim= c(0.40, 0.80))
dev.off()
png("toyglmbcf_2.png", width = 7, height = 7, units = 'in', res = 100)
acf(fitbcf4$sigma, main = "")
dev.off()

png("toyrandbcf_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(fitbcf5$sigma, type = "l",
     ylab = expression(sigma),
     xlab = "Draw",
     ylim= c(0.40, 0.80))
dev.off()
png("toyrandbcf_2.png", width = 7, height = 7, units = 'in', res = 100)
acf(fitbcf5$sigma, main = "")
dev.off()



par(cex=0.5)

#################
#####Vanilla#####
#################

png("toy.vanilla.icex1_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x1,
     ylim = c(-1,3),
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


png("toy.vanilla.icex1_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x1.0025,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("toy.vanilla.icex1_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x1.0975,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("toy.vanilla.icex1_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x1,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[1]))
dev.off()

png("toy.vanilla.icex1_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.vanilla.x1),
     ylab = "Derivative Y",
     xlab = expression(x[1]))
dev.off()


png("toy.vanilla.icex2_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x2,
     ylim = c(-1,3),
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


png("toy.vanilla.icex2_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x2.0025,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("toy.vanilla.icex2_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x2.0975,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("toy.vanilla.icex2_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x2,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[2]))
dev.off()

png("toy.vanilla.icex2_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.vanilla.x2),
     ylab = "Derivative Y",
     xlab = expression(x[2]))
dev.off()


png("toy.vanilla.icex3_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x3,
     ylim = c(-1,3),
     ylab = "Partial Y",
     xlab = expression(x[3]))
lines(as.numeric(names(apply(ice.bart.vanilla.x3.0975$ice_curves,2,mean))),
      apply(ice.bart.vanilla.x3.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.vanilla.x3.0025$ice_curves,2,mean))),
      apply(ice.bart.vanilla.x3.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
{segments(x0 = -3, y0 = 0, x1 = -3/4, y1 = 0, col = "red", lwd = 3)
  segments(x0 = -3/4, y0 = 0, x1 = -3/4, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = -3/4, y0 = 0.5, x1 = 0, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 0, y0 = 0.5, x1 = 0, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 0, y0 = 0.75, x1 = 3/4, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 0.75, x1 = 3/4, y1 = 1, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 1, x1 = 15, y1 = 1, col = "red", lwd = 3)}

dev.off()


png("toy.vanilla.icex3_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x3.0025,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("toy.vanilla.icex3_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x3.0975,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("toy.vanilla.icex3_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.vanilla.x3,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[3]))
dev.off()

png("toy.vanilla.icex3_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.vanilla.x3),
     ylab = "Derivative Y",
     xlab = expression(x[3]))
dev.off()


###################
######Oracle#######
###################

png("toy.oracle.icex1_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x1,
     ylim = c(-1,3),
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


png("toy.oracle.icex1_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x1.0025,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("toy.oracle.icex1_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x1.0975,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("toy.oracle.icex1_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x1,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[1]))
dev.off()

png("toy.oracle.icex1_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.oracle.x1),
     ylab = "Derivative Y",
     xlab = expression(x[1]))
dev.off()


png("toy.oracle.icex2_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x2,
     ylim = c(-1,3),
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


png("toy.oracle.icex2_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x2.0025,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("toy.oracle.icex2_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x2.0975,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("toy.oracle.icex2_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x2,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[2]))
dev.off()

png("toy.oracle.icex2_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.oracle.x2),
     ylab = "Derivative Y",
     xlab = expression(x[2]))
dev.off()


png("toy.oracle.icex3_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x3,
     ylim = c(-1,3),
     ylab = "Partial Y",
     xlab = expression(x[3]))
lines(as.numeric(names(apply(ice.bart.oracle.x3.0975$ice_curves,2,mean))),
      apply(ice.bart.oracle.x3.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.oracle.x3.0025$ice_curves,2,mean))),
      apply(ice.bart.oracle.x3.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
{segments(x0 = -3, y0 = 0, x1 = -3/4, y1 = 0, col = "red", lwd = 3)
  segments(x0 = -3/4, y0 = 0, x1 = -3/4, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = -3/4, y0 = 0.5, x1 = 0, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 0, y0 = 0.5, x1 = 0, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 0, y0 = 0.75, x1 = 3/4, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 0.75, x1 = 3/4, y1 = 1, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 1, x1 = 15, y1 = 1, col = "red", lwd = 3)}

dev.off()


png("toy.oracle.icex3_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x3.0025,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("toy.oracle.icex3_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x3.0975,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("toy.oracle.icex3_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.x3,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[3]))
dev.off()

png("toy.oracle.icex3_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.oracle.x3),
     ylab = "Derivative Y",
     xlab = expression(x[3]))
dev.off()

png("toy.oracle.icepihat_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.pihat,
     ylim = c(-1,3),
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

png("toy.oracle.icepihat_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.pihat.0025,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("toy.oracle.icepihat_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.pihat.0975,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("toy.oracle.icepihat_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.oracle.pihat,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("toy.oracle.icepihat_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.oracle.pihat),
     ylab = "Derivative Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

###########
#####PS####
###########

png("toy.ps.icex1_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x1,
     ylim = c(-1,3),
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


png("toy.ps.icex1_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x1.0025,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("toy.ps.icex1_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x1.0975,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("toy.ps.icex1_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x1,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[1]))
dev.off()

png("toy.ps.icex1_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps.x1),
     ylab = "Derivative Y",
     xlab = expression(x[1]))
dev.off()


png("toy.ps.icex2_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x2,
     ylim = c(-1,3),
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


png("toy.ps.icex2_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x2.0025,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("toy.ps.icex2_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x2.0975,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("toy.ps.icex2_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x2,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[2]))
dev.off()

png("toy.ps.icex2_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps.x2),
     ylab = "Derivative Y",
     xlab = expression(x[2]))
dev.off()


png("toy.ps.icex3_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x3,
     ylim = c(-1,3),
     ylab = "Partial Y",
     xlab = expression(x[3]))
lines(as.numeric(names(apply(ice.bart.ps.x3.0975$ice_curves,2,mean))),
      apply(ice.bart.ps.x3.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.ps.x3.0025$ice_curves,2,mean))),
      apply(ice.bart.ps.x3.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
{segments(x0 = -3, y0 = 0, x1 = -3/4, y1 = 0, col = "red", lwd = 3)
  segments(x0 = -3/4, y0 = 0, x1 = -3/4, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = -3/4, y0 = 0.5, x1 = 0, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 0, y0 = 0.5, x1 = 0, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 0, y0 = 0.75, x1 = 3/4, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 0.75, x1 = 3/4, y1 = 1, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 1, x1 = 15, y1 = 1, col = "red", lwd = 3)}

dev.off()


png("toy.ps.icex3_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x3.0025,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("toy.ps.icex3_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x3.0975,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("toy.ps.icex3_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.x3,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[3]))
dev.off()

png("toy.ps.icex3_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps.x3),
     ylab = "Derivative Y",
     xlab = expression(x[3]))
dev.off()

png("toy.ps.icepihat_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.pihat,
     ylim = c(-1,3),
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


png("toy.ps.icepihat_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.pihat.0025,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("toy.ps.icepihat_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.pihat.0975,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("toy.ps.icepihat_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps.pihat,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("toy.ps.icepihat_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps.pihat),
     ylab = "Derivative Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

#########
###PS3###
#########

png("toy.ps3.icex1_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x1,
     ylim = c(-1,3),
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


png("toy.ps3.icex1_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x1.0025,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("toy.ps3.icex1_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x1.0975,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("toy.ps3.icex1_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x1,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[1]))
dev.off()

png("toy.ps3.icex1_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps3.x1),
     ylab = "Derivative Y",
     xlab = expression(x[1]))
dev.off()


png("toy.ps3.icex2_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x2,
     ylim = c(-1,3),
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


png("toy.ps3.icex2_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x2.0025,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("toy.ps3.icex2_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x2.0975,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("toy.ps3.icex2_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x2,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[2]))
dev.off()

png("toy.ps3.icex2_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps3.x2),
     ylab = "Derivative Y",
     xlab = expression(x[2]))
dev.off()


png("toy.ps3.icex3_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x3,
     ylim = c(-1,3),
     ylab = "Partial Y",
     xlab = expression(x[3]))
lines(as.numeric(names(apply(ice.bart.ps3.x3.0975$ice_curves,2,mean))),
      apply(ice.bart.ps3.x3.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.ps3.x3.0025$ice_curves,2,mean))),
      apply(ice.bart.ps3.x3.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
{segments(x0 = -3, y0 = 0, x1 = -3/4, y1 = 0, col = "red", lwd = 3)
  segments(x0 = -3/4, y0 = 0, x1 = -3/4, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = -3/4, y0 = 0.5, x1 = 0, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 0, y0 = 0.5, x1 = 0, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 0, y0 = 0.75, x1 = 3/4, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 0.75, x1 = 3/4, y1 = 1, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 1, x1 = 15, y1 = 1, col = "red", lwd = 3)}

dev.off()


png("toy.ps3.icex3_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x3.0025,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("toy.ps3.icex3_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x3.0975,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("toy.ps3.icex3_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.x3,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[3]))
dev.off()

png("toy.ps3.icex3_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps3.x3),
     ylab = "Derivative Y",
     xlab = expression(x[3]))
dev.off()

png("toy.ps3.icepihat_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.pihat,
     ylim = c(-1,3),
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


png("toy.ps3.icepihat_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.pihat.0025,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("toy.ps3.icepihat_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.pihat.0975,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("toy.ps3.icepihat_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps3.pihat,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("toy.ps3.icepihat_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps3.pihat),
     ylab = "Derivative Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

#######
##PS4##
#######

png("toy.ps4.icex1_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x1,
     ylim = c(-1,3),
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


png("toy.ps4.icex1_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x1.0025,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("toy.ps4.icex1_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x1.0975,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("toy.ps4.icex1_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x1,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[1]))
dev.off()

png("toy.ps4.icex1_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps4.x1),
     ylab = "Derivative Y",
     xlab = expression(x[1]))
dev.off()


png("toy.ps4.icex2_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x2,
     ylim = c(-1,3),
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


png("toy.ps4.icex2_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x2.0025,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("toy.ps4.icex2_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x2.0975,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("toy.ps4.icex2_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x2,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[2]))
dev.off()

png("toy.ps4.icex2_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps4.x2),
     ylab = "Derivative Y",
     xlab = expression(x[2]))
dev.off()


png("toy.ps4.icex3_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x3,
     ylim = c(-1,3),
     ylab = "Partial Y",
     xlab = expression(x[3]))
lines(as.numeric(names(apply(ice.bart.ps4.x3.0975$ice_curves,2,mean))),
      apply(ice.bart.ps4.x3.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bart.ps4.x3.0025$ice_curves,2,mean))),
      apply(ice.bart.ps4.x3.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
{segments(x0 = -3, y0 = 0, x1 = -3/4, y1 = 0, col = "red", lwd = 3)
  segments(x0 = -3/4, y0 = 0, x1 = -3/4, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = -3/4, y0 = 0.5, x1 = 0, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 0, y0 = 0.5, x1 = 0, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 0, y0 = 0.75, x1 = 3/4, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 0.75, x1 = 3/4, y1 = 1, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 1, x1 = 15, y1 = 1, col = "red", lwd = 3)}

dev.off()


png("toy.ps4.icex3_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x3.0025,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("toy.ps4.icex3_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x3.0975,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("toy.ps4.icex3_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.x3,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[3]))
dev.off()

png("toy.ps4.icex3_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps4.x3),
     ylab = "Derivative Y",
     xlab = expression(x[3]))
dev.off()

png("toy.ps4.icepihat_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.pihat,
     ylim = c(-1,3),
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


png("toy.ps4.icepihat_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.pihat.0025,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("toy.ps4.icepihat_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.pihat.0975,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("toy.ps4.icepihat_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bart.ps4.pihat,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("toy.ps4.icepihat_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bart.ps4.pihat),
     ylab = "Derivative Y",
     xlab = expression(hat(Pi)(x)))
dev.off()


#####################################
############## BCF ##################
#####################################

###############
####Oracle#####
###############


png("toy.bcf.oracle.icex1_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x1,
     ylim = c(-1,3),
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


png("toy.bcf.oracle.icex1_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x1.0025,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("toy.bcf.oracle.icex1_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x1.0975,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("toy.bcf.oracle.icex1_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x1,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[1]))
dev.off()

png("toy.bcf.oracle.icex1_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.oracle.x1),
     ylab = "Derivative Y",
     xlab = expression(x[1]))
dev.off()


png("toy.bcf.oracle.icex2_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x2,
     ylim = c(-1,3),
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


png("toy.bcf.oracle.icex2_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x2.0025,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("toy.bcf.oracle.icex2_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x2.0975,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("toy.bcf.oracle.icex2_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x2,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[2]))
dev.off()

png("toy.bcf.oracle.icex2_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.oracle.x2),
     ylab = "Derivative Y",
     xlab = expression(x[2]))
dev.off()


png("toy.bcf.oracle.icex3_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x3,
     ylim = c(-1,3),
     ylab = "Partial Y",
     xlab = expression(x[3]))
lines(as.numeric(names(apply(ice.bcf.oracle.x3.0975$ice_curves,2,mean))),
      apply(ice.bcf.oracle.x3.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.oracle.x3.0025$ice_curves,2,mean))),
      apply(ice.bcf.oracle.x3.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
{segments(x0 = -3, y0 = 0, x1 = -3/4, y1 = 0, col = "red", lwd = 3)
  segments(x0 = -3/4, y0 = 0, x1 = -3/4, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = -3/4, y0 = 0.5, x1 = 0, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 0, y0 = 0.5, x1 = 0, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 0, y0 = 0.75, x1 = 3/4, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 0.75, x1 = 3/4, y1 = 1, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 1, x1 = 15, y1 = 1, col = "red", lwd = 3)}

dev.off()


png("toy.bcf.oracle.icex3_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x3.0025,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("toy.bcf.oracle.icex3_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x3.0975,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("toy.bcf.oracle.icex3_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.x3,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[3]))
dev.off()

png("toy.bcf.oracle.icex3_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.oracle.x3),
     ylab = "Derivative Y",
     xlab = expression(x[3]))
dev.off()

png("toy.bcf.oracle.icepihat_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.pihat,
     ylim = c(-1,3),
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


png("toy.bcf.oracle.icepihat_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.pihat.0025,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("toy.bcf.oracle.icepihat_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.pihat.0975,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("toy.bcf.oracle.icepihat_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.oracle.pihat,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("toy.bcf.oracle.icepihat_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.oracle.pihat),
     ylab = "Derivative Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

#############
#####Bart####
#############

png("toy.bcf.bart.icex1_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x1,
     ylim = c(-1,3),
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


png("toy.bcf.bart.icex1_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x1.0025,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("toy.bcf.bart.icex1_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x1.0975,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("toy.bcf.bart.icex1_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x1,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[1]))
dev.off()

png("toy.bcf.bart.icex1_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.bart.x1),
     ylab = "Derivative Y",
     xlab = expression(x[1]))
dev.off()


png("toy.bcf.bart.icex2_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x2,
     ylim = c(-1,3),
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


png("toy.bcf.bart.icex2_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x2.0025,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("toy.bcf.bart.icex2_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x2.0975,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("toy.bcf.bart.icex2_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x2,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[2]))
dev.off()

png("toy.bcf.bart.icex2_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.bart.x2),
     ylab = "Derivative Y",
     xlab = expression(x[2]))
dev.off()


png("toy.bcf.bart.icex3_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x3,
     ylim = c(-1,3),
     ylab = "Partial Y",
     xlab = expression(x[3]))
lines(as.numeric(names(apply(ice.bcf.bart.x3.0975$ice_curves,2,mean))),
      apply(ice.bcf.bart.x3.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.bart.x3.0025$ice_curves,2,mean))),
      apply(ice.bcf.bart.x3.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
{segments(x0 = -3, y0 = 0, x1 = -3/4, y1 = 0, col = "red", lwd = 3)
  segments(x0 = -3/4, y0 = 0, x1 = -3/4, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = -3/4, y0 = 0.5, x1 = 0, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 0, y0 = 0.5, x1 = 0, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 0, y0 = 0.75, x1 = 3/4, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 0.75, x1 = 3/4, y1 = 1, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 1, x1 = 15, y1 = 1, col = "red", lwd = 3)}

dev.off()


png("toy.bcf.bart.icex3_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x3.0025,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("toy.bcf.bart.icex3_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x3.0975,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("toy.bcf.bart.icex3_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.x3,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[3]))
dev.off()

png("toy.bcf.bart.icex3_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.bart.x3),
     ylab = "Derivative Y",
     xlab = expression(x[3]))
dev.off()

png("toy.bcf.bart.icepihat_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.pihat,
     ylim = c(-1,3),
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


png("toy.bcf.bart.icepihat_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.pihat.0025,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("toy.bcf.bart.icepihat_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.pihat.0975,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("toy.bcf.bart.icepihat_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.pihat,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("toy.bcf.bart.icepihat_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.bart.pihat),
     ylab = "Derivative Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

#########
###glm###
#########

png("toy.bcf.glm.icex1_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x1,
     ylim = c(-1,3),
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


png("toy.bcf.glm.icex1_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x1.0025,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("toy.bcf.glm.icex1_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x1.0975,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("toy.bcf.glm.icex1_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x1,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[1]))
dev.off()

png("toy.bcf.glm.icex1_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.glm.x1),
     ylab = "Derivative Y",
     xlab = expression(x[1]))
dev.off()


png("toy.bcf.glm.icex2_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x2,
     ylim = c(-1,3),
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


png("toy.bcf.glm.icex2_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x2.0025,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("toy.bcf.glm.icex2_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x2.0975,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("toy.bcf.glm.icex2_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x2,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[2]))
dev.off()

png("toy.bcf.glm.icex2_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.glm.x2),
     ylab = "Derivative Y",
     xlab = expression(x[2]))
dev.off()


png("toy.bcf.glm.icex3_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x3,
     ylim = c(-1,3),
     ylab = "Partial Y",
     xlab = expression(x[3]))
lines(as.numeric(names(apply(ice.bcf.glm.x3.0975$ice_curves,2,mean))),
      apply(ice.bcf.glm.x3.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.glm.x3.0025$ice_curves,2,mean))),
      apply(ice.bcf.glm.x3.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
{segments(x0 = -3, y0 = 0, x1 = -3/4, y1 = 0, col = "red", lwd = 3)
  segments(x0 = -3/4, y0 = 0, x1 = -3/4, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = -3/4, y0 = 0.5, x1 = 0, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 0, y0 = 0.5, x1 = 0, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 0, y0 = 0.75, x1 = 3/4, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 0.75, x1 = 3/4, y1 = 1, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 1, x1 = 15, y1 = 1, col = "red", lwd = 3)}

dev.off()


png("toy.bcf.glm.icex3_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x3.0025,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("toy.bcf.glm.icex3_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x3.0975,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("toy.bcf.glm.icex3_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.x3,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[3]))
dev.off()

png("toy.bcf.glm.icex3_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.glm.x3),
     ylab = "Derivative Y",
     xlab = expression(x[3]))
dev.off()

png("toy.bcf.glm.icepihat_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.pihat,
     ylim = c(-1,3),
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


png("toy.bcf.glm.icepihat_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.pihat.0025,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("toy.bcf.glm.icepihat_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.pihat.0975,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("toy.bcf.glm.icepihat_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.glm.pihat,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("toy.bcf.glm.icepihat_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.glm.pihat),
     ylab = "Derivative Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

#######
##rand##
#######

png("toy.bcf.rand.icex1_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x1,
     ylim = c(-1,3),
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


png("toy.bcf.rand.icex1_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x1.0025,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("toy.bcf.rand.icex1_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x1.0975,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("toy.bcf.rand.icex1_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x1,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[1]))
dev.off()

png("toy.bcf.rand.icex1_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.rand.x1),
     ylab = "Derivative Y",
     xlab = expression(x[1]))
dev.off()


png("toy.bcf.rand.icex2_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x2,
     ylim = c(-1,3),
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


png("toy.bcf.rand.icex2_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x2.0025,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("toy.bcf.rand.icex2_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x2.0975,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("toy.bcf.rand.icex2_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x2,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[2]))
dev.off()

png("toy.bcf.rand.icex2_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.rand.x2),
     ylab = "Derivative Y",
     xlab = expression(x[2]))
dev.off()


png("toy.bcf.rand.icex3_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x3,
     ylim = c(-1,3),
     ylab = "Partial Y",
     xlab = expression(x[3]))
lines(as.numeric(names(apply(ice.bcf.rand.x3.0975$ice_curves,2,mean))),
      apply(ice.bcf.rand.x3.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.rand.x3.0025$ice_curves,2,mean))),
      apply(ice.bcf.rand.x3.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
{segments(x0 = -3, y0 = 0, x1 = -3/4, y1 = 0, col = "red", lwd = 3)
  segments(x0 = -3/4, y0 = 0, x1 = -3/4, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = -3/4, y0 = 0.5, x1 = 0, y1 = 0.5, col = "red", lwd = 3)
  segments(x0 = 0, y0 = 0.5, x1 = 0, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 0, y0 = 0.75, x1 = 3/4, y1 = 0.75, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 0.75, x1 = 3/4, y1 = 1, col = "red", lwd = 3)
  segments(x0 = 3/4, y0 = 1, x1 = 15, y1 = 1, col = "red", lwd = 3)}
dev.off()


png("toy.bcf.rand.icex3_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x3.0025,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("toy.bcf.rand.icex3_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x3.0975,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("toy.bcf.rand.icex3_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.x3,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[3]))
dev.off()

png("toy.bcf.rand.icex3_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.rand.x3),
     ylab = "Derivative Y",
     xlab = expression(x[3]))
dev.off()

png("toy.bcf.rand.icepihat_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.pihat,
     ylim = c(-1,3),
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


png("toy.bcf.rand.icepihat_0025.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.pihat.0025,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("toy.bcf.rand.icepihat_0975.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.pihat.0975,
     ylab = "Partial Y",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("toy.bcf.rand.icepihat_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.rand.pihat,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(hat(Pi)(x)))
dev.off()

png("toy.bcf.rand.icepihat_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.rand.pihat),
     ylab = "Derivative Y",
     xlab = expression(hat(Pi)(x)))
dev.off()
