######################################
##Pedro Henrique Filipini dos Santos##
###University of São Paulo - Brazil###
######################################


######################################
##Real Data Analysis - Smoke Dataset##
######################################

#Libraries
library(BART)
library(bcf)
library(dplyr)
library(ICEbox)
library(caret)

seedaux <- 99
set.seed(seedaux)

#Setting working directory
#setwd("~/Datasets/pscore")


#Parameters to BART functions
burn <- 15000
npost<- 2000
pthin <- 5000
thin <- 100

#Loading data
data<-read.csv("nmesdata.txt", header=T)

#Using their variables
data<-subset(data, packyears>0)
data<-subset(data, TOTALEXP>0)
data<-subset(data, LASTAGE>27)
data<-subset(data, select=c(packyears, LASTAGE, AGESMOKE,
                            MALE, RACE3, marital, educate,
                            SREGION, POVSTALB, beltuse,
                            yearsince, TOTALEXP))

## complete case analysis
data<-na.omit(data)

#Total of elements
n <- nrow(data)

## factor variables
data$RACE3<-factor(data$RACE3, ordered=F)
data$marital<-factor(data$marital, ordered=F)
data$SREGION<-factor(data$SREGION, ordered=F)
data$educate<-factor(data$SREGION, ordered=F)
data$beltuse<-factor(data$beltuse, ordered=F)
data$POVSTALB<-factor(data$POVSTALB, ordered=F)

summary(data)

#Setting the matrix for the propensity score
X <- subset(data, select = -c(TOTALEXP, packyears))
Z <- as.numeric(data$packyears>17)
length(Z)

##Rearranging data
race <- predict(dummyVars(~ RACE3, data = X), newdata = X)
marital <- predict(dummyVars(~ marital, data = X), newdata = X)
educate <- predict(dummyVars(~ educate, data = X), newdata = X)
SREGION <- predict(dummyVars(~ SREGION, data = X), newdata = X)
POVSTALB <- predict(dummyVars(~ POVSTALB, data = X), newdata = X)
beltuse <- predict(dummyVars(~ beltuse, data = X), newdata = X)

X3 <- subset(X, select=c(LASTAGE, AGESMOKE, MALE, yearsince))
X3 <- cbind(X3,race,marital,educate,SREGION,POVSTALB,beltuse)

X3 <- as.matrix(X3)

k <- ncol(X3)

#Estimating the Propensity Score
#set.seed(seedaux)
ps = mc.pbart(X3, Z, nskip = burn,ndpost = npost, keepevery = pthin,
              mc.cores = 3, seed = seedaux)

pihat = ps$prob.train.mean

pdf("smoke_pbart0")
plot(sort(pihat), 
     ylab = expression(hat(pi)(x[i])),
     xlab = "Individuals")
dev.off()

#################################
##Assuring Convergence of pbart##
#################################

i <- floor(seq(1, n, length.out=20))

auto.corr <- acf(ps$yhat.train[ , i], plot=FALSE)
max.lag <- max(auto.corr$lag[ , 1, 1])

j <- seq(-0.5, 0.4, length.out=10)
pdf("smoke_pbart1.pdf")
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

pdf("smoke_pbart2.pdf")
for(j in 1:10) {
  if(j==1)
    plot(pnorm(ps$yhat.train[ , i[j]]),
         type='l', ylim=c(0, 1),
         sub=paste0('N:', n, ', k:', k, ', thin:', 5000),
         ylab=expression(Phi(f(x))), xlab='m')
  else
    lines(pnorm(ps$yhat.train[ , i[j]]),
          type='l', col=j)
}
dev.off()

geweke <- gewekediag(ps$yhat.train)

j <- -10^(log10(n)-1)
pdf("smoke_pbart3.pdf")
plot(geweke$z, pch='.', cex=2, ylab='z', xlab='i',
     sub=paste0('N:', n, ', k:', k, ', thin:', 5000),
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
#Setting the data and test matrix
Y = log(data$TOTALEXP)
X1 = cbind(X3,pihat,Z)
names(X1)

#Test Matrix
Xtest = cbind(rbind(cbind(X3,pihat),cbind(X3,pihat)),c(rep(1,n),rep(0,n)))
colnames(Xtest)[ncol(Xtest)]="Z"
colnames(Xtest)

###################################################################

#Running the model
set.seed(seedaux)
psbart = wbart(X1, Y, Xtest, nskip = burn,
               ndpost = npost, keepevery = thin
                  )


#Calculating estimated ATE
ate_est_psbart = rowMeans(psbart$yhat.test[,1:n] - psbart$yhat.test[,(n+1):(2*n)])

#Calculating estimated ITE
ite_est_psbart = psbart$yhat.test[,1:n] - psbart$yhat.test[,(n+1):(2*n)]



#Running the model with Dirichlet
set.seed(seedaux)
psbart2 = wbart(X1, Y, Xtest, sparse = T, 
                nskip = burn,ndpost = npost, keepevery = thin
                  )

pdf("smoke_dart_1.pdf")
plot(psbart2$varprob.mean, ylim = c(0,1),
     ylab = "DART - Dirichlet Hyperprior Posterior Draws")
for(i in 1:length(psbart2$varprob.mean)){
  segments(x0 = i, y0 = as.numeric(psbart2$varprob.mean)[i],
           x1 = i, y1 = as.numeric(apply(psbart2$varprob,2,quantile, probs = 0.025))[i],
           col = "red", lwd = 2)
}
for(i in 1:length(psbart2$varprob.mean)){
  segments(x0 = i, y0 = as.numeric(psbart2$varprob.mean)[i],
           x1 = i, y1 = as.numeric(apply(psbart2$varprob,2,quantile, probs = 0.975))[i],
           col = "red", lwd = 2)
}
dev.off()

pdf("smoke_dart_2.pdf")
plot(colMeans(psbart2$varcount > 0), ylim=c(0,1),
     ylab = "DART - Posterior Inclusion Probability")
abline(h=0.5, col = "red", lwd = 2)
dev.off()


#Calculating estimated ATE
ate_est_psbart2 = rowMeans(psbart2$yhat.test[,1:n] - psbart2$yhat.test[,(n+1):(2*n)])

#Calculating estimated ITE
ite_est_psbart2 = psbart2$yhat.test[,1:n] - psbart2$yhat.test[,(n+1):(2*n)]


#Looking at the individual treatment effects
#PSBART
require(dplyr)
df_psbart <- data.frame(
  cm= apply(ite_est_psbart,2,median),
  cl= apply(ite_est_psbart,2,quantile,prob=0.025),
  cu= apply(ite_est_psbart,2,quantile,prob=0.975)
)
df_psbart <- df_psbart %>% arrange(cm)
df_psbart <- data.frame(x=1:n,df_psbart)

require(ggplot2)
pdf("smoke_dart_3.pdf")
ggplot(df_psbart, aes(x=x, y=cm)) + 
  geom_point(size = 0.05) +
  geom_point(aes(x=x, y=cu),size = 0.05) +
  geom_point(aes(x=x, y=cl),size = 0.05) +
  geom_line(aes(x=x,y=0), size = 2, linetype ="dotted", colour = "red") +
  ylim(-0.5, 0.6) +
  ylab("Treatment Effect")+
  xlab("Index")+
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.background = element_blank())
dev.off()



#Comparing with the bcf
#BCF with BART (pihat2)
set.seed(seedaux)
fitbcf = bcf(Y, Z, X3, X3, pihat, burn, npost, nthin = thin, include_pi="both")

#Estimated Treatment Effects
ate_est_bartbcf = rowMeans(fitbcf$tau)
ite_est_bartbcf = colMeans(fitbcf$tau)


#BCF
df_bcf <- data.frame(
  cm= apply(fitbcf$tau,2,median),
  cl= apply(fitbcf$tau,2,quantile,prob=0.025),
  cu= apply(fitbcf$tau,2,quantile,prob=0.975)
)
df_bcf <- df_bcf %>% arrange(cm)
df_bcf <- data.frame(x=1:n,df_bcf)

require(ggplot2)
pdf("smoke_bcf_1.pdf")
ggplot(df_bcf, aes(x=x, y=cm)) + 
  geom_point(size = 0.05) +
  geom_point(aes(x=x, y=cu),size = 0.05) +
  geom_point(aes(x=x, y=cl),size = 0.05) +
  geom_line(aes(x=x,y=0), size = 2, linetype ="dotted", colour = "red") +
  ylim(-0.5, 0.6) +
  ylab("Treatment Effect")+
  xlab("Index")+
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.background = element_blank())
dev.off()

pdf("smokecate.pdf")
boxplot(ate_est_psbart2,ate_est_bartbcf,
        names = c("DART", "BCF"),
        ylim = c(-0.05,0.35),
        ylab = "CATE")
dev.off()

pdf("smokedartconv_1.pdf")
plot(psbart2$sigma[seq(15000,215000, by = 100)], type = "l",
     ylab = expression(sigma),
     xlab = "Draw")
dev.off()
pdf("smokedartconv_2.pdf")
acf(psbart2$sigma[seq(15000,215000, by = 100)], main = "")
dev.off()

pdf("smokebcfconv_1.pdf")
plot(fitbcf$sigma, type = "l",
     ylab = expression(sigma),
     xlab = "Draw")
dev.off()
pdf("smokebcfconv_2.pdf")
acf(fitbcf$sigma, main = "")
dev.off()


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
fitbcfaux <- psbart
fitbcfaux$treedraws <- fitbcf$treedraws
fitbcfaux$mu <- mean(Y)
fitbcfaux$vars <- (predict(fitbcfaux, newdata = cbind(X3,pihat))[,1]-mean(Y))/(fitbcf$tau[,1])

#LASTAGE#

Xice = cbind(X3,pihat)
ice.bcf.bart.LASTAGE <- ice(fitbcfaux, Xice, predictor = "LASTAGE", 
                       predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.LASTAGE, file="ice.bcf.bart.LASTAGE.RData")
plot(ice.bcf.bart.LASTAGE)

ice.bcf.bart.LASTAGE.0025 <- ice(fitbcfaux, Xice, predictor = "LASTAGE", 
                       predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.LASTAGE.0025, file="ice.bcf.bart.LASTAGE.0025.RData")
plot(ice.bcf.bart.LASTAGE.0025)

ice.bcf.bart.LASTAGE.0975 <- ice(fitbcfaux, Xice, predictor = "LASTAGE", 
                            predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.LASTAGE.0975, file="ice.bcf.bart.LASTAGE.0975.RData")
plot(ice.bcf.bart.LASTAGE.0975)


#AGESMOKE#

Xice = cbind(X3,pihat)
ice.bcf.bart.AGESMOKE <- ice(fitbcfaux, Xice, predictor = "AGESMOKE", 
                            predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.AGESMOKE, file="ice.bcf.bart.AGESMOKE.RData")
plot(ice.bcf.bart.AGESMOKE)

ice.bcf.bart.AGESMOKE.0025 <- ice(fitbcfaux, Xice, predictor = "AGESMOKE", 
                                 predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.AGESMOKE.0025, file="ice.bcf.bart.AGESMOKE.0025.RData")
plot(ice.bcf.bart.AGESMOKE.0025)

ice.bcf.bart.AGESMOKE.0975 <- ice(fitbcfaux, Xice, predictor = "AGESMOKE", 
                                 predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.AGESMOKE.0975, file="ice.bcf.bart.AGESMOKE.0975.RData")
plot(ice.bcf.bart.AGESMOKE.0975)

#MALE#

Xice = cbind(X3,pihat)
ice.bcf.bart.MALE <- ice(fitbcfaux, Xice, predictor = "MALE", 
                             predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.MALE, file="ice.bcf.bart.MALE.RData")
plot(ice.bcf.bart.MALE)

ice.bcf.bart.MALE.0025 <- ice(fitbcfaux, Xice, predictor = "MALE", 
                                  predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.MALE.0025, file="ice.bcf.bart.MALE.0025.RData")
plot(ice.bcf.bart.MALE.0025)

ice.bcf.bart.MALE.0975 <- ice(fitbcfaux, Xice, predictor = "MALE", 
                                  predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.MALE.0975, file="ice.bcf.bart.MALE.0975.RData")
plot(ice.bcf.bart.MALE.0975)

#yearsince#

Xice = cbind(X3,pihat)
ice.bcf.bart.yearsince <- ice(fitbcfaux, Xice, predictor = "yearsince", 
                             predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.yearsince, file="ice.bcf.bart.yearsince.RData")
plot(ice.bcf.bart.yearsince)

ice.bcf.bart.yearsince.0025 <- ice(fitbcfaux, Xice, predictor = "yearsince", 
                                  predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.yearsince.0025, file="ice.bcf.bart.yearsince.0025.RData")
plot(ice.bcf.bart.yearsince.0025)

ice.bcf.bart.yearsince.0975 <- ice(fitbcfaux, Xice, predictor = "yearsince", 
                                  predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.yearsince.0975, file="ice.bcf.bart.yearsince.0975.RData")
plot(ice.bcf.bart.yearsince.0975)

#RACE3.1#

Xice = cbind(X3,pihat)
ice.bcf.bart.RACE3.1 <- ice(fitbcfaux, Xice, predictor = "RACE3.1", 
                              predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.RACE3.1, file="ice.bcf.bart.RACE3.1.RData")
plot(ice.bcf.bart.RACE3.1)

ice.bcf.bart.RACE3.1.0025 <- ice(fitbcfaux, Xice, predictor = "RACE3.1", 
                                   predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.RACE3.1.0025, file="ice.bcf.bart.RACE3.1.0025.RData")
plot(ice.bcf.bart.RACE3.1.0025)

ice.bcf.bart.RACE3.1.0975 <- ice(fitbcfaux, Xice, predictor = "RACE3.1", 
                                   predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.RACE3.1.0975, file="ice.bcf.bart.RACE3.1.0975.RData")
plot(ice.bcf.bart.RACE3.1.0975)

#RACE3.2#

Xice = cbind(X3,pihat)
ice.bcf.bart.RACE3.2 <- ice(fitbcfaux, Xice, predictor = "RACE3.2", 
                            predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.RACE3.2, file="ice.bcf.bart.RACE3.2.RData")
plot(ice.bcf.bart.RACE3.2)

ice.bcf.bart.RACE3.2.0025 <- ice(fitbcfaux, Xice, predictor = "RACE3.2", 
                                 predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.RACE3.2.0025, file="ice.bcf.bart.RACE3.2.0025.RData")
plot(ice.bcf.bart.RACE3.2.0025)

ice.bcf.bart.RACE3.2.0975 <- ice(fitbcfaux, Xice, predictor = "RACE3.2", 
                                 predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.RACE3.2.0975, file="ice.bcf.bart.RACE3.2.0975.RData")
plot(ice.bcf.bart.RACE3.2.0975)

#RACE3.3#

Xice = cbind(X3,pihat)
ice.bcf.bart.RACE3.3 <- ice(fitbcfaux, Xice, predictor = "RACE3.3", 
                            predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.RACE3.3, file="ice.bcf.bart.RACE3.3.RData")
plot(ice.bcf.bart.RACE3.3)

ice.bcf.bart.RACE3.3.0025 <- ice(fitbcfaux, Xice, predictor = "RACE3.3", 
                                 predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.RACE3.3.0025, file="ice.bcf.bart.RACE3.3.0025.RData")
plot(ice.bcf.bart.RACE3.3.0025)

ice.bcf.bart.RACE3.3.0975 <- ice(fitbcfaux, Xice, predictor = "RACE3.3", 
                                 predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.RACE3.3.0975, file="ice.bcf.bart.RACE3.3.0975.RData")
plot(ice.bcf.bart.RACE3.3.0975)


#marital.1#

Xice = cbind(X3,pihat)
ice.bcf.bart.marital.1 <- ice(fitbcfaux, Xice, predictor = "marital.1", 
                            predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.marital.1, file="ice.bcf.bart.marital.1.RData")
plot(ice.bcf.bart.marital.1)

ice.bcf.bart.marital.1.0025 <- ice(fitbcfaux, Xice, predictor = "marital.1", 
                                 predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.marital.1.0025, file="ice.bcf.bart.marital.1.0025.RData")
plot(ice.bcf.bart.marital.1.0025)

ice.bcf.bart.marital.1.0975 <- ice(fitbcfaux, Xice, predictor = "marital.1", 
                                 predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.marital.1.0975, file="ice.bcf.bart.marital.1.0975.RData")
plot(ice.bcf.bart.marital.1.0975)

#marital.2#

Xice = cbind(X3,pihat)
ice.bcf.bart.marital.2 <- ice(fitbcfaux, Xice, predictor = "marital.2", 
                            predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.marital.2, file="ice.bcf.bart.marital.2.RData")
plot(ice.bcf.bart.marital.2)

ice.bcf.bart.marital.2.0025 <- ice(fitbcfaux, Xice, predictor = "marital.2", 
                                 predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.marital.2.0025, file="ice.bcf.bart.marital.2.0025.RData")
plot(ice.bcf.bart.marital.2.0025)

ice.bcf.bart.marital.2.0975 <- ice(fitbcfaux, Xice, predictor = "marital.2", 
                                 predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.marital.2.0975, file="ice.bcf.bart.marital.2.0975.RData")
plot(ice.bcf.bart.marital.2.0975)

#marital.3#

Xice = cbind(X3,pihat)
ice.bcf.bart.marital.3 <- ice(fitbcfaux, Xice, predictor = "marital.3", 
                            predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.marital.3, file="ice.bcf.bart.marital.3.RData")
plot(ice.bcf.bart.marital.3)

ice.bcf.bart.marital.3.0025 <- ice(fitbcfaux, Xice, predictor = "marital.3", 
                                 predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.marital.3.0025, file="ice.bcf.bart.marital.3.0025.RData")
plot(ice.bcf.bart.marital.3.0025)

ice.bcf.bart.marital.3.0975 <- ice(fitbcfaux, Xice, predictor = "marital.3", 
                                 predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.marital.3.0975, file="ice.bcf.bart.marital.3.0975.RData")
plot(ice.bcf.bart.marital.3.0975)

#marital.4#

Xice = cbind(X3,pihat)
ice.bcf.bart.marital.4 <- ice(fitbcfaux, Xice, predictor = "marital.4", 
                              predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.marital.4, file="ice.bcf.bart.marital.4.RData")
plot(ice.bcf.bart.marital.4)

ice.bcf.bart.marital.4.0025 <- ice(fitbcfaux, Xice, predictor = "marital.4", 
                                   predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.marital.4.0025, file="ice.bcf.bart.marital.4.0025.RData")
plot(ice.bcf.bart.marital.4.0025)

ice.bcf.bart.marital.4.0975 <- ice(fitbcfaux, Xice, predictor = "marital.4", 
                                   predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.marital.4.0975, file="ice.bcf.bart.marital.4.0975.RData")
plot(ice.bcf.bart.marital.4.0975)

#marital.5#

Xice = cbind(X3,pihat)
ice.bcf.bart.marital.5 <- ice(fitbcfaux, Xice, predictor = "marital.5", 
                              predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.marital.5, file="ice.bcf.bart.marital.5.RData")
plot(ice.bcf.bart.marital.5)

ice.bcf.bart.marital.5.0025 <- ice(fitbcfaux, Xice, predictor = "marital.5", 
                                   predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.marital.5.0025, file="ice.bcf.bart.marital.5.0025.RData")
plot(ice.bcf.bart.marital.5.0025)

ice.bcf.bart.marital.5.0975 <- ice(fitbcfaux, Xice, predictor = "marital.5", 
                                   predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.marital.5.0975, file="ice.bcf.bart.marital.5.0975.RData")
plot(ice.bcf.bart.marital.5.0975)



#educate.1#

Xice = cbind(X3,pihat)
ice.bcf.bart.educate.1 <- ice(fitbcfaux, Xice, predictor = "educate.1", 
                              predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.educate.1, file="ice.bcf.bart.educate.1.RData")
plot(ice.bcf.bart.educate.1)

ice.bcf.bart.educate.1.0025 <- ice(fitbcfaux, Xice, predictor = "educate.1", 
                                   predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.educate.1.0025, file="ice.bcf.bart.educate.1.0025.RData")
plot(ice.bcf.bart.educate.1.0025)

ice.bcf.bart.educate.1.0975 <- ice(fitbcfaux, Xice, predictor = "educate.1", 
                                   predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.educate.1.0975, file="ice.bcf.bart.educate.1.0975.RData")
plot(ice.bcf.bart.educate.1.0975)

#educate.2#

Xice = cbind(X3,pihat)
ice.bcf.bart.educate.2 <- ice(fitbcfaux, Xice, predictor = "educate.2", 
                              predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.educate.2, file="ice.bcf.bart.educate.2.RData")
plot(ice.bcf.bart.educate.2)

ice.bcf.bart.educate.2.0025 <- ice(fitbcfaux, Xice, predictor = "educate.2", 
                                   predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.educate.2.0025, file="ice.bcf.bart.educate.2.0025.RData")
plot(ice.bcf.bart.educate.2.0025)

ice.bcf.bart.educate.2.0975 <- ice(fitbcfaux, Xice, predictor = "educate.2", 
                                   predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.educate.2.0975, file="ice.bcf.bart.educate.2.0975.RData")
plot(ice.bcf.bart.educate.2.0975)

#educate.3#

Xice = cbind(X3,pihat)
ice.bcf.bart.educate.3 <- ice(fitbcfaux, Xice, predictor = "educate.3", 
                              predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.educate.3, file="ice.bcf.bart.educate.3.RData")
plot(ice.bcf.bart.educate.3)

ice.bcf.bart.educate.3.0025 <- ice(fitbcfaux, Xice, predictor = "educate.3", 
                                   predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.educate.3.0025, file="ice.bcf.bart.educate.3.0025.RData")
plot(ice.bcf.bart.educate.3.0025)

ice.bcf.bart.educate.3.0975 <- ice(fitbcfaux, Xice, predictor = "educate.3", 
                                   predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.educate.3.0975, file="ice.bcf.bart.educate.3.0975.RData")
plot(ice.bcf.bart.educate.3.0975)

#educate.4#

Xice = cbind(X3,pihat)
ice.bcf.bart.educate.4 <- ice(fitbcfaux, Xice, predictor = "educate.4", 
                              predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.educate.4, file="ice.bcf.bart.educate.4.RData")
plot(ice.bcf.bart.educate.4)

ice.bcf.bart.educate.4.0025 <- ice(fitbcfaux, Xice, predictor = "educate.4", 
                                   predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.educate.4.0025, file="ice.bcf.bart.educate.4.0025.RData")
plot(ice.bcf.bart.educate.4.0025)

ice.bcf.bart.educate.4.0975 <- ice(fitbcfaux, Xice, predictor = "educate.4", 
                                   predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.educate.4.0975, file="ice.bcf.bart.educate.4.0975.RData")
plot(ice.bcf.bart.educate.4.0975)


#SREGION.1#

Xice = cbind(X3,pihat)
ice.bcf.bart.SREGION.1 <- ice(fitbcfaux, Xice, predictor = "SREGION.1", 
                              predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.SREGION.1, file="ice.bcf.bart.SREGION.1.RData")
plot(ice.bcf.bart.SREGION.1)

ice.bcf.bart.SREGION.1.0025 <- ice(fitbcfaux, Xice, predictor = "SREGION.1", 
                                   predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.SREGION.1.0025, file="ice.bcf.bart.SREGION.1.0025.RData")
plot(ice.bcf.bart.SREGION.1.0025)

ice.bcf.bart.SREGION.1.0975 <- ice(fitbcfaux, Xice, predictor = "SREGION.1", 
                                   predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.SREGION.1.0975, file="ice.bcf.bart.SREGION.1.0975.RData")
plot(ice.bcf.bart.SREGION.1.0975)

#SREGION.2#

Xice = cbind(X3,pihat)
ice.bcf.bart.SREGION.2 <- ice(fitbcfaux, Xice, predictor = "SREGION.2", 
                              predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.SREGION.2, file="ice.bcf.bart.SREGION.2.RData")
plot(ice.bcf.bart.SREGION.2)

ice.bcf.bart.SREGION.2.0025 <- ice(fitbcfaux, Xice, predictor = "SREGION.2", 
                                   predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.SREGION.2.0025, file="ice.bcf.bart.SREGION.2.0025.RData")
plot(ice.bcf.bart.SREGION.2.0025)

ice.bcf.bart.SREGION.2.0975 <- ice(fitbcfaux, Xice, predictor = "SREGION.2", 
                                   predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.SREGION.2.0975, file="ice.bcf.bart.SREGION.2.0975.RData")
plot(ice.bcf.bart.SREGION.2.0975)

#SREGION.3#

Xice = cbind(X3,pihat)
ice.bcf.bart.SREGION.3 <- ice(fitbcfaux, Xice, predictor = "SREGION.3", 
                              predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.SREGION.3, file="ice.bcf.bart.SREGION.3.RData")
plot(ice.bcf.bart.SREGION.3)

ice.bcf.bart.SREGION.3.0025 <- ice(fitbcfaux, Xice, predictor = "SREGION.3", 
                                   predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.SREGION.3.0025, file="ice.bcf.bart.SREGION.3.0025.RData")
plot(ice.bcf.bart.SREGION.3.0025)

ice.bcf.bart.SREGION.3.0975 <- ice(fitbcfaux, Xice, predictor = "SREGION.3", 
                                   predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.SREGION.3.0975, file="ice.bcf.bart.SREGION.3.0975.RData")
plot(ice.bcf.bart.SREGION.3.0975)

#SREGION.4#

Xice = cbind(X3,pihat)
ice.bcf.bart.SREGION.4 <- ice(fitbcfaux, Xice, predictor = "SREGION.4", 
                              predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.SREGION.4, file="ice.bcf.bart.SREGION.4.RData")
plot(ice.bcf.bart.SREGION.4)

ice.bcf.bart.SREGION.4.0025 <- ice(fitbcfaux, Xice, predictor = "SREGION.4", 
                                   predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.SREGION.4.0025, file="ice.bcf.bart.SREGION.4.0025.RData")
plot(ice.bcf.bart.SREGION.4.0025)

ice.bcf.bart.SREGION.4.0975 <- ice(fitbcfaux, Xice, predictor = "SREGION.4", 
                                   predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.SREGION.4.0975, file="ice.bcf.bart.SREGION.4.0975.RData")
plot(ice.bcf.bart.SREGION.4.0975)


#POVSTALB.1#

Xice = cbind(X3,pihat)
ice.bcf.bart.POVSTALB.1 <- ice(fitbcfaux, Xice, predictor = "POVSTALB.1", 
                              predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.POVSTALB.1, file="ice.bcf.bart.POVSTALB.1.RData")
plot(ice.bcf.bart.POVSTALB.1)

ice.bcf.bart.POVSTALB.1.0025 <- ice(fitbcfaux, Xice, predictor = "POVSTALB.1", 
                                   predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.POVSTALB.1.0025, file="ice.bcf.bart.POVSTALB.1.0025.RData")
plot(ice.bcf.bart.POVSTALB.1.0025)

ice.bcf.bart.POVSTALB.1.0975 <- ice(fitbcfaux, Xice, predictor = "POVSTALB.1", 
                                   predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.POVSTALB.1.0975, file="ice.bcf.bart.POVSTALB.1.0975.RData")
plot(ice.bcf.bart.POVSTALB.1.0975)

#POVSTALB.2#

Xice = cbind(X3,pihat)
ice.bcf.bart.POVSTALB.2 <- ice(fitbcfaux, Xice, predictor = "POVSTALB.2", 
                              predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.POVSTALB.2, file="ice.bcf.bart.POVSTALB.2.RData")
plot(ice.bcf.bart.POVSTALB.2)

ice.bcf.bart.POVSTALB.2.0025 <- ice(fitbcfaux, Xice, predictor = "POVSTALB.2", 
                                   predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.POVSTALB.2.0025, file="ice.bcf.bart.POVSTALB.2.0025.RData")
plot(ice.bcf.bart.POVSTALB.2.0025)

ice.bcf.bart.POVSTALB.2.0975 <- ice(fitbcfaux, Xice, predictor = "POVSTALB.2", 
                                   predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.POVSTALB.2.0975, file="ice.bcf.bart.POVSTALB.2.0975.RData")
plot(ice.bcf.bart.POVSTALB.2.0975)

#POVSTALB.3#

Xice = cbind(X3,pihat)
ice.bcf.bart.POVSTALB.3 <- ice(fitbcfaux, Xice, predictor = "POVSTALB.3", 
                              predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.POVSTALB.3, file="ice.bcf.bart.POVSTALB.3.RData")
plot(ice.bcf.bart.POVSTALB.3)

ice.bcf.bart.POVSTALB.3.0025 <- ice(fitbcfaux, Xice, predictor = "POVSTALB.3", 
                                   predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.POVSTALB.3.0025, file="ice.bcf.bart.POVSTALB.3.0025.RData")
plot(ice.bcf.bart.POVSTALB.3.0025)

ice.bcf.bart.POVSTALB.3.0975 <- ice(fitbcfaux, Xice, predictor = "POVSTALB.3", 
                                   predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.POVSTALB.3.0975, file="ice.bcf.bart.POVSTALB.3.0975.RData")
plot(ice.bcf.bart.POVSTALB.3.0975)

#POVSTALB.4#

Xice = cbind(X3,pihat)
ice.bcf.bart.POVSTALB.4 <- ice(fitbcfaux, Xice, predictor = "POVSTALB.4", 
                              predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.POVSTALB.4, file="ice.bcf.bart.POVSTALB.4.RData")
plot(ice.bcf.bart.POVSTALB.4)

ice.bcf.bart.POVSTALB.4.0025 <- ice(fitbcfaux, Xice, predictor = "POVSTALB.4", 
                                   predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.POVSTALB.4.0025, file="ice.bcf.bart.POVSTALB.4.0025.RData")
plot(ice.bcf.bart.POVSTALB.4.0025)

ice.bcf.bart.POVSTALB.4.0975 <- ice(fitbcfaux, Xice, predictor = "POVSTALB.4", 
                                   predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.POVSTALB.4.0975, file="ice.bcf.bart.POVSTALB.4.0975.RData")
plot(ice.bcf.bart.POVSTALB.4.0975)

#POVSTALB.5#

Xice = cbind(X3,pihat)
ice.bcf.bart.POVSTALB.5 <- ice(fitbcfaux, Xice, predictor = "POVSTALB.5", 
                              predictfcn = pred.BCF.ATE)
save(ice.bcf.bart.POVSTALB.5, file="ice.bcf.bart.POVSTALB.5.RData")
plot(ice.bcf.bart.POVSTALB.5)

ice.bcf.bart.POVSTALB.5.0025 <- ice(fitbcfaux, Xice, predictor = "POVSTALB.5", 
                                   predictfcn = pred.BCF.ATE.0025)
save(ice.bcf.bart.POVSTALB.5.0025, file="ice.bcf.bart.POVSTALB.5.0025.RData")
plot(ice.bcf.bart.POVSTALB.5.0025)

ice.bcf.bart.POVSTALB.5.0975 <- ice(fitbcfaux, Xice, predictor = "POVSTALB.5", 
                                   predictfcn = pred.BCF.ATE.0975)
save(ice.bcf.bart.POVSTALB.5.0975, file="ice.bcf.bart.POVSTALB.5.0975.RData")
plot(ice.bcf.bart.POVSTALB.5.0975)



#########################################
#########################################
#########################################
#########################################
#########################################


###########################
#########ICE Plots#########
###########################



par(cex=0.5)

png("smoke_ice_MALE.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.MALE,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y")
lines(as.numeric(names(apply(ice.bcf.bart.MALE.0975$ice_curves,2,mean))),
      apply(ice.bcf.bart.MALE.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.bart.MALE.0025$ice_curves,2,mean))),
      apply(ice.bcf.bart.MALE.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
dev.off()

png("smoke_cice_MALE.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.MALE,
     centered = T,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y - Centered")
dev.off()


png("smoke_ice_LASTAGE.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.LASTAGE,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y")
lines(as.numeric(names(apply(ice.bcf.bart.LASTAGE.0975$ice_curves,2,mean))),
      apply(ice.bcf.bart.LASTAGE.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.bart.LASTAGE.0025$ice_curves,2,mean))),
      apply(ice.bcf.bart.LASTAGE.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
dev.off()

png("smoke_cice_LASTAGE.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.LASTAGE,
     centered = T,
     ylab = "Partial Y - Centered")
dev.off()

png("smoke_dice_LASTAGE.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(ice.bcf.bart.LASTAGE),
     ylab = "Derivative Y")
dev.off()

png("smoke_ice_AGESMOKE.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.AGESMOKE,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y")
lines(as.numeric(names(apply(ice.bcf.bart.AGESMOKE.0975$ice_curves,2,mean))),
      apply(ice.bcf.bart.AGESMOKE.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.bart.AGESMOKE.0025$ice_curves,2,mean))),
      apply(ice.bcf.bart.AGESMOKE.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
dev.off()


png("smoke_cice_AGESMOKE.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.AGESMOKE,
     centered = T,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y - Centered")
dev.off()



png("smoke_ice_yearsince.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.yearsince,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y")
lines(as.numeric(names(apply(ice.bcf.bart.yearsince.0975$ice_curves,2,mean))),
      apply(ice.bcf.bart.yearsince.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.bart.yearsince.0025$ice_curves,2,mean))),
      apply(ice.bcf.bart.yearsince.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
dev.off()

png("smoke_cice_yearsince.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.yearsince,
     centered = T,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y - Centered")
dev.off()


png("smoke_ice_RACE3.1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.RACE3.1,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y")
lines(as.numeric(names(apply(ice.bcf.bart.RACE3.1.0975$ice_curves,2,mean))),
      apply(ice.bcf.bart.RACE3.1.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.bart.RACE3.1.0025$ice_curves,2,mean))),
      apply(ice.bcf.bart.RACE3.1.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
dev.off()

png("smoke_cice_RACE3.1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.RACE3.1,
     centered = T,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y - Centered")
dev.off()

png("smoke_ice_RACE3.2.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.RACE3.2,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y")
lines(as.numeric(names(apply(ice.bcf.bart.RACE3.2.0975$ice_curves,2,mean))),
      apply(ice.bcf.bart.RACE3.2.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.bart.RACE3.2.0025$ice_curves,2,mean))),
      apply(ice.bcf.bart.RACE3.2.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
dev.off()

png("smoke_cice_RACE3.2.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.RACE3.2,
     centered = T,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y - Centered")
dev.off()

png("smoke_ice_RACE3.3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.RACE3.3,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y")
lines(as.numeric(names(apply(ice.bcf.bart.RACE3.3.0975$ice_curves,2,mean))),
      apply(ice.bcf.bart.RACE3.3.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.bart.RACE3.3.0025$ice_curves,2,mean))),
      apply(ice.bcf.bart.RACE3.3.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
dev.off()

png("smoke_cice_RACE3.3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.RACE3.3,
     centered = T,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y - Centered")
dev.off()

png("smoke_ice_marital.1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.marital.1,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y")
lines(as.numeric(names(apply(ice.bcf.bart.marital.1.0975$ice_curves,2,mean))),
      apply(ice.bcf.bart.marital.1.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.bart.marital.1.0025$ice_curves,2,mean))),
      apply(ice.bcf.bart.marital.1.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
dev.off()

png("smoke_cice_marital.1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.marital.1,
     centered = T,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y - Centered")
dev.off()

png("smoke_ice_marital.2.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.marital.2,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y")
lines(as.numeric(names(apply(ice.bcf.bart.marital.2.0975$ice_curves,2,mean))),
      apply(ice.bcf.bart.marital.2.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.bart.marital.2.0025$ice_curves,2,mean))),
      apply(ice.bcf.bart.marital.2.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
dev.off()

png("smoke_cice_marital.2.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.marital.2,
     centered = T,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y - Centered")
dev.off()

png("smoke_ice_marital.3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.marital.3,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y")
lines(as.numeric(names(apply(ice.bcf.bart.marital.3.0975$ice_curves,2,mean))),
      apply(ice.bcf.bart.marital.3.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.bart.marital.3.0025$ice_curves,2,mean))),
      apply(ice.bcf.bart.marital.3.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
dev.off()

png("smoke_cice_marital.3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.marital.3,
     centered = T,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y - Centered")
dev.off()

png("smoke_ice_marital.4.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.marital.4,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y")
lines(as.numeric(names(apply(ice.bcf.bart.marital.4.0975$ice_curves,2,mean))),
      apply(ice.bcf.bart.marital.4.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.bart.marital.4.0025$ice_curves,2,mean))),
      apply(ice.bcf.bart.marital.4.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
dev.off()

png("smoke_cice_marital.4.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.marital.4,
     centered = T,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y - Centered")
dev.off()

png("smoke_ice_marital.5.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.marital.5,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y")
lines(as.numeric(names(apply(ice.bcf.bart.marital.5.0975$ice_curves,2,mean))),
      apply(ice.bcf.bart.marital.5.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.bart.marital.5.0025$ice_curves,2,mean))),
      apply(ice.bcf.bart.marital.5.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
dev.off()

png("smoke_cice_marital.5.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.marital.5,
     centered = T,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y - Centered")
dev.off()

png("smoke_ice_POVSTALB.1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.POVSTALB.1,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y")
lines(as.numeric(names(apply(ice.bcf.bart.POVSTALB.1.0975$ice_curves,2,mean))),
      apply(ice.bcf.bart.POVSTALB.1.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.bart.POVSTALB.1.0025$ice_curves,2,mean))),
      apply(ice.bcf.bart.POVSTALB.1.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
dev.off()

png("smoke_cice_POVSTALB.1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.POVSTALB.1,
     centered = T,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y - Centered")
dev.off()

png("smoke_ice_POVSTALB.2.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.POVSTALB.2,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y")
lines(as.numeric(names(apply(ice.bcf.bart.POVSTALB.2.0975$ice_curves,2,mean))),
      apply(ice.bcf.bart.POVSTALB.2.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.bart.POVSTALB.2.0025$ice_curves,2,mean))),
      apply(ice.bcf.bart.POVSTALB.2.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
dev.off()

png("smoke_cice_POVSTALB.2.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.POVSTALB.2,
     centered = T,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y - Centered")
dev.off()

png("smoke_ice_POVSTALB.3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.POVSTALB.3,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y")
lines(as.numeric(names(apply(ice.bcf.bart.POVSTALB.3.0975$ice_curves,2,mean))),
      apply(ice.bcf.bart.POVSTALB.3.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.bart.POVSTALB.3.0025$ice_curves,2,mean))),
      apply(ice.bcf.bart.POVSTALB.3.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
dev.off()

png("smoke_cice_POVSTALB.3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.POVSTALB.3,
     centered = T,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y - Centered")
dev.off()

png("smoke_ice_POVSTALB.4.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.POVSTALB.4,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y")
lines(as.numeric(names(apply(ice.bcf.bart.POVSTALB.4.0975$ice_curves,2,mean))),
      apply(ice.bcf.bart.POVSTALB.4.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.bart.POVSTALB.4.0025$ice_curves,2,mean))),
      apply(ice.bcf.bart.POVSTALB.4.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
dev.off()

png("smoke_cice_POVSTALB.4.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.POVSTALB.4,
     centered = T,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y - Centered")
dev.off()

png("smoke_ice_POVSTALB.5.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.POVSTALB.5,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y")
lines(as.numeric(names(apply(ice.bcf.bart.POVSTALB.5.0975$ice_curves,2,mean))),
      apply(ice.bcf.bart.POVSTALB.5.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.bart.POVSTALB.5.0025$ice_curves,2,mean))),
      apply(ice.bcf.bart.POVSTALB.5.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
dev.off()

png("smoke_cice_POVSTALB.5.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.POVSTALB.5,
     centered = T,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y - Centered")
dev.off()

png("smoke_ice_educate.1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.educate.1,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y")
lines(as.numeric(names(apply(ice.bcf.bart.educate.1.0975$ice_curves,2,mean))),
      apply(ice.bcf.bart.educate.1.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.bart.educate.1.0025$ice_curves,2,mean))),
      apply(ice.bcf.bart.educate.1.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
dev.off()

png("smoke_cice_educate.1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.educate.1,
     centered = T,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y - Centered")
dev.off()

png("smoke_ice_educate.2.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.educate.2,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y")
lines(as.numeric(names(apply(ice.bcf.bart.educate.2.0975$ice_curves,2,mean))),
      apply(ice.bcf.bart.educate.2.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.bart.educate.2.0025$ice_curves,2,mean))),
      apply(ice.bcf.bart.educate.2.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
dev.off()

png("smoke_cice_educate.2.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.educate.2,
     centered = T,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y - Centered")
dev.off()

png("smoke_ice_educate.3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.educate.3,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y")
lines(as.numeric(names(apply(ice.bcf.bart.educate.3.0975$ice_curves,2,mean))),
      apply(ice.bcf.bart.educate.3.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.bart.educate.3.0025$ice_curves,2,mean))),
      apply(ice.bcf.bart.educate.3.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
dev.off()

png("smoke_cice_educate.3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.educate.3,
     centered = T,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y - Centered")
dev.off()

png("smoke_ice_educate.4.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.educate.4,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y")
lines(as.numeric(names(apply(ice.bcf.bart.educate.4.0975$ice_curves,2,mean))),
      apply(ice.bcf.bart.educate.4.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.bart.educate.4.0025$ice_curves,2,mean))),
      apply(ice.bcf.bart.educate.4.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
dev.off()

png("smoke_cice_educate.4.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.educate.4,
     centered = T,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y - Centered")
dev.off()

png("smoke_ice_SREGION.1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.SREGION.1,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y")
lines(as.numeric(names(apply(ice.bcf.bart.SREGION.1.0975$ice_curves,2,mean))),
      apply(ice.bcf.bart.SREGION.1.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.bart.SREGION.1.0025$ice_curves,2,mean))),
      apply(ice.bcf.bart.SREGION.1.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
dev.off()

png("smoke_cice_SREGION.1.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.SREGION.1,
     centered = T,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y - Centered")
dev.off()

png("smoke_ice_SREGION.2.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.SREGION.2,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y")
lines(as.numeric(names(apply(ice.bcf.bart.SREGION.2.0975$ice_curves,2,mean))),
      apply(ice.bcf.bart.SREGION.2.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.bart.SREGION.2.0025$ice_curves,2,mean))),
      apply(ice.bcf.bart.SREGION.2.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
dev.off()

png("smoke_cice_SREGION.2.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.SREGION.2,
     centered = T,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y - Centered")
dev.off()

png("smoke_ice_SREGION.3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.SREGION.3,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y")
lines(as.numeric(names(apply(ice.bcf.bart.SREGION.3.0975$ice_curves,2,mean))),
      apply(ice.bcf.bart.SREGION.3.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.bart.SREGION.3.0025$ice_curves,2,mean))),
      apply(ice.bcf.bart.SREGION.3.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
dev.off()

png("smoke_cice_SREGION.3.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.SREGION.3,
     centered = T,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y - Centered")
dev.off()

png("smoke_ice_SREGION.4.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.SREGION.4,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y")
lines(as.numeric(names(apply(ice.bcf.bart.SREGION.4.0975$ice_curves,2,mean))),
      apply(ice.bcf.bart.SREGION.4.0975$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
lines(as.numeric(names(apply(ice.bcf.bart.SREGION.4.0025$ice_curves,2,mean))),
      apply(ice.bcf.bart.SREGION.4.0025$ice_curves,2,mean),
      lwd = 4, lty = 2, col = "blue")
dev.off()

png("smoke_cice_SREGION.4.png", width = 7, height = 7, units = 'in', res = 100)
plot(ice.bcf.bart.SREGION.4,
     centered = T,
     ylim=c(-0.6,0.6),
     ylab = "Partial Y - Centered")
dev.off()

plot(ice.bcf.bart.LASTAGE, centered = T, color_by = 3)