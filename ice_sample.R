library(ICEbox)
library(BART)

burn <- 2000
npost <- 1000
thin <- 100

set.seed(99)
n <- 1000
x1 <- runif(n,-1,1)
x2 <- runif(n,-1,1)
x3 <- runif(n,-1,1)
e <- rnorm(n)

y <- 0.2*x1-5*x2+10*x2*(x3>=0)+e

set.seed(99)
icebart = wbart(cbind(x1,x2,x3), y, nskip = burn,
               ndpost = npost, keepevery = thin)

save(icebart, file = "ice_sample_bart.RData")

pred.BART <- function(object, newdata){
  pred.aux <- predict(object, newdata)
  return(colMeans(pred.aux))
}

icex1 <- ice(icebart, cbind(x1,x2,x3), predictor = "x1",
             predictfcn = pred.BART)

icex2 <- ice(icebart, cbind(x1,x2,x3), predictor = "x2",
             predictfcn = pred.BART)

icex3 <- ice(icebart, cbind(x1,x2,x3), predictor = "x3",
             predictfcn = pred.BART)

save(icex1, file = "ice_sample_x1.RData")
save(icex2, file = "ice_sample_x2.RData")
save(icex3, file = "ice_sample_x3.RData")

#par(mfrow=c(2,2))
par(cex=0.5)

png("sampleicex1_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(x1,y,
     xlab = expression(x[1]),
     ylab = "Y",
     pch = 16,
     cex = 0.5)
dev.off()

png("sampleicex1_2.png", width = 7, height = 7, units = 'in', res = 100)
plot(icex1,
     ylab = "Partial Y",
     xlab = expression(x[1]))
dev.off()

png("sampleicex1_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(icex1,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[1]))
dev.off()

png("sampleicex1_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(icex1),
     ylab = "Derivative Y",
     xlab = expression(x[1]))
dev.off()


png("sampleicex2_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(x2[x3>=0],y[x3>=0],
     xlab = expression(x[2]),
     ylab = "Y",
     pch = 16,
     col = "red")
points(x2[x3<0],y[x3<=0],
       pch = 17,
       col = "blue")
dev.off()

png("sampleicex2_2.png", width = 7, height = 7, units = 'in', res = 100)
plot(icex2,
     ylab = "Partial Y",
     xlab = expression(x[2]))
dev.off()

png("sampleicex2_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(icex2,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[2]))
dev.off()

png("sampleicex2_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(icex2),
     ylab = "Derivative Y",
     xlab = expression(x[2]))
dev.off()

png("sampleicex3_1.png", width = 7, height = 7, units = 'in', res = 100)
plot(x3,y,
     xlab = expression(x[3]),
     ylab = "Y",
     pch = 16)
dev.off()

png("sampleicex3_2.png", width = 7, height = 7, units = 'in', res = 100)
plot(icex3,
     ylab = "Partial Y",
     xlab = expression(x[3]))
dev.off()

png("sampleicex3_3.png", width = 7, height = 7, units = 'in', res = 100)
plot(icex3,
     centered = T,
     ylab = "Partial Y - Centered",
     xlab = expression(x[3]))
dev.off()

png("sampleicex3_4.png", width = 7, height = 7, units = 'in', res = 100)
plot(dice(icex3),
     ylab = "Derivative Y",
     xlab = expression(x[3]))
dev.off()

png("samplebart1.png", width = 7, height = 7, units = 'in', res = 100)
plot(icebart$sigma[seq(2000,105000,by=100)], type = "l",
     ylab = expression(sigma),
     xlab = "Draw",
     ylim= c(0.9, 1.05))
dev.off()
png("samplebart2.png", width = 7, height = 7, units = 'in', res = 100)
acf(icebart$sigma[seq(2000,105000,by=100)], main = "")
dev.off()
