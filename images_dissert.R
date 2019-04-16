set.seed(12345)
par(mfrow=c(1,2))

n <- 100

x <- rnorm(n,10,2)
x <- sort(x)

y <- x + rnorm(n)
intervals <- predict(lm(y~x),newdata=as.data.frame(x), interval = c("confidence"))
bart1 <- wbart(x, y, nskip = 1000, ndpost = 2000)

pdf("intro1.pdf")
plot(x,y,
     xlab = expression(x[1])
     )
abline(lm(y~x), col = "red", lwd=2)
lines(x, intervals[,2], col = "red", lwd=2, lty = 2)
lines(x, intervals[,3], col = "red", lwd=2, lty = 2)
dev.off()

pdf("intro2.pdf")
plot(x,y,
     xlab = expression(x[1])
)
lines(x, bart1$yhat.train.mean, col = "red", lwd=2, lty = 1)
lines(x, apply(bart1$yhat.train,2,quantile, probs = 0.975), col = "red", lwd=2, lty = 2)
lines(x, apply(bart1$yhat.train,2,quantile, probs = 0.025), col = "red", lwd=2, lty = 2)
dev.off()


y <- x/3*sin(x)+rnorm(n)
intervals <- predict(lm(y~x),newdata=as.data.frame(x), interval = c("confidence"))
pdf("intro3.pdf")
plot(x,y,
     xlab = expression(x[1])
     )
abline(lm(y~x), col = "red", lwd=2)
lines(x, intervals[,2], col = "red", lwd=2, lty = 2)
lines(x, intervals[,3], col = "red", lwd=2, lty = 2)
dev.off()


library(BART)
bart2 <- wbart(x, y, nskip = 1000, ndpost = 2000)
pdf("intro4.pdf")
plot(x,y,
     xlab = expression(x[1])
)
lines(x, bart2$yhat.train.mean, col = "red", lwd=2, lty = 1)
lines(x, apply(bart2$yhat.train,2,quantile, probs = 0.975), col = "red", lwd=2, lty = 2)
lines(x, apply(bart2$yhat.train,2,quantile, probs = 0.025), col = "red", lwd=2, lty = 2)
dev.off()

