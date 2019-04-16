par(mfrow=c(1,3))

plot.default(2,xlim = c(0,1), ylim = c(0,1), xlab = expression(x[1]), ylab = expression(x[2]), frame = F)
box()
abline(v = 0.7)
segments(x0 = 0.7, y0 = 0.6, x1 = 1.05, y1 = 0.6)

text(list(x=0.85,y=0.3),expression(mu[21]), cex=1)
text(list(x=0.85,y=0.8),expression(mu[31]), cex=1)
text(list(x=0.35,y=0.5),expression(mu[11]), cex=1)


plot.default(2,xlim = c(0,1), ylim = c(0,1), xlab = expression(x[1]), ylab = expression(x[2]), frame = F)
box()
abline(v = 0.3)
segments(x0 = -0.05, y0 = 0.4, x1 = 0.3, y1 = 0.4)

text(list(x=0.15,y=0.2),expression(mu[12]), cex=1)
text(list(x=0.15,y=0.7),expression(mu[22]), cex=1)
text(list(x=0.7,y=0.5),expression(mu[32]), cex=1)

plot.default(2,xlim = c(0,1), ylim = c(0,1), xlab = expression(x[1]), ylab = expression(x[2]), frame = F)
box()
abline(v = 0.3)
segments(x0 = -0.05, y0 = 0.4, x1 = 0.3, y1 = 0.4)
abline(v = 0.7)
segments(x0 = 0.7, y0 = 0.6, x1 = 1.05, y1 = 0.6)


text(list(x=0.13,y=0.2),expression(mu[11]+mu[12]), cex=1)
text(list(x=0.13,y=0.7),expression(mu[11]+mu[22]), cex=1)
text(list(x=0.5,y=0.5),expression(mu[11]+mu[32]), cex=1)
text(list(x=0.87,y=0.3),expression(mu[21]+mu[32]), cex=1)
text(list(x=0.87,y=0.8),expression(mu[31]+mu[32]), cex=1)



library(DescTools)

plot.new()
segments(x0 = 0.5, y0 = 0.9, x1 = 0.25, y1 = 0.6)
DrawRegPolygon(x = 0.375, y = 0.75, radius.x = 0.05, radius.y = 0.05, rot = 0, nv = 4)
segments(x0 = 0.5, y0 = 0.9, x1 = 0.75, y1 = 0.6)
DrawRegPolygon(x = 0.625, y = 0.75, radius.x = 0.05, radius.y = 0.05, rot = 0, nv = 4)
segments(x0 = 0.75, y0 = 0.6, x1 = 0.625, y1 = 0.3)
DrawRegPolygon(x = 0.6875, y = 0.45, radius.x = 0.05, radius.y = 0.05, rot = 0, nv = 4)
segments(x0 = 0.75, y0 = 0.6, x1 = 0.875, y1 = 0.3)
DrawRegPolygon(x = 0.8125, y = 0.45, radius.x = 0.05, radius.y = 0.05, rot = 0, nv = 4)

DrawCircle(x=0.5, y=0.9, r.out=0.075, col = c("white"))
DrawCircle(x=0.25, y=0.6, r.out=0.05, col = c("white"))
DrawCircle(x=0.75, y=0.6, r.out=0.075, col = c("white"))
DrawCircle(x=0.625, y=0.3, r.out=0.05, col = c("white"))
DrawCircle(x=0.875, y=0.3, r.out=0.05, col = c("white"))


text(list(x=0.5,y=0.9),expression(paste(x[1], "< 0.4")))
text(list(x = 0.375, y = 0.75),"yes")
text(list(x=0.25,y=0.6),expression(mu[1]))
text(list(x = 0.625, y = 0.75),"no")
text(list(x=0.75,y=0.6),expression(paste(x[2], "< 0.6")))
text(list(x = 0.6875, y = 0.45),"yes")
text(list(x=0.625,y=0.3),expression(mu[2]))
text(list(x = 0.8125, y = 0.45),"no")
text(list(x=0.875,y=0.3),expression(mu[3]))

