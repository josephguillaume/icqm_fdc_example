leg.height=1.3
png("fig05 observed fdc.png",width=3.54,height=3+leg.height,units="in",res=70)
##win.metafile("fig05 observed fdc.wmf",width=3.54,height=3+leg.height)
layout(matrix(c(1,2), nrow = 2), heights = c(leg.height/(3+leg.height), 3/(3+leg.height)))

##################################################################
## Legend
##################################################################
df.leg=list(
  legend=as.expression(c(
    "No flood occurs", 
    "Flood occurs",
    "Flood threshold (normative boundary)",
    "Empirical flow duration curve",
    "Best-fit log-normal approximation",
    bquote(q[1]~". 100 yr flood"),
    bquote(q[2]~". Monthly recurring runoff"),
    "Annually recurring runoff"
  )),
  fill=c("white","grey",rep(NA,6)),
  density=c(NA,25,rep(NA,6)),
  angle=c(NA,45,rep(NA,6)),
  border=c("grey","grey",rep(NA,6)),
  col=c(NA,NA,"grey",rep("black",2),rep("grey50",2),"black"),
  lwd=c(NA,NA,2,NA,1,1,1,2),
  lty=c(NA,NA,1,NA,1,2,3,1),
  pch=c(NA,NA,NA,1,rep(NA,4))
)
par(mar=c(0,1.1,0,1.1),cex=0.8)
plot.new()
do.call(legend,modifyList(as.list(df.leg),list(x="center",ncol=1,bty="n")))

##################################################################
## Log-normal plot of flow duration curve
##################################################################
par(cex=0.8,mar=c(5.1,5.1,1.1,1.1))
plot(xs,log10(fdc.q),
     yaxt="n",ylim=c(-0.65,log10(15)),
     ylab=NA,
     xlab="qnorm(Exceedance probability)")
axis(2,at=log10(c(0.5,1.0,2.0,5.0,10.0)),labels=c("0.5","1.0","2.0","5.0","10.0"),cex.axis=0.8)
title(ylab="Runoff (mm/day)\non log10 scale)",line=2.5)
polygon(y=log10(c(7.6,100,100,7.6)),x=c(-100,-100,100,100),col="grey",density=8,angle=45)
abline(h=log10(7.6),lwd=2,col="grey")
points(xs,log10(fdc.q))
abline(mm)
abline(v=qnorm(1/(c(100)*365)),col="grey50",lty="dashed")
abline(v=qnorm(1/(c(1/12)*365)),col="grey50",lty="dotted")
abline(v=qnorm(1/365),col="black",lty="solid",lwd=2)


dev.off()