## T7 MCMC

leg.height<-1
#win.metafile("fig11 mcmc.wmf",width=7,height=3+leg.height)
png("fig11 mcmc.png",width=7,height=5+leg.height,units="in",res=70)

layout(matrix(c(1,2,3,4), nrow = 2), heights = c(leg.height/(3+leg.height), 3/(3+leg.height)))

##################################################################
## Legend for posterior density plot
##################################################################
df.leg=data.frame(
  legend=c("Flood threshold\n  (normative boundary)\n",
           "Confidence interval\n  (epistemic boundary)\n"
  ),
  lwd=c(2,1),
  lty=c(1,2),
  col=c("grey","blue"),
  stringsAsFactors=FALSE
)
par(mar=c(0,1.1,0,1.1),cex=0.8)
plot.new()
do.call(legend,modifyList(as.list(df.leg),list(x="bottom",ncol=1,bty="n")))

##################################################################
## Posterior density plot
##################################################################

par(cex=0.8,mar=c(5.1,4.1,1.1,1.1))
plot(density(ys),xlab="Annual flow (mm)\nN=3335, Bandwidth=0.22",main="")
points(ys,runif(length(ys),min=-0.005,max=0.005),
       col=ifelse(keep,"grey50","red"),
       pch=ifelse(keep,19,1)
       )
polygon(x=c(7.6,100,100,7.6),y=c(-100,-100,100,100),col="grey",density=2,angle=45)
abline(v=quantile(ys,c((1 - level)/2,1 - (1 - level)/2)),lty=2,col="blue")
abline(v=7.6,lwd=2,col="grey")
mtext(" a) ",side=3,adj=0,line=-1,cex=0.8)
points(lik.ht.pars.ys,c(0,0),pch=3,col="blue",cex=1.5)

##################################################################
## Legend for set diagram
##################################################################

df.leg=data.frame(
  legend=c("Not physically possible",
           "No flood occurs", 
           "Flood occurs",
           "Points in critical region",
           "Points in confidence interval",
           "Critical values"
  ),
  fill=c("grey90","white","grey",NA,NA,NA),
  density=c(NA,NA,25,NA,NA,NA),
  angle=c(NA,NA,45,NA,NA,NA),
  border=c(NA,"grey","grey",NA,NA,NA),
  col=c(NA,NA,NA,"red","grey50","blue"),
  pch=c(rep(NA,3),1,19,3), #3 should be +
  pt.cex=c(rep(NA,3),1,1,1.5),
  stringsAsFactors=FALSE
)
par(mar=c(0,1.1,0,1.1),cex=0.8)
plot.new()
do.call(legend,modifyList(as.list(df.leg),list(x="bottom",ncol=1,bty="n")))

##################################################################
## Normative and epistemic set diagrams in parameter space
##################################################################

par(cex=0.8,mar=c(5.1,4.1,1.1,1.1))
plot.good.flood()
points(psets2,cex=0.4,##cex=1,pch=".",
       col=ifelse(keep,"grey50","red"),
       pch=ifelse(keep,19,1)
)
points(lik.ht.pars,pch=3,col="blue",cex=1.5)
mtext(" b) ",side=3,adj=0,line=-1,cex=0.8)
dev.off()