## Combined hypothesis testing plot
leg.height=2.4
##win.metafile("fig11 hyp test.wmf",width=7,height=5+leg.height)
png("fig11 hyp test.png",width=7,height=5+leg.height,units="in",res=70)
layout(matrix(c(1,2,5,1,3,6,1,4,7), nrow = 3), 
       heights = c(leg.height/(5+leg.height), 2.7/(2.5+leg.height),2.3/(2.5+leg.height)))

##################################################################
## Legend
##################################################################
df.leg=data.frame(
  legend=c("Not physically possible",
           "No flood occurs", 
           "Flood occurs",
           "Normative boundary",
           "Threshold of acceptable performance\n  (epistemic boundary)\n",
           "T8a, T8b Selected model instances at limits\n  of epistemic boundary\n",
           "T8c Selected model instances where user has\n to check whether epistemic boundary is
  breached, i.e. performance is unacceptable",
           "Sampled points with each method"
  ),
  fill=c("grey90","white","grey",rep(NA,5)),
  density=c(NA,NA,25,rep(NA,5)),
  angle=c(NA,NA,45,rep(NA,5)),
  border=c(NA,"grey","grey",rep(NA,5)),
  col=c(NA,NA,NA,"grey","blue","blue","blue","black"),
  pch=c(rep(NA,3),rep(NA,2),3,4,19), #3 should be +
  pt.cex=c(rep(NA,3),rep(NA,2),2,2,0.7),
  lwd=c(rep(NA,3),2,1,rep(NA,3)),
  lty=c(rep(NA,3),1,2,rep(NA,3)),
  stringsAsFactors=FALSE
)
par(mar=c(0,1.1,0,1.1),cex=0.8)
plot.new()
do.call(legend,modifyList(as.list(df.leg),list(x="bottom",ncol=2,trace=TRUE,bty="n")))

##################################################################
## Normative and epistemic set diagrams in parameter space
##################################################################
par(cex=0.8,mar=c(5.1,4.1,3.1,1.1))

## Set membership
plot.good.flood()
points(pars.sm,pch=19,col="black",cex=0.7)
points(pars.sm[w.extreme.sm,],pch=3,col="blue",cex=2)
mtext(" a) ",side=3,adj=0,line=-1.5,cex=0.8)
title("T8a Optimisation - extreme\nvalues with set membership",cex.main=0.9)

## Likelihood
plot.good.flood()
points(pars.lik,pch=19,col="black",cex=0.7)
points(pars.lik[w.extreme.lik,],pch=3,col="blue",cex=2)
mtext(" b) ",side=3,adj=0,line=-1.5,cex=0.8)
title("T8b Optimisation - extreme\nvalues with likelihood",cex.main=0.9)

## Fit
plot.good.flood()
points(pars.fit,pch=19,col="black",cex=0.7)
points(pars.fit[w.extreme.fit,],pch=4,col="blue",cex=2)
mtext(" c) ",side=3,adj=0,line=-1.5,cex=0.8)
title("T8c Optimisation - \nfitting threshold",cex.main=0.9)

##################################################################
## Tradeoffs between fit and prediction
##################################################################
par(cex=0.8,mar=c(5.1,4.1,0,1.1))

## Set Membership
plot(pred.sm,obj.sm,xlim=c(0,30),col="black",pch=19,cex=0.7,
     xlab="Annually recurring runoff\n(mm/day)",ylab="Max absolute error (mm/day)")
polygon(x=c(7.6,100,100,7.6),y=c(-100,-100,100,100),col="grey",density=2,angle=45)
abline(v=7.6,lwd=2,col="grey")
abline(h=0.18,col="blue",lty="dashed")
mtext(" d) ",side=3,adj=0,line=-1.5,cex=0.8)
points(pred.sm[w.extreme.sm],
       obj.sm[w.extreme.sm],
       pch=3,col="blue",cex=2)

## Likelihood
plot(pred.lik,obj.lik,xlim=c(0,30),col="black",pch=19,cex=0.7,
     xlab="Annually recurring runoff\n(mm/day)",ylab="Log Likelihood"
)
polygon(x=c(7.6,100,100,7.6),y=c(-100,-100,100,100),col="grey",density=2,angle=45)
abline(v=7.6,lwd=2,col="grey")
segments(x0=-10,x1=7,y0=thres.lik.lo,y1=thres.lik.lo,
         col="grey",lty="dashed")
abline(h=thres.lik.hi,col="blue",lty="dashed")
mtext(" e) ",side=3,adj=0,line=-2,cex=0.8)
points(pred.lik[w.extreme.lik],
       obj.lik[w.extreme.lik],
       pch=3,col="blue",cex=2)

## Fit
plot(pred.fit,obj.fit,pch=19,cex=0.7,
     xlab="Difference from threshold\n runoff (mm/day)",
     ylab="Fit to observations (Log Lik.)")
polygon(x=c(-100,0,0,-100),y=c(-100,-100,100,100),col="grey",density=2,angle=45)
abline(v=0,lwd=2,col="grey")
mtext("   f) ",side=3,adj=0,line=-1.5,cex=0.8)
points(pred.fit[w.extreme.fit],
       obj.fit[w.extreme.fit],
       pch=4,col="blue",cex=2)


dev.off()


