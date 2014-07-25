## Figure 9 - Parameter bounds and set membership

leg.height<-1.5
#win.metafile("fig09 pars sm.wmf",width=3.54,height=3+leg.height)
png("fig09 pars sm.png",width=3.54,height=3+leg.height,units="in",res=300)
#postscript("fig09 pars sm.eps",width=3.54,height=3+leg.height)
layout(matrix(c(1,2), nrow = 2), heights = c(leg.height/(3+leg.height), 3/(3+leg.height)))

##################################################################
## Legend
##################################################################

df.leg=data.frame(
  legend=c("Not physically\n  possible",
           "No flood occurs", 
           "Flood occurs",
           "Normative bound.",
           "T5 Parameter\n     bounds\n(epistemic bound.)\n",
           "T6 Samples from\n    feasible set\n(epistemic bound.)\n"
  ),
  fill=c("grey90","white","grey",rep(NA,3)),
  density=c(NA,NA,25,rep(NA,3)),
  angle=c(NA,NA,45,rep(NA,3)),
  border=c(NA,"grey","grey",rep(NA,3)),
  col=c(NA,NA,NA,"grey","black","black"),
  lwd=c(NA,NA,NA,2,1,NA),
  lty=c(rep(NA,3),1,1,NA),
  pch=c(rep(NA,5),"."),
  pt.cex=c(rep(NA,5),3),
  stringsAsFactors=FALSE
)
par(mar=c(0,1.1,0,1.1),cex=0.8)
plot.new()
do.call(legend,modifyList(as.list(df.leg),list(x="bottom",ncol=2,bty="n")))

##################################################################
## Normative and epistemic set diagram in parameter space
##################################################################
par(cex=0.8,mar=c(5.1,4.1,1.1,1.1))
## Normative boundary
plot.good.flood()
## Set membership
points(vv[vv$in.feasible.set,1:2],col="black",pch=".",cex=0.5)
# Parameter bounds
rect(xleft=13.5,xright=14.5,ybottom=3.5,ytop=4.5)
dev.off()