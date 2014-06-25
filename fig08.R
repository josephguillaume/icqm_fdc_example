## Plot scenarios (M4)

#win.metafile("fig08 scen.wmf",width=7,height=6)
png("fig08 scen.png",width=7,height=6,units="in",res=70)
layout(matrix(c(2,4,3,1), nrow = 2), heights = c(3, 3))

##################################################################
## Legend
##################################################################

df.leg=data.frame(
  legend=c("Not physically possible",
           "No flood occurs", 
           "Flood occurs",
           "Normative boundary",
           "Status quo",
           "New solution to provide flooding",
           "Possible failure of new solution",
           "Extreme scenarios of trade-off \nunder water scarcity"
  ),
  fill=c("grey90","white","grey",rep(NA,5)),
  density=c(NA,NA,25,rep(NA,5)),
  angle=c(NA,NA,45,rep(NA,5)),
  border=c(NA,"grey","grey",rep(NA,5)),
  col=c(NA,NA,NA,"grey","black","black","black","black"),
  lwd=c(NA,NA,NA,2,rep(NA,4)),
  lty=c(rep(NA,3),1,rep(NA,4)),
  pch=c(rep(NA,4),1,13,4,19),
  stringsAsFactors=FALSE
)
par(mar=c(0,1.1,0,1.1),cex=0.8)
plot.new()
do.call(legend,modifyList(as.list(df.leg),list(x="center",ncol=1,bty="n")))

##################################################################
## Normative and epistemic set diagrams in parameter space
##################################################################

par(cex=0.8,mar=c(5.1,4.1,1.1,1.1))

## Engineering works
plot.good.flood()
mtext(" a) ",side=3,adj=0,line=-1.5,cex=0.8)
Qthres <- 6.6
bb <- partition(1/365,Qthres)
lines(x=seq(6.6,50,length.out=100),y=bb(seq(6.6,50,length.out=100)),col="grey",lwd=2)
arrows(x0=12.5,y0=partition(1/365,7.6)(12.5),x1=10,y1=bb(10),col="grey10",length=0.1,angle=25)
text(13,5.3,labels="Engineering works\nreduce required\nflow for flooding",pos=4)
points(t(start.pos),pch=1)
text(t(start.pos),labels="Status quo",adj=c(1.2,0.2))
points(t(c(14.5,3.6)),pch=4)
text(t(c(14.5,3.6)),labels="Engineering works\ncause unexpected\nchanges in flow\nregime",adj=c(0.9,1.1))

## Releases
plot.good.flood()
mtext(" b) ",side=3,adj=0,line=-1.5,cex=0.8)
points(t(start.pos),pch=1)
text(t(start.pos),labels="Status quo",adj=c(1.1,-0.2))
points(t(c(17,4.3)),pch=13)
text(t(c(17,4.3)),labels="Environmental\nflow release\nsucceeds",adj=c(0,-0.5))
points(t(c(16.5,4.2)),pch=4)
text(t(c(16.5,4.2)),labels="Flow losses\ncause release\nto fail",adj=c(0.1,1.2))

## Marginal scenarios
plot.good.flood()
mtext(" c) ",side=3,adj=0,line=-1.5,cex=0.8)
points(t(c(flood = 16.2020202020202, monthly = 4.38383838383838)),pch=19)
text(t(c(flood = 16.2020202020202, monthly = 4.38383838383838)),labels="Scarce water used\nfor environmental\nflow release. Farmers\nreceive government\nsupport.",adj=c(0.2,-0.1))
points(t(c(11.13131,3.777778)),pch=19)
text(t(c(11.13131,3.777778)),labels="Scarce water used for\nirrigation. Groundwater\nsustains ecosystem.",adj=c(0.5,1.1))

dev.off()
