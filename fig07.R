
#win.metafile("fig07 explore par space.wmf",width=7.48,height=3)
png("fig07 explore par space.png",width=7.48,height=3,units="in",res=72)
layout(matrix(c(1,2), nrow = 1), widths = c(3.54/7.48, 3.45/7.48))
par(cex=0.8,mar=c(5.1,4.1,1.1,1.1))
plot.good.flood()
bb <- partition(1/365,7.6)
curve(bb,0,40,add=T,lwd=2,lty=2)
pbox(sdobj = myboxes.good[[1]]$lbox, xdim = 1, ydim = 2,
     boxnum = NA, fromtype = "oldbox", lwd = 2, gborder = "blue",
     mdborder = "red", col = NA)
pbox(sdobj = myboxes.good2[[1]]$lbox, xdim = 1, ydim = 2,
     boxnum = NA, fromtype = "oldbox", lwd = 2, gborder = "blue",
     mdborder = "red", col = NA)
## From pomore
points(start.pos[1],start.pos[2],col="black",pch="*",cex=1.4)
points(pom.good2[,1],pom.good2[,2],pch="+",col="red",cex=1.4)
## Legend
df.leg=data.frame(
  legend=c("Not physically possible",
           "No flood occurs", 
           "Flood occurs",
           "T1 Analytical solution of normative\n  boundary",
           "T2 PRIM boxes",
           "T3 Reference best guess model\n  scenario",
           "T3 POMORE break-even points"
  ),
  fill=c("grey90","white","grey",NA,NA,NA,NA),
  density=c(NA,NA,25,NA,NA,NA,NA),
  angle=c(NA,NA,45,NA,NA,NA,NA),
  border=c(NA,"grey","grey",NA,NA,NA,NA),
  col=c(NA,NA,NA,"black","blue","black","red"),
  lwd=c(NA,NA,NA,2,2,NA,NA),
  lty=c(rep(NA,3),2,1,rep(NA,2)),
  pch=c(rep(NA,5),"*","+"),
  stringsAsFactors=FALSE
)
par(mar=c(0,0,1.1,0))
plot.new()
do.call(legend,modifyList(as.list(df.leg),list(x="topleft")))
dev.off()



leg.height=1.4
#win.metafile("fig07 explore par space vertical.wmf",width=3.54,height=3+leg.height)
png("fig07 explore par space vertical.png",width=3.54,height=3+leg.height,res=70,units="in")
layout(matrix(c(1,2), nrow = 2), heights = c(leg.height/(3+leg.height), 3/(3+leg.height)))
## Legend
df.leg=data.frame(
  legend=c("Not physically possible",
           "No flood occurs", 
           "Flood occurs",
           "T1 Analytical solution of normative boundary",
           "T2 PRIM boxes",
           "T3 Reference best guess model scenario",
           "T3 POMORE break-even points"
  ),
  fill=c("grey90","white","grey",NA,NA,NA,NA),
  density=c(NA,NA,25,NA,NA,NA,NA),
  angle=c(NA,NA,45,NA,NA,NA,NA),
  border=c(NA,"grey","grey",NA,NA,NA,NA),
  col=c(NA,NA,NA,"black","blue","black","red"),
  lwd=c(NA,NA,NA,2,2,NA,NA),
  lty=c(rep(NA,3),2,1,rep(NA,2)),
  pch=c(rep(NA,5),"*","+"),
  stringsAsFactors=FALSE
)
par(mar=c(0,1.1,0,1.1),cex=0.8)
plot.new()
do.call(legend,modifyList(as.list(df.leg),list(x="bottom",ncol=1)))
par(cex=0.8,mar=c(5.1,4.1,1.1,1.1))
plot.good.flood()
bb <- partition(1/365,7.6)
curve(bb,0,40,add=T,lwd=2,lty=2)
pbox(sdobj = myboxes.good[[1]]$lbox, xdim = 1, ydim = 2,
     boxnum = NA, fromtype = "oldbox", lwd = 2, gborder = "blue",
     mdborder = "red", col = NA)
pbox(sdobj = myboxes.good2[[1]]$lbox, xdim = 1, ydim = 2,
     boxnum = NA, fromtype = "oldbox", lwd = 2, gborder = "blue",
     mdborder = "red", col = NA)
## From pomore
points(start.pos[1],start.pos[2],col="black",pch="*",cex=1.4)
points(pom.good2[,1],pom.good2[,2],pch="+",col="red",cex=1.4)
dev.off()