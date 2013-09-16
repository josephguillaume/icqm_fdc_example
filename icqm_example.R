## R Script to apply ICQM using 8 methods
##  on a simple flow-duration curve example
##
## To run in R:
## install.packages("mco","sdtoolkit")
## install.packages(c("zoo", "latticeExtra", "polynom", "car", "Hmisc","reshape"))
## install.packages("hydromad", repos="http://hydromad.catchment.org")
## install.packages("dream", repos="http://hydromad.catchment.org")
## devAskNewPage(T)
## source("icqm_example.R")
##
################################################################################
##_* Load required packages
library(hydromad) #data
library(dream)  #MCMC method (M7)
library(mco) #NSGA2 multi-objective optimisation (M3,M8)
library(sdtoolkit) #Scenario discovery (M2)
################################################################################
##_* Define functions

## Analytical solution to normatively meaningful boundary
##  Values of coefficients at which flow Qthres is exceeded
##  with probability pthres 
partition <- function(pthres,Qthres){
  p2=1/(1/12*365)
  p1=1/(100*365)
  Qlogthres=log10(Qthres)
    function(Q1) exp(-(qnorm(pthres)*log(Q1))/(qnorm(p1)-qnorm(pthres))+(qnorm(p2)*log(Q1))/(qnorm(p1)-qnorm(pthres))-(log(10)*qnorm(p2)*Qlogthres)/(qnorm(p1)-qnorm(pthres))+(log(10)*qnorm(p1)*Qlogthres)/(qnorm(p1)-qnorm(pthres)))
}

## Plot of analytical solution for a 'good flood',
##   Qthres=7.6 exceeded annually
## Requires global variable vv defining sample of whole parameter space (just to set correct scale)
plot.good.flood <- function(){
  plot(vv[,1:2],col=NA,
       xlab=bquote(q[1]~". 100 yr flood (mm)"),
       ylab=bquote(q[2]~". Monthly flow (mm)")
       )
  bb <- partition(1/365,7.6)
  polygon(x=c(0,seq(7.6,50,length.out=100),100,0),y=c(0,bb(c(seq(7.6,50,length.out=100),100)),0),
          col="blue",density=2,angle=-45)
  polygon(x=c(7.6,100,100,seq(50,7.6,length.out=100)),
          y=c(bb(7.6),bb(7.6),0,bb(seq(50,7.6,length.out=100))),
          col="grey",density=2,angle=45)
}

## pts are delta,log10(Q)
fdc.pred <- function(pts,probs=eval.prob){
    pts[1]=log10(10^pts[1]+10^pts[2])
    xs=qnorm(1/(c(100,1/12)*365))
    10^(diff(pts)/diff(xs)*(qnorm(probs)-xs[1])+pts[1])
}

## Max absolute error using parameters pts compared to data fdc.q
##  evaluated at equivalent values xs and eval.prob
max.abs.err <- function(pts){
  pts[1]=log10(10^pts[1]+10^pts[2])
  xs=qnorm(1/(c(100,1/12)*365))
  max(abs(diff(pts)/diff(xs)*(qnorm(eval.prob)-
                              xs[1])+pts[1]-log10(fdc.q)))
}


## POMORE
## Find closest parameter set on boundary
## Pareto efficient _parameters_ rather than _objectives_
## Greater spread than euclidean distance
POMOREc <- function(FUN.still.better,start.pos,pars.min,pars.max,
                   max.gen=1000,find.max=FALSE,...){

  stopifnot(length(pars.min)==length(pars.max))
  nvars <- length(pars.min)
  
  to.min <- function(x) {
    if (FUN.still.better(x,start.pos)) return(c(Inf,Inf))
    else return(abs(x-start.pos))
  }

  ## Function to be minimised which describes location of furthest point _on boundary_
  ## TODO: improvement?
  ## May be difficult to find
  to.max <- function(x) {
    if (FUN.still.better(x,start.pos)) return(c(Inf,Inf))
    else return(1/abs(x-start.pos))
  }

  if (find.max) {
    nn<-nsga2(to.max,idim=nvars,odim=nvars,
              lower.bounds=pars.min,upper.bounds=pars.max,...
              )
  } else {
    nn<-nsga2(to.min,idim=nvars,odim=nvars,
              lower.bounds=pars.min,upper.bounds=pars.max,...
              )
    
  } ##if
  return(nn)
}   ##POMORE

POMORE <- function(...){
  nn <- POMOREc(...)
  return(unique(nn$par[nn$pareto.optimal,]))
}


## pts are Q
is.good.flood <- function(pts,thres=7.6){
  pts <- log10(pts)
  xs=qnorm(1/(c(100,1/12)*365))
  10^(diff(pts)/diff(xs)*(qnorm(1/365)-xs[1])+pts[1])>=thres
}

in.feasible.set <- function(pts){
  pts <- log10(pts)
  xs=qnorm(1/(c(100,1/12)*365))
  max(abs(diff(pts)/diff(xs)*(qnorm(eval.prob)-
                              xs[1])+pts[1]-log10(fdc.q)))<0.18
}


################################################################################
##_* Prepare flow duration curve data
## Load streamflow from hydromad package
data(Murrindindi)
Q <- coredata(Murrindindi$Q)

## Sample the flow duration curve (FDC)
##  Equi-spaced points in normal distribution
## Evaluate normal distribution quantiles for the number of points in the data
pp <- qnorm(ppoints(NROW(Q))) 
## Difference between points is difference between two most extreme values
space <- diff(pp[1:2])
## Calculate probabilities as which to evaluate FDC
eval.prob <- pnorm(seq(min(pp), max(pp), by = space))
## Calculate FDC
fdc.q <- quantile(Q, 1 - eval.prob)
## Calculate quantiles of those evaluation probabilities
##  (forms x axis for linear representation)
xs <- qnorm(eval.prob)

################################################################################
##_* Fit linear log-normal model to flow duration curve data
mm <- lm(log10(fdc.q)~xs)

## Convert coefficients of the best-fit model
##  i.e. log of 100 year and monthly flow levels
pts <- c(coef(mm)[[1]]+coef(mm)[[2]]*qnorm((1/(100*365))),
           coef(mm)[[1]]+coef(mm)[[2]]*qnorm((1/(1/12*365)))
         )
## calculate actual 100 year and monthly flow levels
start.pos <- 10^pts ## Q values
## Convert for use with fdc.pred
##  i.e. log of monthly flow level
##   and log of difference between monthly and annual
pts[1]=log10(10^pts[1]-10^pts[2]) ## delta log10(Q) and log10(Q)

################################################################################
##_* Plot observations and linear model fit

##win.metafile("fig06 observed fdc.wmf",width=3.54,height=3)
par(cex=0.8,mar=c(5.1,4.1,1.1,1.1))
plot(xs,fdc.q,log="y",
     ylab="Flow (mm)",
     xlab="qnorm(Exceedance probability)")
abline(mm)
points(qnorm(1/(100*365)),quantile(Q,1-1/(100*365)),pch="+")
points(qnorm(1/(1/12*365)),quantile(Q,1-1/(1/12*365)),pch="+")
abline(v=qnorm(1/(c(100,1/12)*365)),col="grey50",lty="dashed")
##dev.off()
##_* Sample whole parameter space
vv <- expand.grid(
                  flood.diff=seq(0,26,length.out=100),
                  monthly=seq(2,6,length.out=100)
                  )
vv$in.feasible.set <- apply(vv[,1:2],1,in.feasible.set)
##_* Explore alternatives (M1,M2,M3)
##_ , Break even analysis - POMORE (M3)
pom.good <- POMORE(FUN.still.better=function(p,start.pos) !is.good.flood(p),
                   pars.min=c(1,1),pars.max=c(30,30),##find.max=T,
                   start.pos=start.pos,
                   generations=200,popsize=200,mprob=1/2
                   )
## Thin result to only 8 points
pom.good2 <- pom.good[order(pom.good[,1]),][seq(1,nrow(pom.good),length.out=8),]
##_  . Plot pomore normalised distances

##win.metafile("fig08 pomore.wmf",width=3.54,height=3)
par(cex=0.7)
par(mar=c(5.1,4.1,1.1,1.1))
boxplot(cbind(abs(pom.good[,1]-start.pos[1])/start.pos[1],
           abs(pom.good[,2]-start.pos[2])/start.pos[2]),
        names=c(expression(q[1]),expression(q[2])),
        xlab="Parameter",ylab="Normalised distance"
        )
##dev.off()
################################################################################
##_ ,Scenario discovery (M2)
## Interactive, so use saved result
load("prim.Rdata")
if(!exists("myboxes.good")){
  myboxes.good <- sdprim(x=vv[,1:2], y=1*vv$will.good.flood)
  ## n - Just pick boxes
  ## Click on box with highest density
  ## Enter the box index
  ## n - don't show scatterplot
  ## n - don't remove variables
  ## n - don't continue covering
  ## n - don't print in csv format
  seqinfo(myboxes.good)
  ##dimplot(myboxes.good, 1)
  ##scatterbox(myboxes.good,xdim=1,ydim=2)

  ## Same
  myboxes.good2 <- sdprim(x=vv[,1:2], y=1-1*vv$will.good.flood)
  seqinfo(myboxes.good2)
  scatterbox(myboxes.good2,xdim=1,ydim=2)
  dimplot(myboxes.good,1)

  save(myboxes.good,myboxes.good2,file="prim.Rdata")
}
##_ , Plot all explore alternatives

##win.metafile("fig07 explore par space.wmf",width=3.54,height=3)
par(cex=0.8,mar=c(5.1,4.1,1.1,1.1))
plot.good.flood()
bb <- partition(1/365,7.6)
curve(bb,0,40,add=T,lwd=2)
pbox(sdobj = myboxes.good[[1]]$lbox, xdim = 1, ydim = 2,
     boxnum = NA, fromtype = "oldbox", lwd = 2, gborder = "blue",
     mdborder = "red", col = NA)
pbox(sdobj = myboxes.good2[[1]]$lbox, xdim = 1, ydim = 2,
     boxnum = NA, fromtype = "oldbox", lwd = 2, gborder = "blue",
     mdborder = "red", col = NA)
## From pomore
points(start.pos[1],start.pos[2],col="black",pch="*",cex=1.4)
points(pom.good2[,1],pom.good2[,2],pch="+",col="red",cex=1.4)
##dev.off()
##_* Plot scenarios (M4)

##win.metafile("fig08 scen.wmf",width=3.5,height=3)
par(cex=0.8,mar=c(5.1,4.1,1.1,1.1))
## Plot as full areas
plot(vv[,1:2],col=NA,
     xlab=bquote(q[1]~". 100 yr flood (mm)"),
     ylab=bquote(q[2]~". Monthly flow (mm)")
     )
Qthres <- 6.6
bb <- partition(1/365,Qthres)
## Red
polygon(x=c(0,seq(Qthres,50,length.out=100),100,0),
        y=c(0,bb(c(seq(Qthres,50,length.out=100),100)),0),
        col="red",density=2,angle=-45)
## Blue
qq1 <- 6.6
qq2 <- 7.6
polygon(x=c(seq(qq1,50,length.out=100),seq(50,qq2,length.out=100),qq1),
        y=c(partition(1/365,qq1)(seq(qq1,50,length.out=100)),
          partition(1/365,qq2)(c(seq(50,qq2,length.out=100),qq1))),
        col="purple",density=2,angle=0)
## Green
Qthres <- 7.6
bb <- partition(1/365,Qthres)
polygon(x=c(Qthres,100,100,seq(50,Qthres,length.out=100)),
        y=c(bb(Qthres),bb(Qthres),bb(100),bb(seq(50,Qthres,length.out=100))),
        col="green",density=2,angle=45)
####
points(t(start.pos),pch="*")
points(t(c(17,4.3)),pch="r")
points(t(c(16.5,4.2)),pch="L")
points(t(c(14.5,3.6)),pch="W")
points(t(c(flood = 16.2020202020202, monthly = 4.38383838383838)),pch="e")
points(t(c(11.13131,3.777778)),pch="i")
##dev.off()
################################################################################
##_* Plot par bounds (M5),set membership (M6)

##win.metafile("fig09 pars sm.wmf",width=3.54,height=3)
par(cex=0.8,mar=c(5.1,4.1,1.1,1.1))
## Partition space for good flood + parameter bounds + set membership
##plot(vv[,1:2],col=ifelse(vv$will.good.flood,"grey","blue"))
plot(vv[,1:2],col=NA,
     xlab=bquote(q[1]~". 100 yr flood (mm)"),
     ylab=bquote(q[2]~". Monthly flow (mm)")
     )
bb <- partition(1/365,7.6)
polygon(x=c(0,seq(7.6,50,length.out=100),100,0),y=c(0,bb(c(seq(7.6,50,length.out=100),100)),0),
        col="blue",density=2,angle=-45)
polygon(x=c(7.6,100,100,seq(50,7.6,length.out=100)),
        y=c(bb(7.6),bb(7.6),0,bb(seq(50,7.6,length.out=100))),
        col="grey",density=2,angle=45)
points(vv[vv$in.feasible.set,1:2],col="black",pch=".",cex=0.5)
rect(xleft=13.5,xright=14.5,ybottom=3.5,ytop=4.5)
##dev.off()
################################################################################
##_* MCMC (M7)
## Initial sample
dd <- dreamCalibrate(function(pts) log10(fdc.pred(pts)),
                     list(flood.diff=c(0,2),monthly=c(0,1)),
                     log10(fdc.q),
                     lik.fun=calc.rmse,
                     control=list(nseq=5))

## Continue sample
dd$control$thin.t <- 1 ## Workaround because simulate thins sequences but not hist.log
dd2 <- simulate(dd,nsim=2e4)

## Thin to effective size
mcmc <- window(dd2,fraction=1)
logp <- as.mcmc(window(dd2$hist.logp,start=start(mcmc)))
effsz <- effectiveSize(mcmc)
thin <- 2 * ceiling(max(nrow(mcmc[[1]])/effsz))
mcmc <- window(mcmc, thin = thin)
logp <- window(logp, thin = thin)
psets <- as.matrix(mcmc)
objseq <- as.vector(logp)

##_ , Confidence intervals on annual
## Calculate log flow
ys <- apply(psets,1,fdc.pred,probs=1/365)
level=0.95
keep <- ys>=quantile(ys,(1-level)/2) & ys<=quantile(ys,1-(1-level)/2)
## Calculate pareto front
paret.max <- hydromad::paretoFilter(cbind(-ys,-objseq)) ##slow
paret.min <- hydromad::paretoFilter(cbind(ys,-objseq))  ##slow
pareto.front <- paret.max|paret.min
## Likelihood threshold for hyp testing
thres.lik.lo <- min(objseq[which(keep & pareto.front)[order(ys[keep & pareto.front],decreasing=F)][1]]) ## Lik of lowest y value on pareto front
thres.lik.hi <- min(objseq[which(keep & pareto.front)[order(ys[keep & pareto.front],decreasing=T)][1]]) ## Lik of highest y value on pareto front
################################################################################
##_ , Plot of density, confidence interval/critical region and corresponding models
##win.metafile("fig11 mcmc.wmf",width=7,height=3)
lik.ht.pars <- matrix(c(8.44598275902987, 21.4774398003501, 2.92690003299676, 4.99546011684257),ncol=2) # hard-coded result from M8b because they are mutually dependent
par(cex=0.8,mar=c(5.1,4.1,1.1,1.1))
par(mfrow=c(1,2))
plot(density(ys),xlab="Annual flow (mm)\nN=3335, Bandwidth=0.22",main="")
abline(v=quantile(ys,c((1 - level)/2,1 - (1 - level)/2)))
abline(v=7.6,lwd=3)
mtext(" a) ",side=3,adj=0,line=-1,cex=0.8)
keep <- ys>=quantile(ys,(1-level)/2) & ys<=quantile(ys,1-(1-level)/2)
psets2 <- cbind(rowSums(10^psets),10^psets[,2])
plot(psets2,cex=0.4,##cex=1,pch=".",
     xlab=bquote(q[1]~". 100 yr flood (mm)"),
     ylab=bquote(q[2]~". Monthly flow (mm)"),
     col=ifelse(keep,"grey50","red"),
     pch=ifelse(keep,19,1)
     )
points(lik.ht.pars,pch="+",col="black",cex=1.5)
bb <- partition(1/365,7.6)
polygon(x=c(0,seq(7.6,50,length.out=100),100,0),y=c(0,bb(c(seq(7.6,50,length.out=100),100)),0),
        col="blue",density=2,angle=-45)
polygon(x=c(7.6,100,100,seq(50,7.6,length.out=100)),
        y=c(bb(7.6),bb(7.6),0,bb(seq(50,7.6,length.out=100))),
        col="grey",density=2,angle=45)
mtext(" b) ",side=3,adj=0,line=-1,cex=0.8)
##dev.off()
################################################################################
##_* Hypothesis testing results (M8)
## Using log10(Q), without difference
##_ , With value directly - set membership (M8a)
trad.minSM <- nsga2(function(p) c(max.abs.err(p),fdc.pred(p,probs=1/365)),
                  idim=2,odim=2,
                  lower.bounds=c(0,0),upper.bounds=c(2,1),
                  generations=200,popsize=200,mprob=1/2
                  )
trad.maxSM <- nsga2(function(p) c(max.abs.err(p),-fdc.pred(p,probs=1/365)),
                  idim=2,odim=2,
                  lower.bounds=c(0,0),upper.bounds=c(2,1),
                  generations=200,popsize=200,mprob=1/2
                  )
##_ , With value directly - likelihood (M8b)
trad.minLik <- nsga2(function(p) c(-calc.rmse(log10(fdc.pred(p)),log10(fdc.q)),
                            fdc.pred(p,probs=1/365)),
              idim=2,odim=2,
              lower.bounds=c(0,0),upper.bounds=c(2,1)
              )
trad.maxLik <- nsga2(function(p) c(-calc.rmse(log10(fdc.pred(p)),log10(fdc.q)),
                            -fdc.pred(p,probs=1/365)),
              idim=2,odim=2,
              lower.bounds=c(0,0),upper.bounds=c(2,1)
              )
##_ , Maximising fit to threshold - good flood (M8c)
set.seed(1)
trad.fit <- nsga2(function(p) c(-calc.rmse(log10(fdc.pred(p)),log10(fdc.q)),
                                abs(fdc.pred(p,probs=1/365)-7.6)),
                  idim=2,odim=2,
                  lower.bounds=c(0,0),upper.bounds=c(2,1),
                  generations=500,popsize=400,mprob=1/2
                  )
##_ , Plot of all hypothesis testing results (M8a,b,c)

##win.metafile("fig12 hyp test.wmf",width=7,height=5)
par(cex=0.8,mar=c(5.1,4.1,1.1,1.1))
## Combined hypothesis testing plot
par(mfrow=c(2,3))
## Par space
## SM
plot.good.flood()
pars.sm <-as.matrix(10^rbind(trad.minSM$par,trad.maxSM$par))
pars.sm <- cbind(pars.sm[,1]+pars.sm[,2],pars.sm[,2]) ## convert from delta
obj.sm <- c(trad.minSM$value[,1],trad.maxSM$value[,1])
pred.sm <- c(trad.minSM$value[,2],-trad.maxSM$value[,2])
points(pars.sm,pch=19,col="black",cex=0.7)
ord <- order(pred.sm)
w.extreme.sm <- c(head(ord[ord %in% which(obj.sm<0.18)],1),tail(ord[ord %in% which(obj.sm<0.18)],1))
points(pars.sm[w.extreme.sm,],pch="+",col="blue",cex=2)
mtext(" a) ",side=3,adj=0,line=-1.5,cex=0.8)
## Lik
plot.good.flood()
pars <-as.matrix(10^rbind(trad.minLik$par,trad.maxLik$par))
obj.lik <- c(-trad.minLik$value[,1],-trad.maxLik$value[,1])
pred.lik <- c(trad.minLik$value[,2],-trad.maxLik$value[,2])
pars <- cbind(pars[,1]+pars[,2],pars[,2])
points(pars,pch=19,col="black",cex=0.7)
ord <- order(pred.lik)
w.extreme.lik <- c(head(ord[ord %in% which(obj.lik>-3.05)],1),
               tail(ord[ord %in% which(obj.lik>-2.37)],1))
lik.ht.pars <- pars[w.extreme.lik,]
points(lik.ht.pars,pch="+",col="blue",cex=2)
mtext(" c) ",side=3,adj=0,line=-1.5,cex=0.8)
## Fit
pars <-as.matrix(10^unique(trad.fit$par))
obj.fit <- c(-trad.fit$value[,1])
pars <- cbind(pars[,1]+pars[,2],pars[,2])
plot.good.flood()
points(pars,pch=19,col="black",cex=0.7)
w.extreme.fit <- c(which.min(obj.fit),which.max(obj.fit))
points(pars[w.extreme.fit,],pch="+",col="blue",cex=2)
mtext(" e) ",side=3,adj=0,line=-1.5,cex=0.8)
pred.fit <- trad.fit$value[,2]
## Tradeoffs
## SM
plot(pred.sm,obj.sm,xlim=c(0,30),col="black",pch=19,cex=0.7,
     xlab="Annual flow (mm)",ylab="Max absolute error (mm)")
abline(v=7.6,lwd=3)
abline(h=0.18,col="grey",lty="dashed")
mtext(" b) ",side=3,adj=0,line=-1.5,cex=0.8)
points(pred.sm[w.extreme.sm],
       obj.sm[w.extreme.sm],
       pch="+",col="blue",cex=2)
## Lik
plot(pred.lik,obj.lik,xlim=c(0,30),col="black",pch=19,cex=0.7,
     xlab="Annual flow (mm)",ylab="Log Likelihood"
     )
abline(v=7.6,lwd=3)
segments(x0=-10,x1=7,y0=thres.lik.lo,y1=thres.lik.lo,
         col="grey",lty="dashed")
abline(h=thres.lik.hi,col="grey",lty="dashed")
mtext(" d) ",side=3,adj=0,line=-2,cex=0.8)
points(pred.lik[w.extreme.lik],
       obj.lik[w.extreme.lik],
       pch="+",col="blue",cex=2)
##plot(vv[,1:2],col=ifelse(vv$will.good.flood,"grey","blue"))
##plot(cbind(rowSums(10^psets),10^psets[,2]),col=ifelse(keep,"black","red"),pch=".")
## Fit
plot(pred.fit,obj.fit,pch=19,cex=0.7,
     xlab="Difference from threshold flow (mm)",
     ylab="Fit to observations (Log Likelihood)")
abline(v=0,lwd=3)
mtext("   f) ",side=3,adj=0,line=-1.5,cex=0.8)
points(pred.fit[w.extreme.fit],
       obj.fit[w.extreme.fit],
       pch="+",col="blue",cex=2)
##dev.off()



