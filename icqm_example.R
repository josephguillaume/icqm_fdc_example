## R Script to apply ICQM using 8 methods
##  on a simple flow-duration curve example

source("packages.R")

source("functions.R")

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
## Calculate probabilities at which to evaluate FDC
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
source("fig05.R")

################################################################################
##_* Sample whole parameter space
vv <- expand.grid(
                  flood.diff=seq(0,26,length.out=100),
                  monthly=seq(2,6,length.out=100)
                  )
vv$in.feasible.set <- apply(vv[,1:2],1,in.feasible.set)

################################################################################
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

png("fig07 pomore.png",width=3.54,height=3,units="in",res=300)
##win.metafile("fig07 pomore.wmf",width=3.54,height=3)
##postscript("fig07 pomore.eps",width=3.54,height=3)
par(cex=0.7)
par(mar=c(5.1,4.1,1.1,1.1))
boxplot(cbind(abs(pom.good[,1]-start.pos[1])/start.pos[1],
           abs(pom.good[,2]-start.pos[2])/start.pos[2]),
        names=NA,
        xlab=NA,ylab="Normalised distance"
        )
axis(side=1,at=1:2,padj=0.6,cex.axis=0.8,
     labels=as.expression(c(bquote(atop(q[1]~"100 yr flood","(mm/day)")),
                                   bquote(atop(q[2]~"Monthly recurring","runoff (mm/day)"))
     )))
title(xlab="Parameters",line=3.5)
dev.off()

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
source("fig06.R")

##_* Plot scenarios (M4)

source("fig08.R")

##_* Plot par bounds (M5),set membership (M6)
source("fig09.R")

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

psets2 <- cbind(rowSums(10^psets),10^psets[,2])

## Parameters of lowest and highest y values (i.e. critical values)
lik.crit.pars <- psets[keep & pareto.front,][c(
  which.min(ys[keep & pareto.front]),
  which.max(ys[keep & pareto.front])
),]
lik.crit.pars.ys <- apply(lik.crit.pars,1,fdc.pred,probs=1/365)
lik.crit.pars.orig <- cbind(rowSums(10^lik.crit.pars),10^lik.crit.pars[,2])

##_ , Plot of density, confidence interval/critical region and corresponding models
source("fig10.R")


## Likelihood threshold for hyp testing
## The results with the lowest and highest y values on the pareto 
##   front and within the confidence interval
thres.lik.lo <- objseq[keep & pareto.front][which.min(ys[keep & pareto.front])]
thres.lik.hi <- objseq[keep & pareto.front][which.max(ys[keep & pareto.front])]


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

pars.sm <-as.matrix(10^rbind(trad.minSM$par,trad.maxSM$par))
pars.sm <- cbind(pars.sm[,1]+pars.sm[,2],pars.sm[,2]) ## convert from delta
obj.sm <- c(trad.minSM$value[,1],trad.maxSM$value[,1])
pred.sm <- c(trad.minSM$value[,2],-trad.maxSM$value[,2])

## Select the model scenarios with the min and max annual runoff within feasible set
ord <- order(pred.sm)
w.extreme.sm <- c(head(ord[ord %in% which(obj.sm<0.18)],1),
                  tail(ord[ord %in% which(obj.sm<0.18)],1))
                      

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

pars.lik <-as.matrix(10^rbind(trad.minLik$par,trad.maxLik$par))
pars.lik <- cbind(pars.lik[,1]+pars.lik[,2],pars.lik[,2])
obj.lik <- c(-trad.minLik$value[,1],-trad.maxLik$value[,1])
pred.lik <- c(trad.minLik$value[,2],-trad.maxLik$value[,2])

## Select the model scenarios with the min and max annual runoff within likelihood thresholds from T7
ord <- order(pred.lik)
w.extreme.lik <- c(head(ord[ord %in% which(obj.lik>thres.lik.lo)],1),
                   tail(ord[ord %in% which(obj.lik>thres.lik.hi)],1))

##_ , Maximising fit to threshold - good flood (M8c)
set.seed(1)
trad.fit <- nsga2(function(p) c(-calc.rmse(log10(fdc.pred(p)),log10(fdc.q)),
                                abs(fdc.pred(p,probs=1/365)-7.6)),
                  idim=2,odim=2,
                  lower.bounds=c(0,0),upper.bounds=c(2,1),
                  generations=500,popsize=400,mprob=1/2
                  )

pars.fit <-as.matrix(10^unique(trad.fit$par))
obj.fit <- c(-trad.fit$value[,1])
pars.fit <- cbind(pars.fit[,1]+pars.fit[,2],pars.fit[,2])
pred.fit <- trad.fit$value[,2]

## Select the model scenarios with the maximum and minimum fit
w.extreme.fit <- c(which.min(obj.fit),which.max(obj.fit))


##_ , Plot of all hypothesis testing results (M8a,b,c)
source("fig11.R")
