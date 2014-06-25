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

png("fig06 observed fdc.png",width=3.54,height=3,units="in",res=70)
##win.metafile("fig06 observed fdc.wmf",width=3.54,height=3)
par(cex=0.8,mar=c(5.1,4.1,1.1,1.1))
plot(xs,fdc.q,log="y",
     ylab="Flow (mm)",
     xlab="qnorm(Exceedance probability)")
abline(mm)
points(qnorm(1/(100*365)),quantile(Q,1-1/(100*365)),pch="+")
points(qnorm(1/(1/12*365)),quantile(Q,1-1/(1/12*365)),pch="+")
abline(v=qnorm(1/(c(100,1/12)*365)),col="grey50",lty="dashed")
dev.off()

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

png("fig08 pomore.png",width=3.54,height=3,units="in",res=70)
##win.metafile("fig08 pomore.wmf",width=3.54,height=3)
par(cex=0.7)
par(mar=c(5.1,4.1,1.1,1.1))
boxplot(cbind(abs(pom.good[,1]-start.pos[1])/start.pos[1],
           abs(pom.good[,2]-start.pos[2])/start.pos[2]),
        names=c(expression(q[1]),expression(q[2])),
        xlab="Parameter",ylab="Normalised distance"
        )
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
source("fig07.R")

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

## Likelihood threshold for hyp testing
thres.lik.lo <- min(objseq[which(keep & pareto.front)[order(ys[keep & pareto.front],decreasing=F)][1]]) ## Lik of lowest y value on pareto front
thres.lik.hi <- min(objseq[which(keep & pareto.front)[order(ys[keep & pareto.front],decreasing=T)][1]]) ## Lik of highest y value on pareto front

psets2 <- cbind(rowSums(10^psets),10^psets[,2])

# hard-coded result from M8b because they are mutually dependent
## Warning: result from any single run may differ
lik.ht.pars <- matrix(c(8.44598275902987, 21.4774398003501, 2.92690003299676, 4.99546011684257),ncol=2)
lik.ht.pars.trans <- cbind(log10(lik.ht.pars[,1]- lik.ht.pars[,2]),
                        log10(lik.ht.pars[,2]))
lik.ht.pars.ys <- apply(lik.ht.pars.trans,1,fdc.pred,probs=1/365)


##_ , Plot of density, confidence interval/critical region and corresponding models
source("fig11.R")

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
source("fig12.R")
