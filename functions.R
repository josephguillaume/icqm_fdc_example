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
  polygon(x=c(-5,100,-5),y=c(-5,100,100),col="grey90",border=NA)
  #   polygon(x=c(0,seq(7.6,50,length.out=100),100,0),y=c(0,bb(c(seq(7.6,50,length.out=100),100)),0),
  #           border="grey",lwd=2)
  lines(x=seq(7.6,50,length.out=100),y=bb(seq(7.6,50,length.out=100)),col="grey",lwd=2)
  polygon(x=c(7.6,100,100,seq(50,7.6,length.out=100)),
          y=c(bb(7.6),bb(7.6),0,bb(seq(50,7.6,length.out=100))),
          col="grey",density=2,angle=45)
  box()
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
