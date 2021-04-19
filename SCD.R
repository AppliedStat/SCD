#===================================================================
# Weighted median
#===================================================================
#--------------
# new version on Oct. 31. 2019
# NOTE: a = 10L; b = exp(log(10)) ;  c(a,b); a==b
#  tiny=.Machine$double.eps; abs(a-b)<tiny  # good idea, but does not wrok
#  tiny = sqrt(.Machine$double.eps); abs(a-b)<tiny # It works!
#--------------
weighted.median <- function(x, w, interpolation=0.5) { 
  # Preparation
  if (missing(w)) w = rep(1L,length(x))
  if (length(w) != length(x)) stop("'x' and 'w' must have the same length")

  x = as.double(as.vector(x))
  w = as.double(as.vector(w))
  ok= complete.cases(x,w); x=x[ok]; w=w[ok]

  stopifnot(all(w >= 0))
  if(all(w <= 0)) stop("All weights are zero", call.=FALSE)

  orderx = order(x)
  x = x[orderx]
  w = w[orderx] / sum(w)
  Fn = cumsum(w)
  tiny = sqrt(.Machine$double.eps)

  # Main part
  if ( all( abs(Fn-0.5)>tiny ) ) {  # any values of Fn is not 1/2.
      k = sum( Fn < 0.5 )
      return( x[k+1] )
  } else {
    k = which.min ( signif(abs(Fn-0.5),digits=12) ) # Find k with Fn=0.5
    if (w[k+1] < tiny) {   # check if w[k+1] == 0 
        return( x[k+1] )
    } else {
      return( (1-interpolation)*x[k] + interpolation*x[k+1] )
    }
  }
}

#
#===================================================================
# Distance
#===================================================================
total.dist <-
function(theta=c(0,0), x,y,w, norm=c("Lsq", "L1", "L2") ) {
   if (missing(w)) w = rep(1L,length(x))
   norm = match.arg(norm)
   DIST = switch (norm,
             Lsq= sum(w*((x-theta[1])^2 + (y-theta[2])^2)),
             L1 = sum(w*(abs(x-theta[1]) + abs(y-theta[2]))),
             L2 = sum(w*sqrt(((x-theta[1])^2 + (y-theta[2])^2)) )
   )
   return(DIST)
}
#-----------------
optimal.location <-
function(x,y,w, norm=c("Lsq", "L1", "L2"), startL2, maxitL2=100 ) {
   if (missing(w)) w = rep(1L,length(x))
   w = w / sum(w)
   norm = match.arg(norm)
   switch (norm,
      Lsq = return( c(weighted.mean(x,w),weighted.mean(y,w)) ),
      L1  = return( c(weighted.median(x,w),weighted.median(y,w)) )
   )
   total.L2 <- function(theta,x,y,w) total.dist(theta,x,y,w,norm="L2")
   if (missing(startL2)) {
       startL2 = c(weighted.median(x,w),weighted.median(y,w))
   }
   OPT = optim(par=startL2, fn=total.L2, x=x,y=y,w=w,
      method="L-BFGS-B", lower=c(min(x),min(y)), upper=c(max(x),max(y)),
      control=list(maxit=maxitL2) )
   return(OPT$par)
}
#
#===================================================================
# Beta Distribution
#===================================================================
rbeta4 = function(n, shape1, shape2, min=0,max=1, ncp = 0) {
  span = max-min 
  rbeta(n, shape1=shape1, shape2=shape2, ncp=ncp)*span + min 
}
#
dbeta4 = function(x, shape1, shape2, min=0,max=1, ncp = 0, log = FALSE) {
  span = max-min 
  if (log) { 
     dbeta((x-min)/span, shape1=shape1, shape2=shape2, ncp=ncp, log=log) - log(span)
  } else {
     dbeta((x-min)/span, shape1=shape1, shape2=shape2, ncp=ncp, log=log) / span
  }
}
#
pbeta4 = function (q, shape1, shape2, min=0,max=1, ncp = 0, lower.tail = TRUE, log.p = FALSE) {
  span = max-min
  pbeta( (q-min)/span, shape1=shape1, shape2=shape2, ncp=ncp, lower.tail=lower.tail, log.p=log.p)
}
#
qbeta4 = function (p, shape1, shape2, min=0,max=1, ncp = 0, lower.tail = TRUE, log.p = FALSE) {
  span = max-min
  qbeta( p, shape1=shape1, shape2=shape2, ncp=ncp, lower.tail=lower.tail, log.p=log.p)*span + min
}
#


#===================================================================
# Empirical Breakdown 
#===================================================================
ebreakdown = function(n, w) {
  if (missing(w)) w = rep(1L,n)
  BIG = 1.0E90 * n
  MID = 1.0E10
  TINY = 1.0E-10 / n
  K = ceiling(n/2)

  SUM.CG <- SUM.L1 <- SUM.L2 <- numeric(K)
  w0 = sort(w)

  # Min Breakdown
  w = rev(w0)
  xx = runif(n, -TINY, TINY )
  yy = runif(n, -TINY, TINY )
  for ( i in 1L:K ) {
      xx[1L:i] = BIG
      yy[1L:i] = BIG
      SUM.CG[i] = sum( abs(optimal.location(xx,yy, w=w, norm="Lsq")) )
      SUM.L1[i] = sum( abs(optimal.location(xx,yy, w=w, norm="L1")) )
      SUM.L2[i] = sum( abs(optimal.location(xx,yy, w=w, norm="L2")) )
  }
  BD.CG = sum(SUM.CG < MID)/n*100
  BD.L1 = sum(SUM.L1 < MID)/n*100
  BD.L2 = sum(SUM.L2 < MID)/n*100
  BREAK = c(BD.CG, BD.L1, BD.L2)
  OUT2 = cbind(SUM.CG, SUM.L1, SUM.L2)
  colnames(OUT2) = c("CG", "L1", "L2")
  rownames(OUT2) = c(1L:K)
  cat("------------------------------------------------\n")
  cat(" n =",n, " weights =", w, "\n")
  prmatrix(OUT2)

  cat("\n* Min. Breakdown *\n")
  names(BREAK) = c("CG(%)", "L1(%)", "L2(%)")
  print( round(BREAK,3) )
  cat("\n================================================\n")

  # Max Breakdown
  w = w0
  xx = runif(n, -TINY, TINY )
  yy = runif(n, -TINY, TINY )
  for ( i in 1L:K ) {
      xx[1L:i] = BIG
      yy[1L:i] = BIG
      SUM.CG[i] = sum( abs(optimal.location(xx,yy, w=w, norm="Lsq")) )
      SUM.L1[i] = sum( abs(optimal.location(xx,yy, w=w, norm="L1")) )
      SUM.L2[i] = sum( abs(optimal.location(xx,yy, w=w, norm="L2")) )
  }
  BD.CG = sum(SUM.CG < MID)/n*100
  BD.L1 = sum(SUM.L1 < MID)/n*100
  BD.L2 = sum(SUM.L2 < MID)/n*100
  BREAK = c(BD.CG, BD.L1, BD.L2)
  OUT1 = cbind(SUM.CG, SUM.L1, SUM.L2)
  colnames(OUT1) = c("CG", "L1", "L2")
  rownames(OUT1) = c(1L:K)
  cat("\n================================================\n")
  cat(" n =",n, " weights =", w, "\n")
  prmatrix(OUT1)
  cat("\n* Max. Breakdown *\n")
  names(BREAK) = c("CG(%)", "L1(%)", "L2(%)")
  print( round(BREAK,3) )
}
##################################


