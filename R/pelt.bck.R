## FUNCTIONS FOR RUNNING PELT - INCLUDING UNIVARIATE COST FUNCTIONS ##


#' Detecting Multiple Changepoints using the PELT Algorithm
#'
#' Calculates the optimal positioning and number of changepoints using PELT.
#'
#' This method uses the PELT algorithm to obtain the optimal number and location of changepoints within a univariate time series. This is done via the minimisation of a penalised cost function using dynamic programming. Inequality-based pruning is used to reduce computation whilst retaining optimality. A range of different cost functions and penalty values can be used.
#'
#' @param data A vector of length \code{n} containing the data within which you wish to find a changepoint.
#' @param pen Numeric value of the linear penalty function.  This value is used in the decision for each individual changepoint so that in total the penalty is k*pen where k is the optimal number of changepoints detected.
#A non-negative numeric value representing the penalty incurred to the cost whenever an additional changepoint in introduced in the model.
#' @param min.dist The minimum distance allowed between any two changepoints. Required to have an integer value of at least 1 for changes in mean, or at least 2 for changes in variance.
#' @param cost The function used to calculate the cost of a given segment of data. The choice of this function dictates the assumed distribution and the type(s) of changes being detected (or simply the generic cost of data if non-parametric). See Details for possible choices.
#' @param sum.stat The function used to generate the summary statistics used by \code{cost}. See Details for possible choices.
#' @param initial.likelihood.value The initial value of the likelihood/cost function.
#'
#' @return The vector of changepoint locations detected by PELT.
#' @export
#'
#' @examples
#' # Normal observations, multiple change in mean.
#' data = c( rnorm(100, mean=0, sd=1), rnorm(100, mean=5, sd=1), rnorm(100, mean=-1, sd=1), rnorm(100, mean=2, sd=1) )
#' # plot.ts(data)
#' n = length(data)
#' pelt.results = pelt(data=data, pen=2*log(n), cost=pelt.norm.mean.cost, sum.stat=pelt.norm.sum)
pelt = function(data, pen = 2*log(length(data)), min.dist=2, cost=pelt.norm.meanvar.cost, sum.stat=pelt.norm.sum, initial.likelihood.value=0)
{

  if(is.matrix(data)){ # data is multivariate
    n = nrow(data)
  } else{ # data is univariate
    n <- length(data)
  }

  d = min.dist # minimum distance between changepoints
  PRUNE = T

  #lastchangelike = array(NA,dim = n+1) #array(0,dim = n+1)
  #lastchangecpts = array(NA,dim = n+1) #array(0,dim = n+1)
  #numchangecpts = array(NA,dim = n+1) #array(0,dim = n+1)

  lastchangelike = array(0,dim = n+1)
  lastchangecpts = array(0,dim = n+1)
  numchangecpts = array(0,dim = n+1)

  lastchangelike[0 +1] <- c(-pen + initial.likelihood.value)
  numchangecpts[0 +1] <- 0
  lastchangecpts[0 +1] <- 0

  sumx <- sum.stat(data)

  if(d>1){

    lastchangelike[1:(d-1) +1] <- NA
    numchangecpts[1:(d-1) +1] <- NA
    lastchangecpts[1:(d-1) +1] <- NA

    for(tau in d:(2*d-1)){
      lastchangelike[tau +1] = cost(tau, 0, sumx) + pen # calculate likelihood from 0 to the d:(2*d-1) points
    }
    numchangecpts[d:(2*d-1) +1] = 0
    lastchangecpts[d:(2*d-1) +1] = 0

  } else{
    lastchangelike[1 +1] <- 0
    numchangecpts[1 +1] <- 0
    lastchangecpts[1 +1] <- 0
  }

  checklist <- array() #stores the candidate changepoint positions
  checklist[1:2] <- c(0, d) #0

  for (tstar in (2*d):n) {
    if(lastchangelike[tstar+1] == 0){

      #if(tstar==7 && all(checklist==c(0,4,5))) debug(cost)
      tmplike <- lastchangelike[checklist+1] + cost(tstar, checklist, sumx) + pen
      if(any(is.na(tmplike))) browser()

      #### Store changepoints and cost function for each tstar ###
      lastchangecpts[tstar+1] <- checklist[min(which(tmplike == min(tmplike[!is.na(tmplike)])))]
      lastchangelike[tstar+1] <- min(tmplike[!is.na(tmplike)])
      numchangecpts[tstar +1] <- numchangecpts[lastchangecpts[tstar+1]+1]+1

      if(PRUNE){
        checklist <- checklist[(tmplike - pen) <= lastchangelike[tstar+1]]
      }
    }
    checklist <- c(checklist,tstar-d+1) #c(checklist,tstar-1)
  }
  cp = n
  lastchangecpts2 <- lastchangecpts[-1]
  while(cp[1]>0){
    cp=c(lastchangecpts2[cp[1]],cp)
  }

  #return.list = list(lastchangecpts, cp[-c(1,length(cp))], lastchangelike, numchangecpts)
  #names(return.list) = list("lastchangecpts", "cpts", "lastchangelike", "numchangecpts")
  #return(return.list)

  return(cp[-c(1,length(cp))])

}



#' Calculate Summary Statistics used for Normal likelihood calculation within PELT.
#'
#' Calculates the necessary summary statistics which are to be used as part of the calculation of the Normal likelihood for the PELT algorithm. Not typically used outside of this purpose.
#'
#' An internal function, not designed for use by the end-user.
#'
#' @param data A vector of time-ordered observations to be analysed using PELT.
#'
#' @return A \code{3} x \code{length(data)} matrix. Row 1 contains the cumulative sums of the observations, row 2 contains the cumulative sums of the squared observations, and row 3 contains the cumulative sums of the squared de-meaned observations.
#' @export
#'
#' @examples
#' data = rnorm(100, 10, 3)
#' pelt.norm.sum(data)
pelt.norm.sum = function(data){
  return(matrix(data = c(cumsum(c(0,data)),cumsum(c(0,data^2)), cumsum(c(0,(data-mean(data))^2))),nrow = 3,byrow = T))
}



#' Calculate Normal likelihood of data segment, assuming variable mean and fixed variance.
#'
#' Calculates the Normal likelihood of multiple segments of data, as defined by the vector of possible start-points (\code{tau}) and a single common end-point (\code{R}). This likelihood calculation assumes that the variance is fixed as \code{1}, and sets the mean equal to its maximum likelihood estimate.
#'
#' @param tau A vector of locations indicating the possible beginnings of a segment.
#' @param R A single location representing the end of the possible segment
#' @param sumx The summary statistics for the entire time series being analysed.
#'
#' @return A vector containing the cost of the different possible segments defined by \code{tau} and \code{R}.
#' @export
#'
#' @examples
#' data = rnorm(100, 10, 3)
#' sumx = pelt.norm.sum(data)
#' tau = c(30,40,50) # start-points of possible segments
#' R = 70 # end-point of possible segments
#' pelt.norm.mean.cost(tau,R,sumx) # costs for each segment
pelt.norm.mean.cost = function(tau,R,sumx){

  cost = (tau-R)*log(2*pi) + (sumx[2,tau +1] - sumx[2,R+1])-(tau-R) *( (sumx[1,tau +1] - sumx[1,R+1])/(tau - R) )^2
  return(cost)

}


#' Calculate Normal likelihood of data segment with segment length penalty, assuming variable mean and fixed variance.
#'
#' Calculates the Normal likelihood of multiple segments of data, as defined by the vector of possible start-points (\code{tau}) and a single common end-point (\code{R}). This likelihood is penalised by the segment length. The calculation assumes that the variance is fixed as \code{1}, and sets the mean equal to its maximum likelihood estimate.
#'
#' @param tau A vector of locations indicating the possible beginnings of a segment.
#' @param R A single location representing the end of the possible segment
#' @param sumx The summary statistics for the entire time series being analysed.
#'
#' @return A vector containing the cost of the different possible segments defined by \code{tau} and \code{R}.
#' @export
#'
#' @examples
#' data = rnorm(100, 10, 3)
#' sumx = pelt.norm.sum(data)
#' tau = c(30,40,50) # start-points of possible segments
#' R = 70 # end-point of possible segments
#' pelt.norm.mean.cost.seglen(tau,R,sumx) # costs for each segment
pelt.norm.mean.cost.seglen = function(tau,R,sumx){

  cost = (tau-R)*log(2*pi) + (sumx[2,tau +1] - sumx[2,R+1])-(tau-R) *( (sumx[1,tau +1] - sumx[1,R+1])/(tau - R) )^2 + log(tau-R) # log(n), segment length penalty.
  return(cost)

}


#' Calculate Normal likelihood of data segment, assuming fixed mean and variable variance.
#'
#' Calculates the Normal likelihood of multiple segments of data, as defined by the vector of possible start-points (\code{tau}) and a single common end-point (\code{R}). This likelihood calculation assumes that the mean is fixed as \code{1}, and sets the variance equal to its maximum likelihood estimate.
#'
#' @param tau A vector of locations indicating the possible beginnings of a segment.
#' @param R A single location representing the end of the possible segment
#' @param sumx The summary statistics for the entire time series being analysed.
#'
#' @return A vector containing the cost of the different possible segments defined by \code{tau} and \code{R}.
#' @export
#'
#' @examples
#' data = rnorm(100, 10, 3)
#' sumx = pelt.norm.sum(data)
#' tau = c(30,40,50) # start-points of possible segments
#' R = 70 # end-point of possible segments
#' pelt.norm.var.cost(tau,R,sumx) # costs for each segment
pelt.norm.var.cost = function(tau,R,sumx){

  mll.var <- function(x,n){
    neg = x<=0
    x[neg==TRUE] = 0.00000000001 #10000
    return(n*(log(2*pi) + log(x/n) + 1))
  }
  cost <- mll.var(sumx[3,tau+1] - sumx[3,R+1], tau - R)
  return(cost)

}


#' Calculate Normal likelihood of data segment with segment length penalty, assuming fixed mean and variable variance.
#'
#' Calculates the Normal likelihood of multiple segments of data, as defined by the vector of possible start-points (\code{tau}) and a single common end-point (\code{R}). This likelihood is penalised by the segment length. The calculation assumes that the variance is fixed as \code{1}, and sets the mean equal to its maximum likelihood estimate.
#'
#' @param tau A vector of locations indicating the possible beginnings of a segment.
#' @param R A single location representing the end of the possible segment
#' @param sumx The summary statistics for the entire time series being analysed.
#'
#' @return A vector containing the cost of the different possible segments defined by \code{tau} and \code{R}.
#' @export
#'
#' @examples
#' data = rnorm(100, 10, 3)
#' sumx = pelt.norm.sum(data)
#' tau = c(30,40,50) # start-points of possible segments
#' R = 70 # end-point of possible segments
#' pelt.norm.var.cost.seglen(tau,R,sumx) # costs for each segment
pelt.norm.var.cost.seglen = function(tau,R,sumx){

  mll.var <- function(x,n){
    neg = x<=0
    x[neg==TRUE] = 10000
    return(n*(log(2*pi) + log(x/n) + 1))
  }
  cost <- mll.var(sumx[3,tau+1] - sumx[3,R+1], tau - R) + log(tau-R) # segment length penalty.
  return(cost)

}


#' Calculate Normal likelihood of data segment, assuming variable mean and variable variance.
#'
#' Calculates the Normal likelihood of multiple segments of data, as defined by the vector of possible start-points (\code{tau}) and a single common end-point (\code{R}). This calculation sets both the mean and variance equal to their maximum likelihood estimates.
#'
#' @param tau A vector of locations indicating the possible beginnings of a segment.
#' @param R A single location representing the end of the possible segment
#' @param sumx The summary statistics for the entire time series being analysed.
#'
#' @return A vector containing the cost of the different possible segments defined by \code{tau} and \code{R}.
#' @export
#'
#' @examples
#' data = rnorm(100, 10, 3)
#' sumx = pelt.norm.sum(data)
#' tau = c(30,40,50) # start-points of possible segments
#' R = 70 # end-point of possible segments
#' pelt.norm.meanvar.cost(tau,R,sumx) # costs for each segment
pelt.norm.meanvar.cost = function(tau,R,sumx){

  cost <- array()
  cost[which((tau-R) <= 1)] = 0
  x <- R[which((tau-R) > 1)]
  #cost[which((tau-R) > 1)]= (tau-x) * (log(2*pi) + log((((sumx[2,tau +1] - sumx[2,x+1])-(tau-x) *((sumx[1,tau +1] - sumx[1,x+1])/(tau - (x)))^2))/(tau-x)) + 1)

  nn = tau-x
  sigmasq.mle = (1/nn) * ( (sumx[2,tau +1] - sumx[2,x+1]) - nn*((sumx[1,tau +1] - sumx[1,x+1])/nn)^2 )
  neg = sigmasq.mle <= 0
  sigmasq.mle[neg==TRUE] = 0.00000000001
  cost[which((tau-R) > 1)] = nn * (log(2*pi) + log(sigmasq.mle) + 1)

  return(cost)

}


#' Calculate Normal likelihood of data segment with segment length penalty, assuming variable mean and variable variance.
#'
#' Calculates the Normal likelihood of multiple segments of data, as defined by the vector of possible start-points (\code{tau}) and a single common end-point (\code{R}). This likelihood is penalised by the segment length. The calculation sets both the mean and variance equal to their maximum likelihood estimates.
#'
#' @param tau A vector of locations indicating the possible beginnings of a segment.
#' @param R A single location representing the end of the possible segment
#' @param sumx The summary statistics for the entire time series being analysed.
#'
#' @return A vector containing the cost of the different possible segments defined by \code{tau} and \code{R}.
#' @export
#'
#' @examples
#' data = rnorm(100, 10, 3)
#' sumx = pelt.norm.sum(data)
#' tau = c(30,40,50) # start-points of possible segments
#' R = 70 # end-point of possible segments
#' pelt.norm.meanvar.cost.seglen(tau,R,sumx) # costs for each segment
pelt.norm.meanvar.cost.seglen = function(tau,R,sumx){

  cost <- array()
  cost[which((tau-R) <= 1)] = 0
  x <- R[which((tau-R) > 1)]
  #cost[which((tau-R) > 1)]= (tau-x) * (log(2*pi) + log((((sumx[2,tau +1] - sumx[2,x+1])-(tau-x) *((sumx[1,tau +1] - sumx[1,x+1])/(tau - (x)))^2))/(tau-x)) + 1)

  nn = tau-x
  sigmasq.mle = (1/nn) * ( (sumx[2,tau +1] - sumx[2,x+1]) - nn*((sumx[1,tau +1] - sumx[1,x+1])/nn)^2 )
  neg = sigmasq.mle <= 0
  sigmasq.mle[neg==TRUE] = 0.00000000001
  cost[which((tau-R) > 1)] = nn * (log(2*pi) + log(sigmasq.mle) + 1)

  cost[which((tau-R) > 1)] = cost[which((tau-R) > 1)] + log(nn) # segment length penalty.
  #cost[which((tau-R) > 1)] = cost[which((tau-R) > 1)] + log(tau-x) # segment length penalty.
  return(cost)

}


#' Calculate Summary Statistics used for Gamma likelihood calculation within PELT.
#'
#' Calculates the necessary summary statistics which are to be used as part of the calculation of the Gamma likelihood for the PELT algorithm. Not typically used outside of this purpose.
#'
#' An internal function, not designed for use by the end-user.
#'
#' @param data A vector of time-ordered observations to be analysed using PELT.
#'
#' @return A vector containing the cumulative sums of the observations.
#' @export
#'
#' @examples
#' data = rgamma(n=100, shape=1, scale=1/2)
#' pelt.gamma.sum(data)
pelt.gamma.sum = function(data){

  return( c(0,cumsum(data)) )

}


#' Calculate Gamma likelihood of data segment, assuming fixed shape parameter and variable scale parameter.
#'
#' Calculates the Gamma likelihood of multiple segments of data, as defined by the vector of possible start-points (\code{tau}) and a single common end-point (\code{R}). This calculation assumes the shape is fixed as some specified value and sets the scale equal to its maximum likelihood estimate.
#'
#' @param tau A vector of locations indicating the possible beginnings of a segment.
#' @param R A single location representing the end of the possible segment
#' @param sumx The summary statistics for the entire time series being analysed.
#' @param shape The value of the shape parameter for the Gamma distribution.
#'
#' @return A vector containing the cost of the different possible segments defined by \code{tau} and \code{R}.
#' @export
#'
#' @examples
#' shape = 1 # the shape parameter of the data
#' data = rgamma(n=100, shape=shape, scale=1/2)
#' sumx = pelt.gamma.sum(data)
#' tau = c(30,40,50) # start-points of possible segments
#' R = 70 # end-point of possible segments
#' pelt.gamma.cost(tau, R, sumx, shape=shape) # costs for each segment
pelt.gamma.cost = function(tau, R, sumx, shape=1){

  N = tau - R

  cost = array()
  cost[which(N <= 1)] = 0
  r <- R[which(N > 1)]
  x = sumx[tau+1] - sumx[r+1]
  n = tau - r

  cost[which(N > 1)] = 2*n*shape*log(x) - 2*n*shape*log(n*shape)

  return(cost)

}


#' Title
#'
#' @param tau
#' @param R
#' @param sumx
#' @param shape
#'
#' @return A vector containing the cost of the different possible segments defined by \code{tau} and \code{R}.
#' @export
#'
pelt.gamma.cost.seglen = function(tau, R, sumx, shape=1){

  N = tau - R

  cost = array()
  cost[which(N <= 1)] = 0
  r <- R[which(N > 1)]
  x = sumx[tau+1] - sumx[r+1]
  n = tau - r

  cost[which(N > 1)] = 2*n*shape*log(x) - 2*n*shape*log(n*shape) + log(n)

  return(cost)

}


#' Title
#'
#' @param tau
#' @param R
#' @param sumx
#'
#' @return A vector containing the cost of the different possible segments defined by \code{tau} and \code{R}.
#' @export
#'
pelt.exp.cost = function(tau, R, sumx){

  N = tau - R

  cost = array()
  cost[which(N <= 1)] = 0
  r <- R[which(N > 1)]
  x = sumx[tau+1] - sumx[r+1]
  n = tau - r

  cost[which(N > 1)] = 2*n*(1 - log(n/x))

  return(cost)

}


#' Title
#'
#' @param tau
#' @param R
#' @param sumx
#'
#' @return A vector containing the cost of the different possible segments defined by \code{tau} and \code{R}.
#' @export
#'
pelt.exp.cost.seglen = function(tau, R, sumx){

  N = tau - R

  cost = array()
  cost[which(N <= 1)] = 0
  r <- R[which(N > 1)]
  x = sumx[tau+1] - sumx[r+1]
  n = tau - r

  cost[which(N > 1)] = 2*n*(1 - log(n/x)) + log(n)

  return(cost)

}


# DON'T NEED ANY OF THE COST FUNCTIONS BELOW FOR A-SMOP, SINCE THEY ARE FOR MULTIVARIATE NORMAL DATA WITHIN PELT.
#
# pelt.multivar.norm.sum = function(data){
#
#   if(is.null(dim(data))){
#     data = as.matrix(data)
#   }
#
#   n = nrow(data)
#   p = ncol(data)
#
#   y = rbind(rep(0,p), apply(data,2,cumsum))
#   y2 = rbind(rep(0,p), apply(data^2,2,cumsum))
#   y2.zeromean = rbind(rep(0,p), apply(data,2, function(x) cumsum((x-mean(x))^2) ))
#
#   return( cbind(y,y2,y2.zeromean) )
#
# }
#
# pelt.multivar.norm.mean.cost = function(tau,R,sumx){
#
#   n = tau - R
#   p = ncol(sumx)/3
#
#   #cost = (tau-R)*log(2*pi) + (sumx[2,tau +1] - sumx[2,R+1])-(tau-R) *( (sumx[1,tau +1] - sumx[1,R+1])/(tau - R) )^2
#   cost.mv = n*log(2*pi) + (sumx[rep(tau+1, length(R)), (p+1):(2*p) ,drop=F] - sumx[R+1, (p+1):(2*p) ,drop=F])-n *( (sumx[rep(tau+1, length(R)), 1:p ,drop=F] - sumx[R+1, 1:p ,drop=F])/n )^2
#   cost = apply( cost.mv, 1, sum )
#
#   return(cost)
#
# }
#
# pelt.multivar.norm.mean.cost.mbic = function(tau,R,sumx){
#
#   n = tau - R
#   p = ncol(sumx)/3
#
#   #cost = (tau-R)*log(2*pi) + (sumx[2,tau +1] - sumx[2,R+1])-(tau-R) *( (sumx[1,tau +1] - sumx[1,R+1])/(tau - R) )^2
#   cost.mv = n*log(2*pi) + (sumx[rep(tau+1, length(R)), (p+1):(2*p) ,drop=F] - sumx[R+1, (p+1):(2*p) ,drop=F])-n *( (sumx[rep(tau+1, length(R)), 1:p ,drop=F] - sumx[R+1, 1:p ,drop=F])/n )^2 + log(n)
#   cost = apply( cost.mv, 1, sum )
#
#   return(cost)
#
# }
#
# pelt.multivar.norm.var.cost = function(tau, R, sumx){
#
#   # tau - single changepoint location
#   # R - checklist of multiple possible previous changepoint locations
#   # sumx - summary statistics for the whole time series
#   # p - number of variables
#
#   #mll.var <- function(x,n){
#   #  neg = x<=0
#   #  x[neg==TRUE] = 0.00000000001 #10000
#   #  return(n*(log(2*pi) + log(x/n) + 1))
#   #}
#
#   n = tau - R
#   p = ncol(sumx)/3
#
#   y = matrix( rep( sumx[tau+1,(2*p+1):(3*p)], length(R) ), nrow=length(R) , byrow=T)
#   x = y - sumx[R+1,(2*p+1):(3*p) ,drop=F]
#
#   neg = x<=0
#   x[neg==TRUE] = 0.00000000001 #10000
#   cost.mv = n*(log(2*pi) + log(x/n) + 1)
#   cost = apply( cost.mv, 1, sum )
#
#   #cost = apply( z, 2, mll.var, n=tau-R)
#
#   return(cost)
#
# }
#
# pelt.multivar.norm.var.cost.mbic = function(tau, R, sumx){
#
#   # tau - single changepoint location
#   # R - checklist of multiple possible previous changepoint locations
#   # sumx - summary statistics for the whole time series
#   # p - number of variables
#
#   #mll.var <- function(x,n){
#   #  neg = x<=0
#   #  x[neg==TRUE] = 0.00000000001 #10000
#   #  return(n*(log(2*pi) + log(x/n) + 1))
#   #}
#
#   n = tau - R
#   p = ncol(sumx)/3
#
#   y = matrix( rep( sumx[tau+1,(2*p+1):(3*p)], length(R) ), nrow=length(R) , byrow=T)
#   x = y - sumx[R+1,(2*p+1):(3*p) ,drop=F]
#
#   neg = x<=0
#   x[neg==TRUE] = 0.00000000001 #10000
#   cost.mv = n*(log(2*pi) + log(x/n) + 1) + log(n)
#   cost = apply( cost.mv, 1, sum )
#
#   #cost = apply( z, 2, mll.var, n=tau-R)
#
#   return(cost)
#
# }
