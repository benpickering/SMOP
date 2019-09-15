## GENERIC MULTIVARIATE COST / CALCULATION FUNCTIONS ##


#' Intra-Variable Penalised Cost Calculation
#'
#' Calculates the intra-variable penalised cost (i.e. using only the \code{alpha} penalty) for all variables for all possible previous changepoint vectors.
#'
#' An internal function, not designed for use by the end-user.
#'
#' @param n A matrix containing the lengths between the changepoints in current changepoint vector being considered and the changepoints in the check-list of previous changepoint vectors.
#' @param n.nonzero A matrix containing the logical values indicating which of the corresponding values in \code{n} are non-zero.
#' @param x A \code{p} x \code{length.checklist} matrix containing the cumulative sums of the observations between the current changepoint vector and each corresponding changepoint vector on the check-list.
#' @param x2 The same as \code{x}, except containing the cumulative sums of the squared observations.
#' @param p The number of variables in the series.
#' @param alpha The variable-specific penalty, used to penalise the addition of a given changepoint into a given variable. A non-negative numeric value.
#' @param length.checklist The number of possible previous changepoint vectors on the check-list being considered.
#' @param cost.func The function being used to calculate the multivariate cost. Required to be a function (not a string name).
#'
#' @return A vector of length \code{length.checklist} containing the intra-variable penalised costs for each changepoint vector on the check-list.
calculate.alpha.cost = function(n, n.nonzero, x, x2, p, alpha, length.checklist, cost.func=norm.meanvar)
{

  # IMPORTANT NOTE: the dimensions here are switched from the usual convention, this makes the computation easier due to R's element-wise handling of matrices and vectors.
  # x and x2 are matrices: p rows, each column represents a different cpt vector.

  # n is a matrix containing the differences between the changepoints in current.cpt.vec and the cpts in the checklist matrix
  # n.nonzero is equivalent to the vector version of: rowSums(c.taustar!=c)

  unpenalised.cost = cost.func(n=n, x=x, x2=x2, p=p, length.checklist=length.checklist)

  return( unpenalised.cost + (n.nonzero)*alpha )

}



#' Multivariate Cost Calculation: Normally distributed observations with fixed variance, detecting changes in mean only.
#'
#' Calculates the cost between a sequence of multivariate observations using the provided cumulative sums of data. The observations are assumed to be Normally distributed with fixed variance (\code{= 1}). The mean parameters are set to their maximum likelihood estimates. Used for detecting changes in mean.
#' An internal function, not designed for use by the end-user.
#'
#' @param n A matrix containing the lengths between the changepoints in current changepoint vector being considered and the changepoints in the check-list of previous changepoint vectors.
#' @param x A \code{p} x \code{length.checklist} matrix containing the cumulative sums of the observations between the current changepoint vector and each corresponding changepoint vector on the check-list.
#' @param x2 The same as \code{x}, except containing the cumulative sums of the squared observations.
#' @param p The number of variables in the series.
#' @param length.checklist The number of possible previous changepoint vectors on the check-list being considered.
#'
#' @return A vector of length \code{length.checklist} containing the multivariate costs for each changepoint vector on the check-list being considered.
norm.mean = function(n, x, x2, p, length.checklist)
{

  # x -- differences in cumulative sums.
  # x2 - differences in squared cumulative sums.

  #sigmasq = 1 # variance is assumed 1 for now.
  #cost = n*log(2*pi) + n*log(sigmasq) + (1/sigmasq)*(x2 - (x^2)/n)

  cost = n*log(2*pi) + (x2 - (x^2)/n)

  return( .colSums(cost, na.rm=T, m=p, n=length.checklist ) )

}



#' Multivariate Cost Calculation: Normally distributed observations with fixed variance, detecting changes in mean only, with additional segment length penalty in likelihood.
#'
#' Calculates the cost between a sequence of multivariate observations using the provided cumulative sums of data, with the additional penalty of the segment length (\code{log(n)}). The observations are assumed to be Normally distributed with fixed variance (\code{= 1}). The mean parameters are set to their maximum likelihood estimates. Used for detecting changes in mean only.
#'
#' An internal function, not designed for use by the end-user.
#'
#' @param n A matrix containing the lengths between the changepoints in current changepoint vector being considered and the changepoints in the check-list of previous changepoint vectors.
#' @param x A \code{p} x \code{length.checklist} matrix containing the cumulative sums of the observations between the current changepoint vector and each corresponding changepoint vector on the check-list.
#' @param x2 The same as \code{x}, except containing the cumulative sums of the squared observations.
#' @param p The number of variables in the series.
#' @param length.checklist The number of possible previous changepoint vectors on the check-list being considered.
#'
#' @return A vector of length \code{length.checklist} containing the multivariate costs for each changepoint vector on the check-list being considered.
norm.mean.seglen = function(n, x, x2, p, length.checklist)
{

  # x -- differences in cumulative sums.
  # x2 - differences in squared cumulative sums.

  cost = n*log(2*pi) + (x2 - (x^2)/n) + log(n) # log(n) is the segment length part of the penalty
  return( .colSums(cost, na.rm=T, m=p, n=length.checklist ) )

}



#' Multivariate Cost Calculation: Normally distributed observations with fixed mean, detecting changes in variance only.
#'
#' Calculates the cost between a sequence of multivariate observations using the provided cumulative sums of data. The observations are assumed to be Normally distributed with fixed mean (\code{= 0}). The variance parameters are set to their maximum likelihood estimates. Used for detecting changes in variance only.
#'
#' An internal function, not designed for use by the end-user.
#'
#' @param n A matrix containing the lengths between the changepoints in current changepoint vector being considered and the changepoints in the check-list of previous changepoint vectors.
#' @param x A \code{p} x \code{length.checklist} matrix containing the cumulative sums of the observations between the current changepoint vector and each corresponding changepoint vector on the check-list.
#' @param x2 The same as \code{x}, except containing the cumulative sums of the squared observations.
#' @param p The number of variables in the series.
#' @param length.checklist The number of possible previous changepoint vectors on the check-list being considered.
#'
#' @return A vector of length \code{length.checklist} containing the multivariate costs for each changepoint vector on the check-list being considered.
norm.var = function(n, x, x2, p, length.checklist)
{

  X = x2 - (x^2)/n
  neg = X<=0
  X[neg==TRUE] = 0.0000000000001
  cost = n*(log(2*pi) + log(X/n) + 1)

  return( .colSums(cost, na.rm=T, m=p, n=length.checklist ) )

}



#' Multivariate Cost Calculation: Normally distributed observations with fixed mean, detecting changes in variance only, with additional segment length penalty in likelihood.
#'
#' Calculates the cost between a sequence of multivariate observations using the provided cumulative sums of data, with the additional penalty of the segment length (\code{log(n)}). The observations are assumed to be Normally distributed with fixed mean (\code{= 0}). The variance parameters are set to their maximum likelihood estimates. Used for detecting changes in variance only.
#'
#' An internal function, not designed for use by the end-user.
#'
#' @param n A matrix containing the lengths between the changepoints in current changepoint vector being considered and the changepoints in the check-list of previous changepoint vectors.
#' @param x A \code{p} x \code{length.checklist} matrix containing the cumulative sums of the observations between the current changepoint vector and each corresponding changepoint vector on the check-list.
#' @param x2 The same as \code{x}, except containing the cumulative sums of the squared observations.
#' @param p The number of variables in the series.
#' @param length.checklist The number of possible previous changepoint vectors on the check-list being considered.
#'
#' @return A vector of length \code{length.checklist} containing the multivariate costs for each changepoint vector on the check-list being considered.
norm.var.seglen = function(n, x, x2, p, length.checklist)
{

  X = x2 - (x^2)/n
  neg = X<= 0
  X[neg==TRUE] = 0.000000000001
  cost = n*(log(2*pi) + log(X/n) + 1) + log(n)

  return( .colSums(cost, na.rm=T, m=p, n=length.checklist ) )

}



#' Multivariate Cost Calculation: Normally distributed observations, detecting changes in both mean and variance.
#'
#' Calculates the cost between a sequence of multivariate observations using the provided cumulative sums of data. The observations are assumed to be Normally distributed. The mean and variance parameters are set to their maximum likelihood estimates. Used for detecting changes in both mean and variance.
#'
#' An internal function, not designed for use by the end-user.
#'
#' @param n A matrix containing the lengths between the changepoints in current changepoint vector being considered and the changepoints in the check-list of previous changepoint vectors.
#' @param x A \code{p} x \code{length.checklist} matrix containing the cumulative sums of the observations between the current changepoint vector and each corresponding changepoint vector on the check-list.
#' @param x2 The same as \code{x}, except containing the cumulative sums of the squared observations.
#' @param p The number of variables in the series.
#' @param length.checklist The number of possible previous changepoint vectors on the check-list being considered.
#'
#' @return A vector of length \code{length.checklist} containing the multivariate costs for each changepoint vector on the check-list being considered.
norm.meanvar = function(n, x, x2, p, length.checklist)
{

  # n, x, x2 -- all (p x length(checklist)) matrices

  sigmasq = (1/n)*(x2-(x^2)/n)
  neg = sigmasq <= 0
  sigmasq[neg==TRUE] = 0.00000000001

  return( .colSums(n*(log(2*pi)+log(sigmasq)+1), na.rm=T, m=p, n=length.checklist ) )

}



#' Multivariate Cost Calculation: Normally distributed observations, detecting changes in both mean and variance, with additional segment length penalty in likelihood.
#'
#' Calculates the cost between a sequence of multivariate observations using the provided cumulative sums of data, with the additional penalty of the segment length (\code{log(n)}). The observations are assumed to be Normally distributed. The mean and variance parameters are set to their maximum likelihood estimates. Used for detecting changes in both mean and variance.
#'
#' An internal function, not designed for use by the end-user.
#'
#' @param n A matrix containing the lengths between the changepoints in current changepoint vector being considered and the changepoints in the check-list of previous changepoint vectors.
#' @param x A \code{p} x \code{length.checklist} matrix containing the cumulative sums of the observations between the current changepoint vector and each corresponding changepoint vector on the check-list.
#' @param x2 The same as \code{x}, except containing the cumulative sums of the squared observations.
#' @param p The number of variables in the series.
#' @param length.checklist The number of possible previous changepoint vectors on the check-list being considered.
#'
#' @return A vector of length \code{length.checklist} containing the multivariate costs for each changepoint vector on the check-list being considered.
norm.meanvar.seglen = function(n, x, x2, p, length.checklist)
{

  # n, x, x2 -- all (p x length(checklist)) matrices

  sigmasq = (1/n)*(x2-(x^2)/n)
  neg = sigmasq <= 0
  sigmasq[neg==TRUE] = 0.00000000001

  cost = n*(log(2*pi)+log(sigmasq)+1) + log(n) # log(n) is the segment length part of the penalty

  return( .colSums(cost, na.rm=T, m=p, n=length.checklist ) )

}



#' Multivariate Cost Calculation for tri-variate series only: Normally distributed (variables 1 and 3) and Gamma distributed (variable 2) observations, detecting changes in both mean and variance.
#'
#' Calculates the cost between a sequence of tri-variate observations using the provided cumulative sums of data. The observations are assumed to be Normally distributed in variables 1 and 3, and Gamma distributed in variable 2. The mean and variance parameters are set to their maximum likelihood estimates. Used for detecting changes in both mean and variance.
#'
#' An internal function, not designed for use by the end-user.
#'
#' @param n A matrix containing the lengths between the changepoints in current changepoint vector being considered and the changepoints in the check-list of previous changepoint vectors.
#' @param x A \code{p} x \code{length.checklist} matrix containing the cumulative sums of the observations between the current changepoint vector and each corresponding changepoint vector on the check-list.
#' @param x2 The same as \code{x}, except containing the cumulative sums of the squared observations.
#' @param p The number of variables in the series.
#' @param length.checklist The number of possible previous changepoint vectors on the check-list being considered.
#'
#' @return A vector of length \code{length.checklist} containing the multivariate costs for each changepoint vector on the check-list being considered.
#' @export
norm13.gamma2 = function(n, x, x2, p, length.checklist)
{

  # n, x, x2 -- (p x length(checklist)) matrices

  norm.vars = c(1,3)
  gam.vars = 2

  shape = 1

  n.norm = n[norm.vars, ,drop=F]
  x.norm = x[norm.vars, ,drop=F]
  x2.norm = x2[norm.vars, ,drop=F]

  n.gam = n[gam.vars, ,drop=F]
  x.gam = x[gam.vars, ,drop=F]

  sigmasq = (1/n.norm)*(x2.norm-(x.norm^2)/n.norm)
  neg = sigmasq <= 0
  sigmasq[neg==TRUE] = 0.00000000001

  norm.cost = .colSums(n.norm*(log(2*pi)+log(sigmasq)+1), na.rm=T, m=2, n=length.checklist)
  #norm.cost = .colSums(n.norm*(log(2*pi)+log(sigmasq)+1), na.rm=T, m=length(norm.vars), n=length.checklist)

  gam.cost = .colSums(2*n.gam*shape*(log(x.gam) - log(n.gam*shape)) , na.rm=T, m= 1, n=length.checklist)
  #gam.cost = .colSums(2*n.gam*shape*log(x.gam) - 2*n.gam*shape*log(n.gam*shape) , na.rm=T, m=length(gam.vars), n=length.checklist)

  return( norm.cost + gam.cost )

}


norm13.gamma2.seglen = function(n, x, x2, p, length.checklist)
{

  # n, x, x2 -- (p x length(checklist)) matrices

  norm.vars = c(1,3)
  gam.vars = 2

  shape = 1

  n.norm = n[norm.vars, ,drop=F]
  x.norm = x[norm.vars, ,drop=F]
  x2.norm = x2[norm.vars, ,drop=F]

  n.gam = n[gam.vars, ,drop=F]
  x.gam = x[gam.vars, ,drop=F]

  sigmasq = (1/n.norm)*(x2.norm-(x.norm^2)/n.norm)
  neg = sigmasq <= 0
  sigmasq[neg==TRUE] = 0.00000000001

  norm.cost = .colSums(n.norm*(log(2*pi)+log(sigmasq)+1) + log(n.norm), na.rm=T, m=2, n=length.checklist)
  #norm.cost = .colSums(n.norm*(log(2*pi)+log(sigmasq)+1), na.rm=T, m=length(norm.vars), n=length.checklist)

  gam.cost = .colSums(2*n.gam*shape*(log(x.gam) - log(n.gam*shape)) + log(n.gam), na.rm=T, m=1, n=length.checklist)
  #gam.cost = .colSums(2*n.gam*shape*log(x.gam) - 2*n.gam*shape*log(n.gam*shape) , na.rm=T, m=length(gam.vars), n=length.checklist)

  return( norm.cost + gam.cost )

}


norm13.gamma2.exp4 = function(n, x, x2, p, length.checklist)
{

  # n, x, x2 -- (p x length(checklist)) matrices

  norm.vars = c(1,3)
  gam.vars = 2
  exp.vars = 4

  shape = 1

  n.norm = n[norm.vars, ,drop=F]
  x.norm = x[norm.vars, ,drop=F]
  x2.norm = x2[norm.vars, ,drop=F]

  n.gam = n[gam.vars, ,drop=F]
  x.gam = x[gam.vars, ,drop=F]

  n.exp = n[exp.vars, ,drop=F]
  x.exp = x[exp.vars, ,drop=F]

  sigmasq = (1/n.norm)*(x2.norm-(x.norm^2)/n.norm)
  neg = sigmasq <= 0
  sigmasq[neg==TRUE] = 0.00000000001

  norm.cost = .colSums(n.norm*(log(2*pi)+log(sigmasq)+1), na.rm=T, m=2, n=length.checklist)
  #norm.cost = .colSums(n.norm*(log(2*pi)+log(sigmasq)+1), na.rm=T, m=length(norm.vars), n=length.checklist)

  gam.cost = .colSums(2*n.gam*shape*(log(x.gam) - log(n.gam*shape)), na.rm=T, m=1, n=length.checklist)
  #gam.cost = .colSums(2*n.gam*shape*log(x.gam) - 2*n.gam*shape*log(n.gam*shape) , na.rm=T, m=length(gam.vars), n=length.checklist)

  exp.cost = .colSums(2*n.exp*(1 - log(n.exp/x.exp)), na.rm=T, m=1, n=length.checklist)
  #exp.cost = .colSums(2*n.exp*(1 - log(n.exp/x.exp)), na.rm=T, m=length(exp.vars), n=length.checklist)

  return( norm.cost + gam.cost + exp.cost )

}


norm13.gamma2.exp4.seglen = function(n, x, x2, p, length.checklist)
{

  # n, x, x2 -- (p x length(checklist)) matrices

  norm.vars = c(1,3)
  gam.vars = 2
  exp.vars = 4

  shape = 1

  n.norm = n[norm.vars, ,drop=F]
  x.norm = x[norm.vars, ,drop=F]
  x2.norm = x2[norm.vars, ,drop=F]

  n.gam = n[gam.vars, ,drop=F]
  x.gam = x[gam.vars, ,drop=F]

  n.exp = n[exp.vars, ,drop=F]
  x.exp = x[exp.vars, ,drop=F]

  sigmasq = (1/n.norm)*(x2.norm-(x.norm^2)/n.norm)
  neg = sigmasq <= 0
  sigmasq[neg==TRUE] = 0.00000000001

  norm.cost = .colSums(n.norm*(log(2*pi)+log(sigmasq)+1) + log(n.norm), na.rm=T, m=2, n=length.checklist)
  #norm.cost = .colSums(n.norm*(log(2*pi)+log(sigmasq)+1), na.rm=T, m=length(norm.vars), n=length.checklist)

  gam.cost = .colSums(2*n.gam*shape*(log(x.gam) - log(n.gam*shape)) + log(n.gam), na.rm=T, m=1, n=length.checklist)
  #gam.cost = .colSums(2*n.gam*shape*log(x.gam) - 2*n.gam*shape*log(n.gam*shape) , na.rm=T, m=length(gam.vars), n=length.checklist)

  exp.cost = .colSums(2*n.exp*(1 - log(n.exp/x.exp)) + log(n.exp), na.rm=T, m=1, n=length.checklist)
  #exp.cost = .colSums(2*n.exp*(1 - log(n.exp/x.exp)), na.rm=T, m=length(exp.vars), n=length.checklist)

  return( norm.cost + gam.cost + exp.cost )

}

