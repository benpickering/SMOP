#' Approximate Subset Multivariate Optimal Partitioning
#'
#' A method which implements the Approximate Subset Multivariate Optimal Partitioning (A-SMOP) algorithm. This algorithm is capable of detecting the presence of changepoints in a multivariate time series, and identifies which of the variables are affected for each detected change.
#'
#' This method implements the Approximate Subset Multivariate Optimal Partitioning (A-SMOP) algorithm of B. Pickering [2016]. This algorithm obtains the changepoint locations within a multivariate time series, and identifies the subsets of variables which are affected by each corresponding change. This is done via the minimisation of a penalised cost function using dynamic programming.
#'
#' A range of different cost functions and penalty values can be used. The use of hard restriction provides a faster but more approximate solution, whereas soft restriction is more accurate but requires more computation. Note that using soft restriction can become very slow even for moderate \code{p}, for practical purposes we recommend using hard restriction.
#'
#' Values currently supported for the cost function \code{cost.func} include:
#' \tabular{ll}{
#'    \code{"norm.mean"} \tab Used for detecting changes in mean in Normally-distributed data. Assumes fixed variance parameters (\code{= 1}). The mean parameters are set to their maximum likelihood estimates. \cr
#'    \code{"norm.var"} \tab Used for detecting changes in variance in Normally-distributed data. Assumes fixed mean parameters (\code{= 0}). The variance parameters are set to their maximum likelihood estimates. \cr
#'    \code{"norm.meanvar"} \tab Used for detecting changes in both mean and variance in Normally-distributed data. The mean and variance parameters are set to their maximum likelihood estimates. \cr
#'    \code{"norm.mean.seglen"}, \code{"norm.var.seglen"}, \code{"norm.meanvar.seglen"} \tab Identical to \code{"norm.mean"}, \code{"norm.var"} and \code{"norm.meanvar"}, respectively, except these contain an additional log(segment length) penalty term in the likelihood for each variable. Designed for use when using the modified BIC penalty (Zhang and Siegmund, 2007) to penalise changes. \cr
#' }
#'
#' @aliases ampelt
#'
#' @param data An \code{n} x \code{p} matrix representing a length \code{n} time series containing observations of \code{p} variables.
#' @param alpha The variable-specific penalty, used to penalise the addition of a given changepoint into a given variable. A non-negative numeric value.
#' @param beta The multivariate penalty, used to penalise the addition of a changepoint into the model. A non-negative numeric value.
#' @param min.dist The minimum distance allowed between any two changepoints. Required to have an integer value of at least 2.
#' @param cost.func The name of the (global) cost function used by the method, given as a string. See details for possible values.
#' @param window.size The size of the window considered to the left and right of a given changepoint when performing subset restriction. A non-negative integer.
#' @param hard.restrict Logical. If \code{TRUE} then hard subset restriction is used. If \code{FALSE} then soft subset restriction is used.
#' @param class Logical. If \code{TRUE} then an object of class \code{cptmv} is returned.
#' @param verbose Logical. If \code{TRUE} then information regarding the changepoint vector check-list is printed during the algorithm.
#'
#' @return If \code{class=TRUE} then an object of S4 class \code{cptmv} is returned. The slot \code{cpts} contains the changepoints that are returned. Otherwise, if \code{class=FALSE}, a list is returned containing the following elements:
#'
#'    \item{data.set}{The data set being analysed for changepoints.}
#'    \item{cost.func}{The name of the function used to calculate the cost.}
#'    \item{cpt.type}{The type of changes which are being detected, e.g. mean, mean and variance.}
#'    \item{alpha}{The value of the alpha penalty used.}
#'    \item{beta}{The value of the beta penalty used.}
#'    \item{num.cpt.vecs}{The number of changepoint vectors within the search-space considered.}
#'    \item{cpt.vecs}{A matrix containing the optimal changepoint vectors for the series.}
#'    \item{like}{The value of the likelihood for the optimal set of changepoint vectors.}
#'    \item{cpts}{The optimal changepoint locations in the series.}
#'    \item{subsets}{A logical matrix containing the optimal affected variable subsets for each of the detected changepoints.}
#'    \item{runtime}{The running time of the algorithm, in seconds.}
#'
#' @export
#' @seealso \code{\link{smop}}
#'
#' @examples
#' # Smaller example: Normal data, single change in mean at mid-point in 2/3 variables
#' n = 20; p=3
#' set.seed(100)
#' data = matrix(NA, n, p)
#' data[,1] = c( rnorm(n/2, 0, 1), rnorm(n/2, 10, 1) )
#' data[,2] = c( rnorm(n/2, 0, 1), rnorm(n/2, 10, 1) )
#' data[,3] = rnorm(n, 0, 1)
#' alpha = 2*log(n)
#' beta = 2*log(p)*log(n)
#' cost.func = "norm.mean.seglen"
#' window.size = 2
#' hard.restrict = TRUE
#' asmop.results = asmop(data=data, alpha=alpha, beta=beta, cost.func=cost.func, window.size=window.size, hard.restrict=hard.restrict)
#'
#' # Larger example: Normal data, multiple changes in variance
#' data("var.change.ex")
#' #plot.ts(var.change.ex, nc=1)
#' n = nrow(var.change.ex) # 500
#' p = ncol(var.change.ex) # 6
#' asmop.results.hard = asmop(data=var.change.ex, alpha=2*log(n), beta=2*log(p)*log(n), cost.func="norm.var.seglen", window.size=10, hard.restrict=TRUE)
#' # WARNING: Using soft restriction (below) can take a few minutes.
#' #asmop.results.soft = asmop(data=var.change.ex, alpha=2*log(n), beta=2*log(p)*log(n), cost.func="norm.var.seglen", window.size=10, hard.restrict=FALSE) # Provides a better segmentation compared to hard.
asmop = function(data, alpha=2*log(nrow(data)), beta=2*log(ncol(data))*log(nrow(data)), min.dist=2, cost.func="norm.meanvar.seglen", window.size, hard.restrict=TRUE, class=TRUE, verbose=FALSE){

  soft.restrict = !hard.restrict

  if(soft.restrict){
    warning("Using soft restriction can lead to long computation times, especially for larger p. We recommend using hard restriction for most cases.", immediate.=T)
  }

  window.restrict.subsets.hard = hard.restrict
  window.restrict.subsets.soft = soft.restrict
  window.restrict.subsets.soft.maximal = F

  alpha.beta.pelt = F
  alpha.beta.pruning.windowed = T
  exact.pruning = F

  start.asmop.runtime = proc.time()

  if(is.null(dim(data))){
    data = as.matrix(data)
  }

  n = nrow(data)
  p = ncol(data)

  ## Running univariate PELT on each variable ##
  uv.cpts.alpha = vector("list", p)
  uv.cpts.alphabeta = vector("list", p)

  initial.likelihood.value = 0

  #if(alpha=="mbic"){
  #  alpha = 2*log(n)
  #}

  cost.func.name = cost.func #deparse(substitute(cost.func))

  if(cost.func.name == "norm.mean"){
    pelt.cost.func = pelt.norm.mean.cost
    sum.stat = pelt.norm.sum
    cpt.type = "mean"
  } else if(cost.func.name == "norm.var"){
    pelt.cost.func = pelt.norm.var.cost
    sum.stat = pelt.norm.sum
    cpt.type = "variance"
  } else if(cost.func.name == "norm.meanvar"){
    pelt.cost.func = pelt.norm.meanvar.cost
    sum.stat = pelt.norm.sum
    cpt.type = "mean and variance"
  } else if(cost.func.name == "norm.mean.seglen"){
    pelt.cost.func = pelt.norm.mean.cost.seglen
    sum.stat = pelt.norm.sum
    initial.likelihood.value = -log(n)
    cpt.type = "mean"
  } else if(cost.func.name == "norm.var.seglen"){
    pelt.cost.func = pelt.norm.var.cost.seglen
    sum.stat = pelt.norm.sum
    initial.likelihood.value = -log(n)
    cpt.type = "variance"
  } else if(cost.func.name == "norm.meanvar.seglen"){
    pelt.cost.func = pelt.norm.meanvar.cost.seglen
    sum.stat = pelt.norm.sum
    initial.likelihood.value = -log(n)
    cpt.type = "mean and variance"
  } else if(cost.func.name == "norm13.gamma2"){
    pelt.cost.func = list(pelt.norm.meanvar.cost, pelt.gamma.cost, pelt.norm.meanvar.cost)
    sum.stat = list(pelt.norm.sum, pelt.gamma.sum, pelt.norm.sum)
    cpt.type = "mean and variance"
  } else if(cost.func.name == "norm13.gamma2.seglen"){
    pelt.cost.func = list(pelt.norm.meanvar.cost.seglen, pelt.gamma.cost.seglen, pelt.norm.meanvar.cost.seglen)
    sum.stat = list(pelt.norm.sum, pelt.gamma.sum, pelt.norm.sum)
    initial.likelihood.value = -log(n)
    cpt.type = "mean and variance"
  } else if(cost.func.name == "norm13.gamma2.exp4"){
    pelt.cost.func = list(pelt.norm.meanvar.cost, pelt.gamma.cost, pelt.norm.meanvar.cost, pelt.exp.cost)
    sum.stat = list(pelt.norm.sum, pelt.gamma.sum, pelt.norm.sum, pelt.gamma.sum)
    cpt.type = "mean and variance"
  } else if(cost.func.name == "norm13.gamma2.exp4.seglen"){
    pelt.cost.func = list(pelt.norm.meanvar.cost.seglen, pelt.gamma.cost.seglen, pelt.norm.meanvar.cost.seglen, pelt.exp.cost.seglen)
    sum.stat = list(pelt.norm.sum, pelt.gamma.sum, pelt.norm.sum, pelt.gamma.sum)
    initial.likelihood.value = -log(n)
    cpt.type = "mean and variance"
  } else{
    stop("Specified cost function type not recognised.")
  }

  for(j in 1:p){

    if(verbose) print(paste0("Running alpha-PELT on variable ", j, "..."))
    if(length(pelt.cost.func)==1){
      uv.cpts.alpha[[j]] = pelt(data[,j], pen=alpha, min.dist=min.dist, cost=pelt.cost.func, sum.stat=sum.stat, initial.likelihood.value=initial.likelihood.value)
    } else{
      uv.cpts.alpha[[j]] = pelt(data[,j], pen=alpha, min.dist=min.dist, cost=pelt.cost.func[[j]], sum.stat=sum.stat[[j]], initial.likelihood.value=initial.likelihood.value)
    }


    if(alpha.beta.pelt == TRUE){

      if(verbose) print(paste0("Running (alpha+beta)-PELT on variable ", j, "..."))
      uv.cpts.alphabeta[[j]] = pelt(data[,j], pen=alpha+beta, min.dist=min.dist, cost=pelt.cost.func, sum.stat=sum.stat, initial.likelihood.value=initial.likelihood.value)

    }

  }

  if(alpha.beta.pelt == TRUE){
    uv.cpts = vector("list", p)
    for(j in 1:p){

      uv.cpts[[j]] = sort(unique(c(uv.cpts.alpha[[j]], uv.cpts.alphabeta[[j]])))

    }
  } else{
    uv.cpts = uv.cpts.alpha
  }

  ## Create overall checklist ##
  overall.checklist = sort(unique(unlist(uv.cpts)))

  ## Obtaining 'free' variables based on alpha-PELT results ##
  free.variables = rep(NA, p)
  for(j in 1:p){

    if( length(uv.cpts.alpha[[j]]) == 0 ){
      free.variables[j] = F
    } else{
      free.variables[j] = T
    }

  }

  if(!any(free.variables)){ # i.e. if there are no detected changepoints in any variable, terminate algorithm.

    final.cpt.vectors = matrix(c(rep(0,p), rep(n,p)), nrow=2, ncol=p, byrow=T)
    final.likelihood.value = 9999999 #NA
    final.cpts = n
    final.subsets = matrix(T, nrow=1, ncol=p)
    rownames(final.subsets) = n
    colnames(final.subsets) = paste0("var_", 1:p)
    num.cpt.vecs = nrow(final.cpt.vectors)
    runtime = (proc.time() - start.asmop.runtime)[1]

    if(class==TRUE){ # return S4 class object

      results = new("cptmv", data.set=as.ts(data), cpt.type=cpt.type, cost.func=cost.func.name, alpha=alpha, beta=beta, num.cpt.vecs=num.cpt.vecs, cpt.vecs=final.cpt.vectors, like=final.likelihood.value, cpts=final.cpts, subsets=final.subsets, runtime=runtime)
      return(results)

    } else{ # return usual S3 list
      return.list = list(num.cpt.vecs, final.cpt.vectors, final.likelihood.value, final.cpts, final.subsets, runtime)
      names(return.list) = c("num.cpt.vecs", "cpt.vecs", "like", "cpts", "subsets")
      return(return.list)
    }

  }

  ## Creating affected variables subsets for each detected tau ##

  affected.variable.subsets = vector("list", length(overall.checklist))

  if( window.restrict.subsets.hard == F && window.restrict.subsets.soft == F && window.restrict.subsets.soft.maximal == F ){ # NO window restricting on subsets => all possible subsets for each detected tau.

    free.variables.indices = which(free.variables)
    p.temp = sum(free.variables)
    #if(p.temp==0) stop("No changepoints detected in any variable.") # this shouldn't happen any more.
    affected.variable.subsets.temp = matrix(F, nrow=2^p.temp-1, ncol=p)
    affected.variable.subsets.temp[, free.variables.indices] = matrix(as.logical(as.integer(unlist(strsplit(binary(1:(2^p.temp-1)), split="")))), nrow=2^p.temp-1, ncol=p.temp, byrow=T)

    for(tau.index in 1:length(overall.checklist)){
      affected.variable.subsets[[tau.index]] = affected.variable.subsets.temp
    }

  } else if( window.restrict.subsets.hard == T && window.restrict.subsets.soft == F && window.restrict.subsets.soft.maximal == F ){ # HARD window restricting on subsets => if for some detected cpt tau, another variable contains a cpt within a window around tau, then we say that those variables are definitely included in the affected subset for tau.

    for(tau.index in 1:length(overall.checklist)){

      tau = overall.checklist[tau.index]
      affected.variable.subsets[[tau.index]] = matrix(NA, nrow=1, ncol=p)

      for(j in 1:p){

        if( any( abs(tau - uv.cpts[[j]]) <= window.size ) ){
          affected.variable.subsets[[tau.index]][,j] = T # variable IS allowed to have a changepoint at tau
        } else{
          affected.variable.subsets[[tau.index]][,j] = F # variable is NOT allowed to have a changepoint at tau
        }

      }

    }

  } else if( window.restrict.subsets.hard == F && window.restrict.subsets.soft == F && window.restrict.subsets.soft.maximal == T ){ # SOFT window restricting on subsets => same as hard window restrict, except the 'hard restrict' is the maximal subset, in soft all subsets of that maximal subset are considered.

    allowed.variables = matrix(NA, nrow=length(overall.checklist), ncol=p)

    for(tau.index in 1:length(overall.checklist)){

      tau = overall.checklist[tau.index]

      for(j in 1:p){

        if( any( abs(tau - uv.cpts[[j]]) <= window.size ) ){
          allowed.variables[tau.index,j] = T # variable IS allowed to have a changepoint at tau
        } else{
          allowed.variables[tau.index,j] = F # variable is NOT allowed to have a changepoint at tau
        }

      }

      tau.allowed.variables.indices = which(allowed.variables[tau.index,])

      p.temp = sum(allowed.variables[tau.index,])
      affected.variable.subsets[[tau.index]] = matrix(F, nrow=2^p.temp-1, ncol=p)
      affected.variable.subsets[[tau.index]][, tau.allowed.variables.indices] = matrix(as.logical(as.integer(unlist(strsplit(binary(1:(2^p.temp-1)), split="")))), nrow=2^p.temp-1, ncol=p.temp, byrow=T)

    }

  } else if( window.restrict.subsets.hard == F && window.restrict.subsets.soft == T && window.restrict.subsets.soft.maximal == F ){ # SOFT window restricting on subsets => same as hard window restrict, except the 'hard restrict' is the minimal subset, in soft the permutations of all other subsets are considered on top of the minimal subset. **NOTE THIS IS WHAT WE WILL BE USING.**

    definite.variables = matrix(NA, nrow=length(overall.checklist), ncol=p)

    for(tau.index in 1:length(overall.checklist)){

      tau = overall.checklist[tau.index]

      for(j in 1:p){

        if( any( abs(tau - uv.cpts[[j]]) <= window.size ) ){
          definite.variables[tau.index,j] = T # affected variable subset for tau HAS to include j
        } else{
          definite.variables[tau.index,j] = F # affected variable subset for tau does not necessarily have to include j
        }

      }

      tau.definite.variables.indices = which(definite.variables[tau.index,])
      tau.not.definite.but.free.variables.indices = which(!definite.variables[tau.index,] & free.variables)

      p.temp = length(tau.not.definite.but.free.variables.indices)
      affected.variable.subsets[[tau.index]] = matrix(F, nrow=2^p.temp, ncol=p)
      affected.variable.subsets[[tau.index]][, tau.definite.variables.indices] = T
      #browser()
      affected.variable.subsets[[tau.index]][, tau.not.definite.but.free.variables.indices] = matrix(as.logical(as.integer(unlist(strsplit(binary(0:(2^p.temp-1)), split="")))), nrow=2^p.temp, ncol=p.temp, byrow=T)

    }

  } else{
    stop("Only one (at most) of window.restrict.subsets.hard, window.restrict.subsets.soft and window.restrict.subsets.soft.maximal can be true.")
  }

  ### Running SMOP ###

  if(verbose) print("Beginning SMOP algorithm...")

  if( cost.func.name=="norm.mean" ){
    smop.cost.func = norm.mean
  } else if( cost.func.name=="norm.var" ){
    smop.cost.func = norm.var
  } else if( cost.func.name=="norm.meanvar" ){
    smop.cost.func = norm.meanvar
  } else if( cost.func.name=="norm.mean.seglen" ){
    smop.cost.func = norm.mean.seglen
  } else if( cost.func.name=="norm.var.seglen" ){
    smop.cost.func = norm.var.seglen
  } else if( cost.func.name=="norm.meanvar.seglen" ){
    smop.cost.func = norm.meanvar.seglen
  } else if( cost.func.name=="norm.gamma" ){
    smop.cost.func = norm13.gamma2
  } else if( cost.func.name=="norm13.gamma2.seglen" ){
    smop.cost.func = norm13.gamma2.seglen
  } else if( cost.func.name=="norm13.gamma2.exp4" ){
    smop.cost.func = norm13.gamma2.exp4
  } else if( cost.func.name=="norm13.gamma2.exp4.seglen" ){
    smop.cost.func = norm13.gamma2.exp4.seglen
  } else{
    warning("Current possible cost functions are 'all variables Normal mean and/or variance'. Setting cost function to 'norm.meanvar.seglen'...")
    smop.cost.func = norm.meanvar.seglen
  }

  if(verbose){ print(paste0("SMOP checklist size: ", length(overall.checklist))) }

  min.dist = min.dist
  asmop.smop.results = smop.for.asmop(data=data, possible.cpts=overall.checklist, possible.subsets=affected.variable.subsets, alpha=alpha, beta=beta, min.dist=min.dist, cost.func=smop.cost.func, verbose=verbose, initial.likelihood.value=initial.likelihood.value)
  #asmop.smop.results = smop.for.asmop(data=data, checklist.cpts=overall.checklist, checklist.subsets=affected.variable.subsets, alpha=alpha, beta=beta, min.dist=min.dist, cost.func=smop.cost.func, alpha.beta.pelt=alpha.beta.pelt, alpha.beta.pruning.windowed=alpha.beta.pruning.windowed, window.size=window.size, print=print, uv.cpts.alphabeta=uv.cpts.alphabeta, exact.pruning=exact.pruning, initial.likelihood.value=initial.likelihood.value)

  asmop.smop.results$runtime = (proc.time() - start.asmop.runtime)[1]
  asmop.smop.results$cpt.type = cpt.type
  #asmop.smop.results$data.set = data
  ##asmop.smop.results$data.set = as.ts(data)
  asmop.smop.results$cost.func = cost.func.name
  asmop.smop.results$alpha = alpha
  asmop.smop.results$beta = beta

  if(class==TRUE){ # return S4 class object

    data.set = data
    #data.set = as.ts(data)

    results = new("cptmv", data.set=data.set, cpt.type=cpt.type, cost.func=cost.func.name, alpha=alpha, beta=beta, num.cpt.vecs=asmop.smop.results$num.cpt.vecs, cpt.vecs=asmop.smop.results$cpt.vecs, like=asmop.smop.results$like, cpts=asmop.smop.results$cpts, subsets=asmop.smop.results$subsets, runtime=asmop.smop.results$runtime)
    return(results)

  } else{ # return usual S3 list

    return(asmop.smop.results)

  }

}

#ampelt = asmop
# var.change.ex: scenario16rep6. rep1 is also good.
