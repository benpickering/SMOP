#' Subset Multivariate Optimal Partitioning
#'
#' A method which implements the Subset Multivariate Optimal Partitioning (SMOP) algorithm. This algorithm is capable of detecting the presence of changepoints in a multivariate time series, and identifies which of the variables are affected for each detected change.
#'
#' This method implements the Approximate Subset Multivariate Optimal Partitioning (A-SMOP) algorithm of B. Pickering [2016]. This algorithm obtains the changepoint locations within a multivariate time series, and identifies the subsets of variables which are affected by each corresponding change. This is done via the minimisation of a penalised cost function using dynamic programming. A range of different cost functions and penalty values can be used.
#'
#' Note that since this algorithm uses an exact search method, it can become increasingly slow even for moderate \code{n} and \code{p}. Therefore, we suggest only using this method for smaller-scale time series, specifically those with a value of \code{p} no greater than 3, or a value of \code{n} no greater than 50. For larger series, we recommend the use of \code{\link{asmop}} with hard subset restriction.
#'
#' Values currently supported for the cost function \code{cost.func} include:
#' \tabular{ll}{
#'    \code{"norm.mean"} \tab Used for detecting changes in mean in Normally-distributed data. Assumes fixed variance parameters (\code{= 1}). The mean parameters are set to their maximum likelihood estimates. \cr
#'    \code{"norm.var"} \tab Used for detecting changes in variance in Normally-distributed data. Assumes fixed mean parameters (\code{= 0}). The variance parameters are set to their maximum likelihood estimates. \cr
#'    \code{"norm.meanvar"} \tab Used for detecting changes in both mean and variance in Normally-distributed data. The mean and variance parameters are set to their maximum likelihood estimates. \cr
#'    \code{"norm.mean.seglen"}, \code{"norm.var.seglen"}, \code{"norm.meanvar.seglen"} \tab Identical to \code{"norm.mean"}, \code{"norm.var"} and \code{"norm.meanvar"}, respectively, except these contain an additional log(segment length) penalty term in the likelihood for each variable. Designed for use when using the modified BIC penalty (Zhang and Siegmund, 2007) to penalise changes. \cr
#' }
#'
#' @param data An \code{n} x \code{p} matrix representing a length \code{n} time series containing observations of \code{p} variables.
#' @param alpha The variable-specific penalty, used to penalise the addition of a given changepoint into a given variable. A non-negative numeric value.
#' @param beta The multivariate penalty, used to penalise the addition of a changepoint into the model. A non-negative numeric value.
#' @param min.dist The minimum distance allowed between any two changepoints. Required to have an integer value of at least 2.
#' @param cost.func The name of the cost function used by the method, given as a string. See details for possible values.
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
#' @seealso \code{\link{asmop}}
#'
#' @examples
#' # Normal data, single change in mean at mid-point in 2/3 variables
#' n = 20; p=3
#' set.seed(100)
#' data = matrix(NA, n, p)
#' data[,1] = c( rnorm(n/2, 0, 1), rnorm(n/2, 10, 1) )
#' data[,2] = c( rnorm(n/2, 0, 1), rnorm(n/2, 10, 1) )
#' data[,3] = rnorm(n, 0, 1)
#' alpha = 2*log(n)
#' beta = 2*log(p)*log(n)
#' cost.func = "norm.mean.seglen" # See help(norm.mean.seglen)
#' smop.results = smop(data=data, alpha=alpha, beta=beta, cost.func=cost.func)
#'
smop = function(data, alpha=2*log(nrow(data)), beta=2*log(ncol(data))*log(nrow(data)), min.dist=2, cost.func="norm.meanvar", class=TRUE, verbose=TRUE){

  start.smop.runtime = proc.time()

  if(is.null(dim(data))){
    data = as.matrix(data)
  }

  n = nrow(data)
  p = ncol(data)

  warning("SMOP requires a long computation time in general. For practical purposes, we recommend usage of the A-SMOP method.", immediate.=T)

  if( (n>=50 && p>=3) || (n>=200 && p>=2) ) {
    warning("SMOP can take a **very** long time for even moderately sized n>50 and p>2.", immediate.=T)
  }

  d = min.dist

  #if(alpha=="mbic"){
  #  alpha = 2*log(n)
  #}

  cost.func.name = cost.func

  if(cost.func.name == "norm.mean"){
    cpt.type = "mean"
  } else if(cost.func.name == "norm.var"){
    cpt.type = "variance"
  } else if(cost.func.name == "norm.meanvar"){
    cpt.type = "mean and variance"
  } else if(cost.func.name == "norm.mean.seglen"){
    cpt.type = "mean"
  } else if(cost.func.name == "norm.var.seglen"){
    cpt.type = "variance"
  } else if(cost.func.name == "norm.meanvar.seglen"){
    cpt.type = "mean and variance"
  } else{
    warning("Current possible cost functions are 'all variables Normal mean and/or variance'. Setting cost function to 'norm.meanvar'...")
    cost.func.name = "norm.meanvar"
  }

  cost.func = eval(as.name(cost.func.name))

  y2 = rbind(rep(0,p), apply(data^2,2,cumsum))
  y = rbind(rep(0,p), apply(data,2,cumsum))

  # create all possible cpt vectors
  source.vector = c(0, d:(n-d))
  n.cpt.vec = (n-(2*d-1) +1)^p + 1
  cpt.vec = matrix(NA, nrow=n.cpt.vec, ncol=p)
  for(j in 1:p){ cpt.vec[-n.cpt.vec,j] = rep(rep(source.vector, each=(n-(2*d-1) +1)^(p-j)), times=(n-(2*d-1) +1)^(j-1)) }
  cpt.vec[n.cpt.vec,] = rep(n,p)
  rm(source.vector)

  ## Storing the locations of any changes which have occurred up to but NOT INCLUDING **EACH** cpt.vec ##
  prev.cpt.locs = matrix(FALSE, nrow = n.cpt.vec, ncol = n)

  ## CREATING LAST-CHANGE LIKELIHOOD LIST FOR **ALL** C.TAUSTAR ##
  FF = rep(NA, length=n.cpt.vec)

  # for each cpt.vec up to when 2 cpts can be added we need to initialize them
  # i.e. any vector where all components are less than 2d
  max.index = sum(cpt.vec[,1] < 2*d)
  initial.vec = 2:max.index
  if(p>1){
    for(i in 2:p){
      initial.vec = initial.vec[cpt.vec[initial.vec,i]<(2*d)]
    }
  }

  # putting initial costs in, creating initial prev.cpt.locs and creating initial checklists
  FF[1] = 0
  for(index in initial.vec){

    temp.cpt.vec = cpt.vec[index,]

    n.temp = as.matrix(cpt.vec[index,] - cpt.vec[1,])
    n.nonzero.temp = colSums(n.temp!=0)

    x.temp  = as.matrix( diag(as.matrix(y[cpt.vec[index,]+1,])) - diag(as.matrix(y[(cpt.vec[1,]+1),])) )
    x2.temp = as.matrix( diag(as.matrix(y2[cpt.vec[index,]+1,])) - diag(as.matrix(y2[(cpt.vec[1,]+1),])) )

    FF[index] = calculate.alpha.cost(n=n.temp, n.nonzero=n.nonzero.temp, x=x.temp, x2=x2.temp, p=p, alpha=alpha, length.checklist=1, cost.func=cost.func)

    prev.cpt.locs[index, temp.cpt.vec[temp.cpt.vec!=0]] = TRUE

  }

  ## CREATING LAST-CHANGEPOINT VECTOR FOR **EVERY** cpt.vec ##
  last.cpt.vec = rep(1, len=n.cpt.vec)

  # Cycling over cpt vectors #
  tindex = setdiff(2:n.cpt.vec, initial.vec) # to remove the ones we have already done
  if(verbose){
    print( paste0("Total number of changepoint vectors: ", n.cpt.vec) )
  }
  for(cpt.index in tindex){

    if(verbose && cpt.index %% 10000 == 0){
      print( paste0("Total number of changepoint vectors: ", n.cpt.vec, ". Changepoint vector index = ", cpt.index ) )
      print(cpt.vec[cpt.index,])
    }

    current.cpt.vec = cpt.vec[cpt.index,]

    # 1. get all valid cpt vectors in the checklist

    valid.indices = 1:sum(cpt.vec[,1]<=current.cpt.vec[1])

    for(j in 1:p){
      valid.indices = valid.indices[ cpt.vec[valid.indices,j] != max(current.cpt.vec) & ( cpt.vec[valid.indices,j] <= current.cpt.vec[j]-d | cpt.vec[valid.indices,j] == current.cpt.vec[j] ) ]
    }

    # 2. remove elements that have been pruned previously and valid for this cpt.index

    checklist.temp = valid.indices

    # 3. Create the temporary matrix containing the previous cpt locations for each element of the checklist (then adding current.cpt.vec to the locations)

    prev.cpt.locs.checklist = prev.cpt.locs[checklist.temp, 1:max(current.cpt.vec)]
    pcl.cvp.checklist = prev.cpt.locs[ last.cpt.vec[checklist.temp] , 1:max(current.cpt.vec)]
    # created such that instead of n columns, we only create as many as necessary (depends on the value of current.cpt.vec)

    # 4. calculate cost for each vector in the checklist

    # [NOTE THE SWITCH OF DIMENSION - MAKES IT EASIER TO CALCULATE THE COST]
    length.checklist.temp = length(checklist.temp)
    cpt.distances = matrix(NA, nrow=p, ncol=length.checklist.temp)
    x.temp = matrix(NA, nrow=p, ncol=length.checklist.temp)
    x2.temp = matrix(NA, nrow=p, ncol=length.checklist.temp)

    for(j in 1:p){

      cpt.distances[j,] = current.cpt.vec[j] - cpt.vec[checklist.temp, j]
      x.temp[j,] = y[current.cpt.vec[j]+1,j] - y[cpt.vec[checklist.temp, j]+1,j]
      x2.temp[j,] = y2[current.cpt.vec[j]+1,j] - y2[cpt.vec[checklist.temp, j]+1,j]

    }

    cpt.distances.nonzero = .colSums( cpt.distances!=0, m=p, n=length.checklist.temp ) # remember - switched dimension!

    pcl.checklist.sizes = .rowSums( prev.cpt.locs.checklist , m=length.checklist.temp, n=max(current.cpt.vec) )
    pcl.cvp.checklist.sizes = .rowSums( pcl.cvp.checklist , m=length.checklist.temp, n=max(current.cpt.vec) )
    checklist.num.cpts = pcl.checklist.sizes - pcl.cvp.checklist.sizes

    checklist.costs = FF[checklist.temp] + calculate.alpha.cost(n=cpt.distances, n.nonzero=cpt.distances.nonzero, x=x.temp, x2=x2.temp, p=p, alpha=alpha, length.checklist=length.checklist.temp, cost.func=cost.func) + beta*checklist.num.cpts

    # 5. decide which is the best and store it, i.e. likelihood and cpt

    min.checklist.cost.index = which.min(checklist.costs)
    last.cpt.vec[cpt.index] = checklist.temp[min.checklist.cost.index]
    FF[cpt.index] = checklist.costs[min.checklist.cost.index]

    # 6. Update the previously detected changepoints for this cpt vector (in prev.cpt.locs; including the current.cpt.vec locations)

    prev.cpt.locs[cpt.index, prev.cpt.locs[last.cpt.vec[cpt.index],]] = TRUE
    prev.cpt.locs[cpt.index, current.cpt.vec[current.cpt.vec!=0]] = TRUE

  }

  # 8. reverse engineer the cpts

  final.cpt.vec.indices = n.cpt.vec
  most.recent.cpt.vec.index = n.cpt.vec
  while(most.recent.cpt.vec.index != 1){

    most.recent.cpt.vec.index = last.cpt.vec[most.recent.cpt.vec.index]
    final.cpt.vec.indices = c( most.recent.cpt.vec.index, final.cpt.vec.indices )

  }

  final.cpt.vectors = cpt.vec[final.cpt.vec.indices,]
  final.likelihood.value = FF[n.cpt.vec] - p*alpha - beta

  final.cpts = sort(unique( setdiff(c(final.cpt.vectors), 0) ))
  final.subsets = matrix(FALSE, nrow=length(final.cpts), ncol=p)

  for(i in 1:length(final.cpts)){
    final.subsets[i, which( colSums( as.matrix(final.cpt.vectors == final.cpts[i]) ) > 0 ) ] = TRUE
  }
  rownames(final.subsets) = as.character(final.cpts)
  colnames(final.subsets) = paste0("var_", 1:p)

  # 9. Return final changepoint vectors and likelihood value

  smop.runtime = (proc.time() - start.smop.runtime)[1]

  if(class==TRUE){ # return S4 class object

    data.set = data
    #data.set = as.ts(data)

    results = new("cptmv", data.set=data.set, cpt.type=cpt.type, cost.func=cost.func.name, alpha=alpha, beta=beta, num.cpt.vecs=n.cpt.vec, cpt.vecs=final.cpt.vectors, like=final.likelihood.value, cpts=final.cpts, subsets=final.subsets, runtime=smop.runtime)
    return(results)

  } else{ # return usual S3 list object

    data.set = data #as.ts(data)

    return.list = list(data.set, cost.func.name, cpt.type, alpha, beta, n.cpt.vec, final.cpt.vectors, final.likelihood.value, final.cpts, final.subsets, smop.runtime)
    names(return.list) = c("data.set", "cost.func", "cpt.type", "alpha", "beta","num.cpt.vecs", "cpt.vecs", "like", "cpts", "subsets", "runtime")
    return(return.list)

  }

}
