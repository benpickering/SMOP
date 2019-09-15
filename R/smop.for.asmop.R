################################################################
### Function which performs SMOP within the A-SMOP algorithm ###
################################################################

#' Subset Multivariate Optimal Partitioning for Possible Changepoint Locations and Corresponding Subsets
#'
#' This method performs the Subset Multivariate Optimal Partitioning (SMOP) algorithm for a given set of possible changepoint locations and given sets of possible corresponding affected variable subsets. The method is designed for use as part of the Approximate SMOP algorithm, with the sets of possible changepoint locations being obtained using PELT, and the affected variable subsets being obtained using hard or soft restriction.
#'
#' The method takes a given set of possible changepoint locations and corresponding sets of possible affected variable subsets, and uses these to generate a search-space consisting of all possible changepoint vectors within the time series. An optimal partitioning search is used to obtain the optimal changepoint vector(s) from this search-space which minimise the penalised cost, as defined by \code{cost.func} and the penalties \code{alpha} and \code{beta}.
#'
#' The possible changepoint locations in \code{checklist.cpts} are obtained by applying PELT to each variable (i.e. column) in \code{data}. The corresponding sets of possible affected variable subsets are obtained by applying either hard restriction or soft restriction (specified by the user) to the set of possible changepoint locations.
#'
#' Values currently supported for the cost function \code{cost.func} include:
#' \tabular{ll}{
#'    \code{"norm.mean"} \tab Used for detecting changes in mean in Normally-distributed data. Assumes fixed variance parameters (\code{= 1}). The mean parameters are set to their maximum likelihood estimates. \cr
#'    \code{"norm.var"} \tab Used for detecting changes in variance in Normally-distributed data. Assumes fixed mean parameters (\code{= 0}). The variance parameters are set to their maximum likelihood estimates. \cr
#'    \code{"norm.meanvar"} \tab Used for detecting changes in both mean and variance in Normally-distributed data. The mean and variance parameters are set to their maximum likelihood estimates. \cr
#'    \code{"norm.mean.seglen"}, \code{"norm.var.seglen"}, \code{"norm.meanvar.seglen"} \tab Identical to \code{"norm.mean"}, \code{"norm.var"} and \code{"norm.meanvar"}, respectively, except these contain an additional log(segment length) penalty term in the likelihood for each variable. Designed for use when using the modified BIC penalty (Zhang and Siegmund, 2007) to penalise changes. \cr
#' }
#'
#' This method is designed to be internal, and is not intended for explicit use by the end-user.
#'
#' @param data An \code{n} x \code{p} matrix representing a length \code{n} time series containing observations of \code{p} variables.
#' @param possible.cpts A vector containing all possible changepoint locations, which will be used to generate the search-space of changepoint vectors.
#' @param possible.subsets A list of length \code{length(possible.cpts)}, where each element is a logical matrix containing \code{ncol(data)} columns. The rows of this matrix represent the different possible affected variable subsets for the corresponding element in \code{possible.cpts}.
#' @param alpha The variable-specific penalty, used to penalise the addition of a given changepoint into a given variable. A non-negative numeric value.
#' @param beta The multivariate penalty, used to penalise the addition of a changepoint into the model. A non-negative numeric value.
#' @param min.dist The minimum distance allowed between any two changepoints. Required to have an integer value of at least 2.
#' @param cost.func A function used to calculate the multivariate cost within SMOP. See details for more information.
#' @param initial.likelihood.value Numeric. The initial value of the likelihood.
#' @param verbose Logical. If \code{TRUE} then information regarding the changepoint vector check-list is printed during the algorithm.
#'
#' @return
#'  A list containing the following elements:
#'    \item{num.cpt.vecs}{The number of changepoint vectors within the search-space considered.}
#'    \item{cpt.vecs}{A matrix containing the optimal changepoint vectors for the series.}
#'    \item{like}{The value of the likelihood for the optimal set of changepoint vectors.}
#'    \item{cpts}{The optimal changepoint locations in the series.}
#'    \item{subsets}{A logical matrix containing the optimal affected variable subsets for each of the detected changepoints.}
#'
#'
smop.for.asmop = function(data, possible.cpts, possible.subsets, alpha, beta, min.dist, cost.func=norm.meanvar, initial.likelihood.value=0, verbose=F)
{

  if(is.null(dim(data))){
    data = as.matrix(data)
  }

  n = nrow(data)
  p = ncol(data)

  overall.checklist = possible.cpts #checklist.cpts
  affected.variable.subsets = possible.subsets #checklist.subsets
  d = min.dist # minimum distance between changepoints

  alpha.beta.pelt = F
  alpha.beta.pruning.windowed = T
  window.size = 10
  uv.cpts.alphabeta = NA
  exact.pruning = F

  # Create initial cpt vectors given the overall.checklist and corresponding possible subsets #

  cpt.vec = matrix(0, nrow=1, ncol=p) # going to rbind things to this

  cpt.vec.storage.list = vector("list", length(overall.checklist)+2)
  cpt.vec.storage.list[[1]] = matrix(0, nrow=1, ncol=p)

  total.num.vecs = 1

  for(tau.index in 1:length(overall.checklist)){

    tau = overall.checklist[tau.index]
    tau.subsets = affected.variable.subsets[[tau.index]]

    cpt.vec.temp = do.call(rbind, cpt.vec.storage.list)

    for(i in 1:nrow(tau.subsets)){

      if(all(tau.subsets[i,])){
        cpt.vec.temp.subset = tau * tau.subsets[i,]
      } else{
        cpt.vec.temp.subset = cpt.vec.temp[!duplicated(cpt.vec.temp[,!tau.subsets[i,] ,drop=F]), ,drop=F] # cpt.vec.temp
        cpt.vec.temp.subset[,tau.subsets[i,]] = tau
      }

      cpt.vec.storage.list[[tau.index+1]] = rbind( cpt.vec.storage.list[[tau.index+1]] , cpt.vec.temp.subset, deparse.level=0)

    }

    total.num.vecs = total.num.vecs + nrow(cpt.vec.storage.list[[tau.index+1]])
    if(verbose){ print(c(tau.index, nrow(cpt.vec.storage.list[[tau.index+1]]), total.num.vecs)) }

    if(nrow(cpt.vec)>300000){

      print(paste0("Number of changepoint vectors too large for SMOP (",  nrow(cpt.vec), " vectors). Terminating algorithm..."))

      final.cpt.vectors = matrix(c(rep(0,p), rep(n,p)), nrow=2, ncol=p, byrow=T)
      final.likelihood.value = 9999999 # NA
      final.cpts = n
      final.subsets = matrix(T, nrow=1, ncol=p)
      rownames(final.subsets) = n
      colnames(final.subsets) = paste0("var_", 1:p)

      return.list = list(final.cpt.vectors, final.likelihood.value, final.cpts, final.subsets)
      names(return.list) = c("cpt.vecs", "like", "cpts", "subsets")
      return(return.list)

    }

  }
  cpt.vec.storage.list[[length(overall.checklist)+2]] = rep(n,p)

  cpt.vec = do.call(rbind, cpt.vec.storage.list)

  num.cpt.vecs = nrow(cpt.vec)

  # Creating storage etc for SMOP #
  prev.cpt.locs = matrix(FALSE, nrow=num.cpt.vecs, ncol=n)
  FF = rep(NA, length=num.cpt.vecs)
  if(exact.pruning==T){
    pruned.vec = matrix(FALSE, nrow=num.cpt.vecs, ncol=num.cpt.vecs)
  }
  last.cpt.vec = rep(1, len=num.cpt.vecs)

  # Initialising values for SMOP #
  y2 = rbind(rep(0,p), apply(data^2,2,cumsum))
  y = rbind(rep(0,p), apply(data,2,cumsum))

  max.index = sum(cpt.vec[,1] < 2*d)
  initial.vec = 1:max.index
  if(p>1){
    for(j in 2:p){
      initial.vec = initial.vec[cpt.vec[initial.vec,j]<(2*d)]
    }
  }

  # putting initial costs in, creating initial prev.cpt.locs and creating initial checklists
  FF[1] = initial.likelihood.value # previously 0. could be -p*alpha ?
  for(index in initial.vec){

    temp.cpt.vec = cpt.vec[index,]

    n.temp = as.matrix(cpt.vec[index,] - cpt.vec[1,])
    n.nonzero.temp = colSums(n.temp!=0)

    x.temp  = as.matrix( diag(as.matrix(y[cpt.vec[index,]+1,])) - diag(as.matrix(y[(cpt.vec[1,]+1),])) )
    x2.temp = as.matrix( diag(as.matrix(y2[cpt.vec[index,]+1,])) - diag(as.matrix(y2[(cpt.vec[1,]+1),])) )

    FF[index] = calculate.alpha.cost(n=n.temp, n.nonzero=n.nonzero.temp, x=x.temp, x2=x2.temp, p=p, alpha=alpha, length.checklist=1, cost.func=cost.func)

    prev.cpt.locs[index, temp.cpt.vec[temp.cpt.vec!=0]] = TRUE

  }

  # Storage number of vectors pruned and calculations saved
  num.calcs.saved = 0
  num.pruned.vectors = 0

  # Cycling over cpt vectors #
  tindex = setdiff(2:num.cpt.vecs, initial.vec) # to remove the ones we have already done
  if(verbose) print( paste0("Total number of changepoint vectors: ", num.cpt.vecs) )
  for(cpt.index in tindex){

    if(verbose && cpt.index %% 10000 == 0){
      print( paste0("Total number of changepoint vectors: ", num.cpt.vecs, ". Changepoint vector index = ", cpt.index ) )
      print(cpt.vec[cpt.index,])
    }

    current.cpt.vec = cpt.vec[cpt.index,]

    # 1. get all valid cpt vectors in the checklist (i.e. where the element for each variable is less than cpt.index[j], for j=1:p)

    valid.indices = 1:sum(cpt.vec[,1]<=current.cpt.vec[1])

    for(j in 1:p){
      valid.indices = valid.indices[ cpt.vec[valid.indices,j] != max(current.cpt.vec) & ( cpt.vec[valid.indices,j] <= current.cpt.vec[j]-d | cpt.vec[valid.indices,j] == current.cpt.vec[j] ) ]
    }

    # 2. remove elements that have been pruned previously and valid for this cpt.index

    if(exact.pruning==T){

      num.valid.indices = length(valid.indices)
      previously.pruned.vecs = .colSums( X=pruned.vec[valid.indices, valid.indices], m=num.valid.indices, n=num.valid.indices )
      num.calcs.saved = num.calcs.saved + sum(previously.pruned.vecs != 0)
      checklist.temp = valid.indices[ previously.pruned.vecs == 0 ]

    } else{

      if( alpha.beta.pelt == F ){

        checklist.temp = valid.indices

      } else if( alpha.beta.pelt == T ){ # we only perform (alpha+beta)-pelt so that we can do this pruning.

        if( alpha.beta.pruning.windowed == F ){

          checklist.temp = valid.indices

          for(j in 1:p){

            pruning.possible = F
            for(tau in sort(uv.cpts.alphabeta[[j]], decreasing=T)){

              if(any(cpt.vec[checklist.temp,j]==tau) && (tau < current.cpt.vec[j]-d)){
                pruning.possible = T
                max.prev.cpt.j = tau
                break
              }

            }

            if(pruning.possible){
              checklist.temp = checklist.temp[which(cpt.vec[checklist.temp,j] >= max.prev.cpt.j)]
            }

          }

        } else if( alpha.beta.pruning.windowed == T ){

          checklist.temp = valid.indices

          for(j in 1:p){

            pruning.possible = F
            for(tau in sort(uv.cpts.alphabeta[[j]], decreasing=T)){

              if(any(cpt.vec[checklist.temp,j]==tau) && (tau < current.cpt.vec[j]-d)){
                pruning.possible = T
                max.prev.cpt.j.leftwindow = tau - max(d, window.size)
                break
              }

            }

            if(pruning.possible){
              checklist.temp = checklist.temp[which(cpt.vec[checklist.temp,j] >= max.prev.cpt.j.leftwindow)]
            }

          }

        }

      }

    }

    # 3. Create the temporary matrix containing the previous cpt locations for each element of the checklist (then adding current.cpt.vec to the locations)

    prev.cpt.locs.checklist = prev.cpt.locs[checklist.temp, 1:max(current.cpt.vec)]
    pcl.cvp.checklist = prev.cpt.locs[ last.cpt.vec[checklist.temp] , 1:max(current.cpt.vec)]

    # 4. calculate cost for each vector in the checklist

    # Creating necessary cumsum stuff [NOTE THE SWITCH OF DIMENSION - MAKES IT EASIER TO CALCULATE THE COST]
    length.checklist.temp = length(checklist.temp)
    cpt.distances = matrix(NA, nrow=p, ncol=length.checklist.temp)
    x.temp = matrix(NA, nrow=p, ncol=length.checklist.temp)
    x2.temp = matrix(NA, nrow=p, ncol=length.checklist.temp)

    for(j in 1:p){

      cpt.distances[j,] = current.cpt.vec[j] - cpt.vec[checklist.temp, j]
      x.temp[j,] = y[current.cpt.vec[j]+1,j] - y[cpt.vec[checklist.temp, j]+1,j]
      x2.temp[j,] = y2[current.cpt.vec[j]+1,j] - y2[cpt.vec[checklist.temp, j]+1,j]

    }

    cpt.distances.nonzero = .colSums( cpt.distances != 0, m=p, n=length.checklist.temp ) # remember - switched dimension!
    pcl.checklist.sizes = .rowSums( prev.cpt.locs.checklist , m=length.checklist.temp, n=max(current.cpt.vec) )
    pcl.cvp.checklist.sizes = .rowSums( pcl.cvp.checklist , m=length.checklist.temp, n=max(current.cpt.vec) )
    checklist.num.cpts = pcl.checklist.sizes - pcl.cvp.checklist.sizes

    checklist.costs = FF[checklist.temp] + calculate.alpha.cost(n=cpt.distances, n.nonzero=cpt.distances.nonzero, x=x.temp, x2=x2.temp, p=p, alpha=alpha, length.checklist=length.checklist.temp, cost.func=cost.func) + beta*checklist.num.cpts
    #checklist.costs = FF[checklist.temp] + calculate.alpha.cost.vectorised(n=cpt.distances, n.nonzero=cpt.distances.nonzero, x2=x2.temp, x=x.temp, p=p, alpha=alpha, length.checklist.temp=length.checklist.temp) + beta*checklist.num.cpts

    # 4b. create 'L' vector for each element of the checklist which is used as part of the pruning threshold

    if(exact.pruning == T){

      k = -(alpha+beta)*p # K-(alpha+beta)*p, but K=0.

      L = alpha*cpt.distances.nonzero + checklist.num.cpts*beta - k

    }

    # 5. decide which is the best and store it, i.e. likelihood and cpt

    min.checklist.cost.index = which.min(checklist.costs)
    last.cpt.vec[cpt.index] = checklist.temp[min.checklist.cost.index]
    FF[cpt.index] = checklist.costs[min.checklist.cost.index]

    # 6. Update the previously detected changepoints for this cpt vector (in prev.cpt.locs; including the current.cpt.vec locations)

    prev.cpt.locs[cpt.index, prev.cpt.locs[last.cpt.vec[cpt.index],]] = TRUE
    prev.cpt.locs[cpt.index, current.cpt.vec[current.cpt.vec!=0]] = TRUE

    # 7. prune this checklist and store indices of pruned changepoint vectors

    if(exact.pruning==T){

      pruned.checklist.indices = checklist.temp[checklist.costs > FF[cpt.index] + L]
      # these are the vectors which ARE pruned - i.e. the ones we DON'T keep.
      num.pruned.vectors = num.pruned.vectors + length(pruned.checklist.indices)
      pruned.vec[cpt.index, pruned.checklist.indices] = TRUE

    }

  }

  # 8. reverse engineer the cpts

  final.cpt.vec.indices = num.cpt.vecs
  most.recent.cpt.vec.index = num.cpt.vecs
  while(most.recent.cpt.vec.index != 1){

    most.recent.cpt.vec.index = last.cpt.vec[most.recent.cpt.vec.index]
    final.cpt.vec.indices = c( most.recent.cpt.vec.index, final.cpt.vec.indices )

  }

  final.cpt.vectors = cpt.vec[final.cpt.vec.indices,]
  colnames(final.cpt.vectors) = paste0("var_", 1:p)

  final.likelihood.value = FF[num.cpt.vecs] - p*alpha - beta

  final.cpts = sort(unique( setdiff(c(final.cpt.vectors), 0) ))
  final.subsets = matrix(FALSE, nrow=length(final.cpts), ncol=p)

  for(i in 1:length(final.cpts)){
    final.subsets[i, which( colSums( as.matrix(final.cpt.vectors == final.cpts[i]) ) > 0 ) ] = TRUE
  }
  rownames(final.subsets) = as.character(final.cpts)
  colnames(final.subsets) = paste0("var_", 1:p)

  # 9. Return final changepoint vectors and likelihood value

  return.list = list(num.cpt.vecs, final.cpt.vectors, final.likelihood.value, final.cpts, final.subsets)
  names(return.list) = c("num.cpt.vecs", "cpt.vecs", "like", "cpts", "subsets")

  #return.list = list(num.cpt.vecs, final.cpt.vectors, final.likelihood.value, final.cpts, final.subsets, num.calcs.saved, num.pruned.vectors)
  #names(return.list) = c("smop.checklist.size", "cpt.vecs", "like", "final.cpts", "final.subsets", "num.calcs.saved", "num.pruned.vectors")

  #return.list = list(overall.checklist, affected.variable.subsets, num.cpt.vecs, final.cpt.vectors, final.likelihood.value, final.cpts, final.subsets, num.calcs.saved, num.pruned.vectors)
  #names(return.list) = c("checklist.cpts", "checklist.subsets", "smop.checklist.size", "cpt.vecs", "like", "final.cpts", "final.subsets", "num.calcs.saved", "num.pruned.vectors")

  #return.list = list(overall.checklist, affected.variable.subsets, final.cpt.vectors, final.likelihood.value, final.cpts, final.subsets)
  #names(return.list) = c("checklist.cpts", "checklist.subsets", "cpt.vecs", "like", "final.cpts", "final.subsets")

  return(return.list)

}
