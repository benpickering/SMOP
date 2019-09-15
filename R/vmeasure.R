#' V-Measure
#'
#' Calculates the V-measure of a subset-multivariate changepoint segmentation produced by \code{smop} or \code{asmop}.
#'
#' @param segmentation A \code{cpt.mv} object produced using \code{smop} or \code{asmop}, representing an estimated set of changepoints and affected variable subsets.
#' @param true.cpts The vector of true changepoint locations.
#' @param true.subsets The matrix of true affected variable subsets corresponding to the true changepoints.
#'
#' @return The V-measure of the provided segmentation produced by SMOP/A-SMOP.
#' @export
#'
vmeasure = function(segmentation, true.cpts, true.subsets)
{

  # true.cpts - SHOULDN'T include 0 but SHOULD include n, and likewise for corresponding true.subsets.

  if(class(segmentation) == "cptmv"){

    n = nrow(segmentation@data.set)
    p = ncol(segmentation@data.set)

    if(true.cpts[length(true.cpts)]!=n){
      stop("true.cpts should contain n as the last changepoint.")
    }

    est.cpts <- segmentation@cpts
    est.subsets <- segmentation@subsets

  } else if(class(segmentation)=="list"){

    n <- true.cpts[length(true.cpts)]
    p <- ncol(true.subsets)

    est.cpts <- segmentation$cpts
    est.subsets <- segmentation$subsets

  } else{
    stop("Input should be a cptmv object, or the list output from asmop().")
  }

  if(length(true.cpts)!=nrow(true.subsets)){
    stop("length(true.cpts) should equal nrow(true.subsets).")
  }

  if( !all(true.subsets[nrow(true.subsets),]) ){
    stop("Final row of true.subsets should contain all TRUE elements.")
  }

  # Creating true labels based on true.cpts and true.subsets

#   if(true.cpts[length(true.cpts)]==n){
#     true.cpts = true.cpts[-length(true.cpts)]
#     true.subsets = true.subsets[-nrow(true.subsets),]
#   }

  true.labels = matrix(NA, nrow=n, ncol=p)
  true.cpt.locs = rbind( rep(0,p), true.cpts * true.subsets )
  #true.cpt.locs = rbind( rep(0,p), true.cpts * true.subsets, rep(n,p) )

  for(i in 2:nrow(true.cpt.locs)){
    for(j in 1:p){
      if(true.cpt.locs[i,j]==0){
        true.cpt.locs[i,j] = true.cpt.locs[i-1,j]
      }
      if(true.cpt.locs[i,j] != true.cpt.locs[i-1,j]){ # i.e. if a change has occurred in variable j
        true.labels[(true.cpt.locs[i-1,j]+1):true.cpt.locs[i,j] ,j] = i-1
      }
    }
  }

  # Creating estimated labels based on segmentation@cpts and segmentation@subsets

  seg.labels = matrix(NA, nrow=n, ncol=p)
  seg.cpt.locs = rbind( rep(0,p), est.cpts * est.subsets )
  #seg.cpts.temp = segmentation@cpts[-length(segmentation@cpts)]
  #seg.subsets.temp = segmentation@subsets[-nrow(segmentation@subsets),]
  #seg.cpt.locs = rbind( rep(0,p), seg.cpts.temp * seg.subsets.temp, rep(n,p) )

  for(i in 2:nrow(seg.cpt.locs)){
    for(j in 1:p){
      if(seg.cpt.locs[i,j]==0){
        seg.cpt.locs[i,j] = seg.cpt.locs[i-1,j]
      }
      if(seg.cpt.locs[i,j] != seg.cpt.locs[i-1,j]){ # i.e. if a change has occurred in variable j
        seg.labels[(seg.cpt.locs[i-1,j]+1):seg.cpt.locs[i,j] ,j] = i-1
      }
    }
  }

  # Calculating V-measure
  N = n*p # = length(true.labels)
  TT = table(seg.labels, true.labels)
  CK = - sum(apply(TT, 1, function(x) return(sum(x[which(x>0)]*log(x[which(x>0)]/sum(x))))))/N
  KC = - sum(apply(TT, 2, function(x) return(sum(x[which(x>0)]*log(x[which(x>0)]/sum(x))))))/N
  K = - sum(apply(TT, 1, function(x) return(sum(x)*log(sum(x)/N))))/N
  C = - sum(apply(TT, 2, function(x) return(sum(x)*log(sum(x)/N))))/N
  if(C!=0){
    h = 1 - CK/C
  }
  else{
    h = 0
  }
  if(K!=0){
    cc = 1 - KC/K
  }
  else{
    cc = 0
  }

  vmeasure.value = 2*h*cc/(h+cc+1e-100)

  return(vmeasure.value)

}

