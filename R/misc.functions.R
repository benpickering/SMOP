#' Permutation Creation
#'
#' Creates all permutations of length \code{p} using the \code{n}-length vector of numbers given by \code{source.vector}.
#'
#' @param n The size of the source vector used to generate the permutations.
#' @param p The length of the permutations to be produced.
#' @param source.vector The vectors of numbers used to generate the permutations.
#'
#' @return perms A matrix containing all of the permutations generated.
#' @export
#'
#' @examples
#' create.perms(n=20, p=3) # produces a 20^3 matrix of permutations.
create.perms = function(n=20, p=3, source.vector=1:n){

  perms = matrix(NA, nrow=n^p, ncol=p)
  for(j in 1:p){ perms[,j] = rep(rep(source.vector, each=n^(p-j)), times=n^(j-1)) }
  return(perms)

}


binary = function(x)
{
  y <- as.integer(x)

  isna <- is.na(y)
  Y <- y[!isna]
  ans0 <- character(length(Y))
  z <- NULL
  while (any(Y > 0) | is.null(z)) {
    z <- Y%%2
    Y <- floor(Y/2)
    ans0 <- paste(z, ans0, sep = "")
  }
  ans <- rep(as.character(NA), length(y))
  ans[!isna] <- ans0

  y <- ans
  dim(y) <- dim(x)
  y
}
