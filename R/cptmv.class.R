## Creating cptmv class ##

#' cptmv: An S4 class for multivariate changepoints
#'
#' An S4 class used to return information on detected changepoints in multivariate time series.
#'
#' This S4 class is used to return the results of subset-multivariate changepoint detection methods. The information included relates the detected changepoints, including their locations, the variables which are affected, the likelihood value of the final segmentation. Other information such as the cost function and penalty values used and the running time of the detection procedure are also included.
#'
#' @slot data.set A ts object. The multivariate data set being analysed for changepoints.
#' @slot cost.func A character object. The name of the function used to calculate the (unpenalised) cost.
#' @slot cpt.type A character object. The type of change(s) which are being detected.
#' @slot alpha A numeric object. The value of the alpha penalty used.
#' @slot beta A numeric object. The value of the beta penalty used.
#' @slot num.cpt.vecs A numeric object. The number of changepoint vectors considered by the detection procedure.
#' @slot cpt.vecs A numeric matrix. The set of estimated changepoint vectors.
#' @slot like A numeric object. The likelihood value of the estimated segmentation.
#' @slot cpts A numeric object. The set of estimated changepoint locations.
#' @slot subsets A boolean matrix. The sets of affected variables corresponding to each detected changepoint.
#' @slot runtime A numeric object. The running time of the detection procedure, in seconds.
#'
#' @aliases cptmv-class
#' @exportClass cptmv
#' @author Benjamin Pickering
#'
cptmv = setClass(

  "cptmv",

  slots = c(
    data.set = "matrix",
    #data.set = "ts",
    cost.func = "character",
    cpt.type = "character",
    alpha = "numeric",
    beta = "numeric",
    num.cpt.vecs = "numeric",
    cpt.vecs = "matrix",
    like = "numeric",
    cpts = "numeric",
    subsets = "matrix",
    runtime = "numeric"
  )

)

# creating plot methdod for 'cptmv' class
#' Title
#'
#' @import zoo
#'
#' @param cptmv
#' @param x The x-axis values of the observations.
#'
#' @export
#'
#' @examples a
setMethod(f = "plot", signature = "cptmv",
          definition = function(x, ...){
            object = x
            data.set = as.matrix(object@data.set)
            cpts = object@cpts
            subsets = object@subsets
            cpt.vecs = object@cpt.vecs
            n = nrow(data.set)
            p = ncol(data.set)
            cpt.type = object@cpt.type

            data.rownames = rownames(data.set)

            cpts = cpts[-length(cpts)] # removing n as a cpt
            subsets = subsets[-nrow(subsets), ,drop=F] # removing subset for n

            panel = function(x, y, ...){

              pan.num = parent.frame()$panel.number
              cpts.current.var = cpts[ subsets[,pan.num] ]

              lines(x,y)
              if(is.ts(data.set)){
                abline(v=time(data.set)[cpts.current.var], col="blue", lty="dashed", lwd=1)
              } else{
                abline(v=cpts.current.var, col="blue", lty="dashed", lwd=1)
              }
              cpts.temp = c(0,cpts.current.var,n)
              if(cpt.type=="mean" || cpt.type=="mean and variance"){
                means = numeric(length(cpts.current.var)+1)
                for(i in 2:(length(cpts.temp))){
                  means[i-1] = mean(data[(cpts.temp[i-1]+1):cpts.temp[i],pan.num])
                  segments(y0=means[i-1], x0=(cpts.temp[i-1]), x1=cpts.temp[i], col="red", lty="solid", lwd=1)
                }
              }

              # Putting rownames on x-axis
              if(pan.num==p){
                if(is.ts(data.set) || is.null(data.rownames)){
                  axis(side=1, at=axTicks(1))
                } else if(!is.null(data.rownames)){
                  #axis(side=1, at=1:n, labels=data.rownames)

                  x.ticks = c(1, axTicks(1)[which( axTicks(1)!=0 )], nrow(data.set) ) #axTicks(1)
                  x.tick.labels = rownames(data.set)[x.ticks]
                  axis(1, at=x.ticks, labels=x.tick.labels)

                }
              }

            }

            plot.zoo(data.set, xaxt="n", panel=panel, nc=if(p<=10){1}else{2}, ylim=range(data.set), main="Plot of Detected Changepoints", xlab="Time", ...)

          }
)



