# Subset Multivariate Optimal Partitioning #

This R package implements the Subset Multivariate Optimal Partitioning (SMOP) method for detecting changepoints in multivariate time series. The method is capable of identifying both the time-points at which changes occur and the subset of variables which are affected by these changes. Details of the SMOP algorithm and related work are presented in my PhD thesis [[1](#references)].

## Getting Started ##

The package has two main functions: `smop` and `asmop`.  

* `smop` performs the original SMOP method without any efficiency modifications. This can become very slow for even moderately-sized problems (series with more than 50 observations and 3 variables). It is included mainly for academic purposes, and is not recommended for general use. 
* `asmop` performs Approximate SMOP, which contains a number of modifications which improve the speed of SMOP drastically with only a minor potential reduction in accuracy (see [[1](#references)] for more details on these). We recommend using this function, and controlling its behaviour using its input arguments. 

See the documentation in the package (`?asmop` and `?smop`) for the more details on these function and examples. Note that these functions will likely be combined in the future. 

## Package ##

The package is currently in an experimental state. I will be intermittently updating and improving the code. If you wish to make any improvements or suggestions then please do not hesitate to make a pull request and submit an issue.

## References ##

[1] Benjamin Pickering. *Changepoint detection for acoustic sensing signals.* PhD thesis, Lancaster University, 2016.
