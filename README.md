# Subset Multivariate Optimal Partitioning #

This R package implements the Subset Multivariate Optimal Partitioning (SMOP) method for detecting changepoints in multivariate time series. The method is capable of identifying both the time-points at which changes occur and the subset of variables which are affected by these changes. Details of the SMOP algorithm and related work are presented in my PhD thesis [[1](#references)].

## Getting Started ##

The main function in the package is `asmop`, which calls the Approximate SMOP method - this improves the speed of SMOP drastically with only a minor potential reduction in accuracy. See the documentation in the package (`?asmop`) for the input arguments and examples.

## Package ##

The package is currently in an experimental state. I will be intermittently updating and improving the code. If you wish to make any improvements or suggestions then please do not hesitate to make a pull request and submit an issue.

## References ##

[1] Benjamin Pickering. *Changepoint detection for acoustic sensing signals.* PhD thesis, Lancaster University, 2016.
