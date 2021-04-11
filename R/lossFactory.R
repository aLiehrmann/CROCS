#' lossFactory
#'
#' @description A function factory that returns several weighted loss forms. 
#' @param type A keyword to select the noise model and associated weighted loss among 
#' * "mean" (gaussian): \code{sum((data-theta)^2*weights)}; 
#' * "poisson": \code{sum((theta-data*log(theta))*weights)}; 
#' * "negbin" (negative binomial): \code{sum((-log(theta)-data*log(1-theta))*weights)}.
#' @return A function which takes as arguments the observations \code{data}, associated \code{weights} and \code{theta} the distribution paramter affected by changes. The function returns the loss computed on the obsevations. 
#' @examples
#' loss <- lossFactory(type="mean")
#' @export
lossFactory <- function(type="mean"){
  f <- switch(type,
    mean = function(data, theta, weights) {
      sum((data-theta)^2*weights)
    },
    poisson = function(data, theta, weights) {
      theta[theta == 0] <- 1
      sum((theta-data*log(theta))*weights)
    },
    negbin = function(data, theta, weights) {
      sum((-log(theta)-data*log(1-theta))*weights)
    }
  )
  attr(f, "type") <- type
  f
}