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