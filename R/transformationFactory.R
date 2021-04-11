#' transformationFactory
#'
#' @description A function factory that returns several useful variance stabilizing transformations for data under Poisson or negative binomial noise models. 
#' @param transformation A keyword to select the type of transformation among 
#' * "raw": \code{x}; 
#' * "anscombe_poisson": \code{sqrt(x+3/8)}; 
#' * "anscombe_negbin": \code{sqrt(x+1/(2*phi))}; 
#' * "delta_negbin": \code{log(2*sqrt(phi*x*(phi*x+1))+2*phi*x+1)}; 
#' * "delta_poisson": \code{sqrt(x)}; 
#' * "logp1": \code{log(x+1)}.
#' @return A function which takes as arguments the observations \code{x} (and \code{phi} when needed) and returns the transformed observations.
#' @examples
#' require(data.table)
#' counts <- as.data.table(CROCS::counts)
#' counts <- counts[sample.id == "McGill0036" & chunk==3,]
#' head(coverage <- counts$coverage)
#' tr <- transformationFactory(transformation="anscombe_poisson")
#' head(tr(x=coverage))
#' tr <- transformationFactory(transformation="anscombe_negbin")
#' head(tr(x=coverage, phi=10))
#' @export
transformationFactory <- function(transformation="anscombe_poisson"){
  f <- switch(transformation,
    raw = function(x,phi=1) {
      x/phi
    },
    anscombe_poisson = function(x,phi=1) {
      sqrt(x+3/8)
    },
    anscombe_negbin = function(x,phi=1) {
      sqrt(x+1/(2*phi))
    },
    delta_negbin = function(x,phi=1) {
      log(2*sqrt(phi**(phi*x+1))+2*phi*x+1)
    },
    delta_poisson = function(x,phi=1) {
      sqrt(x)
    },
    logp1 = function(x,phi=1) {
      log(x+1)
    }
  )
  attr(f, 'type') <- transformation
  f
}