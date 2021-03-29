#' @export
transformationFactory <- function(transformation="anscombe_poisson"){
  f <- switch(transformation,
    raw = function(x,phi) {
      x/phi
    },
    anscombe_poisson = function(x,phi) {
      sqrt(x+3/8)
    },
    anscombe_negbin = function(x,phi) {
      sqrt(x+1/(2*phi))
    },
    delta_negbin = function(x,phi) {
      log(2*sqrt(phi*data*(phi*data+1))+2*phi*data+1)
    },
    delta_poisson = function(x,phi) {
      sqrt(x)
    },
    logp1 = function(x,phi) {
      log(x+1)
    }
  )
  attr(f, 'type') <- transformation
  f
}