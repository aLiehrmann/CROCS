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
      log(2*sqrt(phi*data*(phi*data+1))+2*phi*data+1)
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