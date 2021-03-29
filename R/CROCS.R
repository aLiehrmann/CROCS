#' @export
CROCS <- function(data, weights, lower_bound_peak, upper_bound_peak, solver, 
  phi = 1) {
  if (lower_bound_peak > 0) {
    lower_bound_peak <- lower_bound_peak - 1
  }
  upper_bound_peak <- upper_bound_peak + 1
  lower_bound_lambda <- sequentialSearch(
    data = data, 
    weights = weights, 
    target = upper_bound_peak, 
    solver = solver, 
    phi = phi
  )$lambda
  upper_bound_lambda <- sequentialSearch(
    data = data, 
    weights = weights, 
    target = lower_bound_peak, 
    solver = solver, 
    phi = phi
  )$lambda
  CROPS(
    data = data, 
    weights = weights, 
    lower_bound_lambda = lower_bound_lambda, 
    upper_bound_lambda = upper_bound_lambda, 
    solver = solver, 
    phi = phi
  )
}