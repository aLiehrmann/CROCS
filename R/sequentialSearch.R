#' @export
sequentialSearch <- function(data, weights, target, solver, phi = 1) {
  upper_bound_model <- solver(
    data = data, 
    weights = weights, 
    lambda = 0, 
    phi = phi
  )
  lower_bound_model <- solver(
    data = data, 
    weights = weights, 
    lambda = Inf, 
    phi = phi
  )
  while (target != lower_bound_model$peaks & target != upper_bound_model$peaks) {
    new_candidate_lambda <- (upper_bound_model$loss - lower_bound_model$loss) / (lower_bound_model$changepoints - upper_bound_model$changepoints)
    if (new_candidate_lambda > lower_bound_model$lambda | new_candidate_lambda < upper_bound_model$lambda) {
      new_candidate_lambda <- (lower_bound_model$lambda + upper_bound_model$lambda) / 2
    }
    new_candidate_model <- solver(
      data = data, 
      weights = weights, 
      lambda = new_candidate_lambda, 
      phi = phi
    )
    if (new_candidate_model$changepoints == lower_bound_model$changepoints) {
      if (new_candidate_model$loss < lower_bound_model$loss) { # gfpop misses sometimes the optimal model with k changepoints. It really improves the results with the negbin loss (the roots computation is definitely not straightforward).
        lower_bound_model <- new_candidate_model
      } else {
        break
      }
    } else if (new_candidate_model$changepoints == upper_bound_model$changepoints){ 
      if (new_candidate_model$loss < upper_bound_model$loss) { # gfpop misses sometimes the optimal model with k changepoints. It really improves the results with the negbin loss (the roots computation is definitely not straightforward).
        upper_bound_model <- new_candidate_model
      } else {
        break
      }
    } else if (new_candidate_model$peaks == target) {
      lower_bound_model <- new_candidate_model
    } else if (new_candidate_model$peaks > target) {
      upper_bound_model <- new_candidate_model
    } else {
      lower_bound_model <- new_candidate_model
    }
  }
  lower_bound_model
}