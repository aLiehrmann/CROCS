#' @export
CROPS <- function(data, weights, lower_bound_lambda, upper_bound_lambda, 
  solver, phi = 1) {
  list_of_models <- list()
  last_index <- 1
  list_of_models[[last_index]] <- solver(
    data = data, 
    weights = weights, 
    lambda = lower_bound_lambda,
    phi = phi
  )
  last_index <- last_index + 1
  list_of_models[[last_index]] <- solver(
    data = data, 
    weights = weights, 
    lambda = upper_bound_lambda, 
    phi = phi
  )
  q <- dequer::queue()
  dequer::pushback(q, c(1, 2))
  while (length(q) != 0) {
    d <- dequer::pop(q)
    if (list_of_models[[d[[1]]]]$changepoints > list_of_models[[d[[2]]]]$changepoints + 1) {
      new_candidate_lambda <- (list_of_models[[d[[2]]]]$loss - list_of_models[[d[[1]]]]$loss) / (list_of_models[[d[[1]]]]$changepoints - list_of_models[[d[[2]]]]$changepoints)
      if (new_candidate_lambda < list_of_models[[d[[1]]]]$lambda | new_candidate_lambda > list_of_models[[d[[2]]]]$lambda) {
        new_candidate_lambda <- (list_of_models[[d[[1]]]]$lambda + list_of_models[[d[[2]]]]$lambda) / 2
      }
      new_candidate_model <- solver(
        data = data, 
        weights = weights, 
        lambda = new_candidate_lambda, 
        phi = phi
      )
      if (new_candidate_model$changepoints < list_of_models[[d[[1]]]]$changepoints & new_candidate_model$changepoints > list_of_models[[d[[2]]]]$changepoints) {
        last_index <- last_index + 1
        list_of_models[[last_index]] <- new_candidate_model
        dequer::pushback(q, c(d[[1]], last_index))
        dequer::pushback(q, c(last_index, d[[2]]))
      } else if (new_candidate_model$changepoints == list_of_models[[d[[1]]]]$changepoints & new_candidate_model$loss < list_of_models[[d[[1]]]]$loss) { # gfpop misses sometimes the optimal model with k changepoints. Really improves the results with the negbin loss (the roots computation is definitely not straightforward). 
        list_of_models[[d[[1]]]] <- new_candidate_model
        dequer::pushback(q, c(d[[1]], d[[2]]))
      } else if (new_candidate_model$changepoints == list_of_models[[d[[2]]]]$changepoints & new_candidate_model$loss < list_of_models[[d[[2]]]]$loss) { # gfpop misses sometimes the optimal model with k changepoints. Really improves the results with the negbin loss (the roots computation is definitely not straightforward).
        list_of_models[[d[[2]]]] <- new_candidate_model
        dequer::pushback(q, c(d[[1]], d[[2]]))
      }
    }
  }
  list_of_models
}