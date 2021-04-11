#' graphFactory
#'
#' @description A function factory that returns, on the top of the `gfpop` method (Runge et al. 2020), graphs which encode peak shape assumptions. 
#' @param transformation A keyword to select the graph of peak shape assumptions among 
#' * "updown": We add inequality constraints to the successive segment specific parameters theta_1,...,theta_k so that non-decreasing changes in these parameters are always followed by non-increasing changes; 
#' * "updown_negbin": "updown" with reversed constraints. Specific to the negative binomial noise model.
#' * "std": The sequence of changes is not constrained.
#' @return A function which takes as argument the peanlty \code{lambda} parameter and returns the graph of constraints.
#' @examples
#' print(std_g <- graphFactory(graph="std"))
#' @export
graphFactory <- function(graph="std"){
  f <- switch(graph,
    std = function(lambda) {
      gfpop::graph(
      gfpop::Edge(state1 = "S1", state2 = "S1", type = "std", penalty = lambda),
      gfpop::Edge(state1 = "S1", state2 = "S1", type = "null", penalty = 0)
    )},
    updown = function(lambda) {
      gfpop::graph(
      gfpop::Edge(state1 = "BG", state2 = "UP", type = "up", penalty = lambda),
      gfpop::Edge(state1 = "UP", state2 = "BG", type = "down", penalty = lambda),
      gfpop::Edge(state1 = "BG", state2 = "BG", type = "null", penalty = 0),
      gfpop::Edge(state1 = "UP", state2 = "UP", type = "null", penalty = 0),
      gfpop::StartEnd(start = "BG", end = "BG")
    )},
    updown_negbin = function(lambda) {
      gfpop::graph(
      gfpop::Edge(state1 = "BG", state2 = "UP", type = "down", penalty = lambda),
      gfpop::Edge(state1 = "UP", state2 = "BG", type = "up", penalty = lambda),
      gfpop::Edge(state1 = "BG", state2 = "BG", type = "null", penalty = 0),
      gfpop::Edge(state1 = "UP", state2 = "UP", type = "null", penalty = 0),
      gfpop::StartEnd(start = "UP", end = "UP")
    )}
  )
  attr(f, "type") <- graph
  f
}