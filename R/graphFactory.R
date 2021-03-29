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