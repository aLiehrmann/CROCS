#' @export
peakCallerFactory <- function(mygraph_f=graphFactory(), transformation_f=transformation(), 
  loss_f=lossFactory(), postProcessingRule_f=postProcessingRuleFactory()){
  function(data, weights, lambda, phi=1){
    data <- transformation_f(data,phi)
    if (attr(loss_f,"type") == "negbin"){
      data[data==0] <- 10^-12
    }
    res <- gfpop::gfpop(
      mygraph = mygraph_f(lambda),
      data = data,
      weights = weights,
      type = attr(loss_f,'type'),
    )
    loss <- loss_f(
      data = data, 
      theta = rep(res$parameters, times = diff(c(0, res$changepoints))), 
      weights = weights
    )
    if (attr(loss_f,"type") == "negbin"){
      res$parameters <- 1-res$parameters
    }
    peaks_list <- postProcessingRule_f(
      changepoints = res$changepoints[-length(res$changepoints)],
      theta = res$parameters
    )
    list(
      peaks = ifelse(is.na(peaks_list$start)[[1]],0,length(peaks_list$start)),
      changepoints = length(res$changepoints)-1,
      lambda = lambda,
      phi = phi,
      loss = loss,
      start = peaks_list$start,
      end =peaks_list$end,
      changepoints_vec = res$changepoints,
      theta = res$parameters,
      states = res$states,
      forced = res$forced,
      transformation = attr(transformation_f, "type"),
      mygraph = attr(mygraph_f, "type"),
      type = attr(loss_f, "type"),
      post_processing_rule = attr(postProcessingRule_f, "type")
    ) 
  }
}