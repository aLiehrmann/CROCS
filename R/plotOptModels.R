#' plotOptModels
#'
#' @description Visualize all output optimal segmentations obtained with \code{CROCS::CROCS}.
#' @param CROCS_results_p A list of optimal segmentations (output of \code{CROCS::CROCS}).
#' @param data_p A data table objet. \code{data_p$y} are the observations. \code{data_p$x} are the indices of the observations. 
#' @param xlab_p A label for the x-axis.
#' @param xlab_p A label for the y-axis.
#' @param xmin_p A lower limit for the x-axis.
#' @param xmax_p An upper limit for the x-axis.
#' @return A ggplot2 object.
#' @examples
#' require(data.table)
#' counts <- as.data.table(CROCS::counts)
#' counts <- counts[sample.id == "McGill0036" & chunk==3,]
#' weights <- counts$chromEnd - counts$chromStart
#' coverage <- counts$coverage
#' std_g <- graphFactory(graph="std")
#' loss <- lossFactory(type="mean")
#' tr <- transformationFactory(transformation="anscombe_poisson")
#' rule <- postProcessingRuleFactory(rule="largest_peak")
#' my_peak_caller <- peakCallerFactory(
#'   mygraph=std_g, 
#'   transformation_f=tr, 
#'   loss_f=loss,
#'   postProcessingRule_f=rule
#' )
#' fit <- CROCS(
#'   data = coverage, 
#'   weights = weights, 
#'   lower_bound_peak = 1, 
#'   upper_bound_peak = 3, 
#'   solver = my_peak_caller
#' )
#' data <- data.table(x=counts$chromStart/10^3, y=tr(coverage))
#' plotOptModels(
#'   CROCS_results_p = fit,
#'   data_p = data,
#'   xlab_p = "position on chromosome (kb: kilo bases)",
#'   ylab_p = "anscombe transformation: sqrt(aligned sequence reads + 3/8)",
#'   xmin_p = 40035,
#'   xmax_p = 40140
#' ) 
#' @export
plotOptModels <- function(CROCS_results_p, data_p, xlab_p, ylab_p, xmin_p, xmax_p){
  dt_changepoints <- data.table::rbindlist(purrr::map(CROCS_results_p, ~data.table::data.table(
    x=data_p$x[.x$changepoints_vec[-length(.x$changepoints_vec)]],
    label=paste0("changepoints: ",.x$changepoints,"\npeaks: ",.x$peaks),
    item="changepoints"
  )))
  dt_changepoints <- dt_changepoints[!is.na(dt_changepoints$x),]
  dt_segments <- data.table::rbindlist(purrr::map(CROCS_results_p, ~data.table::data.table(
    x=data_p$x[c(1,.x$changepoints_vec[-length(.x$changepoints_vec)])],
    xend=data_p$x[.x$changepoints_vec],
    y=.x$theta,
    label=paste0("changepoints: ",.x$changepoints,"\npeaks: ",.x$peaks),
    item="segment means",
    changepoints = .x$changepoints
  )))
  dt_peaks <- data.table::rbindlist(purrr::map(CROCS_results_p, ~data.table::data.table(
    xmin=data_p$x[.x$start],
    xmax=data_p$x[.x$end],
    label=paste0("changepoints: ",.x$changepoints,"\npeaks: ",.x$peaks),
    item="peaks"
  ))) 
  dt_peaks <- dt_peaks[!is.na(dt_peaks$xmin),]
  ggplot2::ggplot()+
  ggplot2::geom_line(
    data = data_p, 
    ggplot2::aes(
      y=y, 
      x=x
    ),
    color="grey70"
  )+
  ggplot2::geom_vline(
    data = dt_changepoints, 
    ggplot2::aes(
      xintercept=x, 
      color=factor(item)
    ),
    size=0.4
  )+
  ggplot2::geom_segment(
    data = dt_segments, 
    ggplot2::aes(
      x=x,
      xend=xend, 
      y=y,
      yend=y,
      color=factor(item)
    )
  )+
  ggplot2::geom_rect(
    data = dt_peaks,
    ggplot2::aes(
      xmin=xmin,
      xmax=xmax,
      ymin=-Inf,
      ymax=Inf,
      fill=factor(item)
    ),
    color="red",
    alpha=0.2
  )+
  ggplot2::scale_fill_manual("",values=c("gold"))+
  ggplot2::scale_color_manual("", values=c("grey40", "Royalblue"))+
  ggplot2::facet_grid(
    factor(
      x = label, 
      levels = unique(dt_segments$label)[
        order(unique(dt_segments$changepoints))
      ])~.,
    scales="free"
  )+
  ggplot2::theme_bw()+
  ggplot2::xlab(xlab_p)+
  ggplot2::ylab(ylab_p)+
  ggplot2::theme(
    strip.background = ggplot2::element_rect(fill="grey95"), 
    text = ggplot2::element_text(size=15),
    legend.position = "bottom",
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(), 
    axis.line = ggplot2::element_line(colour = "black")
  )+
  ggplot2::coord_cartesian(xlim=c(xmin_p,xmax_p))
}