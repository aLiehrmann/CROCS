#' postProcessingRuleFactory
#'
#' @description A function factory that returns post-processing rules used to select the start and end of peaks among the returned changepoints in respectively each Up* and Dw*. Up* stands for the successive non-decreasing changes and Dw* the successive non-increasing changes. 
#' @param rule A keyword to select the type of post-processing rule among 
#' * "maxjump": we select the up and down change with the largest mean-difference in Up* and Dw*; 
#' * "thinnest_peak": we select the last up change in Up* and the first down change in Dw*; 
#' * "largest_peak": we select the first up change in Up* and the last down change in Dw*.
#' @return The post-processing rule function which takes as argument the \code{changepoints} and the segment specific parameters \code{theta} and returns the start and the end of each peak.
#' @examples
#' changepoints <- c(1,2,3,4,5,6)
#' theta <- c(2,10,9,2,1,5,1)
#' rule <- postProcessingRuleFactory(rule="maxjump")
#' rule(changepoints=changepoints, theta=theta)
#' rule <- postProcessingRuleFactory(rule="largest_peak")
#' rule(changepoints=changepoints, theta=theta)
#' @export
postProcessingRuleFactory <- function(rule="maxjump"){
  f <- switch(rule,  
    maxjump = function(changepoints, theta) {
      trend <- sign(diff(theta))
      changepoints_tmp <- changepoints[trend!=0]
      theta_tmp <- theta[c(trend!=0,TRUE)]
      trend <- trend[trend!=0]
      if (length(trend)>1){
        trend_rle <- rle(trend)
        targets <- which(
          trend_rle$values[1:(length(trend_rle$values)-1)] == 1 &
          trend_rle$values[2:length(trend_rle$values)] == -1
        )
        if (length(targets)){
          offset_start <- purrr::map_dbl(targets, function(x){
            ifelse(x-1, sum(trend_rle$length[1:(x-1)]), 0)
          })
          offset_end <- purrr::map_dbl(targets+1, function(x){
            ifelse(x-1, sum(trend_rle$length[1:(x-1)]), 0)
          })
          list(
            start =  changepoints_tmp[
              purrr::map2_dbl(
                offset_start, 
                targets, 
                function(x,y){
                  successive_up_changes <- theta_tmp[(x+1):(x+1+trend_rle$length[y])]
                  x+which.max(diff(successive_up_changes))
                }
              )
            ],
            end = changepoints_tmp[
              purrr::map2_dbl(
                offset_end, 
                targets+1, 
                function(x,y){
                  successive_down_changes <- theta_tmp[(x+1):(x+1+trend_rle$length[y])]
                  x+which.min(diff(successive_down_changes))
                }
              )
            ]
          ) 
        } else {
          list(start=NA, end=NA)
        }
      } else {
        list(start=NA, end=NA)
      }
    }, 
    thinnest_peak = function(changepoints, theta) {
      trend <- sign(diff(theta))
      changepoints_tmp <- changepoints[trend!=0]
      trend <- trend[trend!=0]
      if (length(trend)>1){
        trend_rle <- rle(trend)
        targets <- which(
          trend_rle$values[1:(length(trend_rle$values)-1)] == 1 &
          trend_rle$values[2:length(trend_rle$values)] == -1
        )
        if (length(targets)){
          list(
            start =  changepoints_tmp[
              purrr::map_dbl(
                targets, 
                ~sum(trend_rle$length[1:.x])
              )
            ],
            end =  changepoints_tmp[
              purrr::map_dbl(
                targets, 
                ~sum(trend_rle$length[1:(.x+1)])-trend_rle$length[.x+1]+1
              )
            ]
          ) 
        } else {
          list(start=NA, end=NA)
        }
      } else {
        list(start=NA, end=NA)
      }
    },
    largest_peak = function(changepoints, theta) {
      trend <- sign(diff(theta))
      changepoints_tmp <- changepoints[trend!=0]
      trend <- trend[trend!=0]
      if (length(trend)>1){
        trend_rle <- rle(trend)
        targets <- which(
          trend_rle$values[1:(length(trend_rle$values)-1)] == 1 &
          trend_rle$values[2:length(trend_rle$values)] == -1
        )
        if (length(targets)){
          list(
            start = changepoints_tmp[
              purrr::map_dbl(
                targets, 
                ~sum(trend_rle$length[1:.x])-trend_rle$length[.x]+1
              )
            ],
            end = changepoints_tmp[
              purrr::map_dbl(
                targets, 
                ~sum(trend_rle$length[1:(.x+1)])
              )
            ]
          ) 
        } else {
          list(start=NA, end=NA)
        }
      } else {
        list(start=NA, end=NA)
      }
    }
  )
  attr(f, "type") <- rule
  f
} 