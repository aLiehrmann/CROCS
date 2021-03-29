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