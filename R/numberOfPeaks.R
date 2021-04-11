numberOfPeaks <- function(theta) {
  if (length(theta) > 2) {
    sign_theta <- sign(diff(theta))
    sign_theta <- sign_theta[sign_theta!=0]
    if (length(sign_theta)>1){
      sum(sign_theta[1:(length(sign_theta)-1)] == 1 & 
        sign_theta[2:length(sign_theta)] == -1)
    } else {
      0
    }
  } else {
    0
  }
}