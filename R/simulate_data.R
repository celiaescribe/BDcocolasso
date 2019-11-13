#' @title Autoregressive covariance matrix
#' @description Create autoregressive covariance matrix
#' @param p Number of features
#' @return Covariance matrix

#' @export

cov_autoregressive <- function(p){
  cov <- matrix(0,p,p)
  for (i in 1:p){
    for (j in 1:p){
      cov[i,j] = 0.5**(abs(i-j))
    }
  }
  cov
}
