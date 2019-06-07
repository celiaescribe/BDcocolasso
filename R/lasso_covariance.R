#' Lasso in covariance form
#'
#' Solve the least squares loss with lasso penalty written in a form with the covariance matrix : \eq{\frac{1}{2} \beta^{'} \Sigma \beta - \rho^{'} \beta + \lambda \|\beta\|_1}
#'
#' @param n Number of samples of the design matrix
#' @param p Number of features of the matrix
#' @param lambda penalty parameter
#' @param control Including control parameters : max of iterations, tolerance for the convergence of the error, zero threshold to put to zero small beta coefficients
#' @param XX Design matrix corresponding to \eq{\frac{1}{n} X'X} or a modified version in the case of CoCoLasso
#' @param Xy Rho parameter corresponding to \eq{\frac{1}{n} X'y} or a modified version in the case of CoCoLasso
#' @param beta.start Initial value of beta
#' 
#' @return list containing \itemize{
#' \item coefficients : Coefficients corresponding to final beta after convergence of the algoritm
#' \item coef.list : Matrix of coefficients for beta for all iterations
#' \item num.it Number of iterations 
#' }
#' 
#' @example 
#' 
#' @export

lasso_covariance <- function(n,
                             p,
                             lambda, 
                             control = list(maxIter = 1000,
                                            optTol = 10^(-5), 
                                            zeroThreshold = 10^(-6)), 
                             XX, 
                             Xy, 
                             beta.start) {
  
  beta <- beta.start
  wp <- beta
  m <- 1
  
  ## compute the product of XX with beta
  s <- XX %*% beta
  
  while (m < control$maxIter) {
    beta_old <- beta
    for (j in 1:p) {
      # Compute the Shoot and Update the variable
      S0 <- s[j] - XX[j,j]*beta_old[j] - Xy[j]
      if (sum(is.na(S0)) >= 1) {
        beta[j] <- 0
        next
      }
      
      if (S0 > lambda){
        beta[j] <- (lambda - S0)/XX[j, j]
        s <- s + XX[,j]*( beta[j] - beta_old[j])
      } 
      if (S0 < -1 * lambda){
        beta[j] <- (-1 * lambda - S0)/XX[j, j]
        s <- s + XX[,j]*( beta[j] - beta_old[j])
      }
      if (abs(S0) <= lambda) 
        beta[j] <- 0
    }
    # Update
    wp <- cbind(wp, beta)
    # Check termination for early stopping
    if (sum(abs(beta - beta_old), na.rm = TRUE) < control$optTol) {
      break
    }
    m <- m + 1
  }
  w <- beta
  #We impose very small coefficients to be equal to zero
  w[abs(w) < control$zeroThreshold] <- 0
  return(list(coefficients = w, coef.list = wp, num.it = m))
}