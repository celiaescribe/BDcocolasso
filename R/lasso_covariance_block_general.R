#' Lasso in covariance form for the BD-CoCoLasso
#'
#' Solve the least squares loss with lasso penalty written in a form with the covariance matrix : \eqn{\frac{1}{2} \beta^{'} \Sigma \beta - \rho^{'} \beta + \lambda \|\beta\|_1}
#'
#' @param n Number of samples of the design matrix
#' @param p1 Number of uncorrupted predictors
#' @param p2 Number of corrupted predictors containing additive error
#' @param p3 Number of corrupted predictors containing missingness
#' @param X1 first block of the design matrix corresponding to uncorrupted features
#' @param Z2 second block of the design matrix corresponding to corrupted features containing additive error
#' @param Z3 third block of the design matrix corresponding to corrupted features containing missingness
#' @param y Response vector
#' @param sigma1 Covariance matrix for X1 : \eqn{\frac{1}{n} X_1'X_1}. This parameter is automatically furnished in \link{blockwise_coordinate_descent}
#' @param sigma2 Modified covariance matrix for Z2 through the CoCoLasso algorithm. This parameter is automatically furnished in \link{blockwise_coordinate_descent}
#' @param sigma3 Modified covariance matrix for Z3 through the CoCoLasso algorithm. This parameter is automatically furnished in \link{blockwise_coordinate_descent}
#' @param lambda Penalty parameter
#' @param ratio_matrix Observation matrix in the missing data setting (NULL otherwise)
#' @param control Including control parameters : max of iterations, tolerance for the convergence of the error, zero threshold to put to zero small beta coefficients
#' @param beta1.start Initial value for the coefficients of uncorrupted features
#' @param beta2.start Initial value for the coefficients of corrupted features containing additive error
#' @param beta3.start Initial value for the coefficients of corrupted features containing missingness
#' @param penalty Type of penalty used : can be lasso penalty or SCAD penalty
#' 
#' @return list containing \itemize{
#' \item coefficients.beta1 : Coefficients corresponding to final beta1 after convergence of the algoritm
#' \item coefficients.beta2 : Coefficients corresponding to final beta2 after convergence of the algoritm
#' \item num.it : Number of iterations of algorithm
#' }
#' 
#' 
#' @export

lasso_covariance_block_general <- function(n,
                                   p1,
                                   p2,
                                   p3,
                                   X1,
                                   Z2,
                                   Z3,
                                   y,
                                   sigma1,
                                   sigma2,
                                   sigma3,
                                   lambda,
                                   ratio_matrix = NULL,
                                   control = list(maxIter = 1000,
                                                  optTol = 10^(-5), 
                                                  zeroThreshold = 10^(-6)),
                                   beta1.start,
                                   beta2.start,
                                   beta3.start,
                                   penalty=c("lasso","SCAD")){
  
  if (p1 != 0) {
    beta1 <- beta1.start
    beta2 <- beta2.start
    beta3 <- beta3.start
    objective_function <- 1000
    m <- 1
    error_beta1 <- matrix(0,control$maxIter,1)
    error_beta2 <- matrix(0,control$maxIter,1)
    error_beta3 <- matrix(0,control$maxIter,1)
    
    rho1 <- 1/n * (t(X1) %*% y)
    rho2 <- 1/n * (t(Z2) %*% y)
    Z3_tilde <- sapply(1:p3,function(j)Z3[,j] / diag(ratio_matrix)[j])
    rho3 <- 1/n * (t(Z3) %*% y) / diag(ratio_matrix)
    
    while (m < control$maxIter){
      beta1.old <- beta1
      beta2.old <- beta2
      beta3.old <- beta3
      objective_function.old <- objective_function
      
      # First step of the block descent : dealing with uncorrupted predictors
      Z3_tilde <- sapply(1:p3,function(j)Z3[,j] / diag(ratio_matrix)[j])
      Xy1 <- 1/n * t(X1) %*% (y - Z2 %*% beta2 - Z3_tilde %*% beta3)
      beta1 <- lasso_covariance(n=n, p=p1, lambda=lambda, XX=sigma1, Xy = Xy1, beta.start = beta1.old, penalty=penalty)$coefficients
      
      
      # Second step of the block descent : dealing with corrupted predictors - additive
      Xy2 <- 1/n * t(Z2) %*% (y - X1 %*% beta1 - Z3_tilde %*% beta3)
      beta2 <- lasso_covariance(n=n, p=p2, lambda=lambda, XX=sigma2, Xy = Xy2, beta.start = beta2.old, penalty=penalty)$coefficients
      
      # Third step of the block descent : dealing with corrupted predictors - missing
      Xy3 <- 1/n * (t(Z3) %*% (y - X1 %*% beta1 - Z2 %*% beta2)) / diag(ratio_matrix)
      beta3 <- lasso_covariance(n=n, p=p3, lambda=lambda, XX=sigma3, Xy = Xy3, beta.start = beta3.old, penalty=penalty)$coefficients
      
      error_beta1[m,1] = sum(abs(beta1 - beta1.old))
      error_beta2[m,1] = sum(abs(beta2 - beta2.old))
      error_beta3[m,1] = sum(abs(beta3 - beta3.old))
      if ((sum(abs(beta1 - beta1.old), na.rm = TRUE) < control$optTol) && (sum(abs(beta2 - beta2.old), na.rm = TRUE) < control$optTol) && (sum(abs(beta3 - beta3.old), na.rm = TRUE) < control$optTol)) {
        break
      }
      
      
      objective_function <- t(beta1) %*% sigma1 %*% beta1 + 
        t(beta2) %*% sigma2 %*% beta2 +
        t(beta3) %*% sigma3 %*% beta3 -
        2 * t(rho1) %*% beta1 - 
        2 * t(rho2) %*% beta2 -
        2 * t(rho3) %*% beta3 + 
        2 * t(beta2) %*% t(Z2) %*% X1 %*% beta1 +
        2 * t(beta3) %*% t(Z3_tilde) %*% X1 %*% beta1 +
        2 * t(beta3) %*% t(Z3_tilde) %*% Z2 %*% beta2
      # if (objective_function > objective_function.old){
      #   print(paste0("Warning at step ",m))
      # }
      m <- m + 1
    }
    #We impose very small coefficients to be equal to zero
    beta1[abs(beta1) < control$zeroThreshold] <- 0
    beta2[abs(beta2) < control$zeroThreshold] <- 0
    beta3[abs(beta3) < control$zeroThreshold] <- 0
    return(list(coefficients.beta1 = beta1, coefficients.beta2 = beta2, coefficients.beta3 = beta3, num.it = m))
  }
  if (p1 == 0) {
    beta2 <- beta2.start
    beta3 <- beta3.start
    objective_function <- 1000
    m <- 1
    error_beta2 <- matrix(0,control$maxIter,1)
    error_beta3 <- matrix(0,control$maxIter,1)
    
    rho2 <- 1/n * (t(Z2) %*% y)
    Z3_tilde <- sapply(1:p3,function(j)Z3[,j] / diag(ratio_matrix)[j])
    rho3 <- 1/n * (t(Z3) %*% y) / diag(ratio_matrix)
    
    while (m < control$maxIter){
      beta2.old <- beta2
      beta3.old <- beta3
      objective_function.old <- objective_function
      
      # First step of the block descent : dealing with uncorrupted predictors
      Z3_tilde <- sapply(1:p3,function(j)Z3[,j] / diag(ratio_matrix)[j])
      
      # Second step of the block descent : dealing with corrupted predictors - additive
      Xy2 <- 1/n * t(Z2) %*% (y - Z3_tilde %*% beta3)
      beta2 <- lasso_covariance(n=n, p=p2, lambda=lambda, XX=sigma2, Xy = Xy2, beta.start = beta2.old, penalty=penalty)$coefficients
      
      # Third step of the block descent : dealing with corrupted predictors - missing
      Xy3 <- 1/n * (t(Z3) %*% (y - Z2 %*% beta2)) / diag(ratio_matrix)
      beta3 <- lasso_covariance(n=n, p=p3, lambda=lambda, XX=sigma3, Xy = Xy3, beta.start = beta3.old, penalty=penalty)$coefficients
      
      error_beta2[m,1] = sum(abs(beta2 - beta2.old))
      error_beta3[m,1] = sum(abs(beta3 - beta3.old))
      if ((sum(abs(beta2 - beta2.old), na.rm = TRUE) < control$optTol) && (sum(abs(beta3 - beta3.old), na.rm = TRUE) < control$optTol)) {
        break
      }
      
      
      objective_function <- t(beta2) %*% sigma2 %*% beta2 +
        t(beta3) %*% sigma3 %*% beta3 -
        2 * t(rho2) %*% beta2 -
        2 * t(rho3) %*% beta3 + 
        2 * t(beta3) %*% t(Z3_tilde) %*% Z2 %*% beta2
      # if (objective_function > objective_function.old){
      #   print(paste0("Warning at step ",m))
      # }
      m <- m + 1
    }
    #We impose very small coefficients to be equal to zero
    beta2[abs(beta2) < control$zeroThreshold] <- 0
    beta3[abs(beta3) < control$zeroThreshold] <- 0
    return(list(coefficients.beta2 = beta2, coefficients.beta3 = beta3, num.it = m))
  }
}
