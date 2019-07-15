#' Lasso in covariance form for the BD-CoCoLasso
#'
#' Solve the least squares loss with lasso penalty written in a form with the covariance matrix : \eqn{\frac{1}{2} \beta^{'} \Sigma \beta - \rho^{'} \beta + \lambda \|\beta\|_1}
#'
#' @param n Number of samples of the design matrix
#' @param p1 Number of uncorrupted predictors
#' @param p2 Number of corrupted predictors
#' @param X1 first block of the design matrix corresponding to uncorrupted features
#' @param Z2 second block of the design matrix corresponding to corrupted features
#' @param y Response vector
#' @param sigma1 Covariance matrix for X1 : \eqn{\frac{1}{n} X_1'X_1}. This parameter is automatically furnished in \link{blockwise_coordinate_descent}
#' @param sigma2 Modified covariance matrix for Z2 through the CoCoLasso algorithm. This parameter is automatically furnished in \link{blockwise_coordinate_descent}
#' @param lambda Penalty parameter
#' @param noise Type of noise for Z2 : additive or missing
#' @param ratio_matrix Observation matrix in the missing data setting (NULL otherwise)
#' @param control Including control parameters : max of iterations, tolerance for the convergence of the error, zero threshold to put to zero small beta coefficients
#' @param beta1.start Initial value for the coefficients of uncorrupted features
#' @param beta2.start Initial value for the coefficients of corrupted features
#' 
#' @return list containing \itemize{
#' \item coefficients.beta1 : Coefficients corresponding to final beta1 after convergence of the algoritm
#' \item coefficients.beta2 : Coefficients corresponding to final beta2 after convergence of the algoritm
#' \item num.it : Number of iterations of algorithm
#' }
#' 
#' 
#' @export

lasso_covariance_block <- function(n,
                                   p1,
                                   p2,
                                   X1,
                                   Z2,
                                   y,
                                   sigma1,
                                   sigma2,
                                   lambda,
                                   noise = c("additive","missing","HM"),
                                   ratio_matrix = NULL,
                                   control = list(maxIter = 1000,
                                                  optTol = 10^(-5), 
                                                  zeroThreshold = 10^(-6)),
                                   beta1.start,
                                   beta2.start){

  
  beta1 <- beta1.start
  beta2 <- beta2.start
  objective_function <- 1000
  m <- 1
  error_beta1 <- matrix(0,control$maxIter,1)
  error_beta2 <- matrix(0,control$maxIter,1)
  if (noise=="additive"){
    Z2_tilde <- Z2
    rho1 <- 1/n * (t(X1) %*% y)
    rho2 <- 1/n * (t(Z2) %*% y) 
  }
  else if (noise=="missing"){
    Z2_tilde <- sapply(1:p2,function(j)Z2[,j] / diag(ratio_matrix)[j])
    rho1 <- 1/n * (t(X1) %*% y)
    rho2 <- 1/n * (t(Z2) %*% y) / diag(ratio_matrix)
  }
  
  
  while (m < control$maxIter){
    beta1.old <- beta1
    beta2.old <- beta2
    objective_function.old <- objective_function
    
    # First step of the block descent : dealing with uncorrupted predictors
    if (noise == "additive"){
      Xy1 <- 1/n * t(X1) %*% (y - Z2 %*% beta2)
    } else if ((noise == "missing") || (noise == "HM")){
      # Xy1 <- 1/n * t(X1) %*% (y - Z2 %*% beta2)
      Z2_tilde <- sapply(1:p2,function(j)Z2[,j] / diag(ratio_matrix)[j])
      Xy1 <- 1/n * t(X1) %*% (y - Z2_tilde %*% beta2 )
    }
    # Xy1 <- 1/n * t(X1) %*% (y - Z2 %*% beta2)
    beta1 <- lasso_covariance(n=n, p=p1, lambda=lambda, XX=sigma1, Xy = Xy1, beta.start = beta1.old, penalty="lasso")$coefficients
    
    
    # Second step of the block descent : dealing with corrupted predictors
    if (noise == "additive"){
      Xy2 <- 1/n * t(Z2) %*% (y - X1 %*% beta1)
    } else if ((noise == "missing") || (noise == "HM")){
      Xy2 <- 1/n * (t(Z2) %*% (y - X1 %*% beta1)) / diag(ratio_matrix)
    }
    beta2 <- lasso_covariance(n=n, p=p2, lambda=lambda, XX=sigma2, Xy = Xy2, beta.start = beta2.old, penalty="lasso")$coefficients
    
    error_beta1[m,1] = sum(abs(beta1 - beta1.old))
    error_beta2[m,1] = sum(abs(beta2 - beta2.old))
    if ((sum(abs(beta1 - beta1.old), na.rm = TRUE) < control$optTol) && (sum(abs(beta2 - beta2.old), na.rm = TRUE) < control$optTol)) {
      break
    }
    
    
    objective_function <- t(beta1) %*% sigma1 %*% beta1 + t(beta2) %*% sigma2 %*% beta2 - 2 * t(rho1) %*% beta1 - 2 * t(rho2) %*% beta2 + 2 * t(beta2) %*% t(Z2_tilde) %*% X1 %*% beta1
    # if (objective_function > objective_function.old){
    #   print(paste0("Warning at step ",m))
    # }
    m <- m + 1
  }
  #We impose very small coefficients to be equal to zero
  beta1[abs(beta1) < control$zeroThreshold] <- 0
  beta2[abs(beta2) < control$zeroThreshold] <- 0
  return(list(coefficients.beta1 = beta1, coefficients.beta2 = beta2, num.it = m))
}