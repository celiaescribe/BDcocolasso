lambda_max.coordinate_descent <- function(Z,
                                          y,
                                          n,
                                          p,
                                          p1,
                                          p2,
                                          ratio_matrix=NULL,
                                          noise=c("additive","missing","HM")){
  start <- p1 + 1
  X1 <- Z[,1:p1]
  Z2 <- Z[,start:p]
  if (noise=="additive"){
    rho_tilde <- 1/n*t(Z)%*%y
    lambda <- max(abs(rho_tilde))
  }else if ((noise=="missing") || (noise=="HM")){
    rho_tilde1 <- 1/n*t(X1)%*%y
    rho_tilde2 <- 1/n*t(Z2)%*%y/ diag(ratio_matrix)
    lambda = max(max(abs(rho_tilde1)),max(abs(rho_tilde2)))
  }
  lambda
}

cross_validation_function.block_descent <- function(k,
                                                    Z,
                                                    y,
                                                    n,
                                                    n_one_fold,
                                                    n_without_fold,
                                                    p,
                                                    p1,
                                                    p2,
                                                    folds,
                                                    lambda,
                                                    list_PSD_lasso,
                                                    list_sigma_lasso,
                                                    list_PSD_error,
                                                    list_sigma_error,
                                                    ratio_matrix=NULL,
                                                    beta1.start,
                                                    beta2.start,
                                                    noise=c("missing","additive")){
  
  ### Calculating the error for the design matrix without the kth fold
  start = p1 + 1
  #Solving the lasso problem without the kth fold
  sigma_corrupted_train <- list_PSD_lasso[[k]]
  sigma_uncorrupted_train <- list_sigma_lasso[[k]] 
  index <- which(folds==k,arr.ind = TRUE)
  X1_cv_train <- Z[-index,1:p1]
  Z2_cv_train <- Z[-index,start:p]
  y_cv_train <- y[-index]
  out = lasso_covariance_block(n=n_without_fold,p1=p1,p2=p2,X1=X1_cv_train,Z2=Z2_cv_train,y=y_cv_train,
                               sigma1=sigma_uncorrupted_train,sigma2=sigma_corrupted_train,lambda=lambda,
                               noise=noise,ratio_matrix = ratio_matrix,beta1.start = beta1.start, beta2.start = beta2.start)
  beta1.lambda <- out$coefficients.beta1
  beta2.lambda <- out$coefficients.beta2
  
  #Calculating the error on the remaining fold
  sigma_corrupted_test <- list_PSD_error[[k]]
  sigma_uncorrupted_test <- list_sigma_error[[k]]
  X1_cv_test <- Z[index,1:p1]
  Z2_cv_test <- Z[index,start:p]
  y_cv_test <- y[index]
  
  rho_1 <- 1/n_one_fold * t(X1_cv_test) %*% y_cv_test
  if (noise == "additive"){
    rho_2 <- 1/ n_one_fold * t(Z2_cv_test) %*% y_cv_test 
    sigma3 <- 1/ n_one_fold * t(Z2_cv_test) %*% X1_cv_test %*% beta1.lambda
  }else if ((noise== "missing") || (noise=="HM")){
    # rho_2 <- 1/ n_one_fold * (t(Z2_cv_test) %*% y_cv_test) / diag(ratio_matrix)
    # sigma3 <- 1/ n_one_fold * (t(Z2_cv_test) %*% X1_cv_test %*% beta1.lambda) / diag(ratio_matrix)
    Z2_tilde <- sapply(1:p2,function(j)Z2_cv_test[,j] / diag(ratio_matrix)[j])
    rho_2 <- 1/ n_one_fold * (t(Z2_tilde) %*% y_cv_test) 
    sigma3 <- 1/ n_one_fold * (t(Z2_tilde) %*% X1_cv_test %*% beta1.lambda) 
  }
  error <- t(beta1.lambda)%*%sigma_uncorrupted_test%*%beta1.lambda + t(beta2.lambda)%*%sigma_corrupted_test%*%beta2.lambda - 2*t(rho_1)%*%beta1.lambda - 2*t(rho_2)%*%beta2.lambda + 2*t(beta2.lambda)%*%sigma3
  
  error
}

#' Blockwise coordinate descent
#'
#' Implement blockwise coordinate descent algorithm. Best penalty value is evaluated through cross-validation.
#'
#' @param Z Corrupted design matrix (with additive error or missing data)
#' @param y Response vector
#' @param n Number of samples of the design matrix
#' @param p Number of features of the matrix
#' @param p1 Number of uncorrupted predictors
#' @param p2 Number of corrupted predictors
#' @param lambda.factor Range of the lambda interval we are going to explore
#' @param step Number of values of lambda in the interval we are going to test
#' @param K Number of folds for the cross-validation
#' @param mu Penalty parameter for the ADMM algorithm
#' @param tau Standard deviation for the additive error matrix in the additive error setting (NULL in the missing data setting)
#' @param ratio_matrix Observation matrix used in the missing data setting
#' @param etol Tolerance parameter for the ADMM algorithm
#' @param optTol Tolerance parameter for the convergence of the error in the pathwise coordinate descent
#' @param earlyStopping_max Number of iterations allowed when error starts increasing
#' @param noise Type of noise (additive or missing)
#' 
#' @return list containing \itemize{
#' \item \code{lambda.opt} optimal value of lambda corresponding to minimum error
#' \item \code{lambda.sd} Value of lambda corresponding to error higher than minimum error by one standard deviation
#' \item \code{beta.opt} Value of beta corresponding to \code{lambda.opt}
#' \item \code{beta.sd} Value of beta corresponding to \code{lambda.sd}
#' \item \code{data_error} Dataframe containing errors and their standard deviation for each iteration of the algorithm
#' \item \code{data_beta} Dataframe containing the values of beta for each iteration of the algorithm
#' \item \code{earlyStopping} Integer containing the value of iteration when early stopping happens
#' }
#' 
#' @example 
#' 
#' @export
#' 

blockwise_coordinate_descent <- function(Z,
                                         y,
                                         n,
                                         p,
                                         p1,
                                         p2,
                                         lambda.factor=ifelse(dim(Z)[1] < dim(Z)[2],0.01,0.001),
                                         step=100,
                                         K=4,
                                         mu=10,
                                         tau = NULL,
                                         ratio_matrix = NULL,
                                         etol= 1e-4,
                                         optTol = 1e-5,
                                         earlyStopping_max = 10,
                                         noise=c("additive","missing","HM")){
  
  

  #General variables we are going to use in the function
  start = 1 + p1
  n_without_fold = n - floor(n/K)
  n_one_fold = floor(n/K)
  earlyStopping = step
  
  lambda_max <- lambda_max.coordinate_descent(Z,y,n,p1,p2,ratio_matrix,noise)
  lambda_min <- lambda.factor*lambda_max
  lambda_list <- emdbook::lseq(lambda_max,lambda_min,step)
  beta1.start <- rep(0,p1)
  beta2.start <- rep(0,p2)
  beta.start <- c(beta1.start,beta2.start)
  best.lambda <- lambda_max
  beta.opt <- beta.start
  best.error <- 10000
  error_list <- matrix(0,step,4)
  error <- 0
  earlyStopping_high = 0
  
  matrix_beta <- matrix(0,step,p)
  
  ### Creating the K matrices we are going to use for cross validation
  output = cv_covariance_matrices_block_descent(K=K, mat=Z, y=y, p=p, p1=p1, p2=p2, mu=mu, 
                                                tau=tau, ratio_matrix = ratio_matrix, etol=etol, noise = noise)
  list_PSD_lasso = output$list_PSD_lasso
  list_PSD_error = output$list_PSD_error
  list_sigma_lasso = output$list_sigma_lasso
  list_sigma_error = output$list_sigma_error
  sigma1 = output$sigma_global_uncorrupted
  sigma2 = output$sigma_global_corrupted
  folds = output$folds
  X1 = Z[,1:p1]
  Z2 = Z[,start:p]
  
  for (i in 1:step){
    
    lambda_step <- lambda_list[i]
    error_old <- error
    
    out = sapply(1:K, function(k)cross_validation_function.block_descent(k,
                                                                         Z,
                                                                         y,
                                                                         n,
                                                                         n_one_fold,
                                                                         n_without_fold,
                                                                         p,
                                                                         p1,
                                                                         p2,
                                                                         folds,
                                                                         lambda_step,
                                                                         list_PSD_lasso,
                                                                         list_sigma_lasso,
                                                                         list_PSD_error,
                                                                         list_sigma_error,
                                                                         ratio_matrix = ratio_matrix,
                                                                         beta1.start,
                                                                         beta2.start,
                                                                         noise=noise))
    error = mean(out)
    sd_low = stats::quantile(out, probs = c(0.1))
    sd_high = stats::quantile(out, probs = c(0.9))
    error_list[i,1] <- error
    error_list[i,2] <- sd_low
    error_list[i,3] <- sd_high
    error_list[i,4] <- stats::sd(out)
    out = lasso_covariance_block(n=n, p1=p1, p2=p2, X1=X1, Z2=Z2, y=y, sigma1=sigma1, sigma2=sigma2, lambda=lambda_step, 
                                 noise=noise, ratio_matrix = ratio_matrix, beta1.start = beta1.start, beta2.start = beta2.start)
    beta1 <- out$coefficients.beta1
    beta2 <- out$coefficients.beta2
    beta <- c(beta1,beta2)
    
    beta1.start <- beta1
    beta2.start <- beta2
    matrix_beta[i,] <- beta
    
    ### Checking for optimal parameters
    if (error <= best.error){
      best.error <- error
      best.lambda <- lambda_step
      beta.opt <- beta
    }
    
    ## Early stopping
    if (abs(error - error_old) < optTol){
      print("Early Stopping because of convergence of the error")
      earlyStopping = i
      break
    }
    
    # Trying to avoid getting to too high lambda values, leading to too time consuming steps
    if (i>= step/2 && error >= error_list[1]){
      print("Value of lambda yielding too high error : we exit the loop")
      earlyStopping = i
      break
    }
    if (error > best.error){
      earlyStopping_high = earlyStopping_high +1
      if (earlyStopping_high >= earlyStopping_max){
        print("Early stopping because of error getting too high")
        earlyStopping = i
        break
      }
    }
  }
  df <- data.frame(lambda=lambda_list[1:earlyStopping], error=error_list[1:earlyStopping,])
  step.min <- which(df[,"error.1"] == best.error)
  sd.best <- df[step.min,"error.4"]
  step.sd <- max(which(df[,"error.1"] > best.error + sd.best & df[,"lambda"] > df[step.min,"lambda"]))
  lambda.sd <- df[step.sd,"lambda"]
  beta.sd <- matrix_beta[step.sd,]
  
  
  data_intermediate <- data.frame(matrix_beta[1:earlyStopping,])
  names(data_intermediate) <- sapply(1:p, function(i)paste0("beta",i))
  
  data_beta <- data.frame(lambda = lambda_list[1:earlyStopping])
  
  data_beta <- cbind(data_beta,data_intermediate)
  return(list(lambda.opt = best.lambda, lambda.sd=lambda.sd, beta.opt = beta.opt, beta.sd=beta.sd, data_error=df, data_beta = data_beta, earlyStopping = earlyStopping))
}