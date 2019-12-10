lambda_max.coordinate_descent.general <- function(Z,
                                          y,
                                          n,
                                          p,
                                          p1,
                                          p2,
                                          p3,
                                          ratio_matrix=NULL){
  start <- p1 + p2 + 1
  X1 <- Z[,1:(p1+p2)]
  Z2 <- Z[,start:p]
  rho_tilde1 <- 1/n*t(X1)%*%y
  rho_tilde2 <- 1/n*t(Z2)%*%y/ diag(ratio_matrix)
  lambda = max(max(abs(rho_tilde1)),max(abs(rho_tilde2)))
  lambda
}


### FUNCTIONS USED TO PERFORM SCALING OF THE CORRUPTED MATRIX WHEN THERE IS MISSING DATA

scale_manual <- function(j,Z){
  sd <- stats::sd(Z[,j])
  if (sd != 0){
    return(Z[,j]/sd)
  }else{
    return (Z[,j])
  }
  
}

rescale_without_NA_block <- function(j,Z,p1){
  m <- mean(Z[which(!is.na(Z[,p1 + j]), arr.ind = TRUE),p1+j])
  Z[,p1 + j] - m
}

mean_without_NA <- function(j,Z){
  m <- mean(Z[which(!is.na(Z[,j]), arr.ind = TRUE),j])
  m
}

sd_without_NA_block <- function(j,Z){
  sd <- stats::sd(Z[which(!is.na(Z[,j]), arr.ind = TRUE),j])
  sd
}

change_NA_value_block <- function(j,Z,p1){
  Z[which(is.na(Z[,p1 + j]), arr.ind = TRUE),j + p1] <- 0
  Z[,p1 + j]
}

scale_manual_with_sd <- function(j,Z,v){
  sd <- v[j]
  if (sd != 0){
    return(Z[,j]/sd)
  }else{
    return (Z[,j])
  }
}

###CROSS VALIDATION FUNCTION USED IN THE BLOCKWISE COORDINATE DESCENT

cross_validation_function.block_descent_general <- function(k,
                                                    Z,
                                                    y,
                                                    n,
                                                    n_one_fold,
                                                    n_without_fold,
                                                    p,
                                                    p1,
                                                    p2,
                                                    p3,
                                                    folds,
                                                    lambda,
                                                    list_PSD_lasso_additive,
                                                    list_PSD_lasso_missing,
                                                    list_sigma_lasso,
                                                    list_PSD_error_additive,
                                                    list_PSD_error_missing,
                                                    list_sigma_error,
                                                    ratio_matrix=NULL,
                                                    beta1.start,
                                                    beta2.start,
                                                    beta3.start,
                                                    penalty=penalty){
  
  ### Calculating the error for the design matrix without the kth fold
  start = p1 + 1
  #Solving the lasso problem without the kth fold
  if (p1 != 0) {
    sigma_uncorrupted_train <- list_sigma_lasso[[k]] 
  }
  sigma_corrupted_train_additive <- list_PSD_lasso_additive[[k]]
  sigma_corrupted_train_missing <- list_PSD_lasso_missing[[k]]
  
  index <- which(folds==k,arr.ind = TRUE)
  X1_cv_train <- Z[-index,1:p1]
  Z2_cv_train <- Z[-index,start:(p1+p2)]
  Z3_cv_train <- Z[-index,(p1+p2+1):p]
  y_cv_train <- y[-index]
  out = lasso_covariance_block_general(n=n_without_fold,p1=p1,p2=p2,p3=p3,X1=X1_cv_train,Z2=Z2_cv_train,Z3=Z3_cv_train,y=y_cv_train,
                               sigma1=sigma_uncorrupted_train,sigma2=sigma_corrupted_train_additive,sigma3=sigma_corrupted_train_missing,lambda=lambda,
                               ratio_matrix = ratio_matrix,beta1.start = beta1.start, beta2.start = beta2.start, beta3.start = beta3.start,
                               penalty=penalty)
  if (p1 != 0) {
    beta1.lambda <- out$coefficients.beta1
  }
  beta2.lambda <- out$coefficients.beta2
  beta3.lambda <- out$coefficients.beta3
  
  #Calculating the error on the remaining fold
  if (p1 != 0) {
    sigma_uncorrupted_test <- list_sigma_error[[k]]
    sigma_corrupted_test_additive <- list_PSD_error_additive[[k]]
    sigma_corrupted_test_missing <- list_PSD_error_missing[[k]]
    X1_cv_test <- Z[index,1:p1]
    Z2_cv_test <- Z[index,start:(p1+p2)]
    Z3_cv_test <- Z[index,(p1+p2+1):p]
    y_cv_test <- y[index]
    
    rho_1 <- 1/n_one_fold * t(X1_cv_test) %*% y_cv_test
    rho_2 <- 1/n_one_fold * t(Z2_cv_test) %*% y_cv_test
    Z3_tilde <- sapply(1:p3,function(j)Z3_cv_test[,j] / diag(ratio_matrix)[j])
    rho_3 <- 1/n_one_fold * t(Z3_tilde) %*% y_cv_test
    
    sigma21 <- 1/n_one_fold * t(Z2_cv_test) %*% X1_cv_test %*% beta1.lambda
    sigma31 <- 1/n_one_fold * t(Z3_tilde) %*% X1_cv_test %*% beta1.lambda
    sigma32 <- 1/n_one_fold * t(Z3_tilde) %*% Z2_cv_test %*% beta2.lambda
    
    error <- t(beta1.lambda)%*%sigma_uncorrupted_test%*%beta1.lambda + 
      t(beta2.lambda)%*%sigma_corrupted_test_additive%*%beta2.lambda +
      t(beta3.lambda)%*%sigma_corrupted_test_missing%*%beta3.lambda - 
      2*t(rho_1)%*%beta1.lambda - 
      2*t(rho_2)%*%beta2.lambda -
      2*t(rho_3)%*%beta3.lambda + 
      2*t(beta2.lambda)%*%sigma21 +
      2*t(beta3.lambda)%*%sigma31 +
      2*t(beta3.lambda)%*%sigma32
    
    return(error)
  }
  if (p1 == 0) {
    sigma_corrupted_test_additive <- list_PSD_error_additive[[k]]
    sigma_corrupted_test_missing <- list_PSD_error_missing[[k]]
    Z2_cv_test <- Z[index,start:(p1+p2)]
    Z3_cv_test <- Z[index,(p1+p2+1):p]
    y_cv_test <- y[index]
    
    rho_2 <- 1/n_one_fold * t(Z2_cv_test) %*% y_cv_test
    Z3_tilde <- sapply(1:p3,function(j)Z3_cv_test[,j] / diag(ratio_matrix)[j])
    rho_3 <- 1/n_one_fold * t(Z3_tilde) %*% y_cv_test
    
    sigma32 <- 1/n_one_fold * t(Z3_tilde) %*% Z2_cv_test %*% beta2.lambda
    
    error <- t(beta2.lambda)%*%sigma_corrupted_test_additive%*%beta2.lambda +
      t(beta3.lambda)%*%sigma_corrupted_test_missing%*%beta3.lambda - 
      2*t(rho_2)%*%beta2.lambda -
      2*t(rho_3)%*%beta3.lambda + 
      2*t(beta3.lambda)%*%sigma32
    
    return(error)
  }
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
#' @param p2 Number of corrupted predictors containing additive error
#' @param p3 Number of corrupted predictors containing missingness
#' @param center.Z If TRUE, centers Z matrix without taking into account NAs values, and then change NAs to 0 value (in the
#' missing data setting).
#' @param scale.Z If TRUE, divides Z columns by their standard deviation
#' @param center.y If TRUE, centers y
#' @param scale.y If TRUE, divides y by its standard deviation
#' @param lambda.factor Range of the lambda interval we are going to explore
#' @param step Number of values of lambda in the interval we are going to test
#' @param K Number of folds for the cross-validation
#' @param mu Penalty parameter for the ADMM algorithm
#' @param tau Standard deviation for the additive error matrix in the additive error setting (NULL in the missing data setting)
#' @param etol Tolerance parameter for the ADMM algorithm
#' @param optTol Tolerance parameter for the convergence of the error in the pathwise coordinate descent
#' @param earlyStopping_max Number of iterations allowed when error starts increasing
#' @param penalty Type of penalty used : can be lasso penalty or SCAD penalty
#' 
#' @return list containing \itemize{
#' \item \code{lambda.opt} optimal value of lambda corresponding to minimum error
#' \item \code{lambda.sd} Value of lambda corresponding to error higher than minimum error by one standard deviation
#' \item \code{beta.opt} Value of beta corresponding to \code{lambda.opt}
#' \item \code{beta.sd} Value of beta corresponding to \code{lambda.sd}
#' \item \code{data_error} Dataframe containing errors and their standard deviation for each iteration of the algorithm
#' \item \code{data_beta} Dataframe containing the values of beta for each iteration of the algorithm
#' \item \code{earlyStopping} Integer containing the value of iteration when early stopping happens
#' \item \code{vnames} Column names of Z matrix
#' \item \code{mean.Z} Mean of Z matrix without the NAs values
#' \item \code{sd.Z} Standard deviation of Z matrix without the NAs values
#' \item \code{mean.y} Mean of y matrix
#' \item \code{sd.y} Standard deviation of y matrix
#' }
#' 
#' @details It is highly recommended to use center.Z = TRUE for the algorithm to work in the case of missing data. Otherwise 
#' the algorithm will not work the way it is intended.
#' It is recommended to use center.Z = TRUE, scale.Z = TRUE, center.y = TRUE and scale.y = TRUE for both convergence
#' and interpretability reasons. The use of center.Z = TRUE in the additive error setting can be subject to discussion,
#' as it may introduce bias in the algorithm.
#' For computing speed reasons, if model is not converging or running slow, consider changing \code{mu}, decreasing
#' \code{etol} or \code{optTol} or decreasing \code{earlyStopping_max}
#' 
#' @seealso \url{https://arxiv.org/pdf/1510.07123.pdf}
#' 
#' @export
#' 

blockwise_coordinate_descent_general <- function(Z,
                                         y,
                                         n,
                                         p,
                                         p1,
                                         p2,
                                         p3,
                                         center.Z = TRUE,
                                         scale.Z = TRUE,
                                         center.y = TRUE,
                                         scale.y = TRUE,
                                         lambda.factor=ifelse(dim(Z)[1] < dim(Z)[2],0.01,0.001),
                                         step=100,
                                         K=4,
                                         mu=10,
                                         tau = NULL,
                                         etol= 1e-4,
                                         optTol = 1e-5,
                                         earlyStopping_max = 10,
                                         penalty=c("lasso","SCAD")){
  
  nrows = nrow(Z)
  ncols = ncol(Z)
  vnames = colnames(Z)
  
  if(!(is.matrix(Z))){
    stop("Z has to be a matrix")
  }
  if(!(is.matrix(y))){
    stop("y has to be a matrix")
  }
  if(!is.null(tau) & !is.numeric(tau)){
    stop("tau must be numeric")
  }
  if(n != nrows){
    stop(paste("Number of rows in Z (", nrows, ") different from n(", n, ")"),sep="")
  }
  if(p != ncols){
    stop(paste("Number of columns in Z (", ncols, ") different from p (", p, ")"),sep="")
  }
  if (nrows != dim(y)[1]){
    stop(paste("Number of rows in Z (", nrows, ") different from number of rows in y (", dim(y)[1], ") "),sep="")
  }
  if (!is.numeric(y)) {
    stop("The response y must be numeric. Factors must be converted to numeric")
  }
  if(lambda.factor >= 1){
    stop("lambda factor should be smaller than 1")
  }
  if(p1 + p2 + p3 != p){
    stop(paste("Sum of p1, p2 and p3 (", p1 + p2 + p3, ") should be equal to p (", p, ")"),sep="")
  }
  if(mu>500 || mu<1){
    warning(paste("Mu value (", mu, ") is not in the usual range (10-500)"))
  }
  if (n %% K != 0){
    stop("K should be a divider of n")
  }
  if(center.Z == FALSE){
    stop("When the data contain missingness, it is required to center matrix Z. Use center.Z=TRUE.")
  }
  if(scale.Z == FALSE){
    warning("When the data contain missingness, it is recommended to use scale.Z equal to TRUE in order 
            to obtain trustworthy results.")
  }
  
  #General variables we are going to use in the function
  start = 1 + p1 + p2
  n_without_fold = n - floor(n/K)
  n_one_fold = floor(n/K)
  earlyStopping = step
  
  ratio_matrix = NULL
  ratio_matrix = matrix(0,p3,p3)
    
  for (i in 1:p3){
    for (j in i:p3){
      n_ij = length(intersect(which(!is.na(Z[,p1 + p2 + i])),which(!is.na(Z[,p1 + p2 + j]))))
      ratio_matrix[i,j] = n_ij
      ratio_matrix[j,i] = n_ij
    }
  }
  ratio_matrix = ratio_matrix/n
  
  
  mean.Z = sapply(1:p, function(j)mean_without_NA(j,Z))
  sd.Z <- sapply(1:p, function(j)sd_without_NA_block(j,Z))
  
  
  
  if (center.Z == TRUE){
    if (scale.Z == TRUE){
      Z[,start:p] = sapply(1:p3, function(j)rescale_without_NA_block(j,Z,p1+p2))
      Z[,start:p] = sapply(1:p3, function(j)change_NA_value_block(j,Z,p1+p2))
      Z[,1:(p1+p2)] = scale(Z[,1:(p1+p2)], center = TRUE, scale = FALSE)
      if (p1 == 0) {
        Z[,c((p1+p2+1):p)] = sapply(1:(p - p2), function(j)scale_manual_with_sd(j,Z[,c((p1+p2+1):p)],sd.Z[c((p1+p2+1):p)]))
      }
      if (p1 != 0) {
        Z[,c(1:p1,(p1+p2+1):p)] = sapply(1:(p - p2), function(j)scale_manual_with_sd(j,Z[,c(1:p1,(p1+p2+1):p)],sd.Z[c(1:p1,(p1+p2+1):p)]))
      }
    }else{
      Z[,start:p] = sapply(1:p3, function(j)rescale_without_NA_block(j,Z,p1+p2))
      Z[,start:p] = sapply(1:p3, function(j)change_NA_value_block(j,Z,p1+p2))
      Z[,1:(p1+p2)] = scale(Z[,1:(p1+p2)], center = TRUE, scale = FALSE)
    }
  }else{
    if (scale.Z == TRUE){
      Z = sapply(1:p, function(j)scale_manual_with_sd(j,Z,sd.Z))
    }
  }
  
  mean.y = mean(y)
  sd.y = stats::sd(y)
  
  if (center.y == TRUE){
    if (scale.y == TRUE){
      y = scale(y, center = TRUE, scale = TRUE)
    }else{
      y = scale(y, center = TRUE, scale = FALSE)
    }
  }
  
  lambda_max <- lambda_max.coordinate_descent.general(Z=Z,y=y,n=n,p=p,p1=p1,p2=p2,p3=p3,ratio_matrix=ratio_matrix)
  lambda_min <- lambda.factor*lambda_max
  lambda_list <- emdbook::lseq(lambda_max,lambda_min,step)
  beta1.start <- rep(0,p1)
  beta2.start <- rep(0,p2)
  beta3.start <- rep(0,p3)
  beta.start <- c(beta1.start,beta2.start,beta3.start)
  best.lambda <- lambda_max
  beta.opt <- beta.start
  best.error <- 10000
  error_list <- matrix(0,step,4)
  error <- 0
  earlyStopping_high = 0
  
  matrix_beta <- matrix(0,step,p)
  
  ### Creating the K matrices we are going to use for cross validation
  output = cv_covariance_matrices_block_descent_general(K=K, mat=Z, y=y, p=p, p1=p1, p2=p2, p3=p3, mu=mu, 
                                                tau=tau, ratio_matrix = ratio_matrix, etol=etol)
  list_PSD_lasso_additive = output$list_PSD_lasso_additive
  list_PSD_error_additive = output$list_PSD_error_additive
  list_PSD_lasso_missing = output$list_PSD_lasso_missing
  list_PSD_error_missing = output$list_PSD_error_missing
  if (p1!=0) {
    list_sigma_lasso = output$list_sigma_lasso
    list_sigma_error = output$list_sigma_error
    sigma1 = output$sigma_global_uncorrupted
  }
  if (p1==0) {
    list_sigma_lasso = NULL
    list_sigma_error = NULL
    sigma1 = NULL
  }
  sigma2 = output$sigma_global_corrupted_additive
  sigma3 = output$sigma_global_corrupted_missing
  folds = output$folds
  X1 = Z[,1:p1]
  Z2 = Z[,(p1+1):(p1+p2)]
  Z3 = Z[,(p1+p2+1):p]
  
  for (i in 1:step){
    print(paste0("iteration: ",i))
    
    lambda_step <- lambda_list[i]
    error_old <- error
    
    out = sapply(1:K, function(k)cross_validation_function.block_descent_general(k,
                                                                         Z,
                                                                         y,
                                                                         n,
                                                                         n_one_fold,
                                                                         n_without_fold,
                                                                         p,
                                                                         p1,
                                                                         p2,
                                                                         p3,
                                                                         folds,
                                                                         lambda_step,
                                                                         list_PSD_lasso_additive,
                                                                         list_PSD_lasso_missing,
                                                                         list_sigma_lasso,
                                                                         list_PSD_error_additive,
                                                                         list_PSD_error_missing,
                                                                         list_sigma_error,
                                                                         ratio_matrix = ratio_matrix,
                                                                         beta1.start,
                                                                         beta2.start,
                                                                         beta3.start,
                                                                         penalty=penalty))
    error = mean(out)
    sd_low = stats::quantile(out, probs = c(0.1))
    sd_high = stats::quantile(out, probs = c(0.9))
    error_list[i,1] <- error
    error_list[i,2] <- sd_low
    error_list[i,3] <- sd_high
    error_list[i,4] <- stats::sd(out)
    out = lasso_covariance_block_general(n=n, p1=p1, p2=p2, p3=p3, X1=X1, Z2=Z2, Z3=Z3, y=y, sigma1=sigma1, sigma2=sigma2, sigma3=sigma3, lambda=lambda_step, 
                                 ratio_matrix = ratio_matrix, beta1.start = beta1.start, beta2.start = beta2.start, beta3.start = beta3.start,
                                 penalty=penalty)
    if (p1 != 0) {
      beta1 <- out$coefficients.beta1
    }
    beta2 <- out$coefficients.beta2
    beta3 <- out$coefficients.beta3
    if (p1 != 0) {
      beta <- c(beta1,beta2,beta3)
    }
    if (p1 == 0) {
      beta <- c(beta2,beta3)
    }
    
    if (p1 != 0) {
      beta1.start <- beta1
    }
    if (p1 == 0) {
      beta1.start <- rep(0,p1)
    }
    beta2.start <- beta2
    beta3.start <- beta3
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
  df <- data.frame(lambda=lambda_list[1:earlyStopping], error=error_list[1:earlyStopping,1],error.inf=error_list[1:earlyStopping,2],error.sup=error_list[1:earlyStopping,3],error.sd=error_list[1:earlyStopping,4])
  step.min <- which(df[,"error"] == best.error)
  sd.best <- df[step.min,"error.sd"]
  step.sd <- max(which(df[,"error"] > best.error + sd.best & df[,"lambda"] > df[step.min,"lambda"]))
  lambda.sd <- df[step.sd,"lambda"]
  beta.sd <- matrix_beta[step.sd,]
  
  
  data_intermediate <- data.frame(matrix_beta[1:earlyStopping,])
  names(data_intermediate) <- sapply(1:p, function(i)paste0("beta",i))
  
  data_beta <- data.frame(lambda = lambda_list[1:earlyStopping])
  
  data_beta <- cbind(data_beta,data_intermediate)
  
  fit = list(
    lambda.opt = best.lambda,
    lambda.sd = lambda.sd,
    beta.opt = beta.opt,
    beta.sd = beta.sd,
    data_error = df,
    data_beta = data_beta,
    earlyStopping = earlyStopping,
    vnames = vnames,
    mean.Z = mean.Z,
    sd.Z = sd.Z,
    mean.y = mean.y,
    sd.y = sd.y                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
  )
  return(fit)
}

