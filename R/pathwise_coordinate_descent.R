lambda_max <- function(Z,y,n,ratio_matrix=NULL,noise=c("additive","missing")){
  
  if (noise == "additive"){
    rho_tilde <- 1/n*t(Z)%*%y
  } else if ((noise=="missing") || (noise=="HM")){
    #rho_tilde <- 1/n*t(Z)%*%y/ (1 - probs)
    rho_tilde <- 1/n*t(Z)%*%y/ diag(ratio_matrix)
  }
  max(abs(rho_tilde))
}

rescale_without_NA <- function(j,Z){
  m <- mean(Z[which(!is.na(Z[,j]), arr.ind = TRUE),j])
  Z[,j] - m
}

mean_without_NA <- function(j,Z){
  m <- mean(Z[which(!is.na(Z[,j]), arr.ind = TRUE),j])
  m
}

sd_without_NA_block <- function(j,Z){
  sd <- stats::sd(Z[which(!is.na(Z[,j]), arr.ind = TRUE),j])
  sd
}


change_NA_value <- function(j,Z){
  Z[which(is.na(Z[,j]), arr.ind = TRUE),j] <- 0
  Z[,j]
}

scale_manual_with_sd <- function(j,Z,v){
  sd <- v[j]
  if (sd != 0){
    return(Z[,j]/sd)
  }else{
    return (Z[,j])
  }
}

cross_validation_function <- function(k,
                                      n,
                                      p,
                                      lambda_step,
                                      list_matrices_lasso,
                                      list_rho_lasso,
                                      list_matrices_error,
                                      list_rho_error,
                                      beta_start,
                                      penalty=c("lasso","SCAD")){
  
  #Solving the lasso problem without the kth fold
  sigma_train <- list_matrices_lasso[[k]]
  rho_train <- list_rho_lasso[[k]] 
  coef_lambda = lasso_covariance(n=n, p=p,lambda=lambda_step, XX=sigma_train, Xy=rho_train, beta.start=beta_start, penalty=penalty)$coefficients
  
  #Calculating the error on the remaining fold
  sigma_test <- list_matrices_error[[k]]
  rho_test <- list_rho_error[[k]]
  error <- t(coef_lambda)%*%sigma_test%*%coef_lambda - 2*t(rho_test)%*%coef_lambda
  
  error
}


#' Pathwise coordinate descent
#'
#' Implement pathwise coordinate descent algorithm. Best penalty value is evaluated through cross-validation.
#'
#' @param Z Corrupted design matrix (with additive error or missing data)
#' @param y Response vector
#' @param n Number of samples of the design matrix
#' @param p Number of features of the matrix
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
#' @param noise Type of noise (additive or missing)
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
#' \item \code{vnames} Names of the features
#' \item \code{mean.Z} Mean of Z matrix without the NAs values
#' \item \code{sd.Z} Standard deviation of Z matrix without the NAs values
#' \item \code{mean.y} Mean of y matrix
#' \item \code{sd.y} Standard deviation of y matrix
#' }
#' 
#' @details It is highly recommended to use center.Z = TRUE for the algorithm to work in the case of missing data. 
#' It is recommended to use center.Z = TRUE, scale.Z = TRUE, center.y = TRUE and scale.y = TRUE for both convergence
#' and interpretability reasons. The use of center.Z = TRUE in the additive error setting can be subject to discussion,
#' as it may introduce bias in the algorithm.
#' For computing speed reasons, if model is not converging or running slow, consider changing \code{mu}, decreasing
#' \code{etol} or \code{optTol} or decreasing \code{earlyStopping_max}
#'   
#' 
#' @seealso \url{https://arxiv.org/pdf/1510.07123.pdf}
#' 
#' @export
#' 

pathwise_coordinate_descent <- function(Z,
                                        y,
                                        n,
                                        p,
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
                                        optTol = 1e-10,
                                        earlyStopping_max = 10,
                                        noise=c("additive","missing"),
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
  if (any(is.na(y))){
    stop("The response contains NA values. Remove NA values before calling the function.")
  }
  if(lambda.factor >= 1){
    stop("lambda factor should be smaller than 1")
  }
  if (n %% K != 0){
    stop("K should be a divider of n")
  }
  if(mu>500 || mu<1){
    warning(paste("Mu value (", mu, ") is not in the usual range (10-500)"))
  }
  if(noise=="missing" && center.Z == FALSE){
    stop("When noise is equal to missing, it is required to center matrix Z. Use center.Z=TRUE.")
  }
  if(scale.Z == FALSE){
    warning("Is it recommended to use scale.Z equal to TRUE in order to obtain trustworthy results.")
  }
  
  ratio_matrix = NULL
  if (noise=="missing"){
    ratio_matrix = matrix(0,p,p)
    
    for (i in 1:p){
      for (j in i:p){
        n_ij = length(intersect(which(!is.na(Z[,i])),which(!is.na(Z[,j]))))
        ratio_matrix[i,j] = n_ij
        ratio_matrix[j,i] = n_ij
      }
    }
    ratio_matrix = ratio_matrix/n
  }
  
  mean.Z = sapply(1:p, function(j)mean_without_NA(j,Z))
  sd.Z <- sapply(1:p, function(j)sd_without_NA_block(j,Z))
  
  if (center.Z == TRUE){
    if (scale.Z == TRUE){
      Z = sapply(1:p, function(j)rescale_without_NA(j,Z))
      Z = sapply(1:p, function(j)change_NA_value(j,Z))
      Z = sapply(1:p, function(j)scale_manual_with_sd(j,Z,sd.Z))
    }else{
      Z = sapply(1:p, function(j)rescale_without_NA(j,Z))
      Z = sapply(1:p, function(j)change_NA_value(j,Z))
    }
  }else{
    if(scale.Z == TRUE){
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
  }else{
    if (scale.y == TRUE){
      y = scale(y, center = FALSE, scale = TRUE)
    }
  }
  

  
  n_without_fold = n - floor(n/K)
  n_one_fold = floor(n/K)
  earlyStopping = step

  lambda_max <- lambda_max(Z=Z,y=y,n=n,ratio_matrix=ratio_matrix,noise=noise)
  lambda_min <- lambda.factor*lambda_max
  lambda_list <- emdbook::lseq(lambda_max,lambda_min,step)
  beta_start <- rep(0,p)
  best.lambda <- lambda_max
  beta.opt <- beta_start
  best.error <- 1000
  error_list <- matrix(0,step,4)
  error <- 1000
  earlyStopping_high = 0
  
  matrix_beta <- matrix(0,step,p)
  
  ### Creating the K matrices we are going to use for cross validation
  output = cv_covariance_matrices(K=K, mat=Z, y=y, p=p, mu=mu, tau=tau, ratio_matrix = ratio_matrix, etol=etol, noise = noise)
  list_matrices_lasso = output$list_matrices_lasso
  list_matrices_error = output$list_matrices_error
  list_rho_lasso = output$list_rho_lasso
  list_rho_error = output$list_rho_error
  ZZ = output$sigma_global
  Zy = output$rho_global
  for (i in 1:step){
    
    lambda_step <- lambda_list[i]
    error_old <- error
    error <- 0
    
    out = sapply(1:K, function(k)cross_validation_function(k,
                                                           n,
                                                           p,
                                                           lambda_step,
                                                           list_matrices_lasso,
                                                           list_rho_lasso,
                                                           list_matrices_error,
                                                           list_rho_error,
                                                           beta_start,
                                                           penalty=penalty))
    error = mean(out)
    sd_low = stats::quantile(out, probs = c(0.1))
    sd_high = stats::quantile(out, probs = c(0.9))
    error_list[i,1] <- error
    error_list[i,2] <- sd_low
    error_list[i,3] <- sd_high
    error_list[i,4] <- stats::sd(out)
    coef_tot = lasso_covariance(n=n, p=p, lambda=lambda_step, XX=ZZ, Xy=Zy, beta.start = beta_start, penalty=penalty)$coefficients
    beta_start <- coef_tot
    matrix_beta[i,] <- beta_start
    
    ### Checking for optimal parameters
    if (error <= best.error){
      best.error <- error
      best.lambda <- lambda_step
      beta.opt <- coef_tot
    }
    
    ## Early stopping
    if (abs(error - error_old) < optTol){
      print("Early Stopping because of convergence of the error")
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
  fit <- list(
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
