#' Projected covariance matrices for BD-CoCoLasso
#'
#' Creating projected nearest positive semi-definite covariance matrices for the cross validation step of the BD-CoCoLasso algorithm. In that case, the design matrix must be organized as follozs : uncorrupted features must be the first block of the matrix, and corrupted features must be the second block of the matrix.
#'
#' @param K Number of folds for the cross validation
#' @param mat Covariance matrix to be projected
#' @param y Response vector
#' @param p Number of predictors
#' @param p1 Number of uncorrupted predictors
#' @param p2 Number of corrupted predictors
#' @param mu Penalty parameter for the ADMM algorithm.
#' @param tau Standard error of the additive error matrix when the chosen setting is the additive error setting
#' @param ratio_matrix Observation matrix used in the missing data setting
#' @param etol Tolerance used in the ADMM algorithm
#' @param noise Type of setting chosen : additive or missing
#' 
#' @return list containing \itemize{
#' \item \code{sigma_global} projected matrix for \code{mat}
#' \item \code{rho_global} rho parameter for \code{mat}
#' \item \code{list_matrices_lasso} list of the projected matrices for \code{mat} deprived of the k-fold during cross validation
#' \item \code{list_matrices_error} list of the projected matrices for the k-fold of \code{mat} 
#' \item \code{list_rho_lasso} list of the modified \code{rho} for \code{mat} deprived of the K-fold during cross validation
#' \item \code{list_rho_error} list of the modified \code{rho} for the k-fold of \code{mat} 
#' }
#' 
#' 
#' @export

cv_covariance_matrices_block_descent <- function(K,
                                                 mat,
                                                 y,
                                                 p,
                                                 p1,
                                                 p2,
                                                 mu,
                                                 tau=NULL,
                                                 ratio_matrix=NULL,
                                                 etol=1e-4,
                                                 noise=c("additive","missing")){

  n = nrow(mat)
  p = ncol(mat)
  start = p1 + 1
  n_without_fold = n - floor(n/K)
  n_one_fold = floor(n/K)
  
  folds = sample(cut(seq(1,n),breaks=K,labels=FALSE))
  list_PSD_lasso <- list()
  list_PSD_error <- list()
  list_sigma_lasso <- list()
  list_sigma_error <- list()
  
  ### Case where there is additive noise
  if (noise == "additive"){
    
    #We calculate the global nearest PSD cov matrix and the global surrogate rho, when we take into account the whole data set
    
    print("Doing the global data")
    mat_corrupted <- mat[,start:p]
    mat_uncorrupted <- mat[,1:p1]
    cov_modified <- 1/n*t(mat_corrupted)%*%mat_corrupted - tau**2*diag(p2)
    sigma_global_corrupted <- ADMM_proj(cov_modified,mu=mu, etol = etol)$mat
    sigma_global_uncorrupted <- 1/n * t(mat_uncorrupted)%*%mat_uncorrupted
    
    for (i in 1:K){
      # We calculate the necessary matrices for the cross validation
      
      print(paste("Doing the",i,"fold"))
      index <- which(folds==i, arr.ind= TRUE)
      
      #Calculating the nearest PSD cov matrix when we remove the kth fold, to resolve lasso problem during cross validation
      mat_train <- mat[-index,start:p]
      cov_modified_train <- 1/n_without_fold*t(mat_train)%*%mat_train - tau**2*diag(p2)
      mat_cov_train <- ADMM_proj(cov_modified_train,mu=mu, etol = etol)$mat
      list_PSD_lasso <- rlist::list.append(list_PSD_lasso,mat_cov_train)
      
      #Calculating the nearest PSD cov matrix for the kth fold, to calculate the error on the problem solved without the kth fold
      mat_test <- mat[index,start:p]
      cov_modified_test <- 1/n_one_fold*t(mat_test)%*%mat_test - tau**2*diag(p2)
      mat_cov_test <- ADMM_proj(cov_modified_test,mu=mu, etol = etol)$mat
      list_PSD_error <- rlist::list.append(list_PSD_error,mat_cov_test)
      
      #Calculating the cov matrix when we remove the kth fold, to resolve lasso problem during cross validation
      mat_train <- mat[-index,1:p1]
      cov_train <- 1/n_without_fold*t(mat_train)%*%mat_train
      list_sigma_lasso <- rlist::list.append(list_sigma_lasso,cov_train)
      
      #Calculating the cov matrix  for the kth fold, to calculate the error on the problem solved without the kth fold
      mat_test <- mat[index,1:p1]
      cov_test <- 1/n_one_fold*t(mat_test)%*%mat_test
      list_sigma_error <- rlist::list.append(list_sigma_error,cov_test)
      
      
      
    }
  }
  
  else if (noise == "missing"){
    #We calculate the global nearest PSD cov matrix and the global surrogate rho, when we take into account the whole data set
    #mat_for_adjustment <- diag(probs*(1-probs),p,p) + (1-probs) %*% t(1 - probs)
    print("Doing the global data")
    mat_corrupted <- mat[,start:p]
    mat_uncorrupted <- mat[,1:p1]
    cov_modified <- 1/n*t(mat_corrupted)%*%mat_corrupted / ratio_matrix
    sigma_global_corrupted <- ADMM_proj(cov_modified,mu=mu, etol = etol)$mat
    sigma_global_uncorrupted <- 1/n * t(mat_uncorrupted)%*%mat_uncorrupted
    
    for (i in 1:K){
      # We calculate the necessary matrices for the cross validation
      
      print(paste("Doing the",i,"fold"))
      index <- which(folds==i, arr.ind= TRUE)
      
      #Calculating the nearest PSD cov matrix when we remove the kth fold, to resolve lasso problem during cross validation
      mat_train <- mat[-index,start:p]
      # cov_modified_train <- 1/n_without_fold*t(train_mat)%*%train_mat / mat_for_adjustment
      cov_modified_train <- 1/n_without_fold*t(mat_train)%*%mat_train / ratio_matrix
      mat_cov_train <- ADMM_proj(cov_modified_train,mu=mu, etol = etol)$mat
      list_PSD_lasso <- rlist::list.append(list_PSD_lasso,mat_cov_train)
      
      #Calculating the nearest PSD cov matrix for the kth fold, to calculate the error on the problem solved without the kth fold
      mat_test <- mat[index,start:p]
      # cov_modified_test <- 1/n_one_fold*t(test_mat)%*%test_mat / mat_for_adjustment
      cov_modified_test <- 1/n_one_fold*t(mat_test)%*%mat_test / ratio_matrix
      mat_cov_test <- ADMM_proj(cov_modified_test,mu=mu, etol = etol)$mat
      list_PSD_error <- rlist::list.append(list_PSD_error,mat_cov_test)
      
      #Calculating the cov matrix when we remove the kth fold, to resolve lasso problem during cross validation
      mat_train <- mat[-index,1:p1]
      cov_train <- 1/n_without_fold*t(mat_train)%*%mat_train
      list_sigma_lasso <- rlist::list.append(list_sigma_lasso,cov_train)
      
      #Calculating the cov matrix  for the kth fold, to calculate the error on the problem solved without the kth fold
      mat_test <- mat[index,1:p1]
      cov_test <- 1/n_one_fold*t(mat_test)%*%mat_test
      list_sigma_error <- rlist::list.append(list_sigma_error,cov_test)
      
    }
  }
  
  
  return(list(sigma_global_corrupted = sigma_global_corrupted, sigma_global_uncorrupted = sigma_global_uncorrupted, 
              list_PSD_lasso = list_PSD_lasso, list_PSD_error = list_PSD_error, 
              list_sigma_lasso=list_sigma_lasso, list_sigma_error=list_sigma_error, folds = folds))
}