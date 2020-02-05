#' Projected covariance matrices for BD-CoCoLasso
#'
#' Creating projected nearest positive semi-definite covariance matrices for the cross validation step of the BD-CoCoLasso algorithm. In that case, the design matrix must be organized as follozs : uncorrupted features must be the first block of the matrix, and corrupted features must be the second block of the matrix.
#'
#' @param K Number of folds for the cross validation
#' @param mat Covariance matrix to be projected
#' @param y Response vector
#' @param p Number of predictors
#' @param p1 Number of uncorrupted predictors
#' @param p2 Number of corrupted predictors containing additive error
#' @param p3 Number of corrupted predictors containing missingness
#' @param mu Penalty parameter for the ADMM algorithm.
#' @param tau Standard error of the additive error matrix when the chosen setting is the additive error setting
#' @param ratio_matrix Observation matrix used in the missing data setting
#' @param etol Tolerance used in the ADMM algorithm
#' @param mode ADMM or HM
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

cv_covariance_matrices_block_descent_general <- function(K,
                                                 mat,
                                                 y,
                                                 p,
                                                 p1,
                                                 p2,
                                                 p3,
                                                 mu,
                                                 tau=NULL,
                                                 ratio_matrix=NULL,
                                                 etol=1e-4,
                                                 mode="ADMM"){
  
  n = nrow(mat)
  p = ncol(mat)
  start = p1 + 1
  n_without_fold = n - floor(n/K)
  n_one_fold = floor(n/K)
  
  folds = sample(cut(seq(1,n),breaks=K,labels=FALSE))
  list_PSD_lasso_additive <- list()
  list_PSD_error_additive <- list()
  list_PSD_lasso_missing <- list()
  list_PSD_error_missing <- list()
  list_sigma_lasso <- list()
  list_sigma_error <- list()
  
  print("Processing the global data")
  if (p1 != 0) {
    mat_uncorrupted <- mat[,1:p1]
    sigma_global_uncorrupted <- 1/n * t(mat_uncorrupted)%*%mat_uncorrupted
  }
  
  mat_corrupted_additive <- mat[,start:(p1+p2)]
  mat_corrupted_missing <- mat[,(p1+p2+1):p]
  cov_modified_additive <- 1/n*t(mat_corrupted_additive)%*%mat_corrupted_additive - tau**2*diag(p2)
  cov_modified_missing <- 1/n*t(mat_corrupted_missing)%*%mat_corrupted_missing / ratio_matrix
  if (mode=="ADMM") {
    sigma_global_corrupted_additive <- ADMM_proj(cov_modified_additive,mu=mu, etol = etol)$mat
    sigma_global_corrupted_missing <- ADMM_proj(cov_modified_missing,mu=mu, etol = etol)$mat
  }
  if (mode=="HM") {
    sigma_global_corrupted_additive <- HM_proj(sigmaHat = cov_modified_additive,R=ratio_matrix[start:(p1+p2),start:(p1+p2)],mu=mu, tolerance  = etol)
    sigma_global_corrupted_missing <- HM_proj(sigmaHat = cov_modified_missing,R=ratio_matrix[(p1+p2+1):p,(p1+p2+1):p],mu=mu, tolerance = etol)
  }
  
  for (i in 1:K){
    # We calculate the necessary matrices for the cross validation
    
    print(paste("Processing the",i,"fold"))
    index <- which(folds==i, arr.ind= TRUE)
    
    #Calculating the nearest PSD cov matrix when we remove the kth fold, to resolve lasso problem during cross validation
    mat_train <- mat[-index,start:(p1+p2)]
    cov_modified_train <- 1/n_without_fold*t(mat_train)%*%mat_train - tau**2*diag(p2)
    if (mode=="ADMM") {
      mat_cov_train <- ADMM_proj(cov_modified_train,mu=mu, etol = etol)$mat
    }
    if (mode=="HM") {
      mat_cov_train <- HM_proj(sigmaHat = cov_modified_train,R=ratio_matrix[start:(p1+p2),start:(p1+p2)],mu = mu, tolerance = etol)
    }
    list_PSD_lasso_additive <- rlist::list.append(list_PSD_lasso_additive,mat_cov_train)
    
    mat_train <- mat[-index,(p1+p2+1):p]
    cov_modified_train <- 1/n_without_fold*t(mat_train)%*%mat_train / ratio_matrix
    if (mode=="ADMM") {
      mat_cov_train <- ADMM_proj(cov_modified_train,mu=mu, etol = etol)$mat
    }
    if (mode=="HM") {
      mat_cov_train <- HM_proj(sigmaHat = cov_modified_train,R=ratio_matrix[(p1+p2+1):p,(p1+p2+1):p],mu = mu, tolerance = etol)
    }
    list_PSD_lasso_missing <- rlist::list.append(list_PSD_lasso_missing,mat_cov_train)
    
    #Calculating the nearest PSD cov matrix for the kth fold, to calculate the error on the problem solved without the kth fold
    mat_test <- mat[index,start:(p1+p2)]
    cov_modified_test <- 1/n_one_fold*t(mat_test)%*%mat_test - tau**2*diag(p2)
    if (mode=="ADMM") {
      mat_cov_test <- ADMM_proj(cov_modified_test,mu=mu, etol = etol)$mat
    }
    if (mode=="HM") {
      mat_cov_test <- HM_proj(sigmaHat = cov_modified_test,R=ratio_matrix[start:(p1+p2),start:(p1+p2)],mu = mu, tolerance = etol)
    }
    list_PSD_error_additive <- rlist::list.append(list_PSD_error_additive,mat_cov_test)
    
    mat_test <- mat[index,(p1+p2+1):p]
    cov_modified_test <- 1/n_one_fold*t(mat_test)%*%mat_test / ratio_matrix
    mat_cov_test <- ADMM_proj(cov_modified_test,mu=mu, etol = etol)$mat
    if (mode=="ADMM") {
      mat_cov_test <- ADMM_proj(cov_modified_test,mu=mu, etol = etol)$mat
    }
    if (mode=="HM") {
      mat_cov_test <- HM_proj(sigmaHat = cov_modified_test,R=ratio_matrix[(p1+p2+1):p,(p1+p2+1):p],mu = mu, tolerance = etol)
    }
    list_PSD_error_missing <- rlist::list.append(list_PSD_error_missing,mat_cov_test)
    
    #Calculating the cov matrix when we remove the kth fold, to resolve lasso problem during cross validation
    mat_train <- mat[-index,1:p1]
    cov_train <- 1/n_without_fold*t(mat_train)%*%mat_train
    list_sigma_lasso <- rlist::list.append(list_sigma_lasso,cov_train)
    
    #Calculating the cov matrix  for the kth fold, to calculate the error on the problem solved without the kth fold
    mat_test <- mat[index,1:p1]
    cov_test <- 1/n_one_fold*t(mat_test)%*%mat_test
    list_sigma_error <- rlist::list.append(list_sigma_error,cov_test)
  }
  if (p1 != 0) {
    return(list(sigma_global_uncorrupted = sigma_global_uncorrupted,
                sigma_global_corrupted_additive = sigma_global_corrupted_additive, 
                sigma_global_corrupted_missing = sigma_global_corrupted_missing,
                list_PSD_lasso_additive = list_PSD_lasso_additive, list_PSD_error_additive = list_PSD_error_additive,
                list_PSD_lasso_missing =list_PSD_lasso_missing, list_PSD_error_missing = list_PSD_error_missing,
                list_sigma_lasso=list_sigma_lasso, list_sigma_error=list_sigma_error, folds = folds))
  }
  if (p1 == 0) {
    return(list(sigma_global_corrupted_additive = sigma_global_corrupted_additive, 
                sigma_global_corrupted_missing = sigma_global_corrupted_missing,
                list_PSD_lasso_additive = list_PSD_lasso_additive, list_PSD_error_additive = list_PSD_error_additive,
                list_PSD_lasso_missing =list_PSD_lasso_missing, list_PSD_error_missing = list_PSD_error_missing,
                folds = folds))
  }
}
