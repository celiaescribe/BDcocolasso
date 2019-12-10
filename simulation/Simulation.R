library(mvtnorm)
library(clusterGeneration)
library(glmnet)

cov_iid <- function(p) {
  cov <- diag(1,p,p)
  cov
}

cov_symmetric <- function(p) {
  cov <- matrix(0.5,p,p)
  diag(cov) <- 1
  cov
}

cov_autoregressive <- function(p){
  cov <- matrix(0,p,p)
  for (i in 1:p){
    for (j in 1:p){
      cov[i,j] = 0.5**(abs(i-j))
    }
  }
  cov
}

cov_unrestricted <- function(p) {
#  cov <- genPositiveDefMat(p,covMethod = "eigen",lambdaLow = 10)$Sigma
  cov <- matrix(0,p,p)
  for (i in 1:p) {
    for (j in 1:p) {
      cov[i,j] <- runif(1,0,0.2)
      cov[j,i] <- cov[i,j]
    }
  }
  diag(cov) <- 1
  cov %*% t(cov)
}

testBias <- function(beta_true, beta_est) {
  sqrt(sum((beta_true - beta_est)^2))
}

testSparsity <- function(beta_true, beta_est) {
  sum(which(beta_true==0) %in% which(beta_est==0)) + sum(which(beta_true!=0) %in% which(beta_est!=0)) -> correctCount
  correctCount / length(beta_est)
}

testFPR <- function(beta_true, beta_est) {
  sum(which(beta_est!=0) %in% which(beta_true==0)) / length(beta_true[beta_true==0])
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

sigma_e <- 2

nFeature <- 200
nObs <- 10000
causal_rate <- 0.1

# iid setting
X <- scale(rmvnorm(nObs, sigma = cov_iid(nFeature)))
X_test <- scale(rmvnorm(nObs, sigma = cov_iid(nFeature)))

# symmetric covariance setting
X <- scale(rmvnorm(nObs, sigma = cov_symmetric(nFeature)))
X_test <- scale(rmvnorm(nObs, sigma = cov_symmetric(nFeature)))

# autoregressive covariance setting
X <- scale(rmvnorm(nObs, sigma = cov_autoregressive(nFeature)))
X_test <- scale(rmvnorm(nObs, sigma = cov_autoregressive(nFeature)))

# unrestricted covariance setting
X <- scale(rmvnorm(nObs, sigma = cov_unrestricted(nFeature)))
X_test <- scale(rmvnorm(nObs, sigma = cov_unrestricted(nFeature)))

# generate response
true_beta <- rep(0, nFeature)
true_beta[sample(1:nFeature,round(nFeature * causal_rate),replace = F)] <- rnorm(nFeature * causal_rate)
true_beta <- round(true_beta,2)
Y <- X %*% true_beta + rnorm(nObs, sd = sigma_e)
Y_test <- X_test %*% true_beta + rnorm(nObs, sd = sigma_e)

# additive error
p_error <- 0.1
tau <- 0.3
A <- rmvnorm(nObs, sigma = diag(tau, nFeature * p_error))
Z <- X
Z[,(nFeature * (1 - p_error) + 1) : nFeature] <- Z[,(nFeature * (1 - p_error) + 1) : nFeature] + A

NLfit <- cv.glmnet(Z, Y, type.measure = "mse", nfolds = 5)
predict(NLfit, type = "coefficients", s = "lambda.min")[-1] -> NLfitCoef
testBias(true_beta, NLfitCoef)
testFPR(true_beta, NLfitCoef)
testSparsity(true_beta, NLfitCoef)
predict(NLfit, newx = X_test, type = "response", s = "lambda.min") -> NLfitPred
cor(Y_test, NLfitPred)^2

cocofit <- coco(Z, Y, nObs, nFeature, p1 = nFeature * (1 - p_error), p2 = nFeature * p_error,
                center.Z = F, scale.Z = F, center.y = T, scale.y = T, K = 5, tau = sqrt(tau), 
                noise = "additive", block = T, penalty = "lasso")
COCOfitCoef <- cocofit$beta.opt
testBias(true_beta, COCOfitCoef * sd(Y))
testFPR(true_beta, COCOfitCoef)
testSparsity(true_beta, COCOfitCoef)
predict(cocofit, newx = X_test, type = "response", s = cocofit$lambda.opt) -> COCOfitPred
cor(Y_test, COCOfitPred)^2

# missingness
p_error <- 0.1
rho <- 0.2
M <- matrix(rbinom(nObs * nFeature * p_error, 1, 1-rho), nrow = nObs, ncol = nFeature * p_error)
M[M==0] <- NA
Z <- X
Z[,(nFeature * (1 - p_error) + 1) : nFeature] <- Z[,(nFeature * (1 - p_error) + 1) : nFeature] * M

Z_for_NL <- Z
sd.Z <- sapply(1:nFeature, function(j)sd_without_NA_block(j,Z))
Z_for_NL[,(nFeature * (1 - p_error) + 1):nFeature] = sapply(1:(nFeature * p_error), function(j)rescale_without_NA_block(j,Z_for_NL,nFeature*(1 - p_error)))
Z_for_NL[,(nFeature * (1 - p_error) + 1):nFeature] = sapply(1:(nFeature * p_error), function(j)change_NA_value_block(j,Z_for_NL,nFeature*(1 - p_error)))
Z_for_NL[,1:(nFeature * (1 - p_error))] = scale(Z_for_NL[,1:(nFeature * (1 - p_error))], center = TRUE, scale = FALSE)
Z_for_NL = sapply(1:nFeature, function(j)scale_manual_with_sd(j,Z_for_NL,sd.Z))

NLfit <- cv.glmnet(Z_for_NL, Y, type.measure = "mse", nfolds = 5)
predict(NLfit, type = "coefficients", s = "lambda.min")[-1] -> NLfitCoef
testBias(true_beta, NLfitCoef)
testFPR(true_beta, NLfitCoef)
testSparsity(true_beta, NLfitCoef)
predict(NLfit, newx = X, type = "response", s = "lambda.min") -> NLfitPred
cor(Y, NLfitPred)^2

cocofit <- coco(Z, Y, nObs, nFeature, p1 = nFeature * (1 - p_error), p2 = nFeature * p_error,
                center.Z = T, scale.Z = T, center.y = T, scale.y = T, K = 5, tau = NULL, 
                noise = "missing", block = T, penalty = "lasso")
COCOfitCoef <- cocofit$beta.opt
testBias(true_beta, COCOfitCoef * sd(Y))
testFPR(true_beta, COCOfitCoef)
testSparsity(true_beta, COCOfitCoef)
predict(cocofit, newx = X_test, type = "response", s = cocofit$lambda.opt) -> COCOfitPred
cor(Y_test, COCOfitPred)^2

# both types of error
p_error_additive <- 0.1
p_error_missing <- 0.1
tau <- 0.2
rho <- 0.2
A <- rmvnorm(nObs, sigma = diag(tau, nFeature * p_error_additive))
M <- matrix(rbinom(nObs * nFeature * p_error_missing, 1, 1-rho), nrow = nObs, ncol = nFeature * p_error_additive)
M[M==0] <- NA
Z <- X
P2 <- nFeature * p_error_additive
P3 <- nFeature * p_error_missing
P1 <- nFeature - P2 - P3

Z[,(P1 + 1) : (P1 + P2)] <- Z[,(P1 + 1) : (P1 + P2)] + A
Z[,(P1 + P2 + 1) : nFeature] <- Z[,(P1 + P2 + 1) : nFeature] * M

Z_for_NL <- Z
sd.Z <- sapply(1:nFeature, function(j)sd_without_NA_block(j,Z))
Z_for_NL[,(P1 + P2 + 1):nFeature] = sapply(1:P3, function(j)rescale_without_NA_block(j,Z_for_NL,P1+P2))
Z_for_NL[,(P1 + P2 + 1):nFeature] = sapply(1:P3, function(j)change_NA_value_block(j,Z_for_NL,P1+P2))
Z_for_NL[,1:(P1+P2)] = scale(Z[,1:(P1+P2)], center = TRUE, scale = FALSE)
if (P1 == 0) {
  Z_for_NL[,c((P1+P2+1):nFeature)] = sapply(1:(nFeature - P2), function(j)scale_manual_with_sd(j,Z_for_NL[,c((P1+P2+1):nFeature)],sd.Z[c((P1+P2+1):nFeature)]))
}
if (P1 != 0) {
  Z_for_NL[,c(1:P1,(P1+P2+1):nFeature)] = sapply(1:(nFeature - P2), function(j)scale_manual_with_sd(j,Z_for_NL[,c(1:P1,(P1+P2+1):nFeature)],sd.Z[c(1:P1,(P1+P2+1):nFeature)]))
}

NLfit <- cv.glmnet(Z_for_NL, Y, type.measure = "mse", nfolds = 5)
predict(NLfit, type = "coefficients", s = "lambda.min")[-1] -> NLfitCoef
testBias(true_beta, NLfitCoef)
testFPR(true_beta, NLfitCoef)
testSparsity(true_beta, NLfitCoef)
predict(NLfit, newx = X_test, type = "response", s = "lambda.min") -> NLfitPred
cor(Y_test, NLfitPred)^2

cocofit <- generalcoco(Z, Y, nObs, p = nFeature, p1 = P1, p2 = P2, p3 = P3, center.Z = T, scale.Z = T,
                       center.y = T, scale.y = T, K = 5, tau = sqrt(tau),
                       penalty = "lasso")
COCOfitCoef <- cocofit$beta.opt
testBias(true_beta, COCOfitCoef * sd(Y))
testFPR(true_beta, COCOfitCoef)
testSparsity(true_beta, COCOfitCoef)
predict(cocofit, newx = X_test, type = "response", s = cocofit$lambda.opt) -> COCOfitPred
cor(Y_test, COCOfitPred)^2



