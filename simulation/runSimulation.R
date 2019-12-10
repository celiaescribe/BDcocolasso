library(mvtnorm)
library(clusterGeneration)
library(glmnet)
library(data.table)
library(ggplot2)

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

sigma_e <- 3

nFeature <- 200
nObs <- 10000
causal_rate <- 0.1

biasNL <- c()
fprNL <- c()
sparsityNL <- c()
varexpNL <- c()

biasCoCo <- c()
fprCoCo <- c()
sparsityCoCo <- c()
varexpCoCo <- c()

for (tau in seq(0, 0.2, 0.05)) {
  cat(paste0("processing: tau = ",tau),collapse = "\n")
  for (ite in 1:2) {
    cat(paste0("processing iteration: ",ite,"/10"),collapse = "\n")
    # iid setting
    X <- scale(rmvnorm(nObs, sigma = cov_symmetric(nFeature)))
    X_test <- scale(rmvnorm(nObs, sigma = cov_symmetric(nFeature)))
    
    # generate response
    true_beta <- rep(0, nFeature)
    true_beta[sample(1:nFeature,round(nFeature * causal_rate),replace = F)] <- rnorm(nFeature * causal_rate)
    true_beta <- round(true_beta,2)
    Y <- X %*% true_beta + rnorm(nObs, sd = sigma_e)
    Y_test <- X_test %*% true_beta + rnorm(nObs, sd = sigma_e)
    
    # additive error
    p_error <- 0.1
#    tau <- 0.2
    A <- rmvnorm(nObs, sigma = diag(tau, nFeature * p_error))
    Z <- X
    Z[,(nFeature * (1 - p_error) + 1) : nFeature] <- Z[,(nFeature * (1 - p_error) + 1) : nFeature] + A
    
    NLfit <- cv.glmnet(Z, Y, type.measure = "mse", nfolds = 5)
    predict(NLfit, type = "coefficients", s = "lambda.min")[-1] -> NLfitCoef
    biasNL <- c(biasNL, testBias(true_beta, NLfitCoef))
    fprNL <- c(fprNL, testFPR(true_beta, NLfitCoef))
    sparsityNL <- c(sparsityNL, testSparsity(true_beta, NLfitCoef))
    predict(NLfit, newx = X_test, type = "response", s = "lambda.min") -> NLfitPred
    varexpNL <- c(varexpNL, cor(Y_test, NLfitPred)^2)
    
    cocofit <- coco(Z, Y, nObs, nFeature, p1 = nFeature * (1 - p_error), p2 = nFeature * p_error,
                    center.Z = F, scale.Z = F, center.y = T, scale.y = T, K = 5, tau = sqrt(tau), 
                    noise = "additive", block = T, penalty = "lasso")
    COCOfitCoef <- cocofit$beta.opt
    biasCoCo <- c(biasCoCo, testBias(true_beta, COCOfitCoef * sd(Y)))
    fprCoCo <- c(fprCoCo, testFPR(true_beta, COCOfitCoef))
    sparsityCoCo <- c(sparsityCoCo, testSparsity(true_beta, COCOfitCoef))
    predict(cocofit, newx = X_test, type = "response", s = cocofit$lambda.opt) -> COCOfitPred
    varexpCoCo <- c(varexpCoCo, cor(Y_test, COCOfitPred)^2)
  }
}

plotdat <- data.table(Bias = c(biasCoCo, biasNL),
                      FPR = c(fprCoCo, fprNL),
                      Sparsity = c(sparsityCoCo, sparsityNL),
                      R2 = c(varexpCoCo, varexpNL),
                      Method = rep(c("BD-CoCoLasso","Naïve Lasso"), each = length(biasCoCo)),
                      tau = rep(rep(seq(0,0.8,0.05),each = 100),2))
plotdat[, lapply(.SD, mean), by = c("Method","tau")] -> meanplotdat
plotdat[, lapply(.SD, sd), by = c("Method","tau")] -> sdplotdat
cbind(meanplotdat,sdplotdat[,-c(1:2)]) -> plotdat
colnames(plotdat)[7:10] <- c("BiasSD","FPRSD","SparsitySD","R2SD")

ggplot(plotdat, aes(x = tau, y = Bias, color = Method)) + geom_point() + theme_classic() + geom_errorbar(aes(ymin = Bias - BiasSD, ymax = Bias + BiasSD),width=0.03) + theme(legend.position = "bottom",
                                                                                                                                                                             axis.title = element_text(size=20),
                                                                                                                                                                             axis.text = element_text(size=15),
                                                                                                                                                                             legend.text = element_text(size=20),
                                                                                                                                                                             legend.title = element_text(size=20)) +
  xlab(expression(tau)) + ylab(expression(paste("||",beta - hat(beta),"||"[2]))) + ylim(0,0.8) + scale_color_manual(values = c("red","darkblue"))












