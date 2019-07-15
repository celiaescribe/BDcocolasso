context("BD-CoCoLasso fit with simulated dataset")


data("simulated_data_missing_block")
data("simulated_data_additive_block")
p1 <- 180
p2 <- 20

y_missing <- simulated_data_missing_block[,1]
Z_missing = simulated_data_missing_block[,2:dim(simulated_data_missing_block)[2]]
Z_missing = as.matrix(Z_missing)
n_missing <- dim(Z_missing)[1]
p_missing <- dim(Z_missing)[2]
y_missing = as.matrix(y_missing)

fit_missing = BDcocolasso::coco(Z=Z_missing,y=y_missing,n=n_missing,p=p_missing,p1=p1,p2=p2,step=100,K=4,mu=10,tau=NULL,noise="missing",block=TRUE, penalty="lasso")
beta_missing <- fit_missing$beta.opt
beta_missing <- as.matrix(beta_missing)
beta.sd_missing <- fit_missing$beta.sd
beta.sd_missing <- as.matrix(beta.sd_missing)

y_additive <- simulated_data_additive_block[,1]
Z_additive = simulated_data_additive_block[,2:dim(simulated_data_additive_block)[2]]
Z_additive = as.matrix(Z_additive)
n_additive <- dim(Z_additive)[1]
p_additive <- dim(Z_additive)[2]
y_additive = as.matrix(y_additive)

fit_additive = BDcocolasso::coco(Z=Z_additive,y=y_additive,n=n_additive,p=p_additive,p1=p1,p2=p2,center.Z = FALSE, step=100,K=4,mu=10,tau=0.3,noise="additive",block=TRUE, penalty="lasso")
beta_additive <- fit_additive$beta.opt
beta_additive <- as.matrix(beta_additive)
beta.sd_additive <- fit_additive$beta.sd
beta.sd_additive <- as.matrix(beta.sd_additive)


test_that("No error in using BD-CoCoLasso",{
  expect_false(inherits(fit_missing, "try-error"))
  expect_false(inherits(fit_additive, "try-error"))
  expect_equal(dim(beta.sd_missing)[1],dim(Z_missing)[2])
  expect_equal(dim(beta_missing)[1],dim(Z_missing)[2])
  expect_equal(dim(beta.sd_additive)[1],dim(Z_additive)[2])
  expect_equal(dim(beta_additive)[1],dim(Z_additive)[2])
  expect_is(fit_missing, "coco")
  expect_is(fit_additive, "coco")
})

test_that("Number of coefs activated with lambda.sd smaller than number of coefs activated with lambda.opt",{
  
  expect_equal(length(which(beta.sd_missing != 0)) < length(which(beta_missing != 0)), TRUE)
  expect_equal(length(which(beta.sd_additive != 0)) < length(which(beta_additive != 0)), TRUE)
  
})