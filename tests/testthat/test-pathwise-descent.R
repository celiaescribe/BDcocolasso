context("CoCoLasso fit with simulated dataset")


data("simulated_data_missing")
data("simulated_data_additive")

y_missing <- simulated_data_missing[,1]
Z_missing = simulated_data_missing[,2:dim(simulated_data_missing)[2]]
Z_missing = as.matrix(Z_missing)
n_missing <- dim(Z_missing)[1]
p_missing <- dim(Z_missing)[2]
y_missing = as.matrix(y_missing)

output_missing = BDcocolasso::pathwise_coordinate_descent(Z=Z_missing,y=y_missing,n=n_missing,p=p_missing,step=100,K=4,mu=10,tau=NULL, etol = 1e-4,noise = "missing")
beta_missing <- output_missing$beta.opt
beta_missing <- as.matrix(beta_missing)
beta.sd_missing <- output_missing$beta.sd
beta.sd_missing <- as.matrix(beta.sd_missing)
lambda.sd_missing <- output_missing$lambda.sd
lambda.opt_missing <- output_missing$lambda.opt

y_additive <- simulated_data_additive[,1]
Z_additive = simulated_data_additive[,2:dim(simulated_data_additive)[2]]
Z_additive = as.matrix(Z_additive)
n_additive <- dim(Z_additive)[1]
p_additive <- dim(Z_additive)[2]
y_additive = as.matrix(y_additive)

output_additive = BDcocolasso::pathwise_coordinate_descent(Z=Z_additive,y=y_additive,n=n_additive,p=p_additive,center.Z = FALSE, step=100,K=4,mu=10,tau=0.3,etol = 1e-4,noise = "additive")
beta_additive <- output_additive$beta.opt
beta_additive <- as.matrix(beta_additive)
beta.sd_additive <- output_additive$beta.sd
beta.sd_additive <- as.matrix(beta.sd_additive)
lambda.sd_additive <- output_additive$lambda.sd
lambda.opt_additive <- output_additive$lambda.opt



test_that("No error in using CoCoLasso",{
  expect_false(inherits(output_missing, "try-error"))
  expect_false(inherits(output_additive, "try-error"))
  expect_equal(dim(beta.sd_missing)[1],dim(Z_missing)[2])
  expect_equal(dim(beta_missing)[1],dim(Z_missing)[2])
  expect_equal(dim(beta.sd_additive)[1],dim(Z_additive)[2])
  expect_equal(dim(beta_additive)[1],dim(Z_additive)[2])
})


test_that("Number of coefs activated with lambda.sd smaller than number of coefs activated with lambda.opt",{

  expect_equal(length(which(beta.sd_missing != 0)) <= length(which(beta_missing != 0)), TRUE)
  expect_equal(length(which(beta.sd_additive != 0)) <= length(which(beta_additive != 0)), TRUE)
  
})

test_that("Lambda.sd is bigger than lambda.opt",{
  
  expect_equal(lambda.sd_missing >= lambda.opt_missing, TRUE)
  expect_equal(lambda.sd_additive >= lambda.opt_additive, TRUE)
  
})