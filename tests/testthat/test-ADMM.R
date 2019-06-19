context("Testing the projection with ADMM algorithm")

rescale_without_NA <- function(j,Z){
  m <- mean(Z[which(!is.na(Z[,j]), arr.ind = TRUE),j])
  Z[,j] - m
}

change_NA_value <- function(j,Z){
  Z[which(is.na(Z[,j]), arr.ind = TRUE),j] <- 0
  Z[,j]
}

data("simulated_data_missing")
Z = simulated_data[,2:dim(simulated_data)[2]]
Z = as.matrix(Z)
n <- dim(Z)[1]
p <- dim(Z)[2]
Z = sapply(1:p, function(j)rescale_without_NA(j,Z))
Z = sapply(1:p, function(j)change_NA_value(j,Z))

Sigma = 1/n * t(Z) %*% Z
mat_proj <- BDcocolasso::ADMM_proj(mat=Sigma)$mat

test_that("Obtaining a positive semi-definite matrix with ADMM algo", {
  expect_equal(mean(eigen(mat_proj)$values >= 0),1)
})