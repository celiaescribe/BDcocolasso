context("Coef plots")

data("simulated_data_missing")
data("simulated_data_additive")

y_missing <- simulated_data_missing[,1]
Z_missing = simulated_data_missing[,2:dim(simulated_data_missing)[2]]
Z_missing = as.matrix(Z_missing)
n_missing <- dim(Z_missing)[1]
p_missing <- dim(Z_missing)[2]
y_missing = as.matrix(y_missing)
fit_missing = BDcocolasso::pathwise_coordinate_descent(Z=Z_missing,y=y_missing,n=n_missing,p=p_missing,step=100,K=4,mu=10,tau=NULL, etol = 1e-4,noise = "missing")

test_that("Plot returns ggplot object",{
  p1 <- plotCoef(fit_missing)
  expect_is(p1,"ggplot")
  
  p2 <- plotError(fit_missing)
  expect_is(p2,"ggplot")
})

test_that("Printing ggplot object actually works",{
  p1 <- plotCoef(fit_missing)
  expect_error(print(p1), NA)
  
  p2 <- plotError(fit_missing)
  expect_error(print(p2), NA)
})