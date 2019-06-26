#' Simulated data missing
#'
#' @description Simulated data y=X * beta + epsilon where beta=c(3,2,0,0,1.5,0,...) and with missing values
#'
#' @format this data frame has 10000 rows and the following 50 columns:
#' \describe{
#'   \item{y}{the response}
#'   \item{V}{Other features}
#' }
"simulated_data_missing"

#' Simulated data additive
#'
#' @description Simulated data y=X * beta + epsilon where beta=c(3,2,0,0,1.5,0,...) and with additive error
#'
#' @format this data frame has 10000 rows and the following 50 columns:
#' \describe{
#'   \item{y}{the response}
#'   \item{V}{Other features}
#' }
"simulated_data_additive"

#' Simulated data missing block
#' 
#' @description Simulated data y=X1 \* beta1 + Z2 \* beta2 + epsilon where beta1=c(3,2,0,0,1.5,0,...) 
#' qnd beta2 = c(0,...,1.5,0,0,2,3). X1 is uncorrupted while Z2 is corrupted with missing data
#'
#' @format this data frame has 10000 rows and the following 200 columns:
#' \describe{
#'   \item{y}{the response}
#'   \item{V}{180 first columns are uncorrupted covariates}
#'   \item{V}{20 last columns are the corrupted covariates with missing data}
#' }
"simulated_data_missing_block"

#' Simulated data additive block
#' 
#' @description Simulated data y=X1 \* beta1 + Z2 \* beta2 + epsilon where beta1=c(3,2,0,0,1.5,0,...) 
#' qnd beta2 = c(0,...,1.5,0,0,2,3). X1 is uncorrupted while Z2 is corrupted with additive error
#'
#' @format this data frame has 10000 rows and the following 200 columns:
#' \describe{
#'   \item{y}{the response}
#'   \item{V}{180 first columns are uncorrupted covariates}
#'   \item{V}{20 last columns are the corrupted covariates with additive error}
#' }
"simulated_data_additive_block"
