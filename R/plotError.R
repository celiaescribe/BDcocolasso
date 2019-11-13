#' Plot solution path
#'
#' Plot evolution of coefficients for the range of lambda values 
#' 
#' @param object A fitted \code{cocolasso} object as produced by \code{pathwise_coordinate_descent()} or 
#' \code{blockwise_coordinate_descent()}
#' @param linetype Linetype for the lines corresponding to the lambda values. Default is "dashed"
#' @param col Colour for the lines corresponding to lambda values. Default is "black"
#' @param colLine Colour of the error line. Default is "red"
#' 
#' @return A plot is produced and nothing is returned
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot geom_line geom_vline xlab ylab geom_errorbar theme
#' 
#' @export

plotError <- function(object, linetype="dashed", col="black", colLine="red"){
  data_error = object$data_error
  data_beta = object$data_beta
  beta <- object$beta.opt
  beta <- as.matrix(beta)
  best.lambda <- object$lambda.opt
  lambda.sd <- object$lambda.sd
  beta.sd <- object$beta.sd
  beta.sd <- as.matrix(beta.sd)
  
  data_beta.bis <- melt(data_beta, id="lambda", value.name="value" ,variable.name="beta")
  ggplot2::ggplot(data=data_error) + ggplot2::geom_line(data=data_error,ggplot2::aes(log(.data$lambda),.data$error),colour=colLine,size=1) + ggplot2::xlab("Log Lambda") + ggplot2::ylab("Error") + ggplot2::geom_errorbar(data=data_error,ggplot2::aes(x=log(.data$lambda), ymin = .data$error.inf, ymax=.data$error.sup), colour=col) + ggplot2::geom_vline(xintercept = log(best.lambda), linetype=linetype) + ggplot2::geom_vline(xintercept = log(lambda.sd), linetype=linetype)
  
}
