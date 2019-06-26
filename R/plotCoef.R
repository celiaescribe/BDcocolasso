#' Plot solution path
#'
#' Plot evolution of coefficients for the range of lambda values 
#' 
#' @param object A fitted \code{cocolasso} object as produced by \code{pathwise_coordinate_descent()} or 
#' \code{blockwise_coordinate_descent()}
#' @param linetype Linetype for the lines corresponding to the lambda values. Default is "dashed"
#' @param col Colour for the lines corresponding to lambda values. Default is "black"
#' 
#' @return A plot is produced and nothing is returned
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot geom_line geom_vline xlab ylab geom_errorbar theme
#' 
#' @export

plotCoef <- function(object, linetype="dashed", col="black"){
  data_error = object$data_error
  data_beta = object$data_beta
  beta <- object$beta.opt
  beta <- as.matrix(beta)
  best.lambda <- object$lambda.opt
  lambda.sd <- object$lambda.sd
  beta.sd <- object$beta.sd
  beta.sd <- as.matrix(beta.sd)
  
  data_beta.bis <- reshape2::melt(data_beta, id="lambda", value.name="value" ,variable.name="beta")
  ggplot2::ggplot(data=data_beta.bis) + ggplot2::geom_line(ggplot2::aes(x=log(lambda),y=value,colour=beta)) + ggplot2::geom_vline(xintercept = log(best.lambda), linetype=linetype, colour=col) + ggplot2::geom_vline(xintercept = log(lambda.sd), linetype=linetype, colour=col) + ggplot2::theme(legend.position = "none")+ ggplot2::xlab("Log Lambda") + ggplot2::ylab("Coefficients")

}

