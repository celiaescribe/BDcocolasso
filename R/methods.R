######################################
#' R Source code file for predict, coef and print methods for the BDcoco package
#' Author: Celia Escribe
#' Created: 2019
#' Updated: June 27, 2019
#####################################

#' @title Print Method for \code{coco} object
#' @description Print a summary of the \code{coco} path at each step along the
#'   path.
#' @param x fitted \code{coco} object
#' @param ... additional print arguments
#' @return OUTPUT_DESCRIPTION

#' @seealso \code{\link{coco}}
#' @rdname print.coco
#' @export
#' 
print.coco <- function(x, ...) {
  cat("\nCall: ", deparse(x$call), "\n\n")
  print(x$data_error[,c("lambda","error")])
}

#' @title Make predictions from a coco object
#' @description Similar to other predict methods, this functions predicts fitted
#'   values, logits, coefficients and more from a fitted \code{coco} object.
#' @param object Fitted \code{coco} model object
#' @param newx matrix of new values for \code{x} at which predictions are to be
#' made. Do not include the intercept (this function takes care of that). Must
#' be a matrix. This argument is not used for \code{type=c("coefficients")}. This 
#' matrix must have the same number of columns originally supplied to the \code{coco}
#' fitting function.
#' @param s Value(s) of the penalty parameter \code{lambda} at which predictions
#' are required. Default is the entire sequence used to create the model.
#' @param lambda.pred Value(s) of the penalty parameter \code{lambda} at which coefficients
#' are extracted to calculate the response based on the matrix of new values \code{newx}.
#' Default is \code{lambda.sd}.
#' @param type Type of prediction required. Type \code{"coefficients"} computes the coefficients 
#' at the requested values for \code{s}. Type \code{response} computes the response based 
#' on the covariates values in \code{newx}, for coefficients corresponding to \code{lambda} value
#' in \code{lambda.pred}
#' @param ... currently ignored
#' 
#' @return The object returned depends on type.
#' 
#' @rdname predict.coco
#' @export

predict.coco <- function(object, newx, s=NULL, lambda.pred=NULL, type=c("response","coefficients"), ...){
  
  if (type == "response"){
    if (missing(newx)){
      stop("newx is missing. Please supply the vector of covariates for which response is required.")
    }
  }
  lambda.seq = object$data_beta[,"lambda"]
  nbeta = object$data_beta
  if(!is.null(object$vnames)){

    colnames(nbeta) <- c("lambda", object$vnames)
  }
  if (!is.null(s)){
    index <- match(s,lambda.seq)
    nbeta <- nbeta[index,]
  }
  if(type == "coefficients"){
    coef = t(nbeta[,-1])
    colnames(coef) <- c("Coefficient")
    return(coef)
  }else{
    #browser()
    mean.Z = object$mean.Z
    mean.y = object$mean.y
    sd.Z = object$sd.Z
    sd.y = object$sd.y
    matcoef = object$beta.sd
    if(!is.null(lambda.pred)){
      index <- match(lambda.pred,lambda.seq)
      matcoef <- t(object$data_beta[index,2:dim(object$data_beta)[2]])
    }
    newx <- (newx - mean.Z)/sd.Z * sd.y
    y <- newx %*% matcoef + mean.y
    return(y)
  }
}

#' @inheritParams predict.coco
#' @rdname predict.coco
#' @export
coef.coco <- function(object, s = NULL, ...) {
  stats::predict(object, s = s, type = "coefficients", ...)
}


