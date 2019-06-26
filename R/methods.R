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
print.coco <- function(x) {
  cat("\nCall: ", deparse(x$call), "\n\n")
  print(df_error = x$data_error)
}