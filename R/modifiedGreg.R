#' @export
modifiedGreg <- function(x, ...)  UseMethod("modifiedGreg")

#' @export
#' @rawNamespace export(modifiedGreg.numeric)
modifiedGreg.numeric <- function(y,
                                 xsample,
                                 xpop,
                                 condition,
                                 pi = NULL,
                                 pi2 = NULL,
                                 datatype = "raw",
                                 model = "linear",
                                 var_est = F,
                                 var_method = "LinHB",
                                 modelselect = FALSE,
                                 lambda = "lambda.min",
                                 domain_col_name = NULL,
                                 estimation_conditions = NULL,
                                 N = NULL,
                                 messages = T) {

  .args <- as.list(environment())
  do.call(validate_modifiedGreg, .args)
  funcCall <- match.call()

}


# modifiedGreg.data.frame <- function(y, ...) {
#
#   modifiedGreg.numeric(y[[1]], ...)
#
# }



