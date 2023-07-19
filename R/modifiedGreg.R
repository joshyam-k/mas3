#' @export
modifiedGreg <- function(x, ...)  UseMethod("modifiedGreg")


modifiedGreg.default <- function(y,
                                 xsample,
                                 xpop,
                                 conditions,
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

  funcCall <- match.call()
  # add validators

}
