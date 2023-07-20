#' @importFrom cli cli_abort
validate_modifiedGreg <- function(y,
                                  xsample,
                                  xpop,
                                  condition,
                                  pi,
                                  pi2,
                                  datatype,
                                  model,
                                  var_est,
                                  var_method,
                                  modelselect,
                                  lambda,
                                  domain_col_name,
                                  estimation_conditions,
                                  N,
                                  messages) {


  if (!is.data.frame(xsample) || !is.data.frame(xpop)) {
    cli_abort(c("x" = "xsample and xpop must be objects of class `data.frame()`"))
  }

  if (length(condition) != nrow(xsample)) {
    cli_abort(c("x" = "condition input must be the same length as xsample."))
  }

  if(!is.logical(var_est)) {
    cli_abort(c("x" = "var_est must be an object of type logical."))
  }

  if(!is.logical(messages)) {
    cli_abort(c("x" = "messages must be an object of type logical."))
  }

  if(!is.element(lambda, c("lambda.1se", "lambda.min"))) {
    cli_abort(c("x" = "lambda input incorrect. It has to be \"lambda.min\" or \"lambda.1se\"."))
  }

  if (datatype != "raw" && !("N" %in% names(xpop))) {
    cli_abort(c("x" = "xpop must contain a column for population size by domain called 'N' when datatype != raw."))
  }

  if (!all(names(xsample) %in% names(xpop))) {
    cli_abort(c("x" = "All of the column names in `xsample` must exist in `xpop`."))
  }

  if (!is.element(var_method, c("LinHB", "LinHH", "LinHTSRS", "LinHT", "bootstrapSRS"))) {
    cli_abort(c("x" = "Variance method input incorrect. It has to be \"LinHB\", \"LinHH\", \"LinHT\", \"LinHTSRS\", or \"bootstrapSRS\"."))
  }

  if (!is.element(model, c("linear","logistic"))) {
    cli_abort(c("x" = "Method input incorrect, has to be either \"linear\" or \"logistic\""))
  }

  if (model == "logistic" && datatype != "raw") {
    cli_abort(c("x" = "Must supply the raw population data to fit the logistic regression estimator."))
  }

  if (!is.element(datatype, c("raw","totals", "means"))) {
    cli_abort(c("x" = "datatype input incorrect, has to be either \"raw\", \"totals\" or \"means\""))
  }

  ncol_diff <- ncol(xpop) - ncol(xsample)
  if (datatype == "raw" && (ncol_diff != 1)) {
    cli_abort(c("x" = "Incorrect number of columns in either xpop or xsample. When datatype = \"raw\" xpop should only contain columns with the same names as xsample and a column with the domains"))
  } else if (is.element(datatype, c("means", "totals")) && (ncol_diff != 2)) {
    cli_abort(c("x" = "Incorrect number of columns in either xpop or xsample. When datatype != \"raw\" xpop should only contain columns with the same names as xsample as well as column with the domains and a column with the population sizes."))
  }

}
