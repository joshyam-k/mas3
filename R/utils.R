#' @importFrom cli cli_abort
validate_modifiedGreg <- function(y,
                                  xsample,
                                  xpop,
                                  domains,
                                  pi,
                                  pi2,
                                  datatype,
                                  model,
                                  var_est,
                                  var_method,
                                  modelselect,
                                  lambda,
                                  domain_col_name,
                                  estimation_domains,
                                  N,
                                  messages) {


  if (!is.data.frame(xsample) || !is.data.frame(xpop)) {
    cli_abort(c("x" = "xsample and xpop must be objects of class `data.frame()`"))
  }

  if (length(domains) != nrow(xsample)) {
    cli_abort(c("x" = "domains input must be the same length as xsample."))
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

  if (is.null(pi)) {
    if (messages) {
      message("Assuming simple random sampling")
    }

  }

}


varMase <- function(y, pi, pi2 = NULL, method = "LinHB", N = NULL){ #, fpc = fpc){

  n <- length(y)

  if(method=="LinHT" & is.null(pi2)){
    cli_abort(c("x" = "For LinHT variance estimator, need to provide second order inclusion probabilities matrix."))
  }

  if(method == "LinHB"){
    a <- n/(n-1)*(1-pi)
    e <- as.vector(pi^(-1)*y) - c(sum(a)^(-1)*(pi^(-1)*a)%*%y)
    varEst <- sum(a*e^2)
  }

  if(method == "LinHH"){
    t <- pi^(-1)%*%y
    varEst <- 1/(n*(n-1))*t(as.vector(n*y*pi^(-1)) - as.numeric(t))%*%(as.vector(n*y*pi^(-1)) - as.numeric(t))
  }

  if(method == "LinHTSRS"){
    if(is.null(N)){
      N <- sum(pi^(-1))
    }
    varEst <- (N-n)*(N/n)*var(y)
  }

  if(method== "LinHT"){
    a <- (pi2 - pi%*%t(pi))*pi2^(-1)*(pi%*%t(pi))^(-1)*y%*%t(y)
    varEst <- sum(a)
  }

  return(varEst)

}


by_domain_linear <- function(domain_id,
                             xsample,
                             xpop_d,
                             domain_col_name,
                             comp1,
                             comp2,
                             betas,
                             common_pred_vars,
                             y,
                             var_est,
                             var_method,
                             weight,
                             pi,
                             pi2) {

  domain_indic_vec <- as.integer(xsample[domain_col_name] == domain_id)

  xpop_domain <- xpop_d[xpop_d[domain_col_name] == domain_id, , drop = FALSE]
  xpop_d_domain <- unlist(xpop_domain[-which(names(xpop_domain) == domain_col_name)])
  xsample_domain <- xsample[xsample[domain_col_name] == domain_id, , drop = FALSE]
  xsample_d_domain <- model.matrix(~., data = data.frame(xsample_domain[common_pred_vars]))
  xsample_dt_domain <- t(xsample_d_domain)
  weights_domain <- weight[which(domain_indic_vec == 1)]

  w <- as.matrix(
    weight*domain_indic_vec + (
      t(as.matrix(xpop_d_domain) - xsample_dt_domain %*% weights_domain) %*%
        comp1
    ) %*%
      comp2
  )

  t <- w %*% y

  domain_N <- unlist(xpop_domain["N"])

  if (var_est == TRUE) {

    if (var_method != "bootstrapSRS") {

      y_hat <- xsample_d_domain %*% betas
      y_domain <- y[which(domain_indic_vec == 1)]
      e <- y_domain - y_hat
      varEst <- varMase(y = e, pi = pi[which(domain_indic_vec == 1)], pi2 = pi2, method = var_method, N = domain_N)

    } else if (var_method == "boostrapSRS"){

      # need to implement

    }

    return(list(
      domain = domain_id,
      domain_total = as.numeric(t),
      domain_mean = as.numeric(t)/as.numeric(domain_N),
      domain_total_var = as.numeric(varEst),
      domain_mean_var = as.numeric(varEst)/as.numeric(domain_N^2),
      weights = w
    ))

  } else {

    return(list(
      domain = domain_id,
      domain_total = as.numeric(t),
      domain_mean = as.numeric(t)/as.numeric(domain_N),
      weights = w
    ))

  }

}

by_domain_logistic <- function(domain_id,
                               xsample,
                               xpop_d,
                               domain_col_name,
                               xpop_sums,
                               mod,
                               y,
                               var_est,
                               var_method,
                               weight,
                               pi,
                               pi2) {

  domain_indic_vec <- as.integer(xsample[domain_col_name] == domain_id)

  xpop_domain <- xpop_d[xpop_d[domain_col_name] == domain_id, , drop = FALSE]
  xsample_domain <- xsample[xsample[domain_col_name] == domain_id, , drop = FALSE]
  domain_N <- xpop_sums[xpop_sums[domain_col_name] == domain_id, ,drop = FALSE][names(xpop_d)[1]]

  y_hats_U <- as.matrix(predict(mod, newdata = xpop_domain[ , -1], type = "response", family = quasibinomial()))
  y_hats_s <- as.matrix(predict(mod, newdata = xsample_domain, type = "response", family = quasibinomial()))

  y_domain <- y[which(domain_indic_vec == 1)]
  weights_domain <- weight[which(domain_indic_vec == 1)]

  t <- t(y_domain - y_hats_s) %*% weights_domain + sum(y_hats_U)

  if (var_est == TRUE) {

    if (var_method != "bootstrapSRS") {

      e <- y_domain - y_hats_s
      varEst <- varMase(y = e, pi = pi[which(domain_indic_vec == 1)], pi2 = pi2, method = var_method, N = domain_N)

    } else if (var_method == "bootstrapSRS") {

      # need to implement

    }

    return(list(
      domain = domain_id,
      domain_total = as.numeric(t),
      domain_mean = as.numeric(t)/as.numeric(domain_N),
      domain_total_var = as.numeric(varEst),
      domain_mean_var = as.numeric(varEst)/as.numeric(domain_N^2)
    ))

  } else {

    return(list(
      domain = domain_id,
      domain_total = as.numeric(t),
      domain_mean = as.numeric(t)/as.numeric(domain_N)
    ))

  }

}


truncateText <- function(x) {
  if (length(x) > 1)
    x <- paste(x, collapse = "")
  w <- options("width")$width
  if (nchar(x) <= w)
    return(x)

  cont <- TRUE
  out <- x
  while (cont) {
    tmp <- out[length(out)]
    tmp2 <- substring(tmp, 1, w)

    spaceIndex <- gregexpr("[[:space:]]", tmp2)[[1]]
    stopIndex <- spaceIndex[length(spaceIndex) - 1] - 1
    tmp <- c(substring(tmp2, 1, stopIndex),
             substring(tmp, stopIndex + 1))
    out <-
      if (length(out) == 1)
        tmp
    else
      c(out[1:(length(x) - 1)], tmp)
    if (all(nchar(out) <= w))
      cont <- FALSE
  }

  paste(out, collapse = "\n")
}
