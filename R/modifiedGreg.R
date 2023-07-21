#' @export modifiedGreg
#' @import survey
#' @import glmnet
#' @import boot
#' @import utils
#' @importFrom stats model.matrix predict quasibinomial var aggregate as.formula
modifiedGreg <- function(y, ...)  UseMethod("modifiedGreg")

#' @export
#' @rawNamespace export(modifiedGreg.numeric)
modifiedGreg.numeric <- function(y,
                                 xsample,
                                 xpop,
                                 domains,
                                 pi = NULL,
                                 pi2 = NULL,
                                 datatype = "raw",
                                 model = "linear",
                                 var_est = F,
                                 var_method = "LinHB",
                                 modelselect = FALSE,
                                 lambda = "lambda.min",
                                 domain_col_name = NULL,
                                 estimation_domains = NULL,
                                 N = NULL,
                                 messages = T,
                                 ...) {


  .args <- c(as.list(environment()), list(...))
  do.call(validate_modifiedGreg, .args)
  funcCall <- match.call()

  if (is.null(N)) {
    if (datatype == "raw") {
      N <- nrow(xpop)
    } else {
      N <- sum(xpop$N)
    }
  }

  if (is.null(domain_col_name)) {
    if (datatype == "raw") {
      domain_col_name <- setdiff(names(xpop), names(xsample))
    } else {
      domain_col_name <- setdiff(names(xpop), c(names(xsample), "N"))
    }
    if (messages) {
      message(paste0("domain_col_name is not directly specified. ", domain_col_name, " is being used."))
    }
  }


  if (is.null(pi)) {
    pi <- rep(length(y)/N, length(y))
  }

  weight <- as.vector(pi^(-1))
  y <- as.vector(y)

  if (is.null(estimation_domains)) {
    estimation_domains <- unique(xpop[[domain_col_name]])
  }

  # creating a vector of common auxiliary variable names
  common_pred_vars <- intersect(names(xsample), names(xpop))

  # creating the design matrix for entire xsample
  xsample_d <- model.matrix(~., data = data.frame(xsample[common_pred_vars]))
  xsample <- cbind(data.frame(xsample_d[,-1, drop = FALSE]), domains)
  names(xsample) <- c(colnames(xsample_d[,-1, drop = FALSE]), domain_col_name)
  xsample_dt <- t(xsample_d)

  # variable selection

  if (modelselect == TRUE) {
    if(model == "linear"){
      fam <- "gaussian"
    } else{
      fam <- "binomial"
    }

    x <- xsample[ , -which(names(xsample) == domain_col_name)]
    cv <- cv.glmnet(x = as.matrix(x), y = y, alpha = 1, weights = weight, nfolds = 10, family = fam, standardize = FALSE)

    if (lambda == "lambda.min") {
      lambda.opt <- cv$lambda.min
    }

    if (lambda == "lambda.1se") {
      lambda.opt <- cv$lambda.1se
    }

    pred_mod <- glmnet(x = as.matrix(x), y = y, alpha = 1, family = fam, standardize = FALSE, weights = weight)
    lasso_coef <- predict(pred_mod, type = "coefficients", s = lambda.opt)[1:dim(xsample_d)[2],]
    coef_select <- names(lasso_coef[lasso_coef != 0])[-1]

    if (length(coef_select) == 0) {
      if (messages) {
        message("No variables selected in the model selection stage.  Fitting a HT estimator.")
      }

      if (var_est == TRUE) {
        HT <- horvitzThompson(y = y, pi = pi, N = N, pi2 = pi2, var_est = TRUE, var_method = var_method)

        return(list(pop_total = HT$pop_total,
                    pop_mean = HT$pop_total/N,
                    pop_total_var = HT$pop_total_var,
                    pop_mean_var = HT$pop_total_var/N^2,
                    weights = as.vector(pi^(-1))))
      } else {
        HT <- horvitzThompson(y = y, pi = pi, N = N, pi2 = pi2, var_est = FALSE)

        return(list(pop_total = HT$pop_total,
                    pop_mean = HT$pop_total/N,
                    weights = as.vector(pi^(-1))))
      }
    } else {
      xsample <- xsample[ , c(coef_select, domain_col_name), drop = FALSE]
      xsample_d <- model.matrix(~., data = xsample[ , coef_select, drop = FALSE])
      xsample_dt <- t(xsample_d)
      xsample <- cbind(data.frame(xsample_d[,-1, drop = FALSE]), domains)
      names(xsample) <- c(colnames(xsample_d[,-1, drop = FALSE]), domain_col_name)
    }
  }


  if (model == "linear") {
    if (datatype == "raw"){
      # only using columns that occur in xsample too
      xpop_subset <- xpop[common_pred_vars]

      # expand factors and interaction terms and then re-add the domain column
      xpop <- cbind(data.frame(model.matrix(~. -1, data = data.frame(xpop_subset))), xpop[domain_col_name])
      xpop$`N` <- 1

      # sum auxiliary pop values by domain
      bydomain_formula <- as.formula(paste0(". ~", domain_col_name))
      xpop_d <- aggregate(bydomain_formula, xpop, FUN = sum)
      # move N column to the front (where it would normally be)
      xpop_d <- xpop_d[, c(ncol(xpop_d), 1:(ncol(xpop_d) - 1))]
    }

    if (datatype == "totals"){
      xpop_d <- xpop[ ,c("N", common_pred_vars, domain_col_name)]
    }

    if (datatype == "means"){
      xpop[common_pred_vars] <- lapply(xpop[common_pred_vars], function(x) x*(xpop$N))
      xpop_d <- cbind(N = xpop$N, xpop[ ,!(names(xpop) %in% "N")])
    }

    # computing the pieces that remain the same across all domains
    constant_component1 <- solve(xsample_dt %*% diag(weight) %*% xsample_d)
    constant_component2 <- t(weight * xsample_d)
    betas <- solve(xsample_dt %*% diag(weight) %*% xsample_d) %*% (xsample_dt) %*% diag(weight) %*% y

    res <- lapply(estimation_domains,
                  FUN = by_domain_linear,
                  xsample,
                  xpop_d,
                  domain_col_name,
                  constant_component1,
                  constant_component2,
                  betas,
                  common_pred_vars,
                  y,
                  var_est,
                  var_method,
                  weight,
                  pi,
                  pi2)

  } else if (model == "logistic") {

    if(length(levels(as.factor(y))) != 2){
      stop("Function can only handle categorical response with two categories.")
    }

    xpop_subset <- xpop[common_pred_vars]
    xpop <- cbind(data.frame(model.matrix(~., data = data.frame(xpop_subset))), xpop[domain_col_name])

    # need N by domain
    bydomain_formula <- as.formula(paste0(names(xpop)[1], "~", domain_col_name))
    xpop_sums <- aggregate(bydomain_formula, xpop, FUN = sum)

    # preparing logistic model
    xsample_preds <- xsample[ , !(names(xsample) %in% domain_col_name)]
    dat <- data.frame(y, weight, xsample_preds)
    colnames(dat) <- c("y", "weight", names(xsample_preds))
    f <- paste(names(dat)[1], "~", paste(names(dat)[-c(1,2)], collapse = " + "))
    s_design <- survey::svydesign(ids = ~1, weights = ~weight, data = dat)
    mod <- survey::svyglm(f, design = s_design, family = quasibinomial())
    betas <- mod$coefficients

    res <- lapply(estimation_domains,
                  FUN = by_domain_logistic,
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
                  pi2)

  }

  pop_res <- do.call(
    rbind, lapply(
      res, FUN = function(x) data.frame(
        pop_total = x$domain_total,
        pop_totalvar = x$domain_total_var
      )
    )
  ) |>
    colSums()

  out <- list(
    domain_level_res = res,
    pop_res = pop_res,
    coefs = betas,
    call = funcCall,
    modeltype = model,
    var_est = var_est,
    var_method = var_method,

  )

  class(out) <- "modifiedGreg"
  out

}

#' @export
print.modifiedGreg <- function(obj, ...) {

  cat("\nCall:\n",
      truncateText(deparse(obj$call, width.cutoff = 500)),
      "\n\n", sep = "")

  cat(paste0("Model Type: ", obj$modeltype))
  cat("\n")

  if (obj$var_est) {
    cat(paste0("Variance Method: ", obj$var_method))
    cat("\n")
  }



}

#' @export
summary.modifiedGreg <- function(obj, ...) {

  x <- obj$domain_level_res
  by_domain_tab <- do.call(
    rbind, lapply(
      x, FUN = function(t) data.frame(
        domain = t$domain,
        domain_total = t$domain_total,
        domain_totalvar = t$domain_total_var
      )
    )
  )

  pop_tab <- obj$pop_res
  call <- obj$call

  out <- list(call = call,
              pop_tab = pop_tab,
              by_domain_table = by_domain_tab)

  class(out) <- "summary.modifiedGreg"
  out

}

#' @export
print.summary.modifiedGreg <- function(x, ...) {

  cat("\nCall:\n",
      truncateText(deparse(x$call, width.cutoff = 500)),
      "\n\n",
      sep = "")

  print(head(x$by_domain_tab))
  cat("\n")
  print(x$pop_tab)
  cat("\n")
  invisible(x)
}





