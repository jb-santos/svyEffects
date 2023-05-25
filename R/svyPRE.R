#' Proportional reduction in error of survey-weighted binary and ordered logit models
#'
#' Calculates proportional reduction in error (PRE) for survey-weighted
#' binary logit models or ordered logit models estimated using \code{survey::svyglm}
#'
#'
#' @param obj Model object of class \code{svyglm} with \code{family = "binomial"}
#'  or \code{svyolr}
#' @param ... Other options (not currently implemented)
#'
#' @return A tibble
#' @export
#'
#' @author John Santos & Dave Armstrong
#'
#' @examples
#' data(ces19)
#' library(survey)
#' ces19_svy <- svydesign(ids = ~1, strata = NULL, weights = ~pesweight,
#'   data = ces19, digits = 3)
#' VOTECON <- svyglm(votecon ~ agegrp + gender + educ + region + marketlib,
#'   design = ces19_svy, family = binomial)
#' svyPRE(VOTECON)
#'
svyPRE <- function(obj, ...) {UseMethod("svyPRE")}




#' Proportional reduction in error of survey-weighted binary logit models
#'
#' Calculates proportional reduction in error (PRE) for a survey-weighted
#' binary logit model estimated using \code{survey::svyglm}.
#'
#'
#' @param obj Model object of class \code{svyglm} with \code{family = "binomial"}
#'  or \code{svyolr}
#' @param ... Other options (not currently implemented)
#'
#' @return A tibble
#' @export
#'
#' @author John Santos & Dave Armstrong
#'
#' @examples
#' data(ces19)
#' library(survey)
#' ces19_svy <- svydesign(ids = ~1, strata = NULL, weights = ~pesweight,
#'   data = ces19, digits = 3)
#' VOTECON <- svyglm(votecon ~ agegrp + gender + educ + region + marketlib,
#'   design = ces19_svy, family = binomial)
#' svyPRE(VOTECON)
#'
svyPRE.svyglm <- function(obj, ...) {

  if(!(family(obj)$link %in% c("logit"))) {
      stop("PRE only calculated for svyglm models of binomial
           family with logit link.\n")
  }

  data <- survey::svydesign(ids = ~1, strata = NULL,
                            weights = ~`(weights)`,
                            data = model.frame(obj))
  data$variables$y <- as.numeric(obj$y)
  data$variables$yhat <- as.numeric(plogis(predict(obj)) > .5)
  data$variables$correct <- as.numeric(data$variables$yhat == data$variables$y)
  pmc <- survey::svymean(~y, data)[1]  # percent in modal category
  pmc <- ifelse(pmc < .5, 1 - pmc, pmc)
  pcp <- survey::svymean(~correct, data)[1]   # percent correctly classified
  pre <- (pcp - pmc) / (1 - pmc)   # proportional reduction in error

  pre_stats <- tibble::tibble(
    Measure = c("Percent in modal category",
                "Percent correctly classified",
                "Percent reduction in error"),
    Value = c(pmc, pcp, pre))
  attributes(pre_stats)$model <- formula(obj)

  # more detailed output (not currently implemented)
  # out <- list(
  #   pre_stats = pre_stats,
  #   y = data$variables$y,
  #   fitted.values = as.data.frame(plogis(predict(obj)))[1],
  #   yhat = data$variables$yhat,
  #   model.frame = data
  # )

  return(pre_stats)
}




#' Proportional reduction in error of survey-weighted ordered logit models
#'
#' @param obj Model object of class \code{svyglm} with \code{family = "binomial"}
#'   or \code{svyolr}
#' @param ... Other options (not currently implemented)
#'
#' @return A tibble
#' @export
#'
#' @author John Santos & Dave Armstrong
#'
#' @examples
#' data(ces19)
#' library(survey)
#' ces19_svy <- svydesign(ids = ~1, strata = NULL, weights = ~pesweight,
#'   data = ces19, digits = 3)
#' CONLDR <- svyolr(ftconldr ~ agegrp + gender + educ + region + marketlib,
#'   design = ces19_svy)
#' svyPRE(CONLDR)
#'
svyPRE.svyolr <- function(obj, ...) {
  data <- survey::svydesign(ids = ~1, strata = NULL,
                            weights = ~`(weights)`,
                            data = model.frame(obj))
  data$variables$y <- model.frame(obj)[[1]]
  fits <- fitted.values(obj)
  data$variables$yhat <- colnames(fits)[apply(fits, 1, which.max)]
  data$variables$correct <- as.numeric(data$variables$yhat == data$variables$y)
  pmc <- max(survey::svymean(~y, data))   # percent in modal category
  pcp <- survey::svymean(~correct, data)[1]   # percent correctly classified
  pre <- (pcp - pmc) / (1 - pmc)   # proportional reduction in error

  pre_stats <- tibble::tibble(
    Measure = c("Percent in modal category",
                "Percent correctly classified",
                "Percent reduction in error"),
    Value = c(pmc, pcp, pre))
  attributes(pre_stats)$model <- formula(obj)

  # more detailed output (not currently implemented)
  # out <- list(
  #   pre_stats = pre_stats,
  #   y = data$variables$y,
  #   fitted.values = fits,
  #   yhat = data$variables$yhat,
  #   model.frame = data
  # )

  return(pre_stats)
}




#' R-squared of survey-weighted linear models
#'
#' @param obj Model object of class \code{svyglm} with \code{family = "gaussian"}
#' @param ... Other options (not currently implemented)
#'
#' @return A scalar denoting the weighted R-squared for the model
#' @export
#'
#' @author John Santos & Dave Armstrong
#'
#' @examples
#' data(ces19)
#' library(survey)
#' ces19_svy <- svydesign(ids = ~1, strata = NULL, weights = ~pesweight, data = ces19, digits = 3)
#' mod <- svyglm(culturetrad ~ agegrp + gender + educ + region + marketlib, design = ces19_svy)
#' svyR2(mod)
#'
svyR2 <- function(obj, ...) {UseMethod("svyR2")}




#' R-squared of survey-weighted linear models
#'
#' @param obj Model object of class \code{svyglm} with \code{family = "gaussian"}
#' @param ... Other options (not currently implemented)
#'
#' @return A scalar denoting the weighted R-squared for the model
#' @export
#'
#' @author John Santos & Dave Armstrong
#'
#' @examples
#' data(ces19)
#' library(survey)
#' ces19_svy <- svydesign(ids = ~1, strata = NULL, weights = ~pesweight, data = ces19, digits = 3)
#' mod <- svyglm(culturetrad ~ agegrp + gender + educ + region + marketlib, design = ces19_svy)
#' svyR2(mod)
#'
svyR2.svyglm <- function(obj, ...) {

  if(!(family(obj)$family %in% c("gaussian"))) {
    stop("R-squared only calculated for svyglm models of gaussian family.\n")
  }

  data <- survey::svydesign(ids = ~1, strata = NULL,
                            weights = ~`(weights)`,
                            data = model.frame(obj))
  data$variables$y <- as.numeric(obj$y)
  data$variables$yhat <- fitted.values(obj)
  data$variables$y_z <- (data$variables$y - svymean(~y, data)[1]) /
                           sqrt(svyvar(~y, data)[1])
  data$variables$yhat_z <- (data$variables$yhat - svymean(~yhat, data)[1]) /
                           sqrt(svyvar(~yhat, data)[1])

  r <- survey::svyglm(y_z ~ yhat_z, data)
  r2 <- as.numeric(coef(r)[2]^2)
  attributes(r2)$model <- formula(obj)

  return(r2)
}
