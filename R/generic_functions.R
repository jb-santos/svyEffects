# PRIMARY GENERIC FUNCTIONS
# This script contains the S3 generic functions. These functions dispatch
# methods, which are in separate scripts, depending on what kind of dependent
# variable they model--i.e. binary (svyglm), ordinal (svyolr), and multinomial
# (svymultinom).






#' @title Average marginal effects for binary, ordinal, and multinomial
#' dependent variable models of survey-weighted data
#'
#' @description Calculates predicted probabilities for a predictor holding all
#' other observations at observed values (i.e. a true average marginal effect).
#' Uses simulations (the parametric bootstrap) to derive 95% confidence intervals.
#' Note: you will need to specify different arguments depending on whether your
#' model has a binary, ordinal, or categorical dependent variable, so consult
#' the documentation for this function's submethods for more details.
#'
#' @param obj A model object. Currently, this package supports models of class
#' \code{svyglm} and \code{svyolr} from the \pkg{survey} package, and models of
#' class \code{svrepstatmisc} from the \pkg{svrepmisc} package (see note below).
#' Models of class \code{glm}, \code{polr} (from \pkg{MASS}), and \code{multinom}
#' (from \pkg{nnet}) that specify the \code{weight=} option also work.
#' @param varname Character string denoting the name of the predictor variable
#' for which effects are to be calculated.
#' @param ... Other arguments (currently not implemented).
#'
#' @export
#'
svyAME <- function(obj,
                   varname,
                   ...) {UseMethod("svyAME")}






#' @title Marginal effects at reasonable values for binary, ordinal, and
#' multinomial dependent variable models of survey-weighted data
#'
#' @description Calculates predicted probabilities for an independent variable,
#' holding all other observations at reasonable/representative/typical values
#' (i.e. for an "average case"). This involves setting all continuous variables
#' to their medians and all categorical variables to their modes. Uses
#' simulations (the parametric bootstrap) to derive 95% confidence intervals.
#' Note: you will need to specify different arguments depending on whether your
#' model has a binary, ordinal, or categorical dependent variable, so consult
#' the documentation for this function's submethods for more details.
#'
#' @param obj A model object. Currently, this package supports models of class
#' \code{svyglm} and \code{svyolr} from the \pkg{survey} package, and models of
#' class \code{svrepstatmisc} from the \pkg{svrepmisc} package (see note below).
#' Models of class \code{glm}, \code{polr} (from \pkg{MASS}), and \code{multinom}
#' (from \pkg{nnet}) that specify the \code{weight=} option also work.
#' @param varname Character string denoting the name of the predictor variable
#' for which effects are to be calculated.
#' @param ... Other arguments used by svyMER methods.
#'
#' @export
#'
svyMER <- function(obj,
                   varname,
                   ...) {UseMethod("svyMER")}





