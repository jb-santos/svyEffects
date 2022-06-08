# PRIMARY GENERIC FUNCTIONS
# This script contains the S3 generic functions. These functions dispatch
# methods, which are in separate scripts, depending on what kind of dependent
# variable they model--i.e. binary (glm), ordinal (olr), and multinomial
# (multinom).






#' S3 generic function for survey-weighted average marginal effects
#'
#' @param obj Model object on which to conduct post-estimation
#' @param ... Other arguments (currently not implemented).
#'
#' @export
#'
svyAME <- function(obj, ...) {UseMethod("svyAME")}






#' S3 generic function for survey-weighted marginal effects at reasonable values
#'
#' @param obj Model object on which to conduct post-estimation
#' @param ... Other arguments (currently not implemented).
#' @export
#'
svyMER <- function(obj, ...) {UseMethod("svyMER")}
