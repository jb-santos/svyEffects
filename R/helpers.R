# HELPER FUNCTIONS (not exported)



#' Lower 95% CI quantile
#'
#' @param x Value for which the 5th percentile should be returned.
#'
#' @return A numeric indicating the 2.5% quantile
#'
low <- function(x) {unname(quantile(x, .025))}



#' Upper 95% CI quantile
#'
#' @param x Value for which the 95th percentile should be returned.
#'
#' @return A numeric indicating the 97.5% quantile
#'
high <- function(x) {unname(quantile(x, .975))}



#' Find central/average value for all variables (for MERs)
#'
#' @param x Termlist.
#' @param design Survey design object.
#' @param ... Other arguments (currently not implemented).
#'
#' @author Dave Armstrong
#'
#' @return A list for input into the svyMER functions.
#'
svymode <- function(x, design, ...) {
  levs <- levels(design$variables[[x]])
  tab <- survey::svytable(reformulate(x), design)
  wmx <- unname(which.max(tab))
  factor(wmx,
         levels = seq_along(levs),
         labels = levs)
}
