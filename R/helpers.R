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




#' Break Apart Model Formula
#'
#' A sort of inverse of the \code{reformulate} function.
#'
#' @param form A formula
#' @param keep_env Logical indicating whether the formula's
#' environment should be returned with the result
#'
#' @description Works as a sort of inverse to \code{reformulate}
#' by breaking apart the formula into response and the term labels.
#' It also returns the variable names of all of the variables
#' implicated in the formula. (Borrowed from DAMisc; not exported).
#'
#' @return A list with \code{termlabels} giving the rhs terms of the
#' model, \code{response} give the lhs of the model, \code{env} optionally
#' giving the environment of the formula and \code{vars} a vector of the
#' variable names implicated in the formula
#'
#' @author Dave Armstrong
#'
unformulate <- function(form, keep_env=FALSE){
  rhs <- attr(terms(form), "term.labels")
  vn <- all.vars(form)
  l <- as.list(form)
  if(length(l) == 2){
    lhs <- NULL
  }
  if(length(l) == 3){
    lhs <- as.character(l[[2]])
  }
  if(!(length(l) %in% 2:3)){
    stop("formula must transform into a two- or three-element list\n")
  }
  res <- list(termlabels = rhs, response = lhs, vars = vn)
  if(keep_env){
    res$env <- environment(form)
  }
  res
}
