#' @title Print Statistically Significant MNL Coefficients
#'
#' @description By default, the summary for objects of class \code{multinom} is not
#' particularly helpful.  It still requires a lot of work on the part of the
#' user to figure out which coefficients are significantly different from zero
#' and which ones are not.  \code{mnlSig} solves this problem by either
#' flagging significant coefficients with an asterisk or only printing
#' significant coefficients, leaving insignificant ones blank.
#'
#' @param obj Model object for which to display model summary
#' @param ... Other arguments (not currently implemented).
#'
#' @importFrom stringr boundary str_extract_all str_split
#'
#' @author Dave Armstrong
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(ces19)
#' library(survey)
#' ces19_svy <- svydesign(ids = ~1, strata = NULL, weights = ~pesweight,
#'   data = ces19, digits = 3)
#' ces19_svy_r <- as.svrepdesign(ces19_svy, type = "JK1")
#' # remotes::install_github("carlganz/svrepmisc") # (if not already installed)
#' library(svrepmisc)
#' VOTE <- svymultinom(vote ~ agegrp + gender + educ + region + marketlib,
#'   design = ces19_svy_r, trace = FALSE)
#' mnlSig(VOTE)
#' }
#'
#'
mnlSig <- function(obj, ...) {UseMethod("mnlSig")}





#' @title mnlSig method for \code{svrepstatmisc} objects
#'
#' @description This is an adaptation of Dave Armstrong's original \code{mnlsig}
#' function that has been adapted for survey-weight multinomial logit model
#' objects of class \code{svrepstatmisc}.
#'
#' Original description from \code{DAMisc}: By default, the summary for objects
#' of class \code{multinom} is not particularly helpful. It still requires a lot
#' of work on the part of the user to figure out which coefficients are
#' significantly different from zero and which ones are not. \code{mnlSig} solves
#' this problem by either flagging significant coefficients with an asterisk or
#' only printing significant coefficients, leaving insignificant ones blank.
#'
#' @param obj A model object of class \code{svrepstatmisc}.
#' @param pval The desired Type I error rate to identify coefficients as
#' statistically significant.
#' @param two.sided Logical indicating whether calculated p-values should be
#' two-sided (if \code{TRUE}) or one-sided (if \code{FALSE}).
#' @param flag.sig Logical indicating whether an asterisk should be placed
#' beside coefficients which are significant at the \code{pval} level.
#' @param insig.blank Logical indicating whether coefficients which are not
#' significant at the \code{pval} level should be blank in the output.
#' @param ... Other arguments (not currently implemented).
#'
#' @return A data frame suitable for printing with the (optionally
#' significance-flagged) coefficients from a multinomial logit model.
#'
#' @author Dave Armstrong & John Santos
#'
#' @importFrom stringr boundary fixed str_extract_all str_split
#'
#' @export
mnlSig.svrepstatmisc <- function (obj,
                                  pval = 0.05,
                                  two.sided = TRUE,
                                  flag.sig = TRUE,
                                  insig.blank = FALSE,
                                  ...) {

  b <- as.data.frame(obj)$Coefficient
  se <- as.data.frame(obj)$SE
  t  =  b / se
  p  =  (2 ^ as.numeric(two.sided)) * pnorm(abs(t), lower.tail = FALSE)
  nYlev <- length(grep("(Intercept)",  rownames(as.data.frame(obj))))
  nterms <- nrow(as.data.frame(obj)) / nYlev
  termnames <- unlist(stringr::str_split(
    names(obj)[seq.int(from = 1, to = length(b), by = nYlev)],
    fixed("."))) [seq.int(from = 2, to = ((length(b) / nYlev) * 2), by = 2)]
  Ylevnames <- unlist(stringr::str_split(
    names(obj)[seq.int(from = 1, to = nYlev)],
    fixed("."))) [seq.int(from = 1, to = (nYlev * 2), by = 2)]

  b_mat <- matrix(sprintf("%.3f", b), ncol = nterms)
  sig.vec <- c(" ", "*")
  sig.obs <- as.numeric(p < pval) + 1
  if (flag.sig) {
    b_mat <- matrix(paste(b_mat, sig.vec[sig.obs], sep = ""), ncol = nterms)
  }
  if (insig.blank) {
    b_mat[which(p > pval, arr.ind = TRUE)] <- ""
  }
  rownames(b_mat) <- Ylevnames
  colnames(b_mat) <- termnames
  mod_summary <- t(b_mat)
  return(as.data.frame(mod_summary))
}





#' @title mnlSig method for \code{multinom} objects
#'
#' @param obj A model object of class \code{multinom}.
#' @param pval The desired Type I error rate to identify coefficients as
#' statistically significant.
#' @param two.sided Logical indicating whether calculated p-values should be
#' two-sided (if \code{TRUE}) or one-sided (if \code{FALSE}).
#' @param flag.sig Logical indicating whether an asterisk should be placed
#' beside coefficients which are significant at the \code{pval} level.
#' @param insig.blank Logical indicating whether coefficients which are not
#' significant at the \code{pval} level should be blank in the output.
#' @return A data frame suitable for printing with the (optionally
#' significance-flagged) coefficients from a multinomial logit model.
#' @param ... Other arguments (not currently implemented).
#'
#' @author Dave Armstrong
#'
#' @export
mnlSig.multinom <- function (obj,
                             pval = 0.05,
                             two.sided = TRUE,
                             flag.sig = TRUE,
                             insig.blank = FALSE,
                             ...) {
  smulti <- summary(obj)
  multi.t <- smulti$coefficients/smulti$standard.errors
  multi.p <- (2^as.numeric(two.sided)) * pnorm(abs(multi.t),
                                               lower.tail = FALSE)
  b <- matrix(sprintf("%.3f", smulti$coefficients), ncol = ncol(multi.t))
  sig.vec <- c(" ", "*")
  sig.obs <- as.numeric(multi.p < pval) + 1
  if (flag.sig) {
    b <- matrix(paste(b, sig.vec[sig.obs], sep = ""), ncol = ncol(multi.t))
  }
  if (insig.blank) {
    b[which(multi.p > pval, arr.ind = TRUE)] <- ""
  }
  rownames(b) <- rownames(multi.t)
  colnames(b) <- colnames(multi.t)
  b <- as.data.frame(b)
  return(b)
}
