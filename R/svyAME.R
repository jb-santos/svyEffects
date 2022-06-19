#' @title Average marginal effects for survey-weighted binary, ordinal, and
#' multinomial logit models
#'
#'
#' @description Calculates predicted probabilities for a predictor holding all
#' other observations at observed values. Uses simulation methods (the parametric
#' bootstrap) to derive 95% confidence intervals. This is a true "average marginal
#' effect" or a population-level estimate of a predictor's effect, as described
#' by Hanmer and Kalkan 2013.
#'
#' Note: for multinomial logit models, you will also need use arguments that
#' specify the survey design and the model formula. Please consult the
#' documentation for the method \code{svyAME.svrepstatmisc} for more details
#' (accessible by running the command \code{?svyAME.svrepstatmisc}).
#'
#'
#' @param obj A model object. Currently, this package supports models of class
#' \code{svyglm} and \code{svyolr} from the \pkg{survey} package, and models of
#' class \code{svrepstatmisc} from the \pkg{svrepmisc} package (see note below).
#' @param varname Character string denoting the name of the predictor variable
#' for which effects are to be calculated.
#' @param weightvar Character string denoting the name of the sampling weight
#' variable.
#' @param nvals Scalar denoting the sequence length spanning the range of a
#' continuous variable for which effects are to be calculated (default: 11).
#' @param diffchange Character string  denoting over what change in x a first
#' difference is to be calculated for a continuous predictor (default: range).
#' @param byvar For interaction effects, a character string denoting the name of
#' the moderator variable.
#' @param bynvals For interaction effects with a numerical moderator, a scalar
#' denoting the sequence length  spanning the range of the moderator variable
#' for which effects are to be calculated (default: 3).
#' @param sims Scalar denoting how many simulations to conduct (default: 1500).
#' @param seed Seed value for random number generator. By default, the function
#' picks a random integer between 1 and 1,000,000. If you save the output of
#' this function, it will save the seed value used for simulations in the slot
#' \code{$seed}.
#' @param ... Other arguments, depending on the type of model with which you're
#' working. For example, for survey-weight multinomial logit models of class
#' \code{svrepstatmisc}, you will need to specify the survey design object and
#' the model formula.
#'
#'
#' @return A list with dataframes:
#'   \describe{
#'     \item{\code{$preds}}{predicted probabilities}
#'     \item{\code{$diffs}}{differences in predicted probabilities}
#'     \item{\code{$seed}}{seed value used for random number generator}
#'     }
#'
#'
#' @author John Santos & Dave Armstrong
#'
#'
#' @references Hanmer, M.J. and K.O. Kalkan. 2013. \sQuote{Behind the Curve:
#' Clarifying the Best Approach to Calculating Predicted Probabilities and
#' Marginal Effects from Limited Dependent Variable Models}. American Journal
#' of Political Science.  57(1): 263-277.
#'
#'
#' @importFrom dplyr arrange as_tibble filter group_by mutate rename relocate select summarize tibble
#' @importFrom magrittr %>%
#' @importFrom MASS mvrnorm
#' @importFrom rlang sym
#' @importFrom stats as.formula coef formula model.frame model.matrix na.omit plogis pnorm quantile reformulate runif terms vcov weighted.mean
#' @importFrom stringr boundary str_extract_all str_split
#' @importFrom survey svydesign svymean svyquantile svytable
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom utils combn relist
#'
#'
#' @export
svyAME <- function(obj,
                   varname,
                   weightvar,
                   nvals = 11,
                   diffchange = c("range", "unit", "sd"),
                   byvar = NULL,
                   bynvals = 3,
                   sims = 1500,
                   seed = NULL,
                   ...) {UseMethod("svyAME")}