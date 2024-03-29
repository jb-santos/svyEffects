#' @title Print method for \pkg{svyEffects} objects
#' @name print
#'
#' @description Prints results for \code{svyEffects} of objects
#'
#' @param x An object of class \pkg{svyEffects} generated by the either of
#' the functions \code{svyAME} or \code{svyMER}.
#' @param digits = Scalar indicating desired precision (in number digits to the right of the decimal point).
#' @param nrow = Scalar indicating the maximum number of rows to display (useful for interactions).
#' @param ... Other arguments passed on the output.
#'
#' @importFrom dplyr arrange as_tibble filter group_by mutate mutate_if rename relocate select summarize tibble
#' @importFrom glue glue
#' @importFrom magrittr %>%
#' @importFrom rlang sym
#'
#' @return A list
#'
#' @author John Santos
#'
#' @export
#'
#'
#' @examples
#' data(ces19)
#' library(survey)
#' library(ggplot2)
#' ces19_svy <- svydesign(ids = ~1, strata = NULL, weights = ~pesweight,
#'   data = ces19, digits = 3)
#' VOTECON <- svyglm(votecon ~ agegrp + gender + educ + region + marketlib,
#'   design = ces19_svy, family = binomial)
#' ame_educ <- svyAME(VOTECON, varname = "educ", weightvar = "pesweight", seed = 2019)
#' print(ame_educ)
print.svyEffects <- function(x,
                             digits = 3,
                             nrow = 40,
                             ...) {

  print(glue("Dependent variable: ", attr(x, "depvar")))
  print(glue("
             "))
  print(glue("Predictor variable: ", attr(x, "predvar")))
  if("depvar" %in% names(attributes(x))) {
    print(glue("Moderator variable: ", attr(x, "byvar")))
    print(glue("
               "))
  }
  print(glue("Seed value: ", x$seed))
  print(glue("
             "))
  print(glue("Number of simulations: ", x$sims))
  print(glue("
             "))
  print(glue("Model formula: "))
  print(x$formula)
  print(glue("

             "))
  print(glue("Predicted probabilities: "))
  x$preds %>%
    dplyr::mutate_if(is.numeric, round, digits = digits) %>%
    print(n = nrow)
  print(glue("

             "))
  if(!is.null(x$diffs)) {
    print(glue("Differences in predicted probabilities: "))
    x$diffs %>%
      mutate_if(is.numeric, round, digits = digits) %>%
      print(n = nrow)
  }
}

