% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/svyMER.R
\name{svyMER}
\alias{svyMER}
\title{Marginal effects at reasonable values for survey-weighted binary,
ordinal, and multinomial logit models}
\usage{
svyMER(
  obj,
  varname,
  weightvar,
  nvals = 11,
  diffchange = c("range", "unit", "sd"),
  byvar = NULL,
  bynvals = 3,
  sims = 2500,
  seed = NULL,
  ...
)
}
\arguments{
\item{obj}{A model object. Currently, this package supports models of class
\code{svyglm} and \code{svyolr} from the \pkg{survey} package, and models of
class \code{svrepstatmisc} from the \pkg{svrepmisc} package (see note below).}

\item{varname}{Character string denoting the name of the predictor variable
for which effects are to be calculated.}

\item{weightvar}{Character string denoting the name of the sampling weight
variable.}

\item{nvals}{Scalar denoting the sequence length spanning the range of a
continuous variable for which effects are to be calculated (default: 11).}

\item{diffchange}{Character string  denoting over what change in x a first
difference is to be calculated for a continuous predictor (default: range).}

\item{byvar}{For interaction effects, a character string denoting the name of
the moderator variable.}

\item{bynvals}{For interaction effects with a numerical moderator, a scalar
denoting the sequence length  spanning the range of the moderator variable
for which effects are to be calculated (default: 3).}

\item{sims}{Scalar denoting how many simulations to conduct (default: 2500).}

\item{seed}{Seed value for random number generator. By default, the function
picks a random integer between 1 and 1,000,000. If you save the output of
this function, it will save the seed value used for simulations in the slot
\code{$seed}.}

\item{...}{Other arguments, depending on the type of model with which you're
working. For example, for survey-weight multinomial logit models of class
\code{svrepstatmisc}, you will need to specify the survey design object and
the model formula. (Use the command \code{?svyMER.svrepstatmisc} for more information)}
}
\value{
A list with dataframes:
\describe{
\item{\code{$preds}}{predicted probabilities}
\item{\code{$diffs}}{differences in predicted probabilities}
\item{\code{$seed}}{seed value used for random number generator}
\item{\code{$typical}}{the values at which variables are set for the simulations}
}
}
\description{
Calculates predicted probabilities for an independent variable,
holding all other observations at reasonable/representative/typical values.
This involves setting all continuous variables to their medians and all
categorical variables to their modes. Uses simulation methods (the parametric
bootstrap) to derive 95\% confidence intervals. This is the effect of a predictor
for "the average case", as opposed to a true "average marginal effect", as
described by Hanmer and Kalkan 2013.

Note: for multinomial logit models, you will also need use arguments that
specify the survey design and the model formula. Please consult the
documentation for the method \code{svyMER.svrepstatmisc} for more details
(accessible by running the command \code{?svyMER.svrepstatmisc}).
}
\examples{
# Binary models:
data(ces19)
library(survey)
ces19_svy <- svydesign(ids = ~1, strata = NULL, weights = ~pesweight,
  data = ces19, digits = 3)
VOTECON <- svyglm(votecon ~ agegrp + gender + educ + region + marketlib,
  design = ces19_svy, family = binomial)
svyMER(VOTECON, varname = "educ", weightvar = "pesweight", seed = 2019)
svyMER(VOTECON, varname = "marketlib", weightvar = "pesweight", seed = 2019)

# Ordinal models:
data(ces19)
library(survey)
ces19_svy <- svydesign(ids = ~1, strata = NULL, weights = ~pesweight,
  data = ces19, digits = 3)
CONLDR <- svyolr(ftconldr ~ agegrp + gender + educ + region + marketlib,
  design = ces19_svy)
svyMER(CONLDR, varname = "region", weightvar = "pesweight", seed = 2019)
svyMER(CONLDR, varname = "marketlib", weightvar = "pesweight", seed = 2019)


# Multinomial models:
data(ces19)
library(survey)
ces19_svy <- svydesign(ids = ~1, strata = NULL, weights = ~pesweight,
  data = ces19, digits = 3)
ces19_svy_r <- as.svrepdesign(ces19_svy, type = "JK1")
# remotes::install_github("carlganz/svrepmisc")  # (if not already installed)
library(svrepmisc)
VOTE <- svymultinom(vote ~ agegrp + gender + educ + region + marketlib,
  design = ces19_svy_r, trace = FALSE)
svyMER(VOTE, varname = "region", weightvar = "pesweight",
  seed = 2019, design = ces19_svy_r,
  modform = "vote ~ agegrp + gender + educ + region + marketlib")
svyMER(VOTE, varname = "marketlib", weightvar = "pesweight",
  seed = 2019, design = ces19_svy_r,
  modform = "vote ~ agegrp + gender + educ + region + marketlib")

}
\references{
Hanmer, M.J. and K.O. Kalkan. 2013. \sQuote{Behind the Curve:
Clarifying the Best Approach to Calculating Predicted Probabilities and
Marginal Effects from Limited Dependent Variable Models}. American Journal
of Political Science. 57(1): 263-277.
}
\author{
John Santos & Dave Armstrong
}
