% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/svyPRE.R
\name{svyPRE.svyglm}
\alias{svyPRE.svyglm}
\title{Proportional reduction in error of survey-weighted binary logit models}
\usage{
\method{svyPRE}{svyglm}(obj, ...)
}
\arguments{
\item{obj}{Model object of class \code{svyglm} with \code{family = "binomial"}}

\item{...}{Other options (not currently implemented)}
}
\value{
A tibble
}
\description{
Calculates proportional reduction in error (PRE) for a survey-weighted
binary logit model estimated using \code{survey::svyglm}.
}
\examples{
data(ces19)
library(survey)
ces19_svy <- svydesign(ids = ~1, strata = NULL, weights = ~pesweight,
  data = ces19, digits = 3)
VOTECON <- svyglm(votecon ~ agegrp + gender + educ + region + marketlib,
  design = ces19_svy, family = binomial)
svyPRE(VOTECON)

}
\author{
John Santos & Dave Armstrong
}
