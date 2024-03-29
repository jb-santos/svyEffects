% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/svyPRE.R
\name{svyR2}
\alias{svyR2}
\title{R-squared of survey-weighted linear models}
\usage{
svyR2(obj, ...)
}
\arguments{
\item{obj}{Model object of class \code{svyglm} with \code{family = "gaussian"}}

\item{...}{Other options (not currently implemented)}
}
\value{
A scalar denoting the weighted R-squared for the model
}
\description{
R-squared of survey-weighted linear models
}
\examples{
data(ces19)
library(survey)
ces19_svy <- svydesign(ids = ~1, strata = NULL, weights = ~pesweight, data = ces19, digits = 3)
mod <- svyglm(culturetrad ~ agegrp + gender + educ + region + marketlib, design = ces19_svy)
svyR2(mod)

}
\author{
John Santos & Dave Armstrong
}
