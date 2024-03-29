% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mnlSig.R
\name{mnlSig}
\alias{mnlSig}
\title{Print Statistically Significant MNL Coefficients}
\usage{
mnlSig(obj, ...)
}
\arguments{
\item{obj}{Model object for which to display model summary}

\item{...}{Other arguments (not currently implemented).}
}
\description{
By default, the summary for objects of class \code{multinom} is not
particularly helpful.  It still requires a lot of work on the part of the
user to figure out which coefficients are significantly different from zero
and which ones are not.  \code{mnlSig} solves this problem by either
flagging significant coefficients with an asterisk or only printing
significant coefficients, leaving insignificant ones blank.
}
\examples{
\dontrun{
data(ces19)
library(survey)
ces19_svy <- svydesign(ids = ~1, strata = NULL, weights = ~pesweight,
  data = ces19, digits = 3)
ces19_svy_r <- as.svrepdesign(ces19_svy, type = "JK1")
# remotes::install_github("carlganz/svrepmisc") # (if not already installed)
library(svrepmisc)
VOTE <- svymultinom(vote ~ agegrp + gender + educ + region + marketlib,
  design = ces19_svy_r, trace = FALSE)
mnlSig(VOTE)
}


}
\author{
Dave Armstrong
}
