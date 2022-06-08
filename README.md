svyEffects
================

## Developer’s note

The primary functions of this package are complete. However, the
functions in this GitHub repository are not. As I convert functions from
their standalone versions to ones that fit within the overall package
framework, they will be added to GitHub and the documentation will be
fleshed out.

## Introduction

`svyEffects` contains post-estimation functions for GLMs (or limited
dependent variable models) estimated on survey-weighted data.

Primarily, it calculates predicted probabilities using either:

-   the “average marginal effects” approach (also known as “marginal
    effects at observed values”, or “adjusted predictions”); and
-   the “marginal effects at reasonable/representative/typical values”
    approach (also known as “marginal effects for the average case”).

After calculating predicted probabilities, it will then calculate
differences in probabilities (also known as “contrasts” or “first
differences”) using:

-   for continuous variables, either a one-unit or
    one-standard-deviation change centred on the mean, or the change
    across the entire range of the variable; or
-   for categorical variables, all pairwise differences.

For both predictions and differences, it uses simulation methods (the
parametric bootstrap) to derive confidence intervals.

It works with the following survey-weighted models or (non-survey)
weighted models (i.e. models estimated with the `weight=` option
enabled):

-   binary dependent variable models
    -   `survey::svyglm`
    -   `glm`
-   ordinal dependent variable models
    -   `MASS::polr`
    -   `survey::svyolr`
-   multinomial dependent variable models
    -   `svrepmisc::svymultinom`
    -   `nnet::multinom`
