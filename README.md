svyEffects
================

- <a href="#how-to-install-svyeffects"
  id="toc-how-to-install-svyeffects">How to install
  <code>svyEffects</code></a>
- <a href="#introduction" id="toc-introduction">Introduction</a>
- <a href="#development-history-and-differences-from-other-packages"
  id="toc-development-history-and-differences-from-other-packages">Development
  history and differences from other packages</a>
- <a href="#binary-dependent-variable-models"
  id="toc-binary-dependent-variable-models">Binary dependent variable
  models</a>
  - <a href="#predictions-on-a-categorical-variable"
    id="toc-predictions-on-a-categorical-variable">Predictions on a
    categorical variable</a>
  - <a href="#plotting" id="toc-plotting">Plotting</a>
  - <a href="#predictions-on-a-continuous-variable"
    id="toc-predictions-on-a-continuous-variable">Predictions on a
    continuous variable</a>
- <a href="#ordinal-dependent-variable-models"
  id="toc-ordinal-dependent-variable-models">Ordinal dependent variable
  models</a>
  - <a href="#predictions-on-a-categorical-variable-1"
    id="toc-predictions-on-a-categorical-variable-1">Predictions on a
    categorical variable</a>
  - <a href="#predictions-on-a-continuous-variable-1"
    id="toc-predictions-on-a-continuous-variable-1">Predictions on a
    continuous variable</a>
- <a href="#multinomial-dependent-variable-models"
  id="toc-multinomial-dependent-variable-models">Multinomial dependent
  variable models</a>
  - <a href="#predictions-on-a-categorical-variable-2"
    id="toc-predictions-on-a-categorical-variable-2">Predictions on a
    categorical variable</a>
  - <a href="#predictions-on-a-continuous-variable-2"
    id="toc-predictions-on-a-continuous-variable-2">Predictions on a
    continuous variable</a>
- <a href="#marginal-effects-at-reasonable-values"
  id="toc-marginal-effects-at-reasonable-values">Marginal effects at
  reasonable values</a>
- <a href="#interaction-effects" id="toc-interaction-effects">Interaction
  effects</a>
- <a href="#measures-of-model-fit" id="toc-measures-of-model-fit">Measures
  of model fit</a>
- <a href="#planned-updates" id="toc-planned-updates">Planned updates</a>
- <a href="#references" id="toc-references">References</a>
- <a href="#appendix-detailed-comparisons-with-stata-results"
  id="toc-appendix-detailed-comparisons-with-stata-results">Appendix:
  Detailed comparisons with <code>Stata</code> results</a>
  - <a href="#binary-logit" id="toc-binary-logit">Binary logit</a>
  - <a href="#ordered-logit" id="toc-ordered-logit">Ordered logit</a>
  - <a href="#mulitnomial-logit" id="toc-mulitnomial-logit">Mulitnomial
    logit</a>

## How to install `svyEffects`

Run the code below in an `R` or `RStudio` session, and it will install
`{svyEffects}`. You’ll need the `{remotes}` package, if you don’t
already have it (which is why it’s included in the code). If you already
have it, then just skip the first line.

    install.packages("remotes")
    library(remotes)
    remotes::install_github("jb-santos/svyEffects, force = TRUE")

*Note: if you have installed a previous version of this package and want
to update, ensure you include the argument `force = TRUE` in the
function call because all versions of a GitHub packages are
v0.0.0.9000.*

## Introduction

An oft-cited reason why `R` is not more widely used in social science
research is its disjointed and incomplete set of tools to deal with
weights. `{svyEffects}` helps address this problem by providing a suite
of post-estimation tools for working with limited dependent variable
models (binary, ordinal, and multinomial logit) estimated on
survey-weighted data.

Its main set of functions calculate predicted probabilities using
either:

- the *average marginal effects* approach (also known as *marginal
  effects at observed values*); or
- the *marginal effects at reasonable/representative/typical values*
  approach (also known as *marginal effects for the average case*).

These approaches are analogous to `Stata`’s commands `margins x` and
`margins x, at`, respectively.

After calculating predicted probabilities, it will then calculate
differences in probabilities (also known as *contrasts*/*pairwise
comparisons* for categorical predictors and *first differences* for
continuous predictors) using:

- for continuous predictors, the change across the entire range of the
  variable (by default), or a one-unit change centred on the mean or a
  one-standard-deviation change centred on the mean; or
- for categorical predictors, all pairwise differences.

For both predictions and differences, it uses simulation methods (the
parametric bootstrap) to derive 95% confidence intervals.

It works with the following survey-weighted model objects:
`survey::svyglm` (binary logit), `survey::svyolr` (ordered logit),
`svrepmisc::svymultinom` (multinomial logit).

Eventually, support for (non-survey) weighted model objects (i.e. models
estimated with the `weight` option) will be added for `glm`,
`MASS::polr`, and `nnet::multinom`.

Also included in the package are:

- A `plot` method that creates a `ggplot` object of predicted
  probabilities or differences in predicted probabilities. This plot can
  be modified by adding further `ggplot` commands, which is shown below.
- A function `svyPRE` that calculates the proportional reductions in
  error for `svyglm` and `svyolr` models (functionality for
  `svymultinom` models in development).
- A function `mnlSig` that displays a concise summary of multinomial
  logit coefficients with statistical significance stars. This has been
  adapted for use on `svymultinom` objects from Dave Armstrong’s
  original function from `{DAMisc}`, which works for `multinom` objects.
- A snippet of the 2019 Canadian Election Study Online Survey for
  testing and demonstration purposes. This can be loaded with the
  command `data(ces19)`.

------------------------------------------------------------------------

# Development history and differences from other packages

This package extends functions originally written by Dave Armstrong,
some of which are in his `{DAMisc}` package
(<https://github.com/davidaarmstrong/damisc>).

The reporting functions and naming conventions are inspired by Daniel
Ludecke’s excellent `{ggeffects}` package
(<https://github.com/strengejacke/ggeffects>). Current users of
`{ggeffects}` will notice similarities between `{svyEffects}` and
`{ggeffects}`. However, while `{ggeffects}` can estimate MER
probabilities (what it calls *adjusted predictions*) with `svyglm`
objects, it is not compatible with either `svyolr` or `svymultinom`
objects. Moreover, `{svyEffects}` estimates true average marginal
effects, which is the estimate of a variable’s effect on a given outcome
at the population level as opposed to a variable’s effect for a
hypothetical “average case” that may or may not exist or even be
theoretically plausible. (A detailed discussion of the difference is in
Hanmer and Kalkan 2013, *AJPS*, the full citation of which can be found
in the reference section.)

Note: because AMEs run simulations on multiple copies of your dataset,
they can take much more time to calculate than MERs, particularly on
large datasets or when using an older computer. Those needing quick
results can calculate MERs (which, in practice *usually* substantively
similar to AMEs) and then decide from there for which variables they
want to calculate AMEs.

------------------------------------------------------------------------

# Binary dependent variable models

To demonstrate how this function works with binary dependent variables,
we’ll model voting for the Conservative Party of Canada versus voting
for any other party.

``` r
library(svyEffects)
data(ces19)

library(survey)
ces19_svy <- svydesign(ids = ~1, strata = NULL, weights = ~pesweight, 
                       data = ces19, digits = 3)

VOTECON <- svyglm(votecon ~ agegrp + gender + educ + region + relig + marketlib + culturetrad, 
                  design = ces19_svy, family = binomial)
summary(VOTECON)
#> 
#> Call:
#> svyglm(formula = votecon ~ agegrp + gender + educ + region + 
#>     relig + marketlib + culturetrad, design = ces19_svy, family = binomial)
#> 
#> Survey design:
#> svydesign(ids = ~1, strata = NULL, weights = ~pesweight, data = ces19, 
#>     digits = 3)
#> 
#> Coefficients:
#>                             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)                 -0.72420    0.42951  -1.686  0.09202 .  
#> agegrp35-54                  0.18969    0.30665   0.619  0.53630    
#> agegrp55+                    0.37571    0.30225   1.243  0.21407    
#> genderWoman/Other           -0.35954    0.18701  -1.923  0.05476 .  
#> educSome PSE                 0.02885    0.24444   0.118  0.90605    
#> educUni degree               0.08152    0.26961   0.302  0.76244    
#> regionWest                   0.64625    0.19603   3.297  0.00101 ** 
#> regionAtlantic               0.32897    0.35624   0.923  0.35596    
#> religCatholic                0.45252    0.26196   1.727  0.08433 .  
#> religNon-Catholic Christian  0.65283    0.22961   2.843  0.00454 ** 
#> religOther                   0.75055    0.39610   1.895  0.05834 .  
#> marketlib                    2.28579    0.30268   7.552 8.16e-14 ***
#> culturetrad                  1.96134    0.23845   8.225 4.76e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for binomial family taken to be 0.9166063)
#> 
#> Number of Fisher Scoring iterations: 5
```

The estimates from `survey::svyglm` closely resemble the ones from
`Stata`’s `logit` command with `[pweight=]` specified (see the
comparison tables in the Appendix).

## Predictions on a categorical variable

Let’s look at the effect of educational attainment (`educ`), a
categorical predictor with three levels: high school or less, some
post-secondary, and a university degree at the bachelor’s level or
higher.

The function `svyAME` will return average marginal effects for
education, or the effect of a change in education, holding all other
variables at observed values. We’ll specify a seed value for
reproducibility purposes.

The function’s output is a list that contains three data frames:

- `$preds`: predicted probabilities
- `$diffs`: differences in predicted probabilities
- `$seed`: the seed value used for the simulations

``` r
library(svyEffects)
VOTECON_educ_ame <- svyEffects::svyAME(VOTECON,
                                       varname = "educ",
                                       seed = 2019)
VOTECON_educ_ame$preds
#> # A tibble: 3 × 5
#>   educ       predicted conf.low conf.high type       
#>   <fct>          <dbl>    <dbl>     <dbl> <chr>      
#> 1 HS or less     0.404    0.340     0.468 Probability
#> 2 Some PSE       0.408    0.369     0.446 Probability
#> 3 Uni degree     0.416    0.374     0.458 Probability
```

Again, the results from `svyAME` match the `Stata` equivalent of
`margins` within 0.001 for the mean predicted probability and within
0.002 for the upper and lower 95% confidence bounds (see the relevant
table in the Appendix).

`svyAME` also calculates the differences in predicted probabilities for
all pairwise comparisons between levels of our predictor variable.

``` r
VOTECON_educ_ame$diffs
#> # A tibble: 3 × 5
#>   educ                    predicted conf.low conf.high type      
#>   <chr>                       <dbl>    <dbl>     <dbl> <chr>     
#> 1 Some PSE - HS or less     0.00397  -0.0704    0.0783 Difference
#> 2 Uni degree - HS or less   0.0125   -0.0723    0.0933 Difference
#> 3 Uni degree - Some PSE     0.00849  -0.0510    0.0682 Difference
```

The differences also compare favourably to the same results from
`Stata`. The `svyEffects` and `Stata` results differ by \<0.001.

    . lincom _b[2.educ] - _b[1.educ]

     ( 1)  - 1bn.educ + 2.educ = 0

    ------------------------------------------------------------------------------
                 |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
    -------------+----------------------------------------------------------------
             (1) |   .0044692   .0378337     0.12   0.906    -.0696834    .0786218
    ------------------------------------------------------------------------------

    . lincom _b[3.educ] - _b[1.educ]

     ( 1)  - 1bn.educ + 3.educ = 0

    ------------------------------------------------------------------------------
                 |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
    -------------+----------------------------------------------------------------
             (1) |   .0126461   .0417288     0.30   0.762    -.0691408     .094433
    ------------------------------------------------------------------------------

    . lincom _b[3.educ] - _b[2.educ]

     ( 1)  - 2.educ + 3.educ = 0

    ------------------------------------------------------------------------------
                 |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
    -------------+----------------------------------------------------------------
             (1) |   .0081769   .0302207     0.27   0.787    -.0510547    .0674084
    ------------------------------------------------------------------------------

## Plotting

The outputs of this function lend themselves well to plotting using
`{ggplot2}`. As an example, let’s plot the predicted probabilities of
voting Conservative across levels of education.

``` r
library(ggplot2)
ggplot(VOTECON_educ_ame$preds) +
  aes(x = educ,
      y = predicted,
      ymin = conf.low,
      ymax = conf.high) +
  geom_pointrange() +
  labs(title = "Probability of voting Conservative by education",
       y = "Predicted probability",
       x = "Education")
```

![](man/figures/unnamed-chunk-6-1.png)<!-- -->

For convenience, `{svyEffects}` also includes a `plot` method, which
uses the `{ggplot2}` engine to visualize either predicted probabilities
or differences in predicted probabilities.

By default, the predicted probabilities are plotted, as shown below.

``` r
plot(VOTECON_educ_ame)
```

![](man/figures/unnamed-chunk-7-1.png)<!-- -->

Note that labelling is minimal on the automatically-generated plots, but
you can add your own customization using `{ggplot2}`’s code conventions.

``` r
plot(VOTECON_educ_ame) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "My title",
       subtitle = "My subtitle",
       x = "My xvar label",
       y = "My yvar label",
       caption = "My caption") +
  theme_classic()
```

![](man/figures/unnamed-chunk-8-1.png)<!-- -->

You can also plot the differences in predicted probabilities between
levels of education by including the option `what = "diffs"` (or simply
`"diffs"`) in the `plot` function call.

Note, to do this in `Stata`, you would have to calculate each pairwise
difference with a separate command. While it is not difficult to write a
loop to do that, you would still need to output the results to a
separate matrix and then generate the plot. The functions in
`svyEffects` do that for you.

``` r
plot(VOTECON_educ_ame, "diffs") +
  geom_hline(yintercept = 0, linetype = "dotted")
```

![](man/figures/unnamed-chunk-9-1.png)<!-- -->

## Predictions on a continuous variable

Now, let’s look at the effect of market liberalism, a continuous
predictor that ranges from -1 (minimal market liberalism, or the most
left-wing position) to +1 (maximal market liberalism, or the most
right-wing position).

``` r
VOTECON_marketlib_ame <- svyAME(VOTECON,
                                varname = "marketlib",
                                seed = 2019)
VOTECON_marketlib_ame$preds
#> # A tibble: 11 × 5
#>    marketlib predicted conf.low conf.high type       
#>        <dbl>     <dbl>    <dbl>     <dbl> <chr>      
#>  1    -1         0.123   0.0749     0.182 Probability
#>  2    -0.8       0.170   0.118      0.232 Probability
#>  3    -0.6       0.230   0.180      0.281 Probability
#>  4    -0.4       0.300   0.260      0.341 Probability
#>  5    -0.2       0.379   0.345      0.412 Probability
#>  6     0         0.465   0.428      0.505 Probability
#>  7     0.200     0.551   0.501      0.605 Probability
#>  8     0.4       0.633   0.564      0.701 Probability
#>  9     0.6       0.709   0.625      0.791 Probability
#> 10     0.8       0.778   0.685      0.858 Probability
#> 11     1         0.831   0.731      0.909 Probability
```

`svyAME` also produces very similar results to `Stata` for the effect of
market liberalism (see Appendix).

``` r
VOTECON_marketlib_ame$diffs
#> # A tibble: 1 × 5
#>   marketlib              predicted conf.low conf.high type      
#>   <chr>                      <dbl>    <dbl>     <dbl> <chr>     
#> 1 Delta (range) : -1 - 1     0.709    0.550     0.827 Difference
plot(VOTECON_marketlib_ame)
```

![](man/figures/unnamed-chunk-11-1.png)<!-- -->

Note that, because the function returns a first difference for
continuous predictors, the graph is not any more illuminating than the
summary statistic.

``` r
plot(VOTECON_marketlib_ame, "diffs")
```

![](man/figures/unnamed-chunk-12-1.png)<!-- -->

------------------------------------------------------------------------

# Ordinal dependent variable models

To demonstrate ordinal dependent variables, we’ll model feeling
thermometer ratings for the leader of the Conservative Party of Canada.
This variable usually ranges from 0 to 100. But, for this example, we’ll
used a collapsed ordinal measure of “cold” (0-39), “lukewarm” (40-59),
and “hot” (60-100).

``` r
data(ces19)

library(survey)
ces19_svy <- svydesign(ids = ~1, strata = NULL, weights = ~pesweight, 
                        data = ces19, digits = 3)

CONLDR <- svyolr(ftconldr ~ agegrp + gender + educ + region + relig + marketlib + culturetrad, 
                 design = ces19_svy)
summary(CONLDR)
#> Call:
#> svyolr(ftconldr ~ agegrp + gender + educ + region + relig + marketlib + 
#>     culturetrad, design = ces19_svy)
#> 
#> Coefficients:
#>                                   Value Std. Error    t value
#> agegrp35-54                  0.17796372  0.2665702  0.6676055
#> agegrp55+                    0.04601115  0.2686809  0.1712483
#> genderWoman/Other           -0.49866393  0.1707825 -2.9198766
#> educSome PSE                -0.20802018  0.2231985 -0.9319963
#> educUni degree              -0.10211120  0.2430419 -0.4201383
#> regionWest                   0.41572273  0.1792462  2.3192838
#> regionAtlantic              -0.09551799  0.2884670 -0.3311228
#> religCatholic                0.46488829  0.2245240  2.0705501
#> religNon-Catholic Christian  0.60848101  0.2202608  2.7625480
#> religOther                   0.79950865  0.3459580  2.3109994
#> marketlib                    1.96343442  0.2572774  7.6315846
#> culturetrad                  1.95794984  0.2263782  8.6490220
#> 
#> Intercepts:
#>               Value   Std. Error t value
#> Cold|Lukewarm -0.4716  0.3817    -1.2356
#> Lukewarm|Hot   0.1141  0.3704     0.3079
```

## Predictions on a categorical variable

Here’s the effect of region on feelings towards the Conservative Party
leader.

For brevity, only the visualizations of the predicted
probabilities/differences are presented, along with comparisons versus
`Stata`. Detailed comparison tables can be seen in the Appendix.

``` r
CONLDR_region_ame <- svyAME(CONLDR,
                            varname = "region",
                            seed = 2019)
plot(CONLDR_region_ame)
```

![](man/figures/unnamed-chunk-14-1.png)<!-- -->

For ordinal and multinomial probabilities, the plot method follows the
conventions used by the `{ggeffects}` package (i.e. facetting by
response level). But, you can re-create the `Stata` default of
colour-coding the response level by writing your own `ggplot` command,
as shown below.

``` r
ggplot(CONLDR_region_ame$preds) + 
  aes(x = region, y = predicted, ymin = conf.low, ymax = conf.high, colour = y) +
  geom_pointrange(position = position_dodge2(.35)) +
  scale_y_continuous(limits = c(.05,.62)) +
  scale_colour_viridis_d() +
  labs(title = "Effect of region on Conservative leader rating (svyEffects)",
       x = "Region",
       y = "Predicted probability",
       colour = "Conservative \nleader rating") +
  theme_bw()
```

![](man/figures/unnamed-chunk-15-1.png)<!-- -->

The predicted probabilities for region are very similar to the `Stata`
results.

![](images/conldr_region.png)

``` r
plot(CONLDR_region_ame, "diffs") +
  geom_hline(yintercept = 0, linetype = "dotted")
```

![](man/figures/unnamed-chunk-16-1.png)<!-- -->

## Predictions on a continuous variable

Here’s the effect of market liberalism:

``` r
CONLDR_marketlib_ame <- svyAME(CONLDR,
                            varname = "marketlib",
                            diffchange = "range",
                            seed = 2019)
plot(CONLDR_marketlib_ame)
```

![](man/figures/unnamed-chunk-17-1.png)<!-- -->

``` r
plot(CONLDR_marketlib_ame, "diffs")
```

![](man/figures/unnamed-chunk-18-1.png)<!-- -->

The graph below shows the results similar to how `Stata` plots them.

``` r
ggplot(CONLDR_marketlib_ame$preds) +
  aes(x = marketlib, y = predicted, ymin = conf.low, ymax = conf.high, colour = y, fill = y) +
  geom_line() +
  geom_ribbon(colour = "transparent", alpha = 0.2) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, length = 6)) +
  scale_x_continuous(breaks = seq(-1, 1, length = 6)) +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  labs(title = "Effect of market liberalism on Conservative leader ratings (svyEffects)",
       x = "Market liberalism (least to most)",
       y = "Predicted probability",
       fill = "Conservative \nleader rating",
       colour = "Conservative \nleader rating") +
  theme_bw()
```

![](man/figures/unnamed-chunk-19-1.png)<!-- -->

The `Stata` results are very similar.

![](images/conldr_marketlib.png)

------------------------------------------------------------------------

# Multinomial dependent variable models

To demonstrate multinomial dependent variables, we’ll model vote choice
in the 2019 Canadian Federal Election. To keep things simple, we’ll
limit our analysis to the three major parties (the Liberals,
Conservatives, and New Democrats) and exclude the province of Quebec
(which has a different party system and patterns of vote choice).

There is no way to directly estimate a multinomial model with the
`{survey}` package in `R`. The package `{svyrepmisc}` generates an
approximation by turning the weighting scheme into replicate weights and
estimating the model with those. It uses the jackknife to calculate
variances.

We’ll go through this process step-by-step. First, we’ll import the
data, do some data cleaning, and then create our usual survey-design
object.

``` r
data(ces19)
library(survey)
ces19_svy <- svydesign(ids = ~1, strata = NULL, weights = ~pesweight, 
                        data = ces19, digits = 3)
```

Now, we’ll use the function `as.svrepdesign()` from `{survey}` to turn
our sampling weights into replicate weights with variances calculated
using the jackknife.

``` r
ces19_svy_r <- as.svrepdesign(ces19_svy, type = "JK1")
```

After our survey design object with replicate weights and jackknife
variances is created, we can use the function `svymultinom` from
`{svyrepmisc}` to run our vote choice model.

Note: use the option `trace = FALSE` in the `svymultinom` function call
to suppress the reporting of each replication (similar to using the
option `quietly` in Stata).

Included with `{svyEffects}` the function `mnlSig`, which displays
coefficients from multinomial logit models and flags statistically
significant ones. `mnlSig` is adapted from Dave Armstrong’s original
function from his `{DAMisc}` package.

``` r
# remotes::install_github("carlganz/svrepmisc")
library(svrepmisc)

VOTE <- svymultinom(vote ~ agegrp + gender + educ + region + relig + marketlib + culturetrad, 
                    design = ces19_svy_r, trace = FALSE)
mnlSig(VOTE)
#>                             Conservative     NDP
#> (Intercept)                      -0.281  -0.764 
#> agegrp35-54                       0.075  -0.286 
#> agegrp55+                         0.052  -0.959*
#> genderWoman/Other                -0.267   0.266 
#> educSome PSE                      0.195   0.486 
#> educUni degree                    0.032  -0.166 
#> regionWest                        0.898*  0.751*
#> regionAtlantic                    0.334   0.017 
#> religCatholic                     0.418  -0.049 
#> religNon-Catholic Christian       0.581* -0.181 
#> religOther                        0.440  -1.083*
#> marketlib                         2.113* -0.516 
#> culturetrad                       1.990*  0.076
```

## Predictions on a categorical variable

For our post-estimation command, we’ll need to specify a few more
options because `svymultinom` does not store them in its output. These
are:

- `design`: the survey design object used to estimate the model; and
- `modform`: the model formula used in the `svymultinom` call (in the
  form `modform = "y ~ x1 + x2 + x3"`).

Here’s the effect of education:

``` r
VOTE_region_ame <- svyAME(
  VOTE,
  varname = "region",
  weightvar = "pesweight",
  seed = 2019,
  design = ces19_svy_r,
  modform = "vote ~ agegrp + gender + educ + region + relig + marketlib + culturetrad")
VOTE_region_ame$preds
#> # A tibble: 9 × 6
#>   y            region   predicted conf.low conf.high type       
#>   <fct>        <fct>        <dbl>    <dbl>     <dbl> <chr>      
#> 1 Liberal      Ontario      0.445    0.397     0.492 Probability
#> 2 Conservative Ontario      0.361    0.320     0.404 Probability
#> 3 NDP          Ontario      0.194    0.156     0.234 Probability
#> 4 Liberal      West         0.285    0.241     0.335 Probability
#> 5 Conservative West         0.461    0.418     0.503 Probability
#> 6 NDP          West         0.254    0.213     0.298 Probability
#> 7 Liberal      Atlantic     0.408    0.291     0.534 Probability
#> 8 Conservative Atlantic     0.410    0.295     0.525 Probability
#> 9 NDP          Atlantic     0.182    0.116     0.266 Probability
```

``` r
plot(VOTE_region_ame)
```

![](man/figures/unnamed-chunk-24-1.png)<!-- -->

The predicted probabilities by region in `Stata` format:

``` r
ggplot(VOTE_region_ame$preds) +
  aes(x = region, y = predicted, ymin = conf.low, ymax = conf.high, colour = y) +
  geom_pointrange(position = position_dodge2(.2)) +
  scale_y_continuous(limits = c(.1,.55), breaks = c(.1,.2,.3,.4,.5)) +
  scale_colour_manual(values =  c("maroon", "navy", "orange")) +
  labs(title = "Effect of region on vote choice (svyEffects)",
       x = "Region",
       y = "Predicted probability",
       colour = "Vote choice") +
  theme_bw()
```

![](man/figures/unnamed-chunk-25-1.png)<!-- -->

And here are the results from `Stata`:

![](images/vote_region.png)

``` r
plot(VOTE_region_ame, "diffs") +
  geom_hline(yintercept = 0, linetype = "dotted")
```

![](man/figures/unnamed-chunk-26-1.png)<!-- -->

## Predictions on a continuous variable

Here’s the effect of market liberalism:

``` r
VOTE_marketlib_ame <- svyAME(
  VOTE,
  varname = "marketlib",
  weightvar = "pesweight",
  seed = 2019,
  diffchange = "range",
  design = ces19_svy_r,
  modform = "vote ~ agegrp + gender + educ + region + relig + marketlib + culturetrad")
VOTE_marketlib_ame$preds
#> # A tibble: 33 × 6
#>    y            marketlib predicted conf.low conf.high type       
#>    <fct>            <dbl>     <dbl>    <dbl>     <dbl> <chr>      
#>  1 Liberal           -1       0.499   0.406      0.598 Probability
#>  2 Conservative      -1       0.121   0.0727     0.186 Probability
#>  3 NDP               -1       0.380   0.283      0.483 Probability
#>  4 Liberal           -0.8     0.490   0.421      0.568 Probability
#>  5 Conservative      -0.8     0.169   0.117      0.232 Probability
#>  6 NDP               -0.8     0.341   0.268      0.417 Probability
#>  7 Liberal           -0.6     0.472   0.420      0.530 Probability
#>  8 Conservative      -0.6     0.229   0.179      0.284 Probability
#>  9 NDP               -0.6     0.300   0.249      0.353 Probability
#> 10 Liberal           -0.4     0.443   0.403      0.487 Probability
#> # … with 23 more rows
plot(VOTE_marketlib_ame)
```

![](man/figures/unnamed-chunk-27-1.png)<!-- -->

The effect of market liberalism graphed in “`Stata` format”:

``` r
ggplot(VOTE_marketlib_ame$preds) +
  aes(x = marketlib, y = predicted, ymin = conf.low, ymax = conf.high, colour = y, fill = y) +
  geom_line() +
  geom_ribbon(colour = "transparent", alpha = 0.2) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, length = 6)) +
  scale_x_continuous(breaks = seq(-1, 1, length = 6)) +
  scale_colour_manual(values = c("maroon", "navy", "orange")) +
  scale_fill_manual(values = c("maroon", "navy", "orange")) +
  labs(title = "Effect of market liberalism on vote choice (svyEffects)",
       x = "Market liberalism (least to most)",
       y = "Predicted probability",
       fill = "Vote choice",
       colour = "Vote choice") +
  theme_bw()
```

![](man/figures/unnamed-chunk-28-1.png)<!-- -->

Here are the results from `Stata`:

![](images/vote_marketlib.png)

Finally, here are the differences:

``` r
plot(VOTE_marketlib_ame, "diffs")
```

![](man/figures/unnamed-chunk-29-1.png)<!-- -->

------------------------------------------------------------------------

# Marginal effects at reasonable values

*(documentation in progress)*

You can choose to calculate predicted probabilities and differences
using the “marginal effects at reasonable/representative values” (MER)
approach by using the `svyMER` function, which uses the same arguments
as `svyAME`. This command would give you the estimated effect of a
variable “for the ‘typical’ case” (which may or may not be typical,
plausible, or even possible in the real world) as opposed to the effect
of a variable across the population (see Hanmer and Kalkan 2013 for an
in-depth discussion).

MER probabilities are *usually* very similar to AME probabilities, but
not always. However, they are *much* faster to calculate using
simulation methods because they do not use all cases in the data set.

------------------------------------------------------------------------

# Interaction effects

*(documentation in progress)*

Both `svyAME` and `svyMER` support calculating predicted probabilities
of combinations of two predictor variables. This can be done by using
the argument `byvar = "x"` in the function call. This works with or
without a product term.

This will not return differences in predicted probabilities. For
limited-dependent variable models, one would need to calculate a second
difference to test for the significance of an interaction between two
variables, either from the inclusion of a product term or through the
compression inherent in these types of models (see Norton, Wang, and Ai
2004).

Dave Armstrong’s `{DAMisc}` package
(<https://github.com/davidaarmstrong/damisc/>) has an `R` port for
Norton, Wang, and Ai’s original `Stata` function, and this will
eventually be ported to `{svyEffects}` and adapted for use in
survey-weighted models.

``` r
VOTECON_educ_region <- svyAME(VOTECON,
                              varname = "educ",
                              byvar = "region",
                              seed = 2019)
VOTECON_educ_region$preds
#> # A tibble: 9 × 6
#>   educ       region   predicted conf.low conf.high type       
#>   <fct>      <fct>        <dbl>    <dbl>     <dbl> <chr>      
#> 1 HS or less Ontario      0.359    0.289     0.429 Probability
#> 2 Some PSE   Ontario      0.363    0.314     0.411 Probability
#> 3 Uni degree Ontario      0.371    0.322     0.422 Probability
#> 4 HS or less West         0.458    0.388     0.527 Probability
#> 5 Some PSE   West         0.463    0.414     0.515 Probability
#> 6 Uni degree West         0.471    0.417     0.525 Probability
#> 7 HS or less Atlantic     0.408    0.286     0.532 Probability
#> 8 Some PSE   Atlantic     0.411    0.308     0.511 Probability
#> 9 Uni degree Atlantic     0.420    0.315     0.526 Probability
plot(VOTECON_educ_region)
```

![](man/figures/unnamed-chunk-30-1.png)<!-- -->

``` r
VOTECON_marketlib_educ <- svyAME(VOTECON,
                                 varname = "marketlib",
                                 byvar = "educ",
                                 seed = 2019)
VOTECON_marketlib_educ$preds
#> # A tibble: 33 × 5
#>    educ       marketlib predicted conf.low conf.high
#>    <fct>          <dbl>     <dbl>    <dbl>     <dbl>
#>  1 HS or less    -1         0.120   0.0705     0.185
#>  2 HS or less    -0.8       0.167   0.110      0.240
#>  3 HS or less    -0.6       0.226   0.163      0.298
#>  4 HS or less    -0.4       0.297   0.227      0.366
#>  5 HS or less    -0.2       0.375   0.301      0.450
#>  6 HS or less     0         0.458   0.378      0.548
#>  7 HS or less     0.200     0.545   0.454      0.643
#>  8 HS or less     0.4       0.627   0.524      0.727
#>  9 HS or less     0.6       0.702   0.587      0.808
#> 10 HS or less     0.8       0.771   0.653      0.874
#> # … with 23 more rows
plot(VOTECON_marketlib_educ)
```

![](man/figures/unnamed-chunk-31-1.png)<!-- -->

``` r
VOTECON2 <- svyglm(votecon ~ agegrp + gender + educ + region + relig + marketlib + culturetrad +
                     marketlib:educ, 
                   design = ces19_svy, family = binomial)
VOTECON2_marketlib_educ <- svyAME(VOTECON2,
                                  varname = "marketlib",
                                  byvar = "educ",
                                  seed = 2019)
VOTECON2_marketlib_educ$preds
#> # A tibble: 33 × 5
#>    educ       marketlib predicted conf.low conf.high
#>    <fct>          <dbl>     <dbl>    <dbl>     <dbl>
#>  1 HS or less    -1         0.145   0.0492     0.295
#>  2 HS or less    -0.8       0.185   0.0868     0.314
#>  3 HS or less    -0.6       0.240   0.146      0.352
#>  4 HS or less    -0.4       0.304   0.226      0.394
#>  5 HS or less    -0.2       0.377   0.307      0.449
#>  6 HS or less     0         0.456   0.365      0.548
#>  7 HS or less     0.200     0.532   0.406      0.659
#>  8 HS or less     0.4       0.608   0.442      0.764
#>  9 HS or less     0.6       0.676   0.473      0.848
#> 10 HS or less     0.8       0.738   0.501      0.910
#> # … with 23 more rows
plot(VOTECON2_marketlib_educ)
```

![](man/figures/unnamed-chunk-32-1.png)<!-- -->

------------------------------------------------------------------------

# Measures of model fit

*(documentation in progress)*

The `svyPRE` function calculates the proportional reduction in error for
binary (`svyglm`) and ordered (`svyolr`) logit models.

Functionality for multinomial (`svyrepstatmisc`) models is currently
under development.

``` r
svyPRE(VOTECON)
#> # A tibble: 3 × 2
#>   Measure                      Value
#>   <chr>                        <dbl>
#> 1 Percent in modal category    0.592
#> 2 Percent correctly classified 0.747
#> 3 Percent reduction in error   0.380
```

``` r
svyPRE(CONLDR)
#> # A tibble: 3 × 2
#>   Measure                      Value
#>   <chr>                        <dbl>
#> 1 Percent in modal category    0.477
#> 2 Percent correctly classified 0.698
#> 3 Percent reduction in error   0.422
```

------------------------------------------------------------------------

# Planned updates

This package is under active development, and updates will include:

1.  Support for (non-survey) weighted models and non-weighted models.
    While there are other packages that do this, some do not return
    confidence intervals for predictions for some model types. And, to
    my knowledge, none use simulation methods to derive confidence
    intervals. *Note: You can actually already do this with
    `{svyEffects}` by creating a survey design object with a weight of
    “1”, but it would be good to avoid having to use that workaround.*
2.  Expand functionality with `svyrepstatmisc` model objects to reduce
    the number of arguments needed to run the functions.
3.  Support for using an alternative variance-covariance matrix using
    `{sandwich}`. This would only be for binary logit models because
    `{sandwich}` does not play nice with ordinal or multinomial models.
    That said, survey-weighted models do adjust the variance-covariance
    matrix (the documentation for `{survey}` does not specify the
    correction method it uses, but it appears to be HC0, based on what
    I’ve seen).
4.  A second differences function to test for the significance of a
    two-way interaction.
5.  (Eventually) Use of the delta method to calculate confidence
    intervals for AMEs probabilities to speed up computational time.

------------------------------------------------------------------------

# References

Armstrong, Dave. 2022. *DAMisc: Dave Armstrong’s Miscellaneous
Functions.* R package version 1.7.2.

Hanmer, M.J. and K.O. Kalkan. 2013. “Behind the Curve: Clarifying the
Best Approach to Calculating Predicted Probabilities and Marginal
Effects from Limited Dependent Variable Models.” *American Journal of
Political Science*. 57(1): 263-277.

Norton, Edward C., Hua Wang and Chunrong Ai. 2004. Computing Interaction
Effects and Standard Errors in Logit and Probit Models. *The Stata
Journal* 4(2): 154-167.

Rainey, Carlisle. 2016. “Compression and Conditional Effects: A Product
Term Is Essential When Using Logistic Regression to Test for
Interaction.” *Political Science Research and Methods* 4(3): 621-639.

Stephenson, Laura B; Harell, Allison; Rubenson, Daniel; Loewen, Peter
John, 2020, “2019 Canadian Election Study - Online Survey,”
<https://doi.org/10.7910/DVN/DUS88V>, Harvard Dataverse, V1.

------------------------------------------------------------------------

# Appendix: Detailed comparisons with `Stata` results

## Binary logit

Parameter estimates and standard errors

<div style="font-size:9pt">


``` r
logit_b
#> # A tibble: 13 × 5
#>    term                            R_b Stata_b  R_se Stata_se
#>    <chr>                         <dbl>   <dbl> <dbl>    <dbl>
#>  1 (Intercept)                 -0.724  -0.724  0.430    0.430
#>  2 agegrp35-54                  0.190   0.190  0.307    0.307
#>  3 agegrp55+                    0.376   0.376  0.302    0.302
#>  4 genderWoman/Other           -0.360  -0.360  0.187    0.187
#>  5 educSome PSE                 0.0288  0.0289 0.244    0.244
#>  6 educUni degree               0.0815  0.0815 0.270    0.270
#>  7 regionWest                   0.646   0.646  0.196    0.196
#>  8 regionAtlantic               0.329   0.329  0.356    0.356
#>  9 religCatholic                0.453   0.453  0.262    0.262
#> 10 religNon-Catholic Christian  0.653   0.653  0.230    0.230
#> 11 religOther                   0.751   0.751  0.396    0.396
#> 12 marketlib                    2.29    2.29   0.303    0.303
#> 13 culturetrad                  1.96    1.96   0.238    0.238
```

/div\>

Predicted probabilities (AMEs)

<div style="font-size:9pt">

``` r
logit_ame
#> # A tibble: 14 × 8
#>    term      y            R_b Stata_b  R_lwr Stata_lwr R_upr Stata_upr
#>    <chr>     <chr>      <dbl>   <dbl>  <dbl>     <dbl> <dbl>     <dbl>
#>  1 educ      HS or less 0.404   0.403 0.340     0.338  0.468     0.467
#>  2 educ      Some PSE   0.408   0.407 0.369     0.368  0.446     0.446
#>  3 educ      Uni degree 0.416   0.415 0.374     0.373  0.458     0.458
#>  4 marketlib -1         0.123   0.118 0.0639    0.0749 0.182     0.172
#>  5 marketlib -0.8       0.170   0.166 0.110     0.118  0.232     0.222
#>  6 marketlib -0.6       0.230   0.226 0.175     0.180  0.281     0.278
#>  7 marketlib -0.4       0.300   0.298 0.256     0.260  0.341     0.340
#>  8 marketlib -0.2       0.379   0.379 0.345     0.345  0.412     0.412
#>  9 marketlib 0          0.465   0.465 0.427     0.428  0.505     0.502
#> 10 marketlib 0.2        0.551   0.551 0.498     0.501  0.605     0.605
#> 11 marketlib 0.4        0.633   0.635 0.565     0.564  0.701     0.706
#> 12 marketlib 0.6        0.709   0.712 0.629     0.625  0.791     0.796
#> 13 marketlib 0.8        0.778   0.780 0.691     0.685  0.858     0.869
#> 14 marketlib 1          0.831   0.837 0.750     0.731  0.909     0.925
```

</div>

## Ordered logit

Parameter estimates and standard errors

<div style="font-size:9pt">


``` r
ologit_b
#> # A tibble: 14 × 5
#>    term                            R_b Stata_b  R_se Stata_se
#>    <chr>                         <dbl>   <dbl> <dbl>    <dbl>
#>  1 agegrp35-54                  0.178   0.178  0.267    0.267
#>  2 agegrp55+                    0.0460  0.0460 0.269    0.269
#>  3 genderWoman/Other           -0.499  -0.499  0.171    0.171
#>  4 educSome PSE                -0.208  -0.208  0.223    0.223
#>  5 educUni degree              -0.102  -0.102  0.243    0.243
#>  6 regionWest                   0.416   0.416  0.179    0.179
#>  7 regionAtlantic              -0.0955 -0.0955 0.288    0.288
#>  8 religCatholic                0.465   0.465  0.225    0.225
#>  9 religNon-Catholic Christian  0.608   0.608  0.220    0.220
#> 10 religOther                   0.800   0.800  0.346    0.346
#> 11 marketlib                    1.96    1.96   0.257    0.257
#> 12 culturetrad                  1.96    1.96   0.226    0.226
#> 13 Cold|Lukewarm               -0.472  -0.472  0.382    0.382
#> 14 Lukewarm|Hot                 0.114   0.114  0.370    0.370
```

/div\>

Predicted probabilities (AMEs)

<div style="font-size:9pt">

``` r
ologit_ame
#> # A tibble: 42 × 9
#>    term      y        x           R_b Stata_b  R_lwr Stata_lwr R_upr Stata_upr
#>    <chr>     <chr>    <chr>     <dbl>   <dbl>  <dbl>     <dbl> <dbl>     <dbl>
#>  1 region    Cold     Ontario  0.502   0.501  0.461     0.460  0.544     0.542
#>  2 region    Cold     West     0.434   0.0956 0.393     0.0748 0.475     0.116
#>  3 region    Cold     Atlantic 0.517   0.403  0.433     0.365  0.603     0.442
#>  4 region    Hot      Ontario  0.403   0.517  0.364     0.431  0.440     0.602
#>  5 region    Hot      West     0.471   0.0954 0.429     0.0748 0.512     0.116
#>  6 region    Hot      Atlantic 0.389   0.388  0.303     0.302  0.474     0.473
#>  7 region    Lukewarm Ontario  0.0950  0.434  0.0741    0.393  0.116     0.475
#>  8 region    Lukewarm West     0.0945  0.0950 0.0737    0.0743 0.115     0.116
#>  9 region    Lukewarm Atlantic 0.0943  0.471  0.0740    0.429  0.115     0.513
#> 10 marketlib Cold     -1       0.758   0.761  0.681     0.689  0.829     0.834
#> # … with 32 more rows
```

</div>

## Mulitnomial logit

For this section, there are three comparisons:

1)  `R`
2)  `Stata` using `svy jackknife`
3)  `Stata` using `pweight`

Parameter estimates and standard errors

<div style="font-size:9pt">


``` r
mlogit_b
#> # A tibble: 26 × 8
#>    term                      y         R_b Stata…¹ Stata…²  R_se Stata…³ Stata…⁴
#>    <chr>                     <chr>   <dbl>   <dbl>   <dbl> <dbl>   <dbl>   <dbl>
#>  1 (Intercept)               CPC   -0.281  -0.281  -0.281  0.529   0.529   0.493
#>  2 agegrp35-54               CPC    0.0752  0.0752  0.0752 0.358   0.358   0.340
#>  3 agegrp55+                 CPC    0.0517  0.0517  0.0517 0.364   0.364   0.343
#>  4 genderWoman/Other         CPC   -0.267  -0.267  -0.267  0.209   0.209   0.200
#>  5 educSome PSE              CPC    0.195   0.195   0.195  0.277   0.277   0.262
#>  6 educUni degree            CPC    0.0319  0.0319  0.0319 0.307   0.307   0.290
#>  7 regionWest                CPC    0.898   0.898   0.898  0.222   0.222   0.213
#>  8 regionAtlantic            CPC    0.334   0.334   0.334  0.422   0.422   0.384
#>  9 religCatholic             CPC    0.418   0.418   0.418  0.298   0.298   0.283
#> 10 religNon-Catholic Christ… CPC    0.581   0.581   0.581  0.255   0.255   0.245
#> # … with 16 more rows, and abbreviated variable names ¹​Statasvyjk_b,
#> #   ²​Statapweight_b, ³​Statasvyjk_se, ⁴​Statapweight_se
```

/div\>

Predicted probabilities (AMEs)

<div style="font-size:9pt">

``` r
mlogit_ame
#> # A tibble: 42 × 12
#>    term   y     x     R_pred svyjk…¹ pweig…² R_lwr svyjk…³ pweig…⁴ R_upr svyjk…⁵
#>    <chr>  <chr> <chr>  <dbl>   <dbl>   <dbl> <dbl>   <dbl>   <dbl> <dbl>   <dbl>
#>  1 marke… Libe… -1     0.499   0.505   0.505 0.406  0.409   0.413  0.598   0.602
#>  2 marke… Libe… -0.8   0.490   0.496   0.496 0.421  0.422   0.425  0.568   0.570
#>  3 marke… Libe… -0.6   0.472   0.477   0.477 0.420  0.421   0.424  0.530   0.532
#>  4 marke… Libe… -0.4   0.443   0.447   0.447 0.403  0.404   0.406  0.487   0.490
#>  5 marke… Libe… -0.2   0.406   0.408   0.408 0.368  0.370   0.372  0.444   0.447
#>  6 marke… Libe… 0      0.362   0.363   0.363 0.318  0.319   0.322  0.406   0.406
#>  7 marke… Libe… 0.2    0.313   0.313   0.313 0.259  0.259   0.262  0.370   0.368
#>  8 marke… Libe… 0.4    0.264   0.262   0.262 0.201  0.196   0.200  0.333   0.329
#>  9 marke… Libe… 0.6    0.216   0.213   0.213 0.147  0.138   0.142  0.294   0.288
#> 10 marke… Libe… 0.8    0.173   0.168   0.168 0.103  0.0879  0.0924 0.257   0.247
#> # … with 32 more rows, 1 more variable: pweight_upr <dbl>, and abbreviated
#> #   variable names ¹​svyjk_pred, ²​pweight_pred, ³​svyjk_lwr, ⁴​pweight_lwr,
#> #   ⁵​svyjk_upr
```

</div>
