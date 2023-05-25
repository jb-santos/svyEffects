svyEffects
================

## How to install `svyEffects`

Run the code below in an `R` or `RStudio` session, and it will install
`{svyEffects}`. You’ll need the `{remotes}` package, if you don’t
already have it (which is why it’s included in the code). If you already
have it, then just skip the first line.

    install.packages("remotes")
    library(remotes)
    remotes::install_github("jb-santos/svyEffects")

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
  effects at observed values*, or *adjusted predictions*); or
- the *marginal effects at reasonable/representative/typical values*
  approach (also known as *marginal effects for the average case*).

These approaches are analogous to Stata’s commands `margins x` and
`margins x, at`, respectively.

After calculating predicted probabilities, it will then calculate
differences in probabilities (also known as *contrasts*/*pairwise
comparisons* for categorical predictors or *first differences* for
continuous predictors) using:

- for continuous predictors, the change across the entire range of the
  variable (by default), or a one-unit or one-standard-deviation change
  centred on the mean; or
- for categorical predictors, all pairwise differences.

For both predictions and differences, it uses simulation methods (the
parametric bootstrap) to derive 95% confidence intervals.

It works with the following survey-weighted model objects:
`survey::svyglm` (binary logit), `survey::svyolr` (ordered logit),
`svrepmisc::svymultinom` (multinomial logit)

Eventually, support for (non-survey) weighted model objects (i.e. models
estimated with the `weight` option) will be added for `glm`,
`MASS::polr`, and `nnet::multinom`.

Also included in the package are:

- A snippet of the 2019 Canadian Election Study online survey for
  testing and demonstration purposes. This can be loaded with the
  command `data(ces19)`.
- A `plot()` method that creates a `ggplot` object of predicted
  probabilities or differences in predicted probabilities. This plot can
  be modified by adding further `ggplot` commands, which is shown below.
- A function `mnlSig` that displays a concise summary of multinomial
  logit coefficients with statistical significance stars. This has been
  adapted for use on `svymultinom` objects from Dave Armstrong’s
  original function from `{DAMisc}`, which works for `multinom` objects.

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
data(ces19)

library(survey)
ces19_svy <- survey::svydesign(ids = ~1, strata = NULL, weights = ~pesweight, 
                                data = ces19, digits = 3)

VOTECON <- survey::svyglm(votecon ~ agegrp + gender + educ + region + relig + marketlib + culturetrad, 
                          design = ces19_svy, family = binomial)
summary(VOTECON)
#> 
#> Call:
#> svyglm(formula = votecon ~ agegrp + gender + educ + region + 
#>     relig + marketlib + culturetrad, design = ces19_svy, family = binomial)
#> 
#> Survey design:
#> survey::svydesign(ids = ~1, strata = NULL, weights = ~pesweight, 
#>     data = ces19, digits = 3)
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
`Stata`’s `ologit` command with `[pweight=]` specified.

    . logit votecon i.agegrp i.gender i.educ i.region i.relig marketlib culturetrad [pweight=pesweight]

    Iteration 0:   log pseudolikelihood = -748.14844  
    Iteration 1:   log pseudolikelihood = -527.72183  
    Iteration 2:   log pseudolikelihood = -521.12786  
    Iteration 3:   log pseudolikelihood =  -521.0995  
    Iteration 4:   log pseudolikelihood =  -521.0995  

    Logistic regression                             Number of obs     =      1,284
                                                    Wald chi2(12)     =     233.71
                                                    Prob > chi2       =     0.0000
    Log pseudolikelihood =  -521.0995               Pseudo R2         =     0.3035

    -----------------------------------------------------------------------------------------
                            |               Robust
                    votecon |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
    ------------------------+----------------------------------------------------------------
                     agegrp |
                     35-54  |    .189691   .3066531     0.62   0.536    -.4113379      .79072
                       55+  |   .3757134   .3022457     1.24   0.214    -.2166774    .9681042
                            |
                     gender |
               Woman/Other  |  -.3595359    .187009    -1.92   0.055    -.7260668    .0069949
                            |
                       educ |
                  Some PSE  |   .0288547   .2444425     0.12   0.906    -.4502439    .5079532
                Uni degree  |   .0815175   .2696145     0.30   0.762    -.4469172    .6099521
                            |
                     region |
                      West  |   .6462519   .1960296     3.30   0.001     .2620408    1.030463
                  Atlantic  |    .328965   .3562407     0.92   0.356    -.3692538    1.027184
                            |
                      relig |
                  Catholic  |   .4525211   .2619628     1.73   0.084    -.0609167    .9659588
    Non-Catholic Christian  |     .65283   .2296122     2.84   0.004     .2027984    1.102862
                     Other  |   .7505504   .3960992     1.89   0.058    -.0257897    1.526891
                            |
                  marketlib |   2.285789   .3026759     7.55   0.000     1.692555    2.879023
                culturetrad |   1.961341   .2384488     8.23   0.000      1.49399    2.428692
                      _cons |  -.7241952   .4295126    -1.69   0.092    -1.566024    .1176339
    -----------------------------------------------------------------------------------------

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
                                       weightvar = "pesweight",
                                       seed = 2019)
VOTECON_educ_ame$preds
#> # A tibble: 3 × 5
#>   educ       predicted conf.low conf.high type       
#>   <fct>          <dbl>    <dbl>     <dbl> <chr>      
#> 1 HS or less     0.404    0.340     0.468 Probability
#> 2 Some PSE       0.408    0.369     0.446 Probability
#> 3 Uni degree     0.416    0.374     0.458 Probability
```

The results from the equivalent `Stata` command are below. We can see
the results produced by `sveAME` are within 0.001 for the mean predicted
probability and within 0.002 for the upper and lower 95% confidence
bounds.

    . margins educ, post

    Predictive margins                              Number of obs     =      1,284
    Model VCE    : Robust

    Expression   : Pr(votecon), predict()

    ------------------------------------------------------------------------------
                 |            Delta-method
                 |     Margin   Std. Err.      z    P>|z|     [95% Conf. Interval]
    -------------+----------------------------------------------------------------
            educ |
     HS or less  |   .4028122    .032937    12.23   0.000      .338257    .4673675
       Some PSE  |   .4072815   .0198093    20.56   0.000     .3684559     .446107
     Uni degree  |   .4154583   .0218131    19.05   0.000     .3727054    .4582112
    ------------------------------------------------------------------------------

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

These differences also compare favourably to the same results from
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

For convenience, `{svyEffects}` also includes a `plot()` method, which
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
`"diffs"`) in the `plot()` function call.

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

Now, let’s look at the effect of market liberalism, a continuous
predictor that ranges from -1 (minimal market liberalism, or the most
left-wing position) to +1 (maximal market liberalism, or the most
right-wing position).

``` r
VOTECON_marketlib_ame <- svyAME(VOTECON,
                                varname = "marketlib",
                                weightvar = "pesweight",
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
market liberalism. The probabilities differ a bit more at the ends of
the range (up to 0.006), but at the median of the predictor, the results
are within 0.001.

    . margins, at(marketlib=(-1(.2)1)) post

    Predictive margins                              Number of obs     =      1,284
    Model VCE    : Robust

    Expression   : Pr(votecon), predict()

    1._at        : marketlib       =          -1

    2._at        : marketlib       =         -.8

    3._at        : marketlib       =         -.6

    4._at        : marketlib       =         -.4

    5._at        : marketlib       =         -.2

    6._at        : marketlib       =           0

    7._at        : marketlib       =          .2

    8._at        : marketlib       =          .4

    9._at        : marketlib       =          .6

    10._at       : marketlib       =          .8

    11._at       : marketlib       =           1

    ------------------------------------------------------------------------------
                 |            Delta-method
                 |     Margin   Std. Err.      z    P>|z|     [95% Conf. Interval]
    -------------+----------------------------------------------------------------
             _at |
              1  |   .1181695   .0276963     4.27   0.000     .0638859    .1724532
              2  |   .1662161   .0284543     5.84   0.000     .1104467    .2219855
              3  |   .2264479   .0262596     8.62   0.000     .1749801    .2779157
              4  |   .2980519   .0215145    13.85   0.000     .2558844    .3402195
              5  |   .3786382   .0172419    21.96   0.000     .3448446    .4124317
              6  |    .464545   .0192988    24.07   0.000     .4267201      .50237
              7  |   .5514882   .0272903    20.21   0.000     .4980002    .6049763
              8  |    .635286   .0359907    17.65   0.000     .5647454    .7058265
              9  |    .712409   .0424183    16.79   0.000     .6292706    .7955475
             10  |   .7802779   .0453287    17.21   0.000     .6914353    .8691204
             11  |   .8373832   .0445384    18.80   0.000     .7500895     .924677
    ------------------------------------------------------------------------------

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

Here’s the effect of education on feelings towards the Conservative
Party leader.

For brevity, only the visualizations of the predicted
probabilities/differences are presented (along with comparisons versus
`Stata`).

``` r
CONLDR_region_ame <- svyAME(CONLDR,
                            varname = "region",
                            weightvar = "pesweight",
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

The predicted probabilities for region are very similar to the Stata
results.

![](images/conldr_region.png)

``` r
plot(CONLDR_region_ame, "diffs") +
  geom_hline(yintercept = 0, linetype = "dotted")
```

![](man/figures/unnamed-chunk-16-1.png)<!-- -->

Here’s the effect of market liberalism.

``` r
CONLDR_marketlib_ame <- svyAME(CONLDR,
                            varname = "marketlib",
                            weightvar = "pesweight",
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
`{survey}` package in R. The package `{svyrepmisc}` generates an
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

Note: use the option `trace = FALSE` in the `svymultinom()` function
call to suppress the reporting of each replication (similar to using the
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

For our post-estimation command, we’ll need to specify a few more
options because `svymultinom` does not store them in its output. These
are:

- `design`: the survey design object used to estimate the model; and
- `modform`: the model formula used in the `svymultinom` call (in the
  form `modform = "y ~ x1 + x2 + x3"`).

Here’s the effect of education.

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

Here’s the effect of market liberalism.

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

Finally, are the differences:

``` r
plot(VOTE_marketlib_ame, "diffs")
```

![](man/figures/unnamed-chunk-29-1.png)<!-- -->

------------------------------------------------------------------------

# Marginal effects at reasonable values

*(documentation in progress)*

You can choose to calculate marginal effects at reasonable values (MER)
probabilities and differences by using the `svyMER` function, which
follows the same arguments as `svyAME`. This command would give you the
estimated effect of a variable “for the ‘typical’ case” (which may or
may not be typical, plausible, or even possible in the real world) as
opposed to the effect of a variable across the population (see Hanmer
and Kalkan 2013 for an in-depth discussion).

MER probabilities are *usually* very similar to AME probabilities, but
not always. However, they are *much* faster to calculate using
simulation methods.

------------------------------------------------------------------------

# Interaction effects

*(documentation in progress)*

Both `svyAME` and `svyMER` support calculating predicted probabilities
of combinations of two predictor variables. This can be done by using
the argument `byvar = "x"` in the function call.

This will not return differences in predicted probabilities. For
limited-dependent variable models, one would need to calculate a second
difference would be needed to test for the significance of an
interaction between two variables, either from the inclusion of a
product term or through the compression inherent in these types of
models (see Norton, Wang, and Ai 2004).

Dave Armstrong’s `{DAMisc}` package has an `R` port for Norton, Wang,
and Ai’s original `Stata` function, and this will eventually be ported
to `{svyEffects}` for use in survey-weighted models.

------------------------------------------------------------------------

# Planned updates

This package is under active development, and updates will include:

1.  More user-friendly error-checking and reporting (e.g. checking at
    the beginning of the function that the variables in the function
    don’t have typos).
2.  Support for (non-survey) weighted models and weighted models. While
    there are other packages that do this, some do not return confidence
    intervals for predictions for some model types. And, to my
    knowledge, none use simulation methods to derive confidence
    intervals. *Note: You can actually already do this with
    `{svyEffects}` by creating a survey design object with a weight of
    “1”, but it would be good to avoid having to use that workaround.*
3.  In-depth comparisons with Stata results.
4.  Functions for calculating weighted model fit measures.
5.  Support for using an alternative variance-covariance matrix using
    `{sandwich}`. This would only be for binary logit models because
    `{sandwich}` does not play nice with ordinal or multinomial models.
    That said, survey-weighted models do adjust the variance-covariance
    matrix (the documentation for `{survey}` does not specify the
    correction method it uses, but it appears to be HC0, based on what
    I’ve seen).
6.  A second differences function to test for the significance of a
    two-way interaction.
7.  (Eventually) Use of the delta method to calculate confidence
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
