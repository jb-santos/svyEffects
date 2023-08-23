#' Marginal Effects At Reasonable Values For Multinomial Logit Models Of Survey-Weighted Data
#'
#'
#' @description Calculates predicted probabilities and differences in probabilities
#' for a predictor holding all other observations at reasonable/representative/typical
#' values (i.e. for an "average case"). This involves setting all continuous
#' variables to their medians and all categorical variables to their modes. Uses
#' simulations (the parametric bootstrap) to derive 95% confidence intervals.
#'
#' READ THIS BEFORE USING: You need the GitHub package \code{svrepmisc} to
#' use this function! Currently, there are no multinomial logit \code{R} functions
#' that work on survey-weighted data. However, Carl Ganz' \code{svrepmisc} package
#' (available here: \url{https://github.com/carlganz/svrepmisc}) uses replicate
#' weights and jackknife standard errors to approximate the same end result.
#' (The output with this approach is very similar to results obtained in Stata
#' using either sampling weights (i.e. "pweights") or replicate weights with
#' jackknife standard errors.)
#'
#' To use this function...
#' 1. Use \code{survey::svydesign} to create a survey design object with your sampling/design weights.
#' 2. Use \code{survey::as.svrepdesign(x, type="JK1")} convert the sampling weights to replicate weights with jackknife standard errors.
#' 3. Estimate your multinomial logit model on the replicate weights design object using \code{svrepmisc::svymultinom}.
#' 4. Now, you're ready to use \code{svyAME} to calculate probabilities.
#'
#'
#' @param obj Model object of class \code{svrepmisc::svrepstatmisc}.
#' @param varname Character string denoting the name of the predictor variable
#' for which effects are to be calculated.
#' @param nvals Scalar denoting the sequence length spanning the range of a
#' continuous variable for which effects are to be calculated (default: 11).
#' @param diffchange Character string  denoting over what change in x a first
#' difference is to be calculated for a continuous predictor (default: "sd").
#' @param byvar For interaction effects, a character string denoting the name of
#' the moderator variable.
#' @param bychange For interaction effects with a numerical moderator, a
#' character string denoting for what values above and below the mean
#' predictions should be calculated (default: "sd").
#' @param bynvals For interaction effects with a numerical moderator, a scalar
#' denoting the sequence length  spanning the range of the moderator variable
#' for which effects are to be calculated (default: 3).
#' @param sims Scalar denoting how many simulations to conduct (default: 2500).
#' @param seed Seed value for random number generator. By default, the function
#' picks a random integer between 1 and 1,000,000. If you save the output of
#' this function, it will save the seed value used for simulations in the slot
#' \code{$seed}.
#' @param ci Scalar indicating confidence level to be used (default: .95).
#' @param ... Other arguments (currently not implemented).
#'
#'
#' @return A list with dataframes:
#'   \describe{
#'     \item{\code{$preds}}{predicted probabilities}
#'     \item{\code{$diffs}}{differences in predicted probabilities}
#'     \item{\code{$seed}}{seed value used for random number generator}
#'     \item{\code{$typical}}{the values at which variables are set for the simulations}
#'     }
#'
#'
#' @author John Santos & Dave Armstrong
#'
#'
#' @export
#'
#'
#' @examples
#' data(ces19)
#' library(survey)
#' ces19_svy <- svydesign(ids = ~1, strata = NULL, weights = ~pesweight,
#'   data = ces19, digits = 3)
#' ces19_svy_r <- as.svrepdesign(ces19_svy, type = "JK1")
#' # remotes::install_github("carlganz/svrepmisc") # (if not already installed)
#' library(svrepmisc)
#' VOTE <- svymultinom(vote ~ agegrp + gender + educ + region + marketlib,
#'   design = ces19_svy_r, trace = FALSE)
#' svyMER(VOTE, varname = "region", weightvar = "pesweight",
#'   seed = 2019, design = ces19_svy_r,
#'   modform = "vote ~ agegrp + gender + educ + region + marketlib")
#' svyMER(VOTE, varname = "marketlib", weightvar = "pesweight",
#'   seed = 2019, design = ces19_svy_r,
#'   modform = "vote ~ agegrp + gender + educ + region + marketlib")
#'
svyMER.svrepstatmisc <- function(obj,
                                 varname,
                                 nvals = 11,
                                 diffchange = c("sd", "range", "unit"),
                                 byvar = NULL,
                                 bychange = c("sd", "range"),
                                 bynvals = 3,
                                 sims = 2500,
                                 seed = NULL,
                                 ci = .95,
                                 ...) {

  #========== SETUP ============================================================

  modform <- attr(obj, "formula")
  modvars <- unformulate(attr(obj, "formula"))$vars
  data <- attr(obj, "svrep.design")$variables
  data$`(weights)` <- as.vector(attr(obj, "svrep.design")$pweights)
  data <- data %>%
    select(all_of(modvars), `(weights)`) %>%
    na.omit()
  svydata <- survey::svydesign(ids = ~1,
                               strata = NULL,
                               weights = data$`(weights)`,
                               data = data)

  Yname <- as.character(modform[[2]])
  Yvar <- dplyr::select(data, all_of(Yname))
  Yvar <- Yvar[[1]]
  Ylevs <- levels(Yvar)
  Ynlev <- as.numeric(nlevels(Yvar))

  # if(class(data[[varname]]) == "factor") {
  # if(isTRUE(inherits(class(data[[varname]]), "factor"))) {
  if(is.factor(data[[varname]])) {
    Xlevs <- levels(data[[varname]])}
  # if(class(data[[varname]]) == "numeric") {
  # if(isTRUE(inherits(class(data[[varname]]), "numeric"))) {
  if(is.numeric(data[[varname]])) {
    Xlevs <- seq(min(data[[varname]]), max(data[[varname]]), length=nvals)}

  # Check arguments up-front to stop execution before running simulations
  if(isFALSE(varname %in% names(data))) {
    stop(print(paste0(
      "varname ", varname, " not found in survey design object. Maybe check your spelling?")))}
  if(!is.null(nvals) & isFALSE(is.numeric(nvals))) {
    stop(print(paste0(
      "Non-numeric value entered for nvals. Please enter a numeric value.")))}
  if(!is.null(sims) & isFALSE(is.numeric(sims))) {
    stop(print(paste0(
      "Non-numeric value entered for sims. Please enter a numeric value.")))}
  if(!is.null(seed) & isFALSE(is.numeric(seed))) {
    stop(print(paste0(
      "Non-numeric value entered for seed. Please enter a numeric value.")))}

  # Set random number generator
  if(is.null(seed)) {
    seed <- round(runif(1, min = 0, max = 1000000))
  } else {
    seed <- seed
  }
  set.seed(seed)

  # high / low helper functions
  tail <- (1 - ci) / 2
  low <- function(x) {unname(quantile(x, tail))}
  high <- function(x) {unname(quantile(x, 1-tail))}



  #========== NO MODERATOR VARIABLE ============================================

  if(is.null(byvar)) {

    # Define variables to vary
    varylist <- as.list(varname)
    varylist <- lapply(1:length(varylist), function(i) {
      # if(class(data[[varylist[[i]]]])=="numeric") {
      # if(isTRUE(inherits(class(data[[varylist[[i]]]]), "numeric"))) {
      if(is.numeric(data[[varylist[[i]]]])) {
        varylist[[i]] <- seq(min(data[[varylist[[i]]]]), max(data[[varylist[[i]]]]), length = nvals)
      } else {
        varylist[[i]] <- levels(data[[varylist[[i]]]])
      }})
    names(varylist) <- varname

    # Create fake dataset for simulations
    trms <- attr(terms(model.frame(as.formula(modform), data)), "dataClasses")
    l <- vector(mode = "list", length=length(trms))
    names(l) <- names(trms)
    if(!is.null(varylist)) {   # varying terms
      for(i in 1:length(varylist)) {
        l[[names(varylist)[i]]] <- varylist[[i]]
      }}
    for(m in 1:length(l)){   # terms fixed at central observations
      if(length(l[[m]]) == 0){
        if(trms[[m]] == "factor") {
          l[[m]] <- svymode(names(trms)[m], svydata)
        } else {
          l[[m]] <- svyquantile(reformulate(names(trms)[m]), svydata, 0.5)[[1]][1]
        }}}
    fake <- do.call(expand.grid, l)   # matrix for generating predictions

    #========== PREDICTED PROBABILITIES ========================================

    # Parametric bootstrap
    B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))
    fakeX <- model.matrix(formula(as.formula(modform)), data = fake)

    # Generate predictions
    res <- array(dim=c(length(Xlevs),Ynlev,sims))
    for(i in 1:nrow(B)){
      tmpB <- cbind(0, matrix(B[i,], ncol = (Ynlev-1), byrow = TRUE))
      exb <- exp(fakeX %*% tmpB)
      p <- prop.table(exb, 1)
      res[,,i] <- p
    }

    # Assemble table of predictions
    preds <- data.frame(x = Xlevs)
    if(is.factor(data[[varname]])) {preds$x <- factor(preds$x, levels = Xlevs)}
    res_m <- as.data.frame(apply(res, c(1,2), mean))
    names(res_m) <- paste0("predicted_", Ylevs)
    res_l <- as.data.frame(apply(res, c(1,2), low))
    names(res_l) <- paste0("conf.low_", Ylevs)
    res_u <- as.data.frame(apply(res, c(1,2), high))
    names(res_u) <- paste0("conf.high_", Ylevs)

    preds <- cbind(preds, res_m, res_l, res_u)
    preds <- preds %>%
      tidyr::pivot_longer(-x,
                          names_pattern = "(.*)_(.*)",
                          names_to = c(".value", "y")) %>%
      dplyr::mutate(type = "Probability") %>%
      dplyr::relocate(y, .before=1)
    preds$y <- factor(preds$y, levels = Ylevs)


    #========== DIFFERENCES IN PREDICTED PROBABILITIES =========================

    ## Differences for categorical variables -----------------------------------

    if(is.factor(data[[varname]])) {
      P <- lapply(1:length(Xlevs), function(i) (t(res[i,,])))
      P_combs <- combn(P, 2, simplify = FALSE)
      D <- lapply(1:length(P_combs), function(i) P_combs[[i]][[2]] - P_combs[[i]][[1]])
      diff_res <- cbind.data.frame(D)

      diffs <- expand.grid(
        y = factor(Ylevs, levels = Ylevs),
        x1 = 1:length(Xlevs),
        x2 = 1:length(Xlevs)) %>%
        dplyr::filter(x1 > x2) %>%
        dplyr::mutate(x = paste0(Xlevs[x1],
                                 " - ",
                                 Xlevs[x2])) %>%
        dplyr::select(-c("x1", "x2"))
      diffs <- diffs %>%
        dplyr::mutate(
          predicted = colMeans(diff_res),
          conf.low = apply(diff_res, 2, low),
          conf.high = apply(diff_res, 2, high),
          type = "Difference")
      diffs

    } else {

      # Differences for continuous variables -----------------------------------

      diffchange <- match.arg(diffchange)
      tmp0 <- switch(diffchange,
                     NULL = min(data[[varname]]),
                     range = min(data[[varname]]),
                     unit = survey::svymean(data[[varname]], svydata) - .5,
                     sd = survey::svymean(data[[varname]], svydata) -
                       (sqrt(survey::svyvar(data[[varname]], svydata)) / 2)
                     )
      tmp1 <- switch(diffchange,
                     NULL = max(data[[varname]]),
                     range = max(data[[varname]]),
                     unit = survey::svymean(data[[varname]], svydata) + .5,
                     sd = survey::svymean(data[[varname]], svydata) +
                       (sqrt(survey::svyvar(data[[varname]], svydata)) / 2)
                     )
      delta <- c(tmp0, tmp1)

      # Define variables to vary
      varylist <- as.list(varname)
      varylist <- lapply(1:length(varylist), function(i) varylist[[i]] <- delta)
      names(varylist) <- varname

      # Create fake dataset for simulations
      trms <- attr(terms(model.frame(as.formula(modform), data)), "dataClasses")
      l <- vector(mode = "list", length = length(trms))
      names(l) <- names(trms)
      if(!is.null(varylist)) {   # varying terms
        for(i in 1:length(varylist)) {
          l[[names(varylist)[i]]] <- varylist[[i]]
        }}
      for(m in 1:length(l)){   # terms fixed at central observations
        if(length(l[[m]]) == 0){
          if(trms[[m]] == "factor") {
            l[[m]] <- svymode(names(trms)[m], svydata)
          } else {
            l[[m]] <- svyquantile(reformulate(names(trms)[m]), svydata, 0.5)[[1]][1]
          }}}
      fake <- do.call(expand.grid, l)   # matrix for generating predictions

      # Parametric bootstrap
      B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))
      fakeX <- model.matrix(formula(as.formula(modform)), data = fake)

      # Generate predictions
      res <- array(dim=c(2,Ynlev,sims))
      for(i in 1:nrow(B)){
        tmpB <- cbind(0, matrix(B[i,], ncol = (Ynlev-1), byrow = TRUE))
        exb <- exp(fakeX %*% tmpB)
        p <- prop.table(exb, 1)
        res[,,i] <- p
      }

      # Calculate differences
      diff_res <- res[2,,] - res[1,,]

      # Assemble table of differences
      diffs <- data.frame(
        y = factor(Ylevs, levels = Ylevs),
        x = paste0("Delta (", diffchange, ") : ",
                   round(tmp0, 3), " - ", round(tmp1, 3)),
        predicted = apply(diff_res, c(1), mean),
        conf.low = apply(diff_res, c(1), low),
        conf.high = apply(diff_res, c(1), high),
        type = "Difference")
    }

    #========== Output =========================================================

    preds <- rename(preds, !!sym(varname) := x)
    diffs <- rename(diffs, !!sym(varname) := x)
    output <- list(
      preds = tibble::as_tibble(preds),
      diffs = tibble::as_tibble(diffs),
      typical = tibble::as_tibble(fake),
      seed = seed,
      sims = sims,
      formula = as.formula(modform))
    class(output) <- "svyEffects"
    attributes(output)$predvar <- varname
    attributes(output)$depvar <- Yname
    attributes(output)$method <- "MER"
    return(output)

  }



  #==========  WITH MODERATOR VARIABLE =========================================
  # Note: Interactive models only output predictions, not differences.
  # (b/c interactive models should use second differences--this is on to-do list)

  if(!is.null(byvar)) {


    # Check arguments up front to stop execution before running simulations
    if(isFALSE(byvar %in% names(data))) {
      stop(print(paste0(
        "byvar ", byvar, " not found in survey design object. Maybe check your spelling?")))}
    if(!is.null(bynvals) & isFALSE(is.numeric(bynvals))) {
      stop(print(paste0(
        "Non-numeric value entered for bynvals. Please enter a numeric value.")))}


    Yname <- as.character(modform[[2]])
    Yvar <- dplyr::select(data, all_of(Yname))
    Yvar <- Yvar[[1]]
    Ylevs <- levels(Yvar)
    Ynlev <- length(Ylevs)

    if(is.factor(data[[varname]])) {
      Xlevs <- levels(data[[varname]])
      }
    if(is.numeric(data[[varname]])) {
      Xlevs <- seq(min(data[[varname]]), max(data[[varname]]), length = nvals)
      }
    bylevs <- 1   # so length(nlev) * length(bylev) still works in array dimension

    if(is.factor(data[[byvar]])) {
      bylevs <- levels(data[[byvar]])
      }
    if(is.numeric(data[[byvar]])) {
      bychange <- match.arg(bychange)
      byvar_min <- switch(bychange,
                          NULL = as.numeric(survey::svymean(data[byvar], svydata)) -
                            as.numeric(sqrt(survey::svyvar(data[byvar], svydata))),
                          sd = as.numeric(survey::svymean(data[byvar], svydata)) -
                            as.numeric(sqrt(survey::svyvar(data[byvar], svydata))),
                          range = min(data[[byvar]]))
      byvar_max <- switch(bychange,
                          NULL = as.numeric(survey::svymean(data[byvar], svydata)) +
                            as.numeric(sqrt(survey::svyvar(data[byvar], svydata))),
                          sd = as.numeric(survey::svymean(data[byvar], svydata)) +
                            as.numeric(sqrt(survey::svyvar(data[byvar], svydata))),
                          range = max(data[[byvar]]))
      bylevs <- seq(byvar_min, byvar_max, length = bynvals)
      }

    # Define variables to vary
    xvarylist <- as.list(varname)
    xvarylist <- lapply(1:length(varname), function(i) {
      if(is.factor(data[[xvarylist[[i]]]])) {
        xvarylist[[i]] <- levels(data[[varname[[i]]]])
      } else {
        xvarylist[[i]] <- seq(min(data[[varname[[i]]]]), max(data[[varname[[i]]]]), length = nvals)
      }})
    names(xvarylist) <- varname
    varylist <- xvarylist

    byvarylist <- as.list(byvar)
    byvarylist <- lapply(1:length(byvar), function(i) {
      if(is.factor(data[[byvarylist[[i]]]])) {
        byvarylist[[i]] <- levels(data[[byvar[[i]]]])
      } else {
        byvarylist[[i]] <- seq(min(data[[byvar[[i]]]]), max(data[[byvar[[i]]]]), length = bynvals)
      }})
    names(byvarylist) <- byvar
    varylist <- c(xvarylist, byvarylist)

    # Create fake dataset for simulations
    trms <- attr(terms(model.frame(as.formula(modform), data)), "dataClasses")
    l <- vector(mode = "list", length = length(trms))
    names(l) <- names(trms)
    if(!is.null(varylist)) {   # varying terms
      for(i in 1:length(varylist)) {
        l[[names(varylist)[i]]] <- varylist[[i]]
      }}
    for(m in 1:length(l)){   # terms fixed at central observations
      if(length(l[[m]]) == 0){
        if(trms[[m]] == "factor") {
          l[[m]] <- svymode(names(trms)[m], svydata)
        } else {
          l[[m]] <- svyquantile(reformulate(names(trms)[m]), svydata, 0.5)[[1]][1]
        }}}
    fake <- do.call(expand.grid, l)   # matrix for generating predictions

    # Simulate using parametric bootstrap
    B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))
    fakeX <- model.matrix(formula(as.formula(modform)), data = fake)

    # Generate predictions
    res <- array(dim=c((length(Xlevs) * length(bylevs)), length(Ylevs), sims))
    for(i in 1:nrow(B)){
      tmpB <- cbind(0, matrix(B[i,], ncol = (length(Ylevs)-1), byrow = TRUE))
      exb <- exp(fakeX %*% tmpB)
      p <- prop.table(exb, 1)
      res[,,i] <- p
    }

    # Calculate mean, low CI, upper CI for predictions
    res_m <- as.data.frame(apply(res, c(1,2), mean))
    names(res_m) <- paste0("predicted_", Ylevs)
    res_l <- as.data.frame(apply(res, c(1,2), low))
    names(res_l) <- paste0("conf.low_", Ylevs)
    res_u <- as.data.frame(apply(res, c(1,2), high))
    names(res_u) <- paste0("conf.high_", Ylevs)

    # Figure label predictions with correct combination of xvar * byvar (b/c order of
    # preds determined by order of terms in model.formula, not specified varname/byvar)
    if(grep(varname, names(fake)) < grep(byvar, names(fake))){
      v1 <- data[varname]} else {v1 <- data[byvar]}
    if(grep(varname, names(fake)) > grep(byvar, names(fake))){
      v2 <- data[varname] } else {v2 <- data[byvar] }

    # So numerical moderators show up as factors (for plotting purposes)
    if (!is.null(byvar) & is.numeric(data[[byvar]])) {
      varylist[[match(byvar, names(varylist))]] <- factor(round(bylevs, 3))
    }

    # Assemble predictions table
    keygrid <- expand.grid(
      x1 = varylist[[match(names(v1), names(varylist))]],
      x2 = varylist[[match(names(v2), names(varylist))]])
    preds <- cbind(keygrid, res_m, res_l, res_u)
    preds <- cbind(keygrid, res_m, res_l, res_u)
    preds <- preds %>%
      pivot_longer(-c(x1,x2),
                   names_pattern="(.*)_(.*)",
                   names_to=c(".value", "y")) %>%
      mutate(type = "Probability") %>%
      relocate(y, .before=1)
    preds <- preds %>%
      rename(!!sym(names(v1)) := x1,
             !!sym(names(v2)) := x2)
    preds$y <- factor(preds$y, levels = Ylevs)

    # Assemble combined output table
    output <- list(
      preds = tibble::as_tibble(preds),
      typical = tibble::as_tibble(fake),
      seed = seed,
      sims = sims,
      formula = as.formula(modform))
    class(output) <- "svyEffects"
    attributes(output)$predvar <- varname
    attributes(output)$byvar <- byvar
    attributes(output)$depvar <- Yname
    attributes(output)$method <- "MER"
    return(output)

  }

}
