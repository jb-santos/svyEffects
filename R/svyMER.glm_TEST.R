#' Marginal Effects At Reasonable Values For Binary Logit Models Of Survey-Weighted Data
#'
#'
#' @description Calculates predicted probabilities and differences in probabilities
#' for a predictor holding all other observations at reasonable/representative/typical
#' values (i.e. for an "average case"). This involves setting all continuous
#' variables to their medians and all categorical variables to their modes. Uses
#' simulations (the parametric bootstrap) to derive 95% confidence intervals.
#'
#'
#' @param obj Model object of class \code{survey::svyglm} or \code{glm}.
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
#' @param ci Scalar indicating confidence level to be used (default: 0.95).
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
#' ces19_svy <- svydesign(ids = ~1, strata = NULL,
#'   weights = ~pesweight, data = ces19, digits = 3)
#' VOTECON <- svyglm(votecon ~ agegrp + gender + educ + region + marketlib,
#'   design = ces19_svy, family = binomial)
#' svyMER(VOTECON, varname = "educ", seed = 2019)
#' svyMER(VOTECON, varname = "marketlib", seed = 2019)
#'
testfun <- function(obj,
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

  # Check if model is binary logit
  if(!(family(obj)$link %in% c("logit"))) {
    stop("Model should be of binomial family with logit link.\n")
  }

  # Get data
  data <- model.frame(obj)
  data$`(weights)` <- weights(obj$survey.design)

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

  svydata <- survey::svydesign(ids = ~1, strata = NULL, weights = data$`(weights)`, data = data)
  if(is.character(data[[varname]])) {data[[varname]] <- as.factor(data[[varname]])}

  # Set random number generator
  if(is.null(seed)){
    seed <- round(runif(1, min=0, max=1000000))
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
      # if(class(data[[varylist[[i]]]]) == "numeric") {
      # if(isTRUE(inherits(class(data[[varylist[[i]]]]), "numeric"))) {
      if(is.numeric(data[[varylist[[i]]]])) {
        varylist[[i]] <- seq(min(data[[varylist[[i]]]]), max(data[[varylist[[i]]]]), length = nvals)
      } else {
        varylist[[i]] <- levels(data[[varylist[[i]]]])
      }})
    names(varylist) <- varname

    # Create fake dataset for simulations
    trms <- attr(terms(obj), "dataClasses")
    trms <- trms[-length(trms)]
    l <- vector(mode="list", length=length(trms))
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
          l[[m]] <- survey::svyquantile(reformulate(names(trms)[m]), svydata, 0.5)[[1]][1]
        }}}

    fake <- do.call(expand.grid, l)

    #========== PREDICTED PROBABILITIES ========================================

    B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))
    fakeX <- model.matrix(formula(obj), data = fake)
    probs <- t(plogis(fakeX %*% t(B)))

    preds <- probs %>%
      as.data.frame() %>%
      tidyr::pivot_longer(everything(), names_to = "ind", values_to = "probs") %>%
      dplyr::group_by(ind) %>%
      dplyr::summarize(predicted = mean(probs),
                conf.low = low(unname(probs)),
                conf.high = high(unname(probs)))
    preds$ind <- as.numeric(preds$ind)
    preds <- dplyr::arrange(preds, ind)
    preds <- dplyr::select(preds, -ind)
    # if(class(data[[varname]]) == "factor") {
    if(is.factor(data[[varname]])) {
      preds$x <- factor(levels(data[[varname]]), levels = levels(data[[varname]]))}
    else{preds$x <- seq(min(data[[varname]]), max(data[[varname]]), length = nvals)}
    preds <- preds %>%
      dplyr::relocate(x, .before=1) %>%
      dplyr::mutate(type = "Probability")

    #========== DIFFERENCES IN PREDICTED PROBABILITIES =========================

    ## Differences for categorical variables -----------------------------------

    # if(class(data[[varname]]) == "factor") {
    if(is.factor(data[[varname]])) {
      levs <- levels(data[[varname]])
      nlev <- length(levs)

      P <- NULL
      P <- lapply(1:ncol(probs), function(i) P[[i]] <- probs[,i])
      P_combs <- combn(P, 2)
      D <- mapply(x = P_combs[2,], y = P_combs[1,],
                  function(x, y) unlist(x) - unlist(y))
      D <- as.data.frame(D)

      diffs <- tibble::tibble(
        varname1 = combn(levs, 2)[1,],
        varname2 = combn(levs, 2)[2,],
        predicted = sapply(D, mean),
        conf.low = sapply(D, low),
        conf.high = sapply(D, high),
        type="Difference")
      diffs <- diffs %>%
        dplyr::mutate(x = paste0(diffs$varname2,
                                 " - ",
                                 diffs$varname1),
                      .before=1) %>%
        dplyr::select(-c(varname1, varname2))

    } else {

      ## Differences for continuous variables ----------------------------------

      diffchange <- match.arg(diffchange)
      diffrange <- switch(diffchange,
                          NULL = survey::svymean(data[[varname]], svydata) +
                            (c(-.5,.5) * (sqrt(survey::svyvar(data[varname], svydata)))),
                          sd = survey::svymean(data[[varname]], svydata) +
                            (c(-.5,.5) * (sqrt(survey::svyvar(data[varname], svydata)))),
                          range = c(min(data[[varname]]), max(data[[varname]])),
                          unit = survey::svymean(data[[varname]], svydata) + c(-.5,.5)
                          )

      varylist <- as.list(varname)
      varylist <- lapply(1:length(varylist), function(i) varylist[[i]] <- diffrange)
      names(varylist) <- varname

      # Create fake dataset for simulations
      trms <- attr(terms(obj), "dataClasses")
      trms <- trms[-length(trms)]
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
            l[[m]] <- survey::svyquantile(reformulate(names(trms)[m]), svydata, 0.5)[[1]][1]
          }}}
      fake <- do.call(expand.grid, l)

      B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))
      fakeX <- model.matrix(formula(obj), data=fake)
      probs <- t(plogis(fakeX %*% t(B)))
      diff_res <- apply(probs, 1, diff)
      diffs <- tibble::tibble(
        x = paste0("Delta (", diffchange, ") : ",
                   round(diffrange[[1]], 3),
                   " - ",
                   round(diffrange[[2]], 3)),
        predicted = mean(diff_res),
        conf.low = unname(low(diff_res)),
        conf.high = unname(high(diff_res)),
        type = "Difference")
    }

    #========== Output =========================================================

    preds <- dplyr::rename(preds, !!sym(varname) := x)
    diffs <- dplyr::rename(diffs, !!sym(varname) := x)
    output <- list(
      preds = tibble::as_tibble(preds),
      diffs = tibble::as_tibble(diffs),
      typical = tibble::as_tibble(fake),
      seed = seed,
      sims = sims,
      formula = formula(obj))
    class(output) <- "svyEffects"
    attributes(output)$predvar <- varname
    attributes(output)$depvar <- colnames(obj$model)[1]
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
    if(!is.null(nvals) & isFALSE(is.numeric(nvals))) {
      stop(print(paste0(
        "Non-numeric value entered for nvals. Please enter a numeric value.")))}

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
        #bylevs <- seq(byvar_min, byvar_max, length = bynvals)
        #byvarylist[[i]] <- bylevs

        byvarylist[[i]] <- seq(byvar_min, byvar_max, length = bynvals)

        # byvarylist[[i]] <- seq(min(data[[byvar[[i]]]]), max(data[[byvar[[i]]]]), length = bynvals)
      }})
    names(byvarylist) <- byvar
    bylevs <- byvarylist[[1]]
    varylist <- c(xvarylist, byvarylist)

    # Create fake dataset for simulations
    trms <- attr(terms(obj), "dataClasses")
    trms <- trms[-length(trms)]
    l <- vector(mode="list", length=length(trms))
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
          l[[m]] <- survey::svyquantile(reformulate(names(trms)[m]), svydata, 0.5)[[1]][1]
        }}}
    fake <- do.call(expand.grid, l)

    # Generate predictions
    B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))
    fakeX <- model.matrix(formula(obj), data=fake)
    probs <- t(plogis(fakeX %*% t(B)))

    # Assemble predictions table
    preds <- probs %>%
      as.data.frame() %>%
      tidyr::pivot_longer(everything(), names_to = "ind", values_to = "probs") %>%
      dplyr::group_by(ind) %>%
      dplyr::summarize(predicted = mean(probs),
                       conf.low = low(probs),
                       conf.high = high(probs))
    preds$ind <- as.numeric(preds$ind)
    preds <- dplyr::arrange(preds, ind)
    preds <- dplyr::select(preds, -ind)
    preds$type <- "Probability"

    # Label predictions with correct combination of xvar * byvar (b/c order of preds
    # determined by order of terms in model.formula, not specified varname/byvar)
    if(grep(varname, names(model.frame(obj))) < grep(byvar, names(model.frame(obj)))){
      v1 <- data[varname]} else {v1 <- data[byvar]}
    if(grep(varname, names(model.frame(obj))) > grep(byvar, names(model.frame(obj)))){
      v2 <- data[varname] } else {v2 <- data[byvar] }

    # So numerical moderators show up as factors (for plotting purposes)
    if (!is.null(byvar) & is.numeric(data[[byvar]]))
      varylist[[match(byvar, names(varylist))]] <- factor(round(bylevs, 3))

    keygrid <- expand.grid(
      x1 = varylist[[match(names(v1), names(varylist))]],
      x2 = varylist[[match(names(v2), names(varylist))]])
    preds <- cbind(keygrid, preds)
    preds <- preds %>%
      dplyr::rename(!!sym(names(v1)) := x1,
                    !!sym(names(v2)) := x2)

    varname_pos <- grep(varname, names(model.frame(obj)))
    byvar_pos <- grep(byvar, names(model.frame(obj)))

    output <- list(
      preds = tibble::as_tibble(preds),
      typical = tibble::as_tibble(fake),
      seed = seed,
      sims = sims,
      formula = formula(obj))
    class(output) <- "svyEffects"
    attributes(output)$predvar <- varname
    attributes(output)$byvar <- byvar
    attributes(output)$depvar <- colnames(obj$model)[1]
    attributes(output)$method <- "MER"
    return(output)

  }

}
