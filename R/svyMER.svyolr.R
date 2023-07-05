#' Marginal Effects At Reasonable Values For Ordered Logit Models Of Survey-Weighted Data
#'
#'
#' @description Calculates predicted probabilities and differences in probabilities
#' for a predictor holding all other observations at reasonable/representative/typical
#' values (i.e. for an "average case"). This involves setting all continuous
#' variables to their medians and all categorical variables to their modes. Uses
#' simulations (the parametric bootstrap) to derive 95% confidence intervals.
#'
#'
#' @param obj Model object of class \code{survey::svyolr}.
#' @param varname Character string denoting the name of the predictor variable
#' for which effects are to be calculated.
#' @param nvals Scalar denoting the sequence length spanning the range of a
#' continuous variable for which effects are to be calculated (default: 10).
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
#' CONLDR <- svyolr(ftconldr ~ agegrp + gender + educ + region + marketlib,
#'   design = ces19_svy)
#' svyMER(CONLDR, varname = "region", seed = 2019)
#' svyMER(CONLDR, varname = "marketlib", seed = 2019)
#'
svyMER.svyolr <- function(obj,
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

  # Get data
  data <- model.frame(obj)

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

  svydata <- survey::svydesign(ids = ~1, strata = NULL, weights = data$`(weights)`, data=data)
  if(is.character(data[[varname]])) {data[[varname]] <- as.factor(data[[varname]])}

  # Set random number generator
  if(is.null(seed)){
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
      if(is.numeric(data[[varname]])) {
        varylist[[i]] <- seq(min(data[[varylist[[i]]]]), max(data[[varylist[[i]]]]), length=nvals)
      } else {
        varylist[[i]] <- levels(data[[varylist[[i]]]])
      }})
    names(varylist) <- varname

    # Create fake dataset for simulations
    trms <- attr(terms(obj), "dataClasses")
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
    fake <- do.call(expand.grid, l)   # matrix for generating predictions

    #========== PREDICTED PROBABILITIES ========================================

    # Parametric bootstrap
    B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))
    tau_ind <- grep("\\|", names(coef(obj)))
    Tau <- B[,tau_ind]
    B <- B[,-tau_ind]

    # Generate predictions
    fakeX <- model.matrix(formula(obj), data=fake)
    CP <- NULL
    for (i in 1:ncol(Tau)) {
      CP[[i]] <- plogis(matrix(Tau[,i],
                               ncol=nrow(Tau),
                               nrow=nrow(fakeX[,-1]),
                               byrow=TRUE) - fakeX[,-1] %*% t(B))
    }
    P <- NULL
    for (i in 1:ncol(Tau)) {
      if(i==1) {P[[i]] <- CP[[i]]}
      else{P[[i]] <- CP[[i]]- CP[[(i-1)]]
      }
    }
    P[[(ncol(Tau)+1)]] <- 1-CP[[(ncol(Tau))]]

    # Assemble table of predictions
    # if(class(data[[varname]])=="factor") {
    # if(isTRUE(inherits(class(data[[varname]]), "factor"))) {
    if(is.factor(data[[varname]])) {
      levs <- levels(data[[varname]])}
    # if(class(data[[varname]])=="numeric") {
    # if(isTRUE(inherits(class(data[[varname]]), "numeric"))) {
    if(is.numeric(data[[varname]])) {
      levs <- seq(min(data[[varname]]), max(data[[varname]]), length=nvals)}
    preds <- expand.grid(
      x = levs,
      y = obj$lev)  %>%
      mutate(predicted = unlist(lapply(1:length(obj$lev), function(i)apply(P[[i]], 1, mean))),
             conf.low = unlist(lapply(1:length(obj$lev), function(i)apply(P[[i]], 1, low))),
             conf.high = unlist(lapply(1:length(obj$lev), function(i)apply(P[[i]], 1, high))),
             type = "Probability")

    #========== DIFFERENCES IN PREDICTED PROBABILITIES =========================

    ## Differences for categorical variables -----------------------------------

    # if(class(data[[varname]])=="factor") {
    # if(isTRUE(inherits(class(data[[varname]]), "factor"))) {
    if(is.factor(data[[varname]])) {
      levs <- levels(data[[varname]])
      nlev <- length(levs)

      # Calculate differences
      # *Note: This is a bit different from the other difference functions because
      # the prediction here are organized by level of response, not level of predictor.
      D <- lapply(1:length(P), function(i)split(P[[i]], seq(1:nlev)))
      combs <- lapply(1:length(D), function(i)combn(D[[i]], 2))
      diff_list <- lapply(1:length(combs), function(i){
        mapply(x = combs[[i]][2,], y = combs[[i]][1,],
               function(x, y){unlist(x) - unlist(y)})
      })

      # Assemble table of differences
      # *Note: Because predictions are organized by response level, we have to
      # assemble the table of differences by response level too, otherwise the
      # order of the legend entries won't match the order of the differences.
      diff_res <- NULL
      diff_res <- lapply(1:length(diff_list), function(i){
        diff_res[[i]] <- expand.grid(
          x2 = 1:nlev,
          x1 = 1:nlev) %>%
          dplyr::filter(x2 > x1) %>%
          dplyr::mutate(x = paste0(levs[x2],
                                   " - ",
                                   levs[x1])) %>%
          dplyr::select(-c("x1", "x2")) %>%
          dplyr::mutate(
            y = obj$lev[i],
            predicted = colMeans(diff_list[[i]]),
            conf.low = apply(diff_list[[i]], 2, low),
            conf.high = apply(diff_list[[i]], 2, high))
      })
      diffs <- do.call("rbind", diff_res)
      diffs <- diffs %>%
        dplyr::mutate(y = factor(y, levels=obj$lev)) %>%
        dplyr::relocate(y, .before=1) %>%
        dplyr::mutate(type = "Difference") %>%
        dplyr::arrange(x)

    } else {

      ## Differences for continuous variables ----------------------------------

      diffchange <- match.arg(diffchange)
      tmp0 <- switch(diffchange,
                     NULL = survey::svymean(data[[varname]], svydata) -
                       (sqrt(survey::svyvar(data[[varname]], svydata))/2),
                     sd = survey::svymean(data[[varname]], svydata) -
                       (sqrt(survey::svyvar(data[[varname]], svydata))/2),
                     range = min(data[[varname]]),
                     unit = survey::svymean(data[[varname]], svydata) - .5
                     )
      tmp1 <- switch(diffchange,
                     NULL = survey::svymean(data[[varname]], svydata) +
                       (sqrt(survey::svyvar(data[[varname]], svydata))/2),
                     sd = survey::svymean(data[[varname]], svydata) +
                       (sqrt(survey::svyvar(data[[varname]], svydata))/2),
                     range = max(data[[varname]]),
                     unit = survey::svymean(data[[varname]], svydata) + .5
                     )
      delta <- c(tmp0, tmp1)

      # Define variables to vary
      varylist <- as.list(varname)
      varylist <- lapply(1:length(varylist), function(i) {varylist[[i]] <- delta})
      names(varylist) <- varname

      # Create fake dataset for simulations
      trms <- attr(terms(obj), "dataClasses")
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
      fake <- do.call(expand.grid, l)   # matrix for generating predictions

      # Parametric bootstrap
      B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))
      tau_ind <- grep("\\|", names(coef(obj)))
      Tau <- B[,tau_ind]
      B <- B[,-tau_ind]
      fakeX <- model.matrix(formula(obj), data=fake)
      fakeX

      # Generate predictions
      CP <- NULL
      for (i in 1:ncol(Tau)) {
        CP[[i]] <- plogis(matrix(Tau[,i],
                                 ncol=nrow(Tau),
                                 nrow=nrow(fakeX[,-1]),
                                 byrow=TRUE) - fakeX[,-1] %*% t(B))
      }
      P <- NULL
      for (i in 1:ncol(Tau)) {
        if(i==1) {
          P[[i]] <- CP[[i]]
        }else{
          P[[i]] <- CP[[i]]- CP[[(i-1)]]
        }}
      P[[(ncol(Tau)+1)]] <- 1-CP[[(ncol(Tau))]]

      # Calculate differences
      D <- lapply(1:length(P), function(i){P[[i]][2,] - P[[i]][1,]})

      # Assemble table of differences
      diffs <- expand.grid(y = obj$lev,
                           x = paste0("Delta (",
                                      diffchange,
                                      ") : ",
                                      round(tmp0, 3),
                                      " - ",
                                      round(tmp1, 3)))
      diffs$predicted <- sapply(1:length(obj$lev), function(i)mean(D[[i]]))
      diffs$conf.low <- sapply(1:length(obj$lev), function(i)low(D[[i]]))
      diffs$conf.high <- sapply(1:length(obj$lev), function(i)high(D[[i]]))
      diffs$type <- "Contrast"
      diffs$y <- factor(diffs$y, levels=obj$lev)
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
      formula = formula(obj))
    class(output) <- "svyEffects"
    attributes(output)$predvar <- varname
    attributes(output)$depvar <- colnames(model.frame(obj))[1]
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

    if(is.factor(data[[varname]])) {
      Xlevs <- levels(data[[varname]])}
    if(is.numeric(data[[varname]])) {
      Xlevs <- seq(min(data[[varname]]), max(data[[varname]]), length=nvals)}
    bylevs <- 1   # so length(nlev) * length(bylev) still works in array dimension

    if(!is.null(byvar)) {
      if(is.factor(data[[byvar]])) {
        bylevs <- levels(data[[byvar]])}
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
    }

    # Define variables to vary
    xvarylist <- as.list(varname)
    xvarylist <- lapply(1:length(varname), function(i) {
      if(is.factor(data[[xvarylist[[i]]]])) {
        xvarylist[[i]] <- levels(data[[varname[[i]]]])
      } else {
        xvarylist[[i]] <- seq(min(data[[varname[[i]]]]), max(data[[varname[[i]]]]), length=nvals)
      }})
    names(xvarylist) <- varname
    varylist <- xvarylist
    # If there is a moderator variable
    byvarylist <- as.list(byvar)
    byvarylist <- lapply(1:length(byvar), function(i) {
      if(is.factor(data[[byvarylist[[i]]]])) {
        byvarylist[[i]] <- levels(data[[byvar[[i]]]])
      } else {
        byvarylist[[i]] <- seq(min(data[[byvar[[i]]]]), max(data[[byvar[[i]]]]), length=bynvals)
      }})
    names(byvarylist) <- byvar
    varylist <- c(xvarylist, byvarylist)

    # Create fake dataset for simulations
    trms <- attr(terms(obj), "dataClasses")
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
          l[[m]] <- svyquantile(reformulate(names(trms)[m]), svydata, 0.5)[[1]][1]
        }}}
    fake <- do.call(expand.grid, l)   # matrix for generating predictions

    # Parametric bootstrap
    B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))
    tau_ind <- grep("\\|", names(coef(obj)))
    Tau <- B[,tau_ind]
    B <- B[,-tau_ind]

    # Generate predictions
    fakeX <- model.matrix(formula(obj), data=fake)
    CP <- NULL
    for (i in 1:ncol(Tau)) {
      CP[[i]] <- plogis(matrix(Tau[,i],
                               ncol=nrow(Tau),
                               nrow=nrow(fakeX[,-1]),
                               byrow=TRUE) - fakeX[,-1] %*% t(B))
    }
    P <- NULL
    for (i in 1:ncol(Tau)) {
      if(i==1) {P[[i]] <- CP[[i]]}
      else{P[[i]] <- CP[[i]]- CP[[(i-1)]]
      }
    }
    P[[(ncol(Tau)+1)]] <- 1-CP[[(ncol(Tau))]]

    # Assemble predictions
    # Label predictions with correct combination of xvar * byvar (b/c order of preds
    # determined by order of terms in model.formula, not specified varname/byvar)
    if(grep(varname, names(model.frame(obj))) < grep(byvar, names(model.frame(obj)))){
      v1 <- data[varname]} else {v1 <- data[byvar]}
    if(grep(varname, names(model.frame(obj))) > grep(byvar, names(model.frame(obj)))){
      v2 <- data[varname]} else {v2 <- data[byvar]}

    # So numerical moderators show up as factors (for plotting purposes)
    if (!is.null(byvar) & is.numeric(data[[byvar]])) {
      varylist[[match(byvar, names(varylist))]] <- factor(round(bylevs, 3))
    }

    # Assemble predictions table
    preds <- expand.grid(
      x1 = varylist[[match(names(v1), names(varylist))]],
      x2 = varylist[[match(names(v2), names(varylist))]],
      y = obj$lev)  %>%
      mutate(predicted = unlist(lapply(1:length(obj$lev), function(i)apply(P[[i]], 1, mean))),
             conf.low = unlist(lapply(1:length(obj$lev), function(i)apply(P[[i]], 1, low))),
             conf.high = unlist(lapply(1:length(obj$lev), function(i)apply(P[[i]], 1, high))),
             type = "Probability")
    preds <- preds %>%
      rename(!!sym(names(v1)) := x1,
             !!sym(names(v2)) := x2) %>%
      relocate(y, .before=1)


    # Assemble combined output table
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
