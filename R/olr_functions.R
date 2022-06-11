# METHODS FOR ORDINAL DEPENDENT VARIABLE MODELS






#' Average marginal effects (a.k.a. marginal effects at observed values) for
#' ordinal dependent-variable models of survey-weighted data (survey::svyolr).
#'
#' @param obj Model object of class \code{survey::svyolr}.
#' @param varname Character string denoting the name of the predictor variable
#' for which effects are to be calculated.
#' @param weightvar Character string denoting the name of the sampling weight variable.
#' @param design A survey design object of class \code{survey::svydesign}.
#' @param nvals Scalar denoting the sequence length spanning the range of a
#' continuous variable for which effects are to be calculated (default: 11).
#' @param diffchange Character string  denoting over what change in x a first
#' difference is to be calculated for a continuous predictor (default: range).
#' @param sims Scalar denoting how many simulations to conduct (default: 1500).
#' @param seed Seed value for random number generator. By default, the function
#' picks a random integer between 1 and 1,000,000. If you save the output of
#' this function, it will save the seed value used for simulations in the slot
#' \code{$seed}.
#' @param ... Other arguments (currently not implemented).
#'
#' @return A list with dataframes:
#'   \describe{
#'     \item{\code{$preds}}{predicted probabilities}
#'     \item{\code{$preds}}{differences in predicted probabilities}
#'     \item{\code{$seed}}{seed value used for random number generator}
#'     }
#'
#' @importFrom dplyr arrange as_tibble filter mutate rename relocate select summarize tibble
#' @importFrom magrittr %>%
#' @importFrom MASS mvrnorm
#' @importFrom rlang sym
#' @importFrom survey svydesign svymean svyquantile svytable
#' @importFrom tidyr pivot_longer pivot_wider
#'
#' @export
svyAME.svyolr <- function(obj,
                          varname,
                          weightvar,
                          design,
                          nvals = 11,
                          diffchange = c("range", "unit", "sd"),
                          sims = 1500,
                          seed = NULL,
                          ...) {

  #========== SETUP ============================================================

  # Get data
  data <- model.frame(obj)
  data <- data %>% dplyr::rename(!!sym(weightvar) := `(weights)`)
  svydata <- survey::svydesign(ids=~1, strata=NULL, weights=data[weightvar], data=data)
  if(is.character(data[[varname]])) {data[[varname]] <- as.factor(data[[varname]])}

  # Set random number generator
  if(is.null(seed)){
    seed <- round(runif(1, min=0, max=1000000))
  } else {
    seed <- seed
  }
  set.seed(seed)

  #========== NUMERIC / CONTINUOUS predictors ==================================

  if(is.numeric(data[[varname]])) {

    ## Predicted probabilities ------------------------------------------------

    varname_seq <- seq(min(data[[varname]]), max(data[[varname]]), length=nvals)
    preds <- NULL
    for(n in seq_along(varname_seq)){
      new <- dplyr::mutate(data, !!sym(varname) := varname_seq[n])
      Xmat <- model.matrix(formula(obj), data=new)[,-1]
      B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))
      tau_ind <- grep("\\|", names(coef(obj)))
      Tau <- B[,tau_ind]
      B <- B[,-tau_ind]

      CP <- NULL
      for (i in 1:ncol(Tau)) {
        CP[[i]] <- plogis(matrix(Tau[,i],
                                 ncol=nrow(Tau),
                                 nrow=nrow(new),
                                 byrow=TRUE) - Xmat %*% t(B))
      }

      P <- NULL
      for (i in 1:ncol(Tau)) {
        if(i==1) {P[[i]] <- CP[[i]]}
        else{P[[i]] <- CP[[i]]- CP[[(i-1)]]}}
      P[[(ncol(Tau)+1)]] <- 1-CP[[ncol(Tau)]]

      WgtP <- lapply(1:length(P), function(i){
        apply(P[[i]], 2, function(x)weighted.mean(x, data[[weightvar]]))
      })

      tmp <- data.frame(
        y = obj$lev,
        x = varname_seq[n],
        predicted = sapply(WgtP, mean),
        conf.low = sapply(WgtP, low),
        conf.high = sapply(WgtP, high),
        type = "Probability")
      preds <- rbind(preds, tmp)
    }
    preds$y <- factor(preds$y, levels=obj$lev)

    ## Differences in predicted probabilities ----------------------------------

    tmp0 <- switch(diffchange,
                   NULL = min(data[varname]),
                   range = min(data[varname]),
                   unit = survey::svymean(data[varname], svydata) - .5,
                   sd = survey::svymean(data[varname], svydata) -
                     (sqrt(survey::svyvar(data[varname], svydata))/2)
    )
    tmp1 <- switch(diffchange,
                   NULL = max(data[varname]),
                   range = max(data[varname]),
                   unit = survey::svymean(data[varname], svydata) + .5,
                   sd = survey::svymean(data[varname], svydata) +
                     (sqrt(survey::svyvar(data[varname], svydata))/2)
    )

    delta <- c(tmp0, tmp1)
    D_WgtP <- NULL

    for(n in seq_along(delta)){
      new <- dplyr::mutate(data, !!sym(varname) := delta[n])
      Xmat <- model.matrix(formula(obj), data=new)[,-1]
      B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))
      tau_ind <- grep("\\|", names(coef(obj)))
      Tau <- B[,tau_ind]
      B <- B[,-tau_ind]

      CP <- NULL
      for (i in 1:ncol(Tau)) {
        CP[[i]] <- plogis(matrix(Tau[,i], ncol=nrow(Tau),
                                 nrow=nrow(new), byrow=TRUE) - Xmat %*% t(B))
      }
      P <- NULL
      for (i in 1:ncol(Tau)) {
        if(i==1) {P[[i]] <- CP[[i]]}
        else{P[[i]] <- CP[[i]]- CP[[(i-1)]]}}
      P[[(ncol(Tau)+1)]] <- 1-CP[[ncol(Tau)]]

      WgtP <- lapply(1:length(P), function(i){
        apply(P[[i]], 2, function(x)weighted.mean(x, data[[weightvar]]))
      })

      D_WgtP[[n]] <- WgtP
    }

    D_WgtP <- mapply(x = D_WgtP[[2]], y = D_WgtP[[1]],
                     function(x, y){unlist(x) - unlist(y)})

    diffs <- expand.grid(y = obj$lev,
                         x = paste0("Delta (",
                                    diffchange,
                                    ") : ",
                                    round(tmp0, 3),
                                    " - ",
                                    round(tmp1, 3)))
    diffs$predicted <- colMeans(D_WgtP)
    diffs$conf.low <- apply(D_WgtP, 2, low)
    diffs$conf.high <- apply(D_WgtP, 2, high)
    diffs$type <- "Difference"

    ## Output ------------------------------------------------------------------

    preds <- rename(preds, !!sym(varname) := x)
    diffs <- rename(diffs, !!sym(varname) := x)
    output <- list(
      preds = dplyr::as_tibble(preds),
      diffs = dplyr::as_tibble(diffs),
      seed = seed)
    class(output) <- "svyEffects"
    attributes(output)$predvar <- varname
    attributes(output)$depvar <- colnames(model.frame(obj))[1]
    return(output)

  } else {

    #========== FACTOR / CATEGORICAL predictors ================================

    ## Predicted probabilities -------------------------------------------------

    nlev <- as.numeric(nlevels(data[[varname]]))
    levs <- levels(data[[varname]])
    nlevC2 <- (factorial(nlev)/(2*factorial(nlev-2)))

    # Create dataframes for simulations
    DFs <- lapply(1:nlev, function(x) {x <- data})
    for (i in 1:nlev) {DFs[[i]] <- DFs[[i]] %>%
      dplyr::mutate(!!sym(varname) := factor(i, levels=1:nlev, labels=levs))}
    X <- NULL
    for (i in 1:nlev) {X[[i]] <- model.matrix(formula(obj), data=DFs[[i]])[,-1]}

    # Create coefficient vectors
    B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))
    tau_ind <- grep("\\|", names(coef(obj)))
    Tau <- B[,tau_ind]
    B <- B[,-tau_ind]

    # Generate predictions
    CP <- NULL
    P <- NULL
    Pr <- NULL
    for (k in 1:nlev) {
      for (i in 1:ncol(Tau)) {
        CP[[i]] <- plogis(matrix(Tau[,i], ncol=nrow(Tau),
                                 nrow=nrow(X[[k]]), byrow=TRUE) - X[[k]] %*% t(B))
      }
      for (j in 1:ncol(Tau)) {
        if(j==1) {
          P[[j]] <- CP[[j]]
        }
        else{
          P[[j]] <- CP[[j]]- CP[[(j-1)]]
        }
        P[[(ncol(Tau)+1)]] <- 1-CP[[(ncol(Tau))]]
      }
      Pr[[k]] <- P
    }

    # Apply weights
    Pr_ul <- unlist(Pr, recursive=FALSE)
    WgtPr <- lapply(1:length(Pr_ul), function(i){
      apply(Pr_ul[[i]], 2, function(x)weighted.mean(x, data[[weightvar]]))
    })

    # Assemble table of predicted probabilities
    preds <- expand.grid(y = obj$lev, x = levs)
    preds$predicted <- sapply(WgtPr, mean)
    preds$conf.low <- sapply(WgtPr, low)
    preds$conf.high <- sapply(WgtPr, high)
    preds$type <- "Probability"

    ## Differences in predicted probabilities ----------------------------------

    WgtPr <- as.data.frame(WgtPr)
    i <- 1
    j <- 1
    P <- NULL
    while (i <= nlev) {
      P[[i]] <- as.data.frame(WgtPr[,j:(j+length(obj$lev)-1)])
      i <- i + 1
      j <- j + length(obj$lev)
    }

    P_combs <- combn(P, 2, simplify=FALSE)
    D <- lapply(1:length(P_combs), function(i)P_combs[[i]][[2]] - P_combs[[i]][[1]])
    diff_res <- cbind.data.frame(D)

    # Assemble table of differences
    diffs <- expand.grid(
      y = obj$lev,
      x2 = 1:nlev,
      x1 = 1:nlev) %>%
      dplyr::filter(x2 > x1) %>%
      dplyr::mutate(x = paste0(levs[x2],
                               " - ",
                               levs[x1])) %>%
      dplyr::select(-c("x1", "x2"))
    diffs <- diffs %>%
      dplyr::mutate(
        predicted = colMeans(diff_res),
        conf.low = apply(diff_res, 2, low),
        conf.high = apply(diff_res, 2, high),
        type = "Difference")

    ## Output ------------------------------------------------------------------

    preds <- rename(preds, !!sym(varname) := x)
    diffs <- rename(diffs, !!sym(varname) := x)
    output <- list(
      preds = dplyr::as_tibble(preds),
      diffs = dplyr::as_tibble(diffs),
      seed = seed)
    class(output) <- "svyEffects"
    attributes(output)$predvar <- varname
    attributes(output)$depvar <- colnames(model.frame(obj))[1]
    return(output)

  }
}






#' Marginal effects at reasonable values for ordinal dependent-variable models
#' ordinal dependent-variable models of survey-weighted data (survey::svyolr).
#'
#' @param obj Model object of class \code{survey::svyolr}.
#' @param varname Character string denoting the name of the predictor variable
#' for which effects are to be calculated.
#' @param weightvar Character string denoting the name of the sampling weight  variable.
#' @param design A survey design object of class \code{survey::svydesign}.
#' @param nvals Scalar denoting the sequence length spanning the range of a
#' continuous variable for which effects are to be calculated (default: 10).
#' @param diffchange Character string  denoting over what change in x a first
#' difference is to be calculated for a continuous predictor (default: range).
#' @param sims Scalar denoting how many simulations to conduct (default: 1500).
#' @param seed Seed value for random number generator. By default, the function
#' picks a random integer between 1 and 1,000,000. If you save the output of
#' this function, it will save the seed value used for simulations in the slot
#' \code{$seed}.
#' @param ... Other arguments (currently not implemented).
#'
#' @return A list with dataframes:
#'   \describe{
#'     \item{\code{$preds}}{predicted probabilities}
#'     \item{\code{$preds}}{differences in predicted probabilities}
#'     \item{\code{$seed}}{seed value used for random number generator}
#'     \item{\code{$typical}}{the values at which variables are set for the simulations}
#'     }
#'
#' @importFrom dplyr arrange as_tibble mutate rename relocate select summarize tibble
#' @importFrom magrittr %>%
#' @importFrom MASS mvrnorm
#' @importFrom rlang sym
#' @importFrom survey svydesign svymean svyquantile svytable
#' @importFrom tidyr pivot_longer pivot_wider
#'
#' @examples
#'
#' @export
svyMER.svyolr <- function(obj,
                          varname,
                          weightvar,
                          design,
                          nvals = 11,
                          diffchange = c("range", "unit", "sd"),
                          sims = 1500,
                          seed = NULL,
                          ...) {

  #========== SETUP ============================================================

  # Get data
  data <- model.frame(obj)
  data <- data %>% dplyr::rename(!!sym(weightvar) := `(weights)`)
  svydata <- survey::svydesign(ids=~1, strata=NULL, weights=data[weightvar], data=data)
  if(is.character(data[[varname]])) {data[[varname]] <- as.factor(data[[varname]])}

  # Set random number generator
  if(is.null(seed)){
    seed <- round(runif(1, min=0, max=1000000))
  } else {
    seed <- seed
  }
  set.seed(seed)

  # Define variables to vary
  varylist <- as.list(varname)
  varylist <- lapply(1:length(varylist), function(i) {
    if(class(data[[varylist[[i]]]])=="numeric") {
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

  #========== PREDICTED PROBABILITIES ==========================================

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
  if(class(data[[varname]])=="factor") {
    levs <- levels(data[[varname]])}
  if(class(data[[varname]])=="numeric") {
    levs <- seq(min(data[[varname]]), max(data[[varname]]), length=nvals)}
  preds <- expand.grid(
    x = levs,
    y = obj$lev)  %>%
    mutate(predicted = unlist(lapply(1:length(obj$lev), function(i)apply(P[[i]], 1, mean))),
           conf.low = unlist(lapply(1:length(obj$lev), function(i)apply(P[[i]], 1, low))),
           conf.high = unlist(lapply(1:length(obj$lev), function(i)apply(P[[i]], 1, high))),
           type = "Probability")

  #========== DIFFERENCES IN PREDICTED PROBABILITIES ===========================

  ## Differences for categorical variables -------------------------------------

  if(class(data[[varname]])=="factor") {
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

    ## Differences for continuous variables ------------------------------------

    tmp0 <- switch(diffchange,
                   NULL = min(data[[varname]]),
                   range = min(data[[varname]]),
                   unit = survey::svymean(data[[varname]], svydata) - .5,
                   sd = survey::svymean(data[[varname]], svydata) -
                     (sqrt(survey::svyvar(data[[varname]], svydata))/2)
    )
    tmp1 <- switch(diffchange,
                   NULL = max(data[[varname]]),
                   range = max(data[[varname]]),
                   unit = survey::svymean(data[[varname]], svydata) + .5,
                   sd = survey::svymean(data[[varname]], svydata) +
                     (sqrt(survey::svyvar(data[[varname]], svydata))/2)
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

  #========== Output ===========================================================
  preds <- rename(preds, !!sym(varname) := x)
  diffs <- rename(diffs, !!sym(varname) := x)
  output <- list(
    preds = dplyr::as_tibble(preds),
    diffs = dplyr::as_tibble(diffs),
    seed = seed,
    typical = dplyr::as_tibble(fake))
  class(output) <- "svyEffects"
  attributes(output)$predvar <- varname
  attributes(output)$depvar <- colnames(obj$model)[1]
  return(output)

}





#' Average marginal effects (a.k.a. marginal effects at observed values) for
#' ordinal dependent-variable models of survey-weighted data.
#'
#' @param obj Model object of class \code{survey::svyolr}.
#' @param varname Character string denoting the name of the predictor variable
#' for which effects are to be calculated.
#' @param weightvar Character string denoting the name of the sampling weight  variable.
#' @param design A survey design object of class \code{survey::svydesign}.
#' @param nvals Scalar denoting the sequence length spanning the range of a
#' continuous variable for which effects are to be calculated (default: 10).
#' @param diffchange Character string  denoting over what change in x a first
#' difference is to be calculated for a continuous predictor (default: range).
#' @param sims Scalar denoting how many simulations to conduct (default: 1500).
#' @param seed Seed value for random number generator. By default, the function
#' picks a random integer between 1 and 1,000,000. If you save the output of
#' this function, it will save the seed value used for simulations in the slot
#' \code{$seed}.
#' @param ... Other arguments (currently not implemented).
#'
#' @return A list with dataframes:
#'   \describe{
#'     \item{\code{$preds}}{predicted probabilities}
#'     \item{\code{$preds}}{differences in predicted probabilities}
#'     \item{\code{$seed}}{seed value used for random number generator}
#'     }
#'
#' @importFrom dplyr arrange as_tibble filter mutate rename relocate select summarize tibble
#' @importFrom magrittr %>%
#' @importFrom MASS mvrnorm
#' @importFrom rlang sym
#' @importFrom survey svydesign svymean svyquantile svytable
#' @importFrom tidyr pivot_longer pivot_wider
#'
#' @export
svyAME.polr <- function(obj,
                        varname,
                        weightvar,
                        design,
                        nvals = 11,
                        diffchange = c("range", "unit", "sd"),
                        sims = 1500,
                        seed = NULL,
                        ...) {

  #========== SETUP ============================================================

  # Get data
  data <- model.frame(obj)
  data <- data %>% dplyr::rename(!!sym(weightvar) := `(weights)`)
  svydata <- survey::svydesign(ids=~1, strata=NULL, weights=data[weightvar], data=data)
  if(is.character(data[[varname]])) {data[[varname]] <- as.factor(data[[varname]])}

  # Set random number generator
  if(is.null(seed)){
    seed <- round(runif(1, min=0, max=1000000))
  } else {
    seed <- seed
  }
  set.seed(seed)

  #========== NUMERIC / CONTINUOUS predictors ==================================

  if(is.numeric(data[[varname]])) {

    ## Predicted probabilities ------------------------------------------------

    varname_seq <- seq(min(data[[varname]]), max(data[[varname]]), length=nvals)
    preds <- NULL
    for(n in seq_along(varname_seq)){
      new <- dplyr::mutate(data, !!sym(varname) := varname_seq[n])
      Xmat <- model.matrix(formula(obj), data=new)[,-1]
      B <- MASS::mvrnorm(sims, c(coef(obj), obj$zeta), vcov(obj))
      tau_ind <- grep("\\|", names(c(coef(obj), obj$zeta)))
      Tau <- B[,tau_ind]
      B <- B[,-tau_ind]

      CP <- NULL
      for (i in 1:ncol(Tau)) {
        CP[[i]] <- plogis(matrix(Tau[,i],
                                 ncol=nrow(Tau),
                                 nrow=nrow(new),
                                 byrow=TRUE) - Xmat %*% t(B))
      }

      P <- NULL
      for (i in 1:ncol(Tau)) {
        if(i==1) {P[[i]] <- CP[[i]]}
        else{P[[i]] <- CP[[i]]- CP[[(i-1)]]}}
      P[[(ncol(Tau)+1)]] <- 1-CP[[ncol(Tau)]]

      WgtP <- lapply(1:length(P), function(i){
        apply(P[[i]], 2, function(x)weighted.mean(x, data[[weightvar]]))
      })

      tmp <- data.frame(
        y = obj$lev,
        x = varname_seq[n],
        predicted = sapply(WgtP, mean),
        conf.low = sapply(WgtP, low),
        conf.high = sapply(WgtP, high),
        type = "Probability")
      preds <- rbind(preds, tmp)
    }
    preds$y <- factor(preds$y, levels=obj$lev)

    ## Differences in predicted probabilities ----------------------------------

    tmp0 <- switch(diffchange,
                   NULL = min(data[varname]),
                   range = min(data[varname]),
                   unit = survey::svymean(data[varname], svydata) - .5,
                   sd = survey::svymean(data[varname], svydata) -
                     (sqrt(survey::svyvar(data[varname], svydata))/2)
    )
    tmp1 <- switch(diffchange,
                   NULL = max(data[varname]),
                   range = max(data[varname]),
                   unit = survey::svymean(data[varname], svydata) + .5,
                   sd = survey::svymean(data[varname], svydata) +
                     (sqrt(survey::svyvar(data[varname], svydata))/2)
    )

    delta <- c(tmp0, tmp1)
    D_WgtP <- NULL

    for(n in seq_along(delta)){
      new <- dplyr::mutate(data, !!sym(varname) := delta[n])
      Xmat <- model.matrix(formula(obj), data=new)[,-1]
      B <- MASS::mvrnorm(sims, c(coef(obj), obj$zeta), vcov(obj))
      tau_ind <- grep("\\|", names(c(coef(obj), obj$zeta)))
      Tau <- B[,tau_ind]
      B <- B[,-tau_ind]

      CP <- NULL
      for (i in 1:ncol(Tau)) {
        CP[[i]] <- plogis(matrix(Tau[,i], ncol=nrow(Tau),
                                 nrow=nrow(new), byrow=TRUE) - Xmat %*% t(B))
      }
      P <- NULL
      for (i in 1:ncol(Tau)) {
        if(i==1) {P[[i]] <- CP[[i]]}
        else{P[[i]] <- CP[[i]]- CP[[(i-1)]]}}
      P[[(ncol(Tau)+1)]] <- 1-CP[[ncol(Tau)]]

      WgtP <- lapply(1:length(P), function(i){
        apply(P[[i]], 2, function(x)weighted.mean(x, data[[weightvar]]))
      })

      D_WgtP[[n]] <- WgtP
    }

    D_WgtP <- mapply(x = D_WgtP[[2]], y = D_WgtP[[1]],
                     function(x, y){unlist(x) - unlist(y)})

    diffs <- expand.grid(y = obj$lev,
                         x = paste0("Delta (",
                                    diffchange,
                                    ") : ",
                                    round(tmp0, 3),
                                    " - ",
                                    round(tmp1, 3)))
    diffs$predicted <- colMeans(D_WgtP)
    diffs$conf.low <- apply(D_WgtP, 2, low)
    diffs$conf.high <- apply(D_WgtP, 2, high)
    diffs$type <- "Difference"

    ## Output ------------------------------------------------------------------

    preds <- rename(preds, !!sym(varname) := x)
    diffs <- rename(diffs, !!sym(varname) := x)
    output <- list(
      preds = dplyr::as_tibble(preds),
      diffs = dplyr::as_tibble(diffs),
      seed = seed)
    class(output) <- "svyEffects"
    attributes(output)$predvar <- varname
    attributes(output)$depvar <- colnames(model.frame(obj))[1]
    return(output)

  } else {

    #========== FACTOR / CATEGORICAL predictors ================================

    ## Predicted probabilities -------------------------------------------------

    nlev <- as.numeric(nlevels(data[[varname]]))
    levs <- levels(data[[varname]])
    nlevC2 <- (factorial(nlev)/(2*factorial(nlev-2)))

    # Create dataframes for simulations
    DFs <- lapply(1:nlev, function(x) {x <- data})
    for (i in 1:nlev) {DFs[[i]] <- DFs[[i]] %>%
      dplyr::mutate(!!sym(varname) := factor(i, levels=1:nlev, labels=levs))}
    X <- NULL
    for (i in 1:nlev) {X[[i]] <- model.matrix(formula(obj), data=DFs[[i]])[,-1]}

    # Create coefficient vectors
    B <- MASS::mvrnorm(sims, c(coef(obj), obj$zeta), vcov(obj))
    tau_ind <- grep("\\|", names(c(coef(obj), obj$zeta)))
    Tau <- B[,tau_ind]
    B <- B[,-tau_ind]

    # Generate predictions
    CP <- NULL
    P <- NULL
    Pr <- NULL
    for (k in 1:nlev) {
      for (i in 1:ncol(Tau)) {
        CP[[i]] <- plogis(matrix(Tau[,i], ncol=nrow(Tau),
                                 nrow=nrow(X[[k]]), byrow=TRUE) - X[[k]] %*% t(B))
      }
      for (j in 1:ncol(Tau)) {
        if(j==1) {
          P[[j]] <- CP[[j]]
        }
        else{
          P[[j]] <- CP[[j]]- CP[[(j-1)]]
        }
        P[[(ncol(Tau)+1)]] <- 1-CP[[(ncol(Tau))]]
      }
      Pr[[k]] <- P
    }

    # Apply weights
    Pr_ul <- unlist(Pr, recursive=FALSE)
    WgtPr <- lapply(1:length(Pr_ul), function(i){
      apply(Pr_ul[[i]], 2, function(x)weighted.mean(x, data[[weightvar]]))
    })

    # Assemble table of predicted probabilities
    preds <- expand.grid(y = obj$lev, x = levs)
    preds$predicted <- sapply(WgtPr, mean)
    preds$conf.low <- sapply(WgtPr, low)
    preds$conf.high <- sapply(WgtPr, high)
    preds$type <- "Probability"

    ## Differences in predicted probabilities ----------------------------------

    WgtPr <- as.data.frame(WgtPr)
    i <- 1
    j <- 1
    P <- NULL
    while (i <= nlev) {
      P[[i]] <- as.data.frame(WgtPr[,j:(j+length(obj$lev)-1)])
      i <- i + 1
      j <- j + length(obj$lev)
    }

    P_combs <- combn(P, 2, simplify=FALSE)
    D <- lapply(1:length(P_combs), function(i)P_combs[[i]][[2]] - P_combs[[i]][[1]])
    diff_res <- cbind.data.frame(D)

    # Assemble table of differences
    diffs <- expand.grid(
      y = obj$lev,
      x2 = 1:nlev,
      x1 = 1:nlev) %>%
      dplyr::filter(x2 > x1) %>%
      dplyr::mutate(x = paste0(levs[x2],
                               " - ",
                               levs[x1])) %>%
      dplyr::select(-c("x1", "x2"))
    diffs <- diffs %>%
      dplyr::mutate(
        predicted = colMeans(diff_res),
        conf.low = apply(diff_res, 2, low),
        conf.high = apply(diff_res, 2, high),
        type = "Difference")

    ## Output ------------------------------------------------------------------

    preds <- rename(preds, !!sym(varname) := x)
    diffs <- rename(diffs, !!sym(varname) := x)
    output <- list(
      preds = dplyr::as_tibble(preds),
      diffs = dplyr::as_tibble(diffs),
      seed = seed)
    class(output) <- "svyEffects"
    attributes(output)$predvar <- varname
    attributes(output)$depvar <- colnames(model.frame(obj))[1]
    return(output)

  }
}






#' Marginal effects at reasonable values for ordinal dependent-variable models
#' of weighted data (using MASS::polr).
#'
#' @param obj Model object of class \code{MASS::polr}.
#' @param varname Character string denoting the name of the predictor variable
#' for which effects are to be calculated.
#' @param weightvar Character string denoting the name of the sampling weight  variable.
#' @param design A survey design object of class \code{survey::svydesign}.
#' @param nvals Scalar denoting the sequence length spanning the range of a
#' continuous variable for which effects are to be calculated (default: 10).
#' @param diffchange Character string  denoting over what change in x a first
#' difference is to be calculated for a continuous predictor (default: range).
#' @param sims Scalar denoting how many simulations to conduct (default: 1500).
#' @param seed Seed value for random number generator. By default, the function
#' picks a random integer between 1 and 1,000,000. If you save the output of
#' this function, it will save the seed value used for simulations in the slot
#' \code{$seed}.
#' @param ... Other arguments (currently not implemented).
#'
#' @return A list with dataframes:
#'   \describe{
#'     \item{\code{$preds}}{predicted probabilities}
#'     \item{\code{$preds}}{differences in predicted probabilities}
#'     \item{\code{$seed}}{seed value used for random number generator}
#'     \item{\code{$typical}}{the values at which variables are set for the simulations}
#'     }
#'
#' @importFrom dplyr arrange as_tibble filter mutate rename relocate select summarize tibble
#' @importFrom magrittr %>%
#' @importFrom MASS mvrnorm
#' @importFrom rlang sym
#' @importFrom survey svydesign svymean svyquantile svytable
#' @importFrom tidyr pivot_longer pivot_wider
#'
#' @export
svyMER.polr <- function(obj,
                        varname,
                        weightvar,
                        design,
                        nvals = 11,
                        diffchange = c("range", "unit", "sd"),
                        sims = 1500,
                        seed = NULL,
                        ...) {

  #========== SETUP ============================================================

  # Get data
  data <- model.frame(obj)
  data <- data %>% dplyr::rename(!!sym(weightvar) := `(weights)`)
  svydata <- survey::svydesign(ids=~1, strata=NULL, weights=data[weightvar], data=data)
  if(is.character(data[[varname]])) {data[[varname]] <- as.factor(data[[varname]])}

  # Set random number generator
  if(is.null(seed)){
    seed <- round(runif(1, min=0, max=1000000))
  } else {
    seed <- seed
  }
  set.seed(seed)

  # Define variables to vary
  varylist <- as.list(varname)
  varylist <- lapply(1:length(varylist), function(i) {
    if(class(data[[varylist[[i]]]])=="numeric") {
      varylist[[i]] <- seq(min(data[[varylist[[i]]]]), max(data[[varylist[[i]]]]), length=nvals)
    } else {
      varylist[[i]] <- levels(data[[varylist[[i]]]])
    }})
  names(varylist) <- varname

  # Create fake dataset for simulations
  trms <- attr(terms(obj), "dataClasses")
  trms <- trms[-length(trms)]   # `polr` class objects have weights in the termlist!
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

  #========== PREDICTED PROBABILITIES ==========================================

  # Parametric bootstrap
  B <- MASS::mvrnorm(sims, c(coef(obj), obj$zeta), vcov(obj))
  tau_ind <- grep("\\|", names(c(coef(obj), obj$zeta)))
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
  if(class(data[[varname]])=="factor") {
    levs <- levels(data[[varname]])}
  if(class(data[[varname]])=="numeric") {
    levs <- seq(min(data[[varname]]), max(data[[varname]]), length=nvals)}
  preds <- expand.grid(
    x = levs,
    y = obj$lev)  %>%
    mutate(predicted = unlist(lapply(1:length(obj$lev), function(i)apply(P[[i]], 1, mean))),
           conf.low = unlist(lapply(1:length(obj$lev), function(i)apply(P[[i]], 1, low))),
           conf.high = unlist(lapply(1:length(obj$lev), function(i)apply(P[[i]], 1, high))),
           type = "Probability")

  #========== DIFFERENCES IN PREDICTED PROBABILITIES ===========================

  ## Differences for categorical variables -------------------------------------

  if(class(data[[varname]])=="factor") {
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
        filter(x2 > x1) %>%
        dplyr::mutate(x = paste0(levs[x2],
                                 " - ",
                                 levs[x1])) %>%
        dplyr::select(-c("x1", "x2")) %>%
        mutate(
          y = obj$lev[i],
          predicted = colMeans(diff_list[[i]]),
          conf.low = apply(diff_list[[i]], 2, low),
          conf.high = apply(diff_list[[i]], 2, high))
    })
    diffs <- do.call("rbind", diff_res)
    diffs <- diffs %>%
      mutate(y = factor(y, levels=obj$lev)) %>%
      relocate(y, .before=1) %>%
      mutate(type = "Difference") %>%
      arrange(x)

  } else {

    ## Differences for continuous variables ------------------------------------

    tmp0 <- switch(diffchange,
                   NULL = min(data[[varname]]),
                   range = min(data[[varname]]),
                   unit = survey::svymean(data[[varname]], svydata) - .5,
                   sd = survey::svymean(data[[varname]], svydata) -
                     (sqrt(survey::svyvar(data[[varname]], svydata))/2)
    )
    tmp1 <- switch(diffchange,
                   NULL = max(data[[varname]]),
                   range = max(data[[varname]]),
                   unit = survey::svymean(data[[varname]], svydata) + .5,
                   sd = survey::svymean(data[[varname]], svydata) +
                     (sqrt(survey::svyvar(data[[varname]], svydata))/2)
    )
    delta <- c(tmp0, tmp1)

    # Define variables to vary
    varylist <- as.list(varname)
    varylist <- lapply(1:length(varylist), function(i) {varylist[[i]] <- delta})
    names(varylist) <- varname

    # Create fake dataset for simulations
    trms <- attr(terms(obj), "dataClasses")
    trms <- trms[-length(trms)]   # `polr` class objects have weights in the termlist!
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
    B <- MASS::mvrnorm(sims, c(coef(obj), obj$zeta), vcov(obj))
    tau_ind <- grep("\\|", names(c(coef(obj), obj$zeta)))
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

  #========== Output ===========================================================
  preds <- rename(preds, !!sym(varname) := x)
  diffs <- rename(diffs, !!sym(varname) := x)
  output <- list(
    preds = as_tibble(preds),
    diffs = as_tibble(diffs),
    seed = seed,
    typical = as_tibble(fake))
  class(output) <- "svyEffects"
  attributes(output)$predvar <- varname
  attributes(output)$depvar <- colnames(model.frame(obj))[1]
  return(output)

}
