# METHODS FOR MULTINOMIAL DEPENDENT VARIABLE MODELS






#' Average marginal effects (a.k.a. marginal effects at observed values) for
#' multinomial dependent-variable models of survey-weighted data
#' (svrepmisc::svrepstatmisc).
#'
#'
#' READ THIS FIRST: YOU NEED THE GITHUB PACKAGE \code{svrepmisc} TO USE THIS
#' FUNCTION!This function is designed to be used on a multinomial dependent
#' variable model of class \code{svrepmisc::svymultinom}) estimated on a survey
#' design object of class \code{survey::svrepstatmisc}). Currently, there are no
#' multinomial logit \code{R} functions that work on survey-weighted data.
#' However, Carl Ganz' \code{svrepmisc} package (available on GitHub,
#' https://github.com/carlganz/svrepmisc) uses replicate weights and jackknife
#' standard errors to approximate the same end result. (The output with this
#' approach is nearly identical to results obtained in Stata using either
#' sampling weights (i.e. "pweights") or replicate weights with jackknife
#' standard errors.)
#'
#' To use this function, you use \code{survey::svydesign} to create a survey
#' design object with your sampling/design weights. Then, you create
#' survey design object with replicate weights. Then, you estimate your
#' multinomial logit model on the second survey design object using
#' \code{svrepmisc::svymultinom} (see the examples below).
#'
#'
#' @param obj Model object of class \code{svrepmisc::svrepstatmisc}.
#' @param varname Character string denoting the name of the predictor variable
#' for which effects are to be calculated.
#' @param weightvar A survey design object of class \code{survey::svydesign}.
#' @param design A survey design object with replicate weights of class
#' \code{survey::svyrep.design}.
#' @param modform Character string denoting the model formula used of the model.
#' Must be in the format of class \code{modform = "y ~ x"} and match the model's
#' formula exactly.
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
#' @importFrom stringr boundary str_extract_all
#' @importFrom survey svydesign svymean svyquantile svytable
#' @importFrom tidyr pivot_longer pivot_wider
#'
#' @export
svyAME.svrepstatmisc <- function(obj,
                                 varname,
                                 weightvar,
                                 design,
                                 modform,
                                 nvals = 11,
                                 diffchange = c("range", "unit", "sd"),
                                 sims = 1500,
                                 seed = NULL,
                                 ...) {
  # Object class check
  if(!inherits(obj, "svrepstatmisc")) {
    stop(print("This function work with only survey-weighted multinomial logit model objects of class svrepstat."))
  }

  #========== SETUP ============================================================

  modform <- as.formula(modform)
  data <- design$variables
  varlist <- unlist(stringr::str_extract_all(modform, boundary("word")))
  varlist <- append(varlist, weightvar)
  data <- dplyr::select(data, all_of(varlist)) %>% na.omit()
  svydata <- survey::svydesign(ids=~1, strata=NULL, weights=data[weightvar], data=data)

  Yname <- as.character(modform[[2]])
  Yvar <- dplyr::select(data, all_of(Yname))
  Yvar <- Yvar[[1]]
  Ylevs <- levels(Yvar)
  Ynlev <- as.numeric(nlevels(Yvar))

  if(is.character(data[[varname]])) {
    data[[varname]] <- as.factor(data[[varname]])}
  diffchange <- match.arg(diffchange)


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

    B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))

    varname_seq <- seq(min(data[[varname]]), max(data[[varname]]), length=nvals)
    res <- array(dim=c(sims, Ynlev, length(varname_seq)))
    for(i in seq_along(varname_seq)){
      new <- dplyr::mutate(data, !!sym(varname) := varname_seq[i])
      Xmat <- model.matrix(modform, data=new)
      res[,,i] <- t(apply(B, 1, function(b)
        apply(
          prop.table(
            exp(Xmat %*% cbind(0, matrix(b, ncol=(Ynlev-1), byrow=TRUE))),
            1),
          2,
          weighted.mean,
          w=data[[weightvar]])))
    }

    preds <- data.frame(x = varname_seq)
    m <- as.data.frame(t(apply(res, c(2,3), mean)))
    names(m) <- paste0("predicted_", levels(droplevels(Yvar)))
    l <- as.data.frame(t(apply(res, c(2,3), low)))
    names(l) <- paste0("conf.low_", levels(droplevels(Yvar)))
    u <- as.data.frame(t(apply(res, c(2,3), high)))
    names(u) <- paste0("conf.high_", levels(droplevels(Yvar)))

    preds <- cbind(preds, m, l, u)
    preds <- preds %>%
      pivot_longer(-x,
                   names_pattern="(.*)_(.*)",
                   names_to=c(".value", "y"))
    preds$type <- "Probability"
    preds$y <- factor(preds$y, levels=Ylevs)
    preds <- preds %>% relocate(y, .before=1)

    ## Differences in predicted probabilities ----------------------------------

    tmp0 <- switch(diffchange,
                   range = min(data[varname]),
                   unit = survey::svymean(data[varname], svydata) - .5,
                   sd = survey::svymean(data[varname], svydata) -
                     (sqrt(survey::svyvar(data[varname], svydata))/2)
    )
    tmp1 <- switch(diffchange,
                   range = max(data[varname]),
                   unit = survey::svymean(data[varname], svydata) + .5,
                   sd = survey::svymean(data[varname], svydata) +
                     (sqrt(survey::svyvar(data[varname], svydata))/2)
    )
    diff_seq <- c(tmp0, tmp1)

    DFs <- lapply(1:2, function(x) {x <- data})
    DFs <- lapply(1:2, function(i) {
      DFs[[i]] <- DFs[[i]] %>% dplyr::mutate(!!sym(varname) := diff_seq[[i]])
    })
    X <- NULL
    X <- lapply(1:2, function(i) {
      X[[i]] <- model.matrix(modform, data=DFs[[i]])
    })

    B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))

    res <- NULL
    for(i in 1:nrow(B)){
      tmpB <- cbind(0, matrix(B[i,], ncol=(Ynlev-1), byrow=TRUE))
      EXB <- lapply(1:2, function(x){exp(X[[x]] %*% tmpB)})
      P <- lapply(1:2, function(x){prop.table(EXB[[x]], 1)})
      res_temp <- sapply(as.data.frame(P), weighted.mean, w=data[[weightvar]])
      res_temp <- t(as.data.frame(res_temp))
      rownames(res_temp) <- NULL
      res <- rbind(res, res_temp)
    }

    P <- NULL
    P <- lapply(1:2, function(i) {
      P[[i]] <- res[,i:(i+Ynlev-1)]
    })

    P_combs <- combn(P, 2)
    Pr_diffs <- as.data.frame(P_combs[2,]) - as.data.frame(P_combs[1,])

    diffs <- expand.grid(
      y = levels(droplevels(Yvar)),
      x1 = 1:2,
      x2 = 1:2) %>%
      dplyr::filter(x1 < x2) %>%
      dplyr::mutate(x = paste0("Delta (", diffchange, ") : ",
                        round(diff_seq[x1], 3),
                        " - ",
                        round(diff_seq[x2], 3))) %>%
      dplyr::select(-c("x1", "x2"))
    diffs <- diffs %>%
      dplyr::mutate(
        predicted = colMeans(Pr_diffs),
        conf.low = apply(Pr_diffs, 2, low),
        conf.high = apply(Pr_diffs, 2, high),
        type = "Difference")
    diffs$y <- factor(diffs$y, levels=Ylevs)
    diffs <- diffs %>% relocate(y, .before=1)

    ## Output ------------------------------------------------------------------

    preds <- rename(preds, !!sym(varname) := x)
    diffs <- rename(diffs, !!sym(varname) := x)
    output <- list(
      preds = dplyr::as_tibble(preds),
      diffs = dplyr::as_tibble(diffs))
    class(output) <- "svyEffects"
    attributes(output)$predvar <- varname
    attributes(output)$depvar <- Yname
    return(output)

  } else {

    #========== FACTOR / CATEGORICAL predictors ================================

    ## Predicted probabilities -------------------------------------------------

    levs <- levels(data[[varname]])
    nlev <- length(levs)
    nlevC2 <- (factorial(nlev)/(2*factorial(nlev-2)))

    # Create dataframes for simulations
    DFs <- lapply(1:nlev, function(i) {i <- data})
    DFs <- lapply(1:nlev, function(i) {
      DFs[[i]] <- DFs[[i]] %>% dplyr::mutate(!!sym(varname) := factor(i, levels=1:nlev, labels=levs))
    })
    X <- NULL
    X <- lapply(1:nlev, function(i) {
      X[[i]] <- model.matrix(modform, data=DFs[[i]])
    })

    # Create coefficient vectors
    B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))

    # Generate predictions
    res <- NULL
    res <- lapply(1:nrow(B), function(i) {
      tmpB <- cbind(0, matrix(B[i,], ncol=(Ynlev-1), byrow=TRUE))
      EXB <- lapply(1:nlev, function(j){exp(X[[j]] %*% tmpB)})
      P <- lapply(1:nlev, function(k){prop.table(EXB[[k]], 1)})
      res_temp <- sapply(as.data.frame(P), weighted.mean, w=data[[weightvar]])
      res_temp <- t(as.data.frame(res_temp))
      rownames(res_temp) <- NULL
      res[[i]] <- res_temp
    })
    res <- t(sapply(res, function(x)unlist(x)))

    # Assemble table of predicted probabilities
    preds <- expand.grid(
      y = levels(droplevels(Yvar)),
      x = levs
    )
    preds$predicted <- colMeans(res)
    preds$conf.low <- apply(res, 2, low)
    preds$conf.high <- apply(res, 2, high)
    preds$type <- "Probability"
    preds$y <- factor(preds$y, levels=Ylevs)
    preds$x <- factor(preds$x, levels=levs)
    preds <- preds %>% relocate(y, .before=1)

    ## Differences in predicted probabilities ----------------------------------

    i <- 1
    j <- 1
    P <- NULL
    while (i <= nlev) {
      P[[i]] <- res[,j:(j+Ynlev-1)]
      i <- i + 1
      j <- j + Ynlev
    }

    P_combs <- combn(P, 2, simplify=FALSE)
    D <- lapply(1:length(P_combs), function(i)P_combs[[i]][[2]] - P_combs[[i]][[1]])
    diff_res <- cbind.data.frame(D)

    # Assemble table of differences

    diffs <- expand.grid(
      y = levels(droplevels(Yvar)),
      x1 = 1:nlev,
      x2 = 1:nlev) %>%
      dplyr::filter(x1 < x2) %>%
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
    diffs$y <- factor(diffs$y, levels=Ylevs)
    diffs <- diffs %>% relocate(y, .before=1)

    ## Output ------------------------------------------------------------------

    preds <- rename(preds, !!sym(varname) := x)
    diffs <- rename(diffs, !!sym(varname) := x)
    output <- list(
      preds = as_tibble(preds),
      diffs = as_tibble(diffs))
    class(output) <- "svyEffects"
    attributes(output)$predvar <- varname
    attributes(output)$depvar <- Yname
    return(output)

  }
}







#' Marginal effects at reasonable values for multinomial dependent-variable
#' models of survey-weighted data (svrepmisc::svrepstatmisc).
#'
#'
#' READ THIS FIRST: YOU NEED THE GITHUB PACKAGE \code{svrepmisc} TO USE THIS
#' FUNCTION!This function is designed to be used on a multinomial dependent
#' variable model of class \code{svrepmisc::svymultinom}) estimated on a survey
#' design object of class \code{survey::svrepstatmisc}). Currently, there are no
#' multinomial logit \code{R} functions that work on survey-weighted data.
#' However, Carl Ganz' \code{svrepmisc} package (available on GitHub,
#' https://github.com/carlganz/svrepmisc) uses replicate weights and jackknife
#' standard errors to approximate the same end result. (The output with this
#' approach is nearly identical to results obtained in Stata using either
#' sampling weights (i.e. "pweights") or replicate weights with jackknife
#' standard errors.)
#'
#' To use this function, you use \code{survey::svydesign} to create a survey
#' design object with your sampling/design weights. Then, you create
#' survey design object with replicate weights. Then, you estimate your
#' multinomial logit model on the second survey design object using
#' \code{svrepmisc::svymultinom} (see the examples below).
#'
#'
#' @param obj Model object of class \code{svrepmisc::svrepstatmisc}.
#' @param varname Character string denoting the name of the predictor variable
#' for which effects are to be calculated.
#' @param weightvar A survey design object of class \code{survey::svydesign}.
#' @param design A survey design object with replicate weights of class
#' \code{survey::svyrep.design}.
#' @param modform Character string denoting the model formula used of the model.
#' Must be in the format of class \code{modform = "y ~ x"} and match the model's
#' formula exactly.
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
#'     \item{\code{$typical}}{the values at which variables are set for the simulations}
#'     }
#'
#' @importFrom dplyr arrange as_tibble filter mutate rename relocate select summarize tibble
#' @importFrom magrittr %>%
#' @importFrom MASS mvrnorm
#' @importFrom rlang sym
#' @importFrom stringr boundary str_extract_all
#' @importFrom survey svydesign svymean svyquantile svytable
#' @importFrom tidyr pivot_longer pivot_wider
#'
#' @export
svyMER.svrepstatmisc <- function(obj,
                                 varname,
                                 weightvar,
                                 design,
                                 modform,
                                 nvals = 11,
                                 diffchange = c("range", "unit", "sd"),
                                 sims = 1500,
                                 seed = NULL,
                                 ...) {
  # Object class check
  if(!inherits(obj, "svrepstatmisc")){
    stop(print("This function work with only survey-weighted multinomial logit model objects of class svrepstat."))
  }

  #========== SETUP ============================================================
  modform <- as.formula(modform)
  data <- design$variables
  varlist <- unlist(stringr::str_extract_all(modform, boundary("word")))
  varlist <- append(varlist, weightvar)
  data <- dplyr::select(data, all_of(varlist)) %>% na.omit()
  svydata <- survey::svydesign(ids=~1, strata=NULL, weights=data[weightvar], data=data)

  Yname <- as.character(modform[[2]])
  Yvar <- dplyr::select(data, all_of(Yname))
  Yvar <- Yvar[[1]]
  Ylevs <- levels(Yvar)
  Ynlev <- as.numeric(nlevels(Yvar))

  if(class(data[[varname]])=="factor") {
    Xlevs <- levels(data[[varname]])}
  if(class(data[[varname]])=="numeric") {
    Xlevs <- seq(min(data[[varname]]), max(data[[varname]]), length=nvals)}

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
  trms <- attr(terms(model.frame(as.formula(modform), data)), "dataClasses")
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

  #========== PREDICTED PROBABILITIES ==========================================

  # Parametric bootstrap
  B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))
  fakeX <- model.matrix(formula(as.formula(modform)), data=fake)

  # Generate predictions
  res <- array(dim=c(length(Xlevs),Ynlev,sims))
  for(i in 1:nrow(B)){
    tmpB <- cbind(0, matrix(B[i,], ncol=(Ynlev-1), byrow=TRUE))
    exb <- exp(fakeX %*% tmpB)
    p <- prop.table(exb, 1)
    res[,,i] <- p
  }

  # Assemble table of predictions
  preds <- data.frame(x = Xlevs)
  if(class(data[[varname]])=="factor") {preds$x <- factor(preds$x, levels=Xlevs)}
  res_m <- as.data.frame(apply(res, c(1,2), mean))
  names(res_m) <- paste0("predicted_", Ylevs)
  res_l <- as.data.frame(apply(res, c(1,2), low))
  names(res_l) <- paste0("conf.low_", Ylevs)
  res_u <- as.data.frame(apply(res, c(1,2), high))
  names(res_u) <- paste0("conf.high_", Ylevs)

  preds <- cbind(preds, res_m, res_l, res_u)
  preds <- preds %>%
    tidyr::pivot_longer(-x,
                 names_pattern="(.*)_(.*)",
                 names_to=c(".value", "y")) %>%
    dplyr::mutate(type = "Probability") %>%
    dplyr::relocate(y, .before=1)
  preds$y <- factor(preds$y, levels=Ylevs)


  #========== DIFFERENCES IN PREDICTED PROBABILITIES ===========================

  ## Differences for categorical variables -------------------------------------

  if(is.factor(data[[varname]])) {
    P <- lapply(1:length(Xlevs), function(i)(t(res[i,,])))
    P_combs <- combn(P, 2, simplify=FALSE)
    D <- lapply(1:length(P_combs), function(i){P_combs[[i]][[2]] - P_combs[[i]][[1]]})
    diff_res <- cbind.data.frame(D)

    diffs <- expand.grid(
      y = factor(Ylevs, levels=Ylevs),
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

    # Differences for continuous variables ------------------------------------

    tmp0 <- switch(diffchange,
                   range = min(data[[varname]]),
                   unit = survey::svymean(data[[varname]], svydata) - .5,
                   sd = survey::svymean(data[[varname]], svydata) -
                     (sqrt(survey::svyvar(data[[varname]], svydata))/2)
    )
    tmp1 <- switch(diffchange,
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
    trms <- attr(terms(model.frame(as.formula(modform), data)), "dataClasses")
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
    fakeX <- model.matrix(formula(as.formula(modform)), data=fake)

    # Generate predictions
    res <- array(dim=c(2,Ynlev,sims))
    for(i in 1:nrow(B)){
      tmpB <- cbind(0, matrix(B[i,], ncol=(Ynlev-1), byrow=TRUE))
      exb <- exp(fakeX %*% tmpB)
      p <- prop.table(exb, 1)
      res[,,i] <- p
    }

    # Calculate differences
    diff_res <- res[2,,] - res[1,,]

    # Assemble table of differences
    diffs <- data.frame(
      y = factor(Ylevs, levels=Ylevs),
      x = paste0("Delta (", diffchange, ") : ",
                 round(tmp0, 3), " - ", round(tmp1, 3)),
      predicted = apply(diff_res, c(1), mean),
      conf.low = apply(diff_res, c(1), low),
      conf.high = apply(diff_res, c(1), high),
      type = "Difference")
  }

  #========== Output ===========================================================

  preds <- rename(preds, !!sym(varname) := x)
  diffs <- rename(diffs, !!sym(varname) := x)
  output <- list(
    preds = dplyr::as_tibble(preds),
    diffs = dplyr::as_tibble(diffs),
    typical = dplyr::as_tibble(fake)
    )
  class(output) <- "svyEffects"
  attributes(output)$predvar <- varname
  attributes(output)$depvar <- Yname
  return(output)

}






#' Average marginal effects (a.k.a. marginal effects at observed values) for
#' (non-survey) weighted multinomial dependent-variable models (nnet::multinom).
#'
#'
#' @param obj Model object of class \code{svrepmisc::svrepstatmisc}.
#' @param varname Character string denoting the name of the predictor variable
#' for which effects are to be calculated.
#' @param weightvar A survey design object of class \code{survey::svydesign}.
#' @param design A survey design object with replicate weights of class
#' \code{survey::svyrep.design}.
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
#' @importFrom stringr boundary str_extract_all
#' @importFrom survey svydesign svymean svyquantile svytable
#' @importFrom tidyr pivot_longer pivot_wider
#'
#' @export
svyAME.multinom <- function(obj,
                            varname,
                            weightvar,
                            design,
                            nvals = 11,
                            diffchange = c("range", "unit", "sd"),
                            sims = 1500,
                            seed = NULL,
                            ...) {
  # Object class check
  if(!inherits(obj, "multinom")){
    stop(print("This function work with only multinomial logit model objects of class multinom."))
  }

  #========== SETUP ============================================================

  modform <- as.formula(modform)
  data <- design$variables
  varlist <- unlist(stringr::str_extract_all(modform, boundary("word")))
  varlist <- append(varlist, weightvar)
  data <- dplyr::select(data, all_of(varlist)) %>% na.omit()
  svydata <- survey::svydesign(ids=~1, strata=NULL, weights=data[weightvar], data=data)

  Yname <- as.character(modform[[2]])
  Yvar <- dplyr::select(data, all_of(Yname))
  Yvar <- Yvar[[1]]
  Ylevs <- levels(Yvar)
  Ynlev <- as.numeric(nlevels(Yvar))

  if(is.character(data[[varname]])) {
    data[[varname]] <- as.factor(data[[varname]])}
  diffchange <- match.arg(diffchange)


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

    B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))

    varname_seq <- seq(min(data[[varname]]), max(data[[varname]]), length=nvals)
    res <- array(dim=c(sims, Ynlev, length(varname_seq)))
    for(i in seq_along(varname_seq)){
      new <- dplyr::mutate(data, !!sym(varname) := varname_seq[i])
      Xmat <- model.matrix(modform, data=new)
      res[,,i] <- t(apply(B, 1, function(b)
        apply(
          prop.table(
            exp(Xmat %*% cbind(0, matrix(b, ncol=(Ynlev-1), byrow=TRUE))),
            1),
          2,
          weighted.mean,
          w=data[[weightvar]])))
    }

    preds <- data.frame(x = varname_seq)
    m <- as.data.frame(t(apply(res, c(2,3), mean)))
    names(m) <- paste0("predicted_", levels(droplevels(Yvar)))
    l <- as.data.frame(t(apply(res, c(2,3), low)))
    names(l) <- paste0("conf.low_", levels(droplevels(Yvar)))
    u <- as.data.frame(t(apply(res, c(2,3), high)))
    names(u) <- paste0("conf.high_", levels(droplevels(Yvar)))

    preds <- cbind(preds, m, l, u)
    preds <- preds %>%
      pivot_longer(-x,
                   names_pattern="(.*)_(.*)",
                   names_to=c(".value", "y"))
    preds$type <- "Probability"
    preds$y <- factor(preds$y, levels=Ylevs)
    preds <- preds %>% relocate(y, .before=1)

    ## Differences in predicted probabilities ----------------------------------

    tmp0 <- switch(diffchange,
                   range = min(data[varname]),
                   unit = survey::svymean(data[varname], svydata) - .5,
                   sd = survey::svymean(data[varname], svydata) -
                     (sqrt(survey::svyvar(data[varname], svydata))/2)
    )
    tmp1 <- switch(diffchange,
                   range = max(data[varname]),
                   unit = survey::svymean(data[varname], svydata) + .5,
                   sd = survey::svymean(data[varname], svydata) +
                     (sqrt(survey::svyvar(data[varname], svydata))/2)
    )
    diff_seq <- c(tmp0, tmp1)

    DFs <- lapply(1:2, function(x) {x <- data})
    DFs <- lapply(1:2, function(i) {
      DFs[[i]] <- DFs[[i]] %>% dplyr::mutate(!!sym(varname) := diff_seq[[i]])
    })
    X <- NULL
    X <- lapply(1:2, function(i) {
      X[[i]] <- model.matrix(modform, data=DFs[[i]])
    })

    B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))

    res <- NULL
    for(i in 1:nrow(B)){
      tmpB <- cbind(0, matrix(B[i,], ncol=(Ynlev-1), byrow=TRUE))
      EXB <- lapply(1:2, function(x){exp(X[[x]] %*% tmpB)})
      P <- lapply(1:2, function(x){prop.table(EXB[[x]], 1)})
      res_temp <- sapply(as.data.frame(P), weighted.mean, w=data[[weightvar]])
      res_temp <- t(as.data.frame(res_temp))
      rownames(res_temp) <- NULL
      res <- rbind(res, res_temp)
    }

    P <- NULL
    P <- lapply(1:2, function(i) {
      P[[i]] <- res[,i:(i+Ynlev-1)]
    })

    P_combs <- combn(P, 2)
    Pr_diffs <- as.data.frame(P_combs[2,]) - as.data.frame(P_combs[1,])

    diffs <- expand.grid(
      y = levels(droplevels(Yvar)),
      x1 = 1:2,
      x2 = 1:2) %>%
      dplyr::filter(x1 < x2) %>%
      dplyr::mutate(x = paste0("Delta (", diffchange, ") : ",
                               round(diff_seq[x1], 3),
                               " - ",
                               round(diff_seq[x2], 3))) %>%
      dplyr::select(-c("x1", "x2"))
    diffs <- diffs %>%
      dplyr::mutate(
        predicted = colMeans(Pr_diffs),
        conf.low = apply(Pr_diffs, 2, low),
        conf.high = apply(Pr_diffs, 2, high),
        type = "Difference")
    diffs$y <- factor(diffs$y, levels=Ylevs)
    diffs <- diffs %>% relocate(y, .before=1)

    ## Output ------------------------------------------------------------------

    preds <- rename(preds, !!sym(varname) := x)
    diffs <- rename(diffs, !!sym(varname) := x)
    output <- list(
      preds = dplyr::as_tibble(preds),
      diffs = dplyr::as_tibble(diffs))
    output

  } else {

    #========== FACTOR / CATEGORICAL predictors ================================

    ## Predicted probabilities -------------------------------------------------

    levs <- levels(data[[varname]])
    nlev <- length(levs)
    nlevC2 <- (factorial(nlev)/(2*factorial(nlev-2)))

    # Create dataframes for simulations
    DFs <- lapply(1:nlev, function(i) {i <- data})
    DFs <- lapply(1:nlev, function(i) {
      DFs[[i]] <- DFs[[i]] %>% dplyr::mutate(!!sym(varname) := factor(i, levels=1:nlev, labels=levs))
    })
    X <- NULL
    X <- lapply(1:nlev, function(i) {
      X[[i]] <- model.matrix(modform, data=DFs[[i]])
    })

    # Create coefficient vectors
    B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))

    # Generate predictions
    res <- NULL
    res <- lapply(1:nrow(B), function(i) {
      tmpB <- cbind(0, matrix(B[i,], ncol=(Ynlev-1), byrow=TRUE))
      EXB <- lapply(1:nlev, function(j){exp(X[[j]] %*% tmpB)})
      P <- lapply(1:nlev, function(k){prop.table(EXB[[k]], 1)})
      res_temp <- sapply(as.data.frame(P), weighted.mean, w=data[[weightvar]])
      res_temp <- t(as.data.frame(res_temp))
      rownames(res_temp) <- NULL
      res[[i]] <- res_temp
    })
    res <- t(sapply(res, function(x)unlist(x)))

    # Assemble table of predicted probabilities
    preds <- expand.grid(
      y = levels(droplevels(Yvar)),
      x = levs
    )
    preds$predicted <- colMeans(res)
    preds$conf.low <- apply(res, 2, low)
    preds$conf.high <- apply(res, 2, high)
    preds$type <- "Probability"
    preds$y <- factor(preds$y, levels=Ylevs)
    preds$x <- factor(preds$x, levels=levs)
    preds <- preds %>% relocate(y, .before=1)

    ## Differences in predicted probabilities ----------------------------------

    i <- 1
    j <- 1
    P <- NULL
    while (i <= nlev) {
      P[[i]] <- res[,j:(j+Ynlev-1)]
      i <- i + 1
      j <- j + Ynlev
    }

    P_combs <- combn(P, 2, simplify=FALSE)
    D <- lapply(1:length(P_combs), function(i)P_combs[[i]][[2]] - P_combs[[i]][[1]])
    diff_res <- cbind.data.frame(D)

    # Assemble table of differences

    diffs <- expand.grid(
      y = levels(droplevels(Yvar)),
      x1 = 1:nlev,
      x2 = 1:nlev) %>%
      dplyr::filter(x1 < x2) %>%
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
    diffs$y <- factor(diffs$y, levels=Ylevs)
    diffs <- diffs %>% relocate(y, .before=1)

    ## Output ------------------------------------------------------------------

    preds <- rename(preds, !!sym(varname) := x)
    diffs <- rename(diffs, !!sym(varname) := x)
    output <- list(
      preds = as_tibble(preds),
      diffs = as_tibble(diffs))
    class(output) <- "svyEffects"
    attributes(output)$predvar <- varname
    attributes(output)$depvar <- Yname
    return(output)

  }
}







#' Marginal effects at reasonable values for multinomial dependent-variable
#' models of (non-survey) weighted data (nnet::multinom).
#'
#'
#' @param obj Model object of class \code{nnet::multinom}.
#' @param varname Character string denoting the name of the predictor variable
#' for which effects are to be calculated.
#' @param weightvar A survey design object of class \code{survey::svydesign}.
#' @param design A survey design object with replicate weights of class
#' \code{survey::svyrep.design}.
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
#'     \item{\code{$typical}}{the values at which variables are set for the simulations}
#'     }
#'
#' @importFrom dplyr arrange as_tibble filter mutate rename relocate select summarize tibble
#' @importFrom magrittr %>%
#' @importFrom MASS mvrnorm
#' @importFrom rlang sym
#' @importFrom stringr boundary str_extract_all
#' @importFrom survey svydesign svymean svyquantile svytable
#' @importFrom tidyr pivot_longer pivot_wider
#'
#' @export
svyMER.multinom <- function(obj,
                            varname,
                            weightvar,
                            design,
                            nvals = 11,
                            diffchange = c("range", "unit", "sd"),
                            sims = 1500,
                            seed = NULL,
                            ...) {
  # Object class check
  if(!inherits(obj, "multinom")){
    stop(print("This function work with only multinomial logit model objects of class multinom."))
  }

  #========== SETUP ============================================================
  data <- model.frame(obj)
  data[weightvar] <- as.numeric(vote$weights)
  svydata <- survey::svydesign(ids=~1, strata=NULL, weights=data[weightvar], data=data)

  Yvar <- data[[1]]
  Ylevs <- levels(Yvar)
  Ynlev <- as.numeric(nlevels(Yvar))

  if(class(data[[varname]])=="factor") {
    Xlevs <- levels(data[[varname]])}
  if(class(data[[varname]])=="numeric") {
    Xlevs <- seq(min(data[[varname]]), max(data[[varname]]), length=nvals)}

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
  trms <- trms[-length(trms)]   # `multinom` class objects have weights in the termlist!
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
  B <- MASS::mvrnorm(sims, as.numeric(coef(obj)), vcov(obj))
  fakeX <- model.matrix(formula(obj), data=fake)

  # Generate predictions
  res <- array(dim=c(length(Xlevs),Ynlev,sims))
  for(i in 1:nrow(B)){
    tmpB <- cbind(0, matrix(B[i,], ncol=(Ynlev-1), byrow=TRUE))
    exb <- exp(fakeX %*% tmpB)
    p <- prop.table(exb, 1)
    res[,,i] <- p
  }

  # Assemble table of predictions
  preds <- data.frame(x = Xlevs)
  if(class(data[[varname]])=="factor") {preds$x <- factor(preds$x, levels=Xlevs)}
  res_m <- as.data.frame(apply(res, c(1,2), mean))
  names(res_m) <- paste0("predicted_", Ylevs)
  res_l <- as.data.frame(apply(res, c(1,2), low))
  names(res_l) <- paste0("conf.low_", Ylevs)
  res_u <- as.data.frame(apply(res, c(1,2), high))
  names(res_u) <- paste0("conf.high_", Ylevs)

  preds <- cbind(preds, res_m, res_l, res_u)
  preds <- preds %>%
    tidyr::pivot_longer(-x,
                        names_pattern="(.*)_(.*)",
                        names_to=c(".value", "y")) %>%
    dplyr::mutate(type = "Probability") %>%
    dplyr::relocate(y, .before=1)
  preds$y <- factor(preds$y, levels=Ylevs)


  #========== DIFFERENCES IN PREDICTED PROBABILITIES ===========================

  ## Differences for categorical variables -------------------------------------

  if(is.factor(data[[varname]])) {
    P <- lapply(1:length(Xlevs), function(i)(t(res[i,,])))
    P_combs <- combn(P, 2, simplify=FALSE)
    D <- lapply(1:length(P_combs), function(i){P_combs[[i]][[2]] - P_combs[[i]][[1]]})
    diff_res <- cbind.data.frame(D)

    diffs <- expand.grid(
      y = factor(Ylevs, levels=Ylevs),
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

    # Differences for continuous variables ------------------------------------

    tmp0 <- switch(diffchange,
                   range = min(data[[varname]]),
                   unit = survey::svymean(data[[varname]], svydata) - .5,
                   sd = survey::svymean(data[[varname]], svydata) -
                     (sqrt(survey::svyvar(data[[varname]], svydata))/2)
    )
    tmp1 <- switch(diffchange,
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
    trms <- trms[-length(trms)]   # `multinom` class objects have weights in the termlist!
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
    B <- MASS::mvrnorm(sims, as.numeric(coef(obj)), vcov(obj))
    fakeX <- model.matrix(formula(obj), data=fake)

    # Generate predictions
    res <- array(dim=c(2,Ynlev,sims))
    for(i in 1:nrow(B)){
      tmpB <- cbind(0, matrix(B[i,], ncol=(Ynlev-1), byrow=TRUE))
      exb <- exp(fakeX %*% tmpB)
      p <- prop.table(exb, 1)
      res[,,i] <- p
    }

    # Calculate differences
    diff_res <- res[2,,] - res[1,,]

    # Assemble table of differences
    diffs <- data.frame(
      y = factor(Ylevs, levels=Ylevs),
      x = paste0("Delta (", diffchange, ") : ",
                 round(tmp0, 3), " - ", round(tmp1, 3)),
      predicted = apply(diff_res, c(1), mean),
      conf.low = apply(diff_res, c(1), low),
      conf.high = apply(diff_res, c(1), high),
      type = "Difference")
  }

  #========== Output ===========================================================

  preds <- rename(preds, !!sym(varname) := x)
  diffs <- rename(diffs, !!sym(varname) := x)
  output <- list(
    preds = dplyr::as_tibble(preds),
    diffs = dplyr::as_tibble(diffs),
    typical = dplyr::as_tibble(fake))
  class(output) <- "svyEffects"
  attributes(output)$predvar <- varname
  attributes(output)$depvar <- colnames(obj$model)[1]
  return(output)

}
