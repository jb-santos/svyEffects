# METHODS FOR BINARY-DEPENDENT VARIABLE MODELS






#' Average marginal effects (a.k.a. marginal effects at observed values) for
#' binary dependent-variable models of survey-weighted data.
#'
#' @param obj Model object of class \code{survey::svyglm} or \code{glm}.
#' @param varname Character string denoting the name of the predictor variable
#' for which effects are to be calculated.
#' @param weightvar Character string denoting the name of the sampling weight
#' variable.
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
#' @importFrom dplyr arrange as_tibble mutate rename relocate select summarize tibble
#' @importFrom magrittr %>%
#' @importFrom MASS mvrnorm
#' @importFrom rlang sym
#' @importFrom survey svydesign svymean svyquantile svytable
#' @importFrom tidyr pivot_longer pivot_wider
#'
#' @examples
#' mod <- svyglm(votecon ~ agegrp + region + educ + market, design, family=binomial)
#' svyAME(mod, varname="educ", weightvar="weight")
#' svyAME(mod, varname="market", weightvar="weight", diffchange="sd")
#'
#' @export
#'
svyAME.glm <- function(obj,
                       varname,
                       weightvar,
                       nvals = 11,
                       diffchange = c("range", "unit", "sd"),
                       sims = 1500,
                       seed = NULL,
                       ...) {

  #========== SETUP ============================================================

  # Get data
  data <- model.frame(obj)
  data <- data %>% rename(!!sym(weightvar) := `(weights)`)
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

    ## Predicted probabilities -------------------------------------------------

    varname_seq <- seq(min(data[[varname]]), max(data[[varname]]), length=nvals)

    res <- sapply(1:length(varname_seq), function(i) {
      tmp_d <- data
      tmp_d[[varname]] <- varname_seq[i]
      B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))
      X <- model.matrix(formula(obj), tmp_d)
      Xb <- X %*% t(B)
      p <- plogis(Xb)
      m <- apply(p, 2, weighted.mean, w=data[[weightvar]])
      m
    })

    preds <- tibble(
      x = varname_seq,
      predicted = colMeans(res),
      conf.low = apply(unname(res), 2, low),
      conf.high = apply(unname(res), 2, high),
      type = "Probability")

    ## Differences in predicted probabilities ----------------------------------

    diffchange <- match.arg(diffchange)
    tmp0 <- tmp1 <- data
    tmp0[varname] <- switch(diffchange,
                            range = min(data[varname]),
                            unit = survey::svymean(data[varname], svydata) - .5,
                            sd = survey::svymean(data[varname], svydata) -
                              (sqrt(survey::svyvar(data[varname], svydata))/2)
    )
    tmp1[varname] <- switch(diffchange,
                            range = max(data[varname]),
                            unit = survey::svymean(data[varname], svydata) + .5,
                            sd = survey::svymean(data[varname], svydata) +
                              (sqrt(survey::svyvar(data[varname], svydata))/2)
    )

    B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))
    X0 <- model.matrix(formula(obj), tmp0)
    X1 <- model.matrix(formula(obj), tmp1)
    Xb0 <- X0 %*% t(B)
    Xb1 <- X1 %*% t(B)
    p0 <- as.data.frame(plogis(Xb0))
    p1 <- as.data.frame(plogis(Xb1))
    diff <- p1 - p0
    mean_diff <- apply(diff, 2, weighted.mean, w=data[[weightvar]])

    diffs <- tibble(
      x = paste0("Delta (",
                 diffchange,
                 ") : ",
                 round(tmp0[[varname]][1], 3),
                 " - ",
                 round(tmp1[[varname]][1], 3)),
      predicted = mean(mean_diff),
      conf.low = unname(quantile(mean_diff, .025)),
      conf.high = unname(quantile(mean_diff, .975)),
      type = "Difference")

    ## Output ------------------------------------------------------------------

    preds <- dplyr::rename(preds, !!sym(varname) := x)
    diffs <- dplyr::rename(diffs, !!sym(varname) := x)
    output <- list(
      preds = preds,
      diffs = diffs,
      seed = seed)
    return(output)

  } else {

    #========== FACTOR / CATEGORICAL predictors ================================

    ## Predicted probabilities -------------------------------------------------

    levs <- levels(data[[varname]])
    nlev <- length(levs)
    nlevC2 <- factorial(nlev)/(2*factorial(nlev-2))
    DFs <- lapply(1:nlev, function(i) i <- data)
    DFs <- lapply(1:nlev, function(i)
      DFs[[i]] %>% dplyr::mutate(!!sym(varname) := factor(i, levels=1:nlev, labels=levs)))
    X <- lapply(1:nlev, function(i) model.matrix(formula(obj), data=DFs[[i]]))

    B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))
    XB <- lapply(1:nlev, function(i) X[[i]] %*% t(B))
    Pr <- lapply(1:nlev, function(i) plogis(XB[[i]]))
    WgtPr <- lapply(1:nlev, function(i)
      apply(Pr[[i]], 2, function(x)weighted.mean(x, data[[weightvar]])))

    preds <- tibble(
      x = factor(levs, levels=levs),
      predicted = sapply(WgtPr, mean),
      conf.low = sapply(WgtPr, low),
      conf.high = sapply(WgtPr, high),
      type="Predicted probability")

    ## Differences in predicted probabilities ----------------------------------

    Pr_combs <- combn(Pr, 2)
    Pr_diffs <- mapply(x = Pr_combs[2,], y = Pr_combs[1,],
                       function(x, y){unlist(x) - unlist(y)})   # calculate first differences
    Pr_diffs <- relist(Pr_diffs, Pr_combs)
    WgtPr_D <- NULL  # create list of weighted first differences
    for (i in 1:nlevC2) {WgtPr_D[[i]] <-
      apply(Pr_diffs[[i]], 2, function(x)weighted.mean(x, data[[weightvar]]))}

    diffs <- tibble(
      varname1 = combn(levs, 2)[1,] ,
      varname2 = combn(levs, 2)[2,] ,
      predicted = sapply(WgtPr_D, mean),
      conf.low = sapply(WgtPr_D, low),
      conf.high = sapply(WgtPr_D, high))
    diffs <- diffs %>%
      mutate(x = paste(varname2, " - ", varname1, sep=""), .before=1) %>%
      mutate(type = "Difference") %>%
      dplyr::select(-c("varname1", "varname2"))

    ## Output ------------------------------------------------------------------

    preds <- rename(preds, !!sym(varname) := x)
    diffs <- rename(diffs, !!sym(varname) := x)
    output <- list(
      preds = preds,
      diffs = diffs,
      seed = seed)
    return(output)
  }
}






#' Marginal effects at reasonable values for binary dependent-variable models of
#' survey-weighted data.
#'
#' @param obj Model object of class \code{survey::svyglm} or \code{glm}.
#' @param varname Character string denoting the name of the predictor variable
#' for which effects are to be calculated.
#' @param weightvar Character string denoting the name of the sampling weight
#' variable.
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
#' @importFrom dplyr arrange as_tibble mutate rename relocate select summarize tibble
#' @importFrom magrittr %>%
#' @importFrom MASS mvrnorm
#' @importFrom rlang sym
#' @importFrom survey svydesign svymean svyquantile svytable
#' @importFrom tidyr pivot_longer pivot_wider
#'
#' @examples
#' mod <- svyglm(votecon ~ agegrp + region + educ + market, design, family=binomial)
#' svyMER(mod, varname="educ", weightvar="weight")
#' svyMER(mod, varname="market", weightvar="weight", diffchange="sd")
#'
#' @export
#'
svyMER.glm <- function(obj,
                       varname,
                       weightvar,
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

  #========== PREDICTED PROBABILITIES ==========================================

  B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))
  fakeX <- model.matrix(formula(obj), data=fake)
  probs <- t(plogis(fakeX %*% t(B)))

  preds <- probs %>%
    as.data.frame() %>%
    tidyr::pivot_longer(everything(), names_to="ind", values_to="probs") %>%
    dplyr::group_by(ind) %>%
    dplyr::summarize(predicted = mean(probs),
              conf.low = low(unname(probs)),
              conf.high = high(unname(probs)))
  preds$ind <- as.numeric(preds$ind)
  preds <- arrange(preds, ind)
  preds <- dplyr::select(preds, -ind)
  if(class(data[[varname]])=="factor") {
    preds$x <- factor(levels(data[[varname]]), levels=levels(data[[varname]]))}
  else{preds$x <- seq(min(data[[varname]]), max(data[[varname]]), length=nvals)}
  preds <- preds %>%
    relocate(x, .before=1) %>%
    mutate(type = "Probability")

  #========== DIFFERENCES IN PREDICTED PROBABILITIES ===========================

  ## Differences for categorical variables -------------------------------------

  if(class(data[[varname]])=="factor") {
    levs <- levels(data[[varname]])
    nlev <- length(levs)

    P <- NULL
    P <- lapply(1:ncol(probs), function(i){P[[i]] <- probs[,i]})
    P_combs <- combn(P, 2)
    D <- mapply(x = P_combs[2,], y = P_combs[1,],
                function(x, y){unlist(x) - unlist(y)})
    D <- as.data.frame(D)

    diffs <- dplyr::tibble(
      varname1 = combn(levs, 2)[1,],
      varname2 = combn(levs, 2)[2,],
      predicted = sapply(D, mean),
      conf.low = sapply(D, low),
      conf.high = sapply(D, high),
      type="Difference")
    diffs <- diffs %>%
      mutate(x = paste0(diffs$varname2,
                        " - ",
                        diffs$varname1),
             .before=1) %>%
      select(-c(varname1, varname2))

  } else {

    ## Differences for continuous variables ------------------------------------
    diffchange <- match.arg(diffchange)
    diffrange <- switch(diffchange,
                        range = c(min(data[[varname]]), max(data[[varname]])),
                        unit = survey::svymean(data[varname], svydata) + c(-.5,.5),
                        sd = survey::svymean(data[varname], svydata) +
                          (c(-.5,.5)*(sqrt(survey::svyvar(data[varname], svydata))))
    )

    varylist <- as.list(varname)
    varylist <- lapply(1:length(varylist), function(i) {varylist[[i]] <- diffrange})
    names(varylist) <- varname
    # Create dataframe with varying variables and variables fixed at central observations
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

    B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))
    fakeX <- model.matrix(formula(obj), data=fake)
    probs <- t(plogis(fakeX %*% t(B)))
    diff_res <- apply(probs, 1, diff)
    diffs <- dplyr::tibble(
      x = paste0("Delta (", diffchange, ") : ",
                 round(diffrange[[1]], 3),
                 " - ",
                 round(diffrange[[2]], 3)),
      predicted = mean(diff_res),
      conf.low = unname(low(diff_res)),
      conf.high = unname(high(diff_res)),
      type = "Difference")
  }

  #========== Output ===========================================================

  preds <- dplyr::rename(preds, !!sym(varname) := x)
  diffs <- dplyr::rename(diffs, !!sym(varname) := x)
  output <- list(
    preds = preds,
    diffs = diffs,
    seed = seed,
    typical = dplyr::as_tibble(fake))
  return(output)

}
