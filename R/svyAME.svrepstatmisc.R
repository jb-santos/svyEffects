#' Average Marginal Effects For Multinomial Logit Models Of Survey-Weighted Data
#'
#'
#' @description Calculates predicted probabilities and differences in probabilities
#' for a predictor holding all other observations at observed values
#' (i.e. a true average marginal effect). Uses simulations (the parametric
#' bootstrap) to derive 95% confidence intervals.
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
#' \dontrun{
#' data(ces19)
#' library(survey)
#' ces19_svy <- svydesign(ids = ~1, strata = NULL, weights = ~pesweight,
#'   data = ces19, digits = 3)
#' ces19_svy_r <- as.svrepdesign(ces19_svy, type = "JK1")
#' # remotes::install_github("carlganz/svrepmisc") # (if not already installed)
#' library(svrepmisc)
#' VOTE <- svymultinom(vote ~ agegrp + gender + educ + region + marketlib,
#'   design = ces19_svy_r, trace = FALSE)
#' svyAME(VOTE, varname = "region", weightvar = "pesweight", seed = 2019,
#'   design = ces19_svy_r,
#'   modform = "vote ~ agegrp + gender + educ + region + marketlib")
#' svyAME(VOTE, varname = "marketlib", weightvar = "pesweight", seed = 2019,
#'   design = ces19_svy_r,
#'   modform = "vote ~ agegrp + gender + educ + region + marketlib")
#' }
#'
svyAME.svrepstatmisc <- function(obj,
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

  # SETUP ======================================================================

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

  if(is.character(data[[varname]])) {
    data[[varname]] <- as.factor(data[[varname]])
  }

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

    # NUMERIC / CONTINUOUS predictors ============================================

    if(is.numeric(data[[varname]])) {

      ## Predicted probabilities ------------------------------------------------

      B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))

      varname_seq <- seq(min(data[[varname]]), max(data[[varname]]), length = nvals)
      res <- array(dim = c(sims, Ynlev, length(varname_seq)))
      for(i in seq_along(varname_seq)){
        new <- dplyr::mutate(data, !!sym(varname) := varname_seq[i])
        Xmat <- model.matrix(as.formula(modform), data = new)
        res[,,i] <- t(apply(B, 1, function(b)
          apply(
            prop.table(
              exp(Xmat %*% cbind(0, matrix(b, ncol = (Ynlev-1), byrow = TRUE))),
              1),
            2,
            weighted.mean,
            w = data$`(weights)`)))
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
      preds <- preds %>% relocate(y, .before = 1)

      # res_pp <- res   # save for building feature to include yhats in output

      ## Differences in predicted probabilities ----------------------------------

      diffchange <- match.arg(diffchange)

      tmp0 <- switch(diffchange,
                     NULL = survey::svymean(data[varname], svydata) -
                       (sqrt(survey::svyvar(data[varname], svydata))/2),
                     sd = survey::svymean(data[varname], svydata) -
                       (sqrt(survey::svyvar(data[varname], svydata))/2),
                     range = min(data[varname]),
                     unit = survey::svymean(data[varname], svydata) - .5
                     )
      tmp1 <- switch(diffchange,
                     NULL = survey::svymean(data[varname], svydata) +
                       (sqrt(survey::svyvar(data[varname], svydata)) / 2),
                     sd = survey::svymean(data[varname], svydata) +
                       (sqrt(survey::svyvar(data[varname], svydata)) / 2),
                     range = max(data[varname]),
                     unit = survey::svymean(data[varname], svydata) + .5
                     )
      diff_seq <- c(tmp0, tmp1)

      DFs <- lapply(1:2, function(x) x <- data)
      DFs <- lapply(1:2, function(i)
        DFs[[i]] <- DFs[[i]] %>%
          dplyr::mutate(!!sym(varname) := diff_seq[[i]]))
      X <- NULL
      X <- lapply(1:2, function(i)
        X[[i]] <- model.matrix(as.formula(modform), data = DFs[[i]]))

      B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))

      res <- NULL
      for(i in 1:nrow(B)){
        tmpB <- cbind(0, matrix(B[i,], ncol=(Ynlev-1), byrow=TRUE))
        EXB <- lapply(1:2, function(x){exp(X[[x]] %*% tmpB)})
        P <- lapply(1:2, function(x){prop.table(EXB[[x]], 1)})
        res_temp <- sapply(as.data.frame(P), weighted.mean, w = data$`(weights)`)
        res_temp <- t(as.data.frame(res_temp))
        rownames(res_temp) <- NULL
        res <- rbind(res, res_temp)
      }

      Pr <- as.list(as.data.frame(res))

      Pr_diffs <- NULL
      i <- Ynlev
      j <- 1
      while(i < Ynlev * 2) {
        Pr_diffs[[j]] <- Pr[[i+1]] - Pr[[j]]
        i <- i + 1
        j <- j + 1
      }

      Pr_diffs <- do.call("rbind", Pr_diffs) %>%
        t() %>%
        as.data.frame()
      names(Pr_diffs) <- Ylevs

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
      diffs$y <- factor(diffs$y, levels = Ylevs)
      diffs <- diffs %>% relocate(y, .before = 1)


      ## Output ------------------------------------------------------------------

      preds <- rename(preds, !!sym(varname) := x)
      diffs <- rename(diffs, !!sym(varname) := x)
      output <- list(
        preds = tibble::as_tibble(preds),
        diffs = tibble::as_tibble(diffs),
        seed = seed,
        sims = sims,
        formula = as.formula(modform)
        # res_pp = res_pp,   # outputs array of PPs
        # res_diffs = res,   # outputs data frame of PPs used to calculate diffs
        # diffs_sims = Pc_diffs   # calculated first differences
        )
      class(output) <- "svyEffects"
      attributes(output)$predvar <- varname
      attributes(output)$depvar <- Yname
      attributes(output)$method <- "AME"
      return(output)

    } else {

      # FACTOR / CATEGORICAL predictors ==========================================

      ## Predicted probabilities -------------------------------------------------

      levs <- levels(data[[varname]])
      nlev <- length(levs)
      nlevC2 <- (factorial(nlev)/(2*factorial(nlev-2)))

      # Create dataframes for simulations
      DFs <- lapply(1:nlev, function(i) i <- data)
      DFs <- lapply(1:nlev, function(i)
        DFs[[i]] <- DFs[[i]] %>%
          dplyr::mutate(!!sym(varname) := factor(i, levels = 1:nlev, labels = levs)))
      X <- NULL
      X <- lapply(1:nlev, function(i)
        X[[i]] <- model.matrix(as.formula(modform), data = DFs[[i]]))

      # Create coefficient vectors
      B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))

      # Generate predictions
      res <- NULL
      res <- lapply(1:nrow(B), function(i) {
        tmpB <- cbind(0, matrix(B[i,],
                                ncol = (Ynlev-1),
                                byrow = TRUE))
        EXB <- lapply(1:nlev, function(j) {exp(X[[j]] %*% tmpB)})
        P <- lapply(1:nlev, function(k) {prop.table(EXB[[k]], 1)})
        res_temp <- sapply(as.data.frame(P), weighted.mean, w = data$`(weights)`)
        res_temp <- t(as.data.frame(res_temp))
        rownames(res_temp) <- NULL
        res[[i]] <- res_temp
      })
      res <- t(sapply(res, function(x) unlist(x)))

      # Assemble table of predicted probabilities
      preds <- expand.grid(
        y = levels(droplevels(Yvar)),
        x = levs)
      preds$predicted <- colMeans(res)
      preds$conf.low <- apply(res, 2, low)
      preds$conf.high <- apply(res, 2, high)
      preds$type <- "Probability"
      preds$y <- factor(preds$y, levels = Ylevs)
      preds$x <- factor(preds$x, levels = levs)
      preds <- preds %>% relocate(y, .before = 1)

      ## Differences in predicted probabilities --------------------------------

      i <- 1
      j <- 1
      P <- NULL
      while (i <= nlev) {
        P[[i]] <- res[, j:(j + Ynlev-1)]
        i <- i + 1
        j <- j + Ynlev
      }

      P_combs <- combn(P, 2, simplify = FALSE)
      D <- lapply(1:length(P_combs), function(i) P_combs[[i]][[2]] - P_combs[[i]][[1]])
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
      diffs$y <- factor(diffs$y, levels = Ylevs)
      diffs <- diffs %>% relocate(y, .before = 1)

      ## Output ----------------------------------------------------------------

      preds <- rename(preds, !!sym(varname) := x)
      diffs <- rename(diffs, !!sym(varname) := x)
      output <- list(
        preds = tibble::as_tibble(preds),
        diffs = tibble::as_tibble(diffs),
        seed = seed,
        sims = sims,
        formula = as.formula(modform))
      class(output) <- "svyEffects"
      attributes(output)$predvar <- varname
      attributes(output)$depvar <- Yname
      attributes(output)$method <- "AME"
      return(output)

    }
  }



  #==========  WITH MODERATOR VARIABLE =========================================
  # Note: Interactive models only output predictions, not differences.
  # (b/c interactive models should use second differences--this is on to-do list)

  if(!is.null(byvar)) {


    # Check arguments up front to stop executionn before running simulations
    if(isFALSE(byvar %in% names(data))) {
      stop(print(paste0(
        "byvar ", byvar, " not found in survey design object. Maybe check your spelling?")))}
    if(!is.null(bynvals) & isFALSE(is.numeric(bynvals))) {
      stop(print(paste0(
        "Non-numeric value entered for bynvals. Please enter a numeric value.")))}


    # CATEGORICAL * CATEGORICAL ================================================
    # (creates as many simulation dataframes as there are combinations of varname and byvar)

    if(is.factor(data[[varname]]) & is.factor(data[[byvar]])) {

      # Define predictor (varname) and moderator (byvar) variables
      levs <- levels(data[[varname]])
      nlev <- length(levs)
      bylevs <- levels(data[[byvar]])

      # Create dataframes for simulations
      df_list <- lapply(1:(length(levs) * length(bylevs)), function(x) x <- data)

      # Set main predictor levels across all dataframes
      # Sequence should be loop through each lev and repeat same number times as bylev
      i <- 1
      j <- 1
      k <- 1
      while(j <= length(bylevs)) {
        if (i <= length(levs)) {
          df_list[[k]] <- df_list[[k]] %>%
            mutate(!!sym(varname) := factor(i, levels = 1:length(levs), labels = levs))
          k <- k + 1
          i <- i + 1
        } else {
          i <- 1
          j <- j+1
        }
      }

      # Set byvar levels across all dataframes
      # Sequence should be set bylev to first level same number of times as levs
      # and then repeat same numer of times as bylev
      i <- 1
      j <- 1
      k <- 1
      while(j <= length(bylevs)) {
        if (i <= length(levs)) {
          df_list[[k]] <- df_list[[k]] %>%
            mutate(!!sym(byvar) := factor(j, levels = 1:length(bylevs), labels = bylevs))
          k <- k + 1
          i <- i + 1
        } else {
          i <- 1
          j <- j+1
        }
      }

      # Create model matrices
      X_list <- lapply(1:length(df_list), function(i)
        model.matrix(as.formula(modform), data = df_list[[i]]))

      # Parametric bootstrap
      B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))

      # Calculate probabilities and weight
      res <- NULL
      for(i in 1:nrow(B)){
        tmpB <- cbind(0, matrix(B[i,], ncol = (Ynlev-1), byrow = TRUE))
        EXB <- lapply(1:(length(levs) * length(bylevs)), function(i) exp(X_list[[i]] %*% tmpB))
        P <- lapply(1:(length(levs) * length(bylevs)), function(i) prop.table(EXB[[i]], 1))
        res_temp <- sapply(as.data.frame(P), weighted.mean, w = data$`(weights)`)
        res_temp <- t(as.data.frame(res_temp))
        rownames(res_temp) <- NULL
        res <- rbind(res, res_temp)
      }

      # Assemble table of predictions
      preds <- expand.grid(
        y = levels(droplevels(Yvar)),   # depvar
        x = levels(data[[varname]]),   # varname
        z = levels(data[[byvar]]))   # byvar
      preds$predicted<- colMeans(res)
      preds$conf.low <- apply(res, 2, low)
      preds$conf.high <- apply(res, 2, high)
      preds$type <- "Probability"
      preds <- preds %>%
        dplyr::rename(
          !!sym(varname) := "x",
          !!sym(byvar) := "z") %>%
        relocate(y, .before=1)

    }

    #========== CONTINUOUS * CATEGORICAL =======================================
    # (this is, basically, a loop of the univariate continuous function)

    if(is.numeric(data[[varname]]) & is.factor(data[[byvar]])) {

      # Define predictor sequence
      varname_seq <- seq(min(data[[varname]]), max(data[[varname]]), length = nvals)

      # Create simulation dataframe for each level of byvar
      bylevs <- levels(data[[byvar]])
      df_list <- lapply(1:length(bylevs), function(i) i <- data)
      df_list <- lapply(1:length(bylevs), function(i)
        df_list[[i]] <- df_list[[i]] %>%
          dplyr::mutate(!!sym(byvar) := factor(i, levels = 1:length(bylevs), labels = bylevs))
      )

      # Parametric bootstrap
      B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))

      # Calculate probabilities and weight
      res_list <- NULL
      for(j in 1:length(df_list)) {
        res <- array(dim = c(sims, length(Ylevs), length(varname_seq)))
        for(i in seq_along(varname_seq)){
          new <- df_list[[j]] %>% mutate(!!sym(varname) := varname_seq[i])
          Xmat <- model.matrix(as.formula(modform), data = new)
          res[,,i] <- t(apply(B, 1, function(b) {
            apply(
              prop.table(
                exp(Xmat %*% cbind(0, matrix(b, ncol = (length(Ylevs)-1), byrow = TRUE))),
                1),
              2,
              weighted.mean,
              w = data$`(weights)`)
          }))
        }
        res_list[[j]] <- res
      }

      # Generate table of probabilities
      preds_list <- NULL
      preds_list <- lapply(1:length(res_list), function(i) {
        preds_tmp <- data.frame(x = varname_seq)
        m <- as.data.frame(t(apply(res_list[[i]], c(2,3), mean)))
        names(m) <- paste0("predicted_", levels(droplevels(Yvar)))
        l <- as.data.frame(t(apply(res_list[[i]], c(2,3), low)))
        names(l) <- paste0("conf.low_", levels(droplevels(Yvar)))
        u <- as.data.frame(t(apply(res_list[[i]], c(2,3), high)))
        names(u) <- paste0("conf.high_", levels(droplevels(Yvar)))
        preds_tmp <- cbind(preds_tmp, m, l, u)
        preds_tmp <- preds_tmp %>%
          pivot_longer(-x,
                       names_pattern = "(.*)_(.*)",
                       names_to = c(".value", "y"))
        preds_tmp$z <- factor(bylevs)[[i]]
        preds_list[[i]] <- preds_tmp
      })
      preds <- do.call(rbind, preds_list)
      preds$type <- "Probability"
      preds$y <- factor(preds$y, levels = Ylevs)
      preds$z <- factor(preds$z, levels = bylevs)
      preds <- preds %>%
        dplyr::rename(
          !!sym(varname) := "x",
          !!sym(byvar) := "z") %>%
        relocate(y, .before=1)

    }

    # CATEGORICAL * CONTINUOUS =================================================
    # (loops the univariate categorical function over multiple levels of moderator)

    if(is.factor(data[[varname]]) & is.numeric(data[[byvar]])) {

      # Define predictor (varname) and moderator (byvar) variables
      levs <- levels(data[[varname]])
      nlev <- length(levs)

      bychange <- match.arg(bychange)
      byvar_min <- switch(bychange,
                          NULL = as.numeric(survey::svymean(data[byvar], svydata)) -
                            as.numeric(sqrt(survey::svyvar(data[byvar], svydata))),
                          sd = as.numeric(survey::svymean(data[byvar], svydata)) -
                            as.numeric(sqrt(survey::svyvar(data[byvar], svydata))),
                          range = min(data[[byvar]])
      )
      byvar_max <- switch(bychange,
                          NULL = as.numeric(survey::svymean(data[byvar], svydata)) +
                            as.numeric(sqrt(survey::svyvar(data[byvar], svydata))),
                          sd = as.numeric(survey::svymean(data[byvar], svydata)) +
                            as.numeric(sqrt(survey::svyvar(data[byvar], svydata))),
                          range = max(data[[byvar]])
      )
      bylevs <- seq(byvar_min, byvar_max, length = bynvals)

      # Create dataframes for simulations
      df_list <- lapply(1:(length(levs) * length(bylevs)), function(x) x <- data)

      # Set main predictor levels across all dataframes
      # Sequence should be loop through each lev and repeat same number times as bylev
      i <- 1
      j <- 1
      k <- 1
      while(j <= length(bylevs)) {
        if (i <= length(levs)) {
          df_list[[k]] <- df_list[[k]] %>%
            mutate(!!sym(varname) := factor(i, levels = 1:length(levs), labels = levs))
          k <- k + 1
          i <- i + 1
        } else {
          i <- 1
          j <- j+1
        }
      }

      # Set byvar levels across all dataframes
      # Sequence should be set bylev to first level same number of times as levs
      # and then repeat same numer of times as bylev
      i <- 1
      j <- 1
      k <- 1
      while(j <= length(bylevs)) {
        if (i <= length(levs)) {
          df_list[[k]] <- df_list[[k]] %>%
            mutate(!!sym(byvar) := bylevs[[j]])
          k <- k + 1
          i <- i + 1
        } else {
          i <- 1
          j <- j+1
        }
      }

      # Create model matrices
      X_list <- lapply(1:length(df_list), function(i)
        model.matrix(as.formula(modform), data = df_list[[i]]))

      # Parametric bootstrap
      B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))

      # Calculate predicted probabilities and weight
      res <- NULL
      for(i in 1:nrow(B)){
        tmpB <- cbind(0, matrix(B[i,], ncol = (Ynlev-1), byrow = TRUE))
        EXB <- lapply(1:(length(levs) * length(bylevs)), function(i) exp(X_list[[i]] %*% tmpB))
        P <- lapply(1:(length(levs) * length(bylevs)), function(i) prop.table(EXB[[i]], 1))
        res_temp <- sapply(as.data.frame(P), weighted.mean, w = data$`(weights)`)
        res_temp <- t(as.data.frame(res_temp))
        rownames(res_temp) <- NULL
        res <- rbind(res, res_temp)
      }

      # Assemble table of probabilities
      preds <- expand.grid(
        y = levels(droplevels(Yvar)),   # depvar
        x = levels(data[[varname]]),   # varname
        z = as.factor(round(bylevs, 3)))   # byvar
      preds$predicted<- colMeans(res)
      preds$conf.low <- apply(res, 2, low)
      preds$conf.high <- apply(res, 2, high)
      preds$type <- "Probability"
      preds$y <- factor(preds$y, levels = Ylevs)
      preds <- preds %>%
        dplyr::rename(
          !!sym(varname) := "x",
          !!sym(byvar) := "z") %>%
        relocate(y, .before=1)

    }

    # CONTINUOUS * CONTINUOUS ==================================================
    # (same as continues * categorical, but first factorizes the moderator)

    if(is.numeric(data[[varname]]) & is.numeric(data[[byvar]])) {

      # Define variables
      varname_seq <- seq(min(data[[varname]]), max(data[[varname]]), length = nvals)

      bychange <- match.arg(bychange)
      byvar_min <- switch(bychange,
                          NULL = as.numeric(survey::svymean(data[byvar], svydata)) -
                            as.numeric(sqrt(survey::svyvar(data[byvar], svydata))),
                          sd = as.numeric(survey::svymean(data[byvar], svydata)) -
                            as.numeric(sqrt(survey::svyvar(data[byvar], svydata))),
                          range = min(data[[byvar]])
      )
      byvar_max <- switch(bychange,
                          NULL = as.numeric(survey::svymean(data[byvar], svydata)) +
                            as.numeric(sqrt(survey::svyvar(data[byvar], svydata))),
                          sd = as.numeric(survey::svymean(data[byvar], svydata)) +
                            as.numeric(sqrt(survey::svyvar(data[byvar], svydata))),
                          range = max(data[[byvar]])
      )
      bylevs <- seq(byvar_min, byvar_max, length = bynvals)

      # Create as many simulation dataframes as there are levels of moderator
      df_list <- lapply(1:length(bylevs), function(i) i <- data)
      df_list <- lapply(1:length(bylevs), function(i)
        df_list[[i]] <- df_list[[i]] %>%
          dplyr::mutate(!!sym(byvar) := bylevs[[i]]))

      # Parametric bootstrap
      B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))

      # Calculate predicted probabilities and weight
      res_list <- NULL
      for(j in 1:length(df_list)) {
        res <- array(dim = c(sims, length(Ylevs), length(varname_seq)))
        for(i in seq_along(varname_seq)){
          new <- df_list[[j]] %>% mutate(!!sym(varname) := varname_seq[i])
          Xmat <- model.matrix(as.formula(modform), data = new)
          res[,,i] <- t(apply(B, 1, function(b) {
            apply(
              prop.table(
                exp(Xmat %*% cbind(0, matrix(b, ncol = (length(Ylevs)-1), byrow = TRUE))),
                1),
              2,
              weighted.mean,
              w = data$`(weights)`)
          }))
        }
        res_list[[j]] <- res
      }

      # Assemble table of probabilities
      preds_list <- NULL
      preds_list <- lapply(1:length(res_list), function(i) {
        preds_tmp <- data.frame(x = varname_seq)
        m <- as.data.frame(t(apply(res_list[[i]], c(2,3), mean)))
        names(m) <- paste0("predicted_", levels(droplevels(Yvar)))
        l <- as.data.frame(t(apply(res_list[[i]], c(2,3), low)))
        names(l) <- paste0("conf.low_", levels(droplevels(Yvar)))
        u <- as.data.frame(t(apply(res_list[[i]], c(2,3), high)))
        names(u) <- paste0("conf.high_", levels(droplevels(Yvar)))
        preds_tmp <- cbind(preds_tmp, m, l, u)
        preds_tmp <- preds_tmp %>%
          pivot_longer(-x,
                       names_pattern = "(.*)_(.*)",
                       names_to = c(".value", "y"))
        preds_tmp$z <- factor(round(bylevs, 3))[[i]]
        preds_list[[i]] <- preds_tmp
      })
      preds <- do.call(rbind, preds_list)
      preds$type <- "Probability"
      preds$y <- factor(preds$y, levels = Ylevs)
      preds <- preds %>%
        dplyr::rename(
          !!sym(varname) := "x",
          !!sym(byvar) := "z") %>%
        relocate(y, .before=1)

    }

    # OUTPUT ===================================================================

    output <- list(
      preds = tibble::as_tibble(preds),
      seed = seed,
      sims = sims,
      formula = as.formula(modform))
    class(output) <- "svyEffects"
    attributes(output)$predvar <- varname
    attributes(output)$byvar <- byvar
    attributes(output)$depvar <- Yname
    attributes(output)$method <- "AME"
    return(output)

  }

}
