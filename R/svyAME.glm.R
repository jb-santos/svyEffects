#' @title Average marginal effects for binary logit models of survey-weighted data
#'
#'
#' @description Calculates predicted probabilities and differences in probabilities
#' for a predictor holding all other observations at observed values
#' (i.e. a true average marginal effect). Uses simulations (the parametric
#' bootstrap) to derive 95% confidence intervals.
#'
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
#' @param byvar For interaction effects, a character string denoting the name of
#' the moderator variable.
#' @param bynvals For interaction effects with a numerical moderator, a scalar
#' denoting the sequence length  spanning the range of the moderator variable
#' for which effects are to be calculated (default: 3).
#' @param sims Scalar denoting how many simulations to conduct (default: 2500).
#' @param seed Seed value for random number generator. By default, the function
#' picks a random integer between 1 and 1,000,000. If you save the output of
#' this function, it will save the seed value used for simulations in the slot
#' \code{$seed}.
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
#' @author John Santos & Dave Armstrong
#'
#'
#' @export
#'
#'
#' @examples
#' data(ces19)
#' library(survey)
#' ces19_svy <- svydesign(ids = ~1, strata = NULL, weights = ~pesweight, data = ces19, digits = 3)
#' VOTECON <- svyglm(votecon ~ agegrp + gender + educ + region + marketlib, design = ces19_svy, family = binomial)
#' svyAME(VOTECON, varname = "educ", weightvar = "pesweight", seed = 2019)
#' svyAME(VOTECON, varname = "marketlib", weightvar = "pesweight", seed = 2019)
#'
svyAME.glm <- function(obj,
                       varname,
                       weightvar,
                       nvals = 11,
                       diffchange = c("range", "unit", "sd"),
                       byvar = NULL,
                       bynvals = 3,
                       sims = 2500,
                       seed = NULL,
                       ...) {

  # SETUP ======================================================================

  # Get data
  data <- model.frame(obj)
  data <- data %>% rename(!!sym(weightvar) := `(weights)`)
  svydata <- survey::svydesign(ids = ~1, strata = NULL, weights = data[weightvar], data = data)
  if(is.character(data[[varname]])) {data[[varname]] <- as.factor(data[[varname]])}

  # Set random number generator
  if(is.null(seed)){
    seed <- round(runif(1, min=0, max=1000000))
  } else {
    seed <- seed
  }
  set.seed(seed)


  # NO MODERATOR VARIABLE ======================================================

  if(is.null(byvar)) {

    ## NUMERIC / CONTINUOUS predictors =========================================

    if(is.numeric(data[[varname]])) {

      ### Predicted probabilities ----------------------------------------------

      varname_seq <- seq(min(data[[varname]]), max(data[[varname]]), length = nvals)

      res <- sapply(1:length(varname_seq), function(i) {
        tmp_d <- data
        tmp_d[[varname]] <- varname_seq[i]
        B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))
        X <- model.matrix(formula(obj), tmp_d)
        Xb <- X %*% t(B)
        p <- plogis(Xb)
        m <- apply(p, 2, weighted.mean, w = data[[weightvar]])
        m
      })

      preds <- tibble(
        x = varname_seq,
        predicted = colMeans(res),
        conf.low = apply(res, 2, low),
        conf.high = apply(res, 2, high),
        type = "Probability")

      ### Differences in predicted probabilities -------------------------------

      diffchange <- match.arg(diffchange)
      tmp0 <- tmp1 <- data
      tmp0[varname] <- switch(diffchange,
                              NULL = min(data[varname]),
                              range = min(data[varname]),
                              unit = survey::svymean(data[varname], svydata) - .5,
                              sd = survey::svymean(data[varname], svydata) -
                                (sqrt(survey::svyvar(data[varname], svydata)) / 2)
                              )
      tmp1[varname] <- switch(diffchange,
                              NULL = max(data[varname]),
                              range = max(data[varname]),
                              unit = survey::svymean(data[varname], svydata) + .5,
                              sd = survey::svymean(data[varname], svydata) +
                                (sqrt(survey::svyvar(data[varname], svydata)) / 2)
                              )

      B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))
      X0 <- model.matrix(formula(obj), tmp0)
      X1 <- model.matrix(formula(obj), tmp1)
      Xb0 <- X0 %*% t(B)
      Xb1 <- X1 %*% t(B)
      p0 <- as.data.frame(plogis(Xb0))
      p1 <- as.data.frame(plogis(Xb1))
      diff <- p1 - p0
      mean_diff <- apply(diff, 2, weighted.mean, w = data[[weightvar]])

      diffs <- tibble(
        x = paste0("Delta (",
                   diffchange,
                   ") : ",
                   round(tmp0[[varname]][1], 3),
                   " - ",
                   round(tmp1[[varname]][1], 3)),
        predicted = mean(mean_diff),
        conf.low = quantile(mean_diff, .025),
        conf.high = quantile(mean_diff, .975),
        type = "Difference")

      ### Output ---------------------------------------------------------------

      preds <- dplyr::rename(preds, !!sym(varname) := x)
      diffs <- dplyr::rename(diffs, !!sym(varname) := x)
      output <- list(
        preds = dplyr::as_tibble(preds),
        diffs = dplyr::as_tibble(diffs),
        seed = seed,
        sims = sims)
      class(output) <- "svyEffects"
      attributes(output)$predvar <- varname
      attributes(output)$depvar <- colnames(obj$model)[1]
      return(output)

    } else {

      ## FACTOR / CATEGORICAL predictors =======================================

      ### Predicted probabilities ----------------------------------------------

      levs <- levels(data[[varname]])
      nlev <- length(levs)
      nlevC2 <- factorial(nlev) / (2 * factorial(nlev-2))
      DFs <- lapply(1:nlev, function(i) i <- data)
      DFs <- lapply(1:nlev, function(i)
        DFs[[i]] %>% dplyr::mutate(!!sym(varname) := factor(i, levels = 1:nlev, labels = levs)))
      X <- lapply(1:nlev, function(i) model.matrix(formula(obj), data = DFs[[i]]))

      B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))
      XB <- lapply(1:nlev, function(i) X[[i]] %*% t(B))
      Pr <- lapply(1:nlev, function(i) plogis(XB[[i]]))
      WgtPr <- lapply(1:nlev, function(i)
        apply(Pr[[i]], 2, function(x) weighted.mean(x, data[[weightvar]])))

      preds <- tibble(
        x = factor(levs, levels = levs),
        predicted = sapply(WgtPr, mean),
        conf.low = sapply(WgtPr, low),
        conf.high = sapply(WgtPr, high),
        type="Probability")

      ### Differences in predicted probabilities --------------------------------

      Pr_combs <- combn(Pr, 2)
      Pr_diffs <- mapply(x = Pr_combs[2,], y = Pr_combs[1,],
                         function(x, y) {unlist(x) - unlist(y)})   # calculate first differences
      Pr_diffs <- relist(Pr_diffs, Pr_combs)
      WgtPr_D <- NULL  # create list of weighted first differences
      for (i in 1:nlevC2) {
        WgtPr_D[[i]] <-
          apply(Pr_diffs[[i]], 2, function(x)weighted.mean(x, data[[weightvar]]))
        }

      diffs <- tibble(
        varname1 = combn(levs, 2)[1,] ,
        varname2 = combn(levs, 2)[2,] ,
        predicted = sapply(WgtPr_D, mean),
        conf.low = sapply(WgtPr_D, low),
        conf.high = sapply(WgtPr_D, high))
      diffs <- diffs %>%
        mutate(x = paste(varname2, " - ", varname1, sep = ""), .before=1) %>%
        mutate(type = "Difference") %>%
        dplyr::select(-c("varname1", "varname2"))

      ### Output ---------------------------------------------------------------

      preds <- rename(preds, !!sym(varname) := x)
      diffs <- rename(diffs, !!sym(varname) := x)
      output <- list(
        preds = dplyr::as_tibble(preds),
        diffs = dplyr::as_tibble(diffs),
        seed = seed,
        sims = sims)
      class(output) <- "svyEffects"
      attributes(output)$predvar <- varname
      attributes(output)$depvar <- colnames(obj$model)[1]
      return(output)
    }

  }


  #==========  WITH MODERATOR VARIABLE =========================================
  # Note: Interactive models only output predictions, not differences.
  # (b/c interactive models should use second differences--this is on to-do list)

  if(!is.null(byvar)) {

    ## CATEGORICAL * CATEGORICAL -----------------------------------------------

    # (creates as many simulation dataframes as there are combinations of varname and byvar)

    if(is.factor(data[[varname]]) & is.factor(data[[byvar]])) {

      # Setup
      levs <- levels(data[[varname]])
      nlev <- as.numeric(nlevels(data[[varname]]))
      nseq <- seq(1:nlev)
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
      # and then repeat same number of times as bylev
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
        model.matrix(formula(obj), data = df_list[[i]]))

      # Parametric bootstrap
      B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))

      # Calculate predicted probabilities
      XB_list <- lapply(1:length(X_list), function(i) X_list[[i]] %*% t(B))
      pr_list <- lapply(1:length(XB_list), function(i) plogis(XB_list[[i]]))

      # Weight probabilities
      wgtpr_list <- lapply(1:length(pr_list), function(i)
        apply(pr_list[[i]], 2, function(x)
          weighted.mean(x, data[[weightvar]])
        )
      )

      # Create table of predictions
      preds <-  expand.grid(levs, bylevs) %>%
        dplyr::mutate(
          predicted = sapply(wgtpr_list, mean),
          conf.low = sapply(wgtpr_list, low),
          conf.high = sapply(wgtpr_list, high),
          type = "Probability") %>%
        dplyr::rename(
          !!sym(varname) := "Var1",
          !!sym(byvar) := "Var2")

    }

    ## CONTINUOUS * CATEGORICAL ------------------------------------------------
    # (this is, basically, a loop of the univariate continuous function)

    if(is.numeric(data[[varname]]) & is.factor(data[[byvar]])) {

      varname_seq <- seq(min(data[[varname]]), max(data[[varname]]), length = nvals)

      bylevs <- levels(data[[byvar]])
      bydata <- lapply(1:length(bylevs), function(i) i <- data)
      bydata <- lapply(1:length(bylevs), function(i)
        bydata[[i]] %>%
          mutate(!!sym(byvar) := factor(i, levels = 1:length(bylevs), labels = bylevs)))

      byres <- NULL
      byres <- lapply(1:length(bylevs), function(j) {
        res <- NULL
        for(i in seq_along(varname_seq)){
          tmp_d <- bydata[[j]]
          tmp_d[[varname]] <- varname_seq[i]
          B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))
          X <- model.matrix(formula(obj), tmp_d)
          Xb <- X %*% t(B)
          p <- plogis(Xb)
          m <- apply(p, 2, weighted.mean, w = data[[weightvar]])
          res <- rbind(res, data.frame(
            z = bylevs[[j]],
            x = varname_seq[i],
            predicted = mean(m),
            conf.low = low(m),
            conf.high = high(m)))
          rownames(res) <- NULL
        }
        byres[[j]] <- as.data.frame(res)
      })

      preds <- do.call(rbind, byres)
      preds <- dplyr::rename(preds, !!sym(varname) := x)
      preds <- dplyr::rename(preds, !!sym(byvar) := z)
      preds[[byvar]] <- factor(preds[[byvar]], levels = bylevs)

    }

    ## CATEGORICAL * CONTINUOUS ------------------------------------------------
    # (loops the univariate categorical function over multiple levels of moderator)

    if(is.factor(data[[varname]]) & is.numeric(data[[byvar]])) {

      levs <- levels(data[[varname]])
      nlev <- length(levs)
      bylevs <- seq(min(data[[byvar]]), max(data[[byvar]]), length = bynvals)

      df_list <- lapply(1:(length(levs) * length(bylevs)), function(x)
        x <- data)

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
          j <- j + 1
        }
      }

      # Parametric bootstrap
      B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))

      # Calculate predicted probabilities
      X_list <- lapply(1:length(df_list), function(i) model.matrix(formula(obj), data = df_list[[i]]))
      XB_list <- lapply(1:length(X_list), function(i) X_list[[i]] %*% t(B))
      pr_list <- lapply(1:length(XB_list), function(i) plogis(XB_list[[i]]))

      # Weight probabilities
      wgtpr_list <- lapply(1:length(pr_list), function(i)
        apply(pr_list[[i]], 2, function(x) weighted.mean(x, data[[weightvar]])))

      # Assemble table of probabilities
      preds <-  expand.grid(
        x = levs,
        z = as.factor(round(bylevs, 3))) %>%
        mutate(
          predicted = sapply(wgtpr_list, mean),
          conf.low = sapply(wgtpr_list, low),
          conf.high = sapply(wgtpr_list, high),
          type = "Predicted probability")
      preds <- dplyr::rename(preds, !!sym(varname) := x)
      preds <- dplyr::rename(preds, !!sym(byvar) := z)

    }

    ## CONTINUOUS * CONTINUOUS -------------------------------------------------
    # (same as continues * categorical, but first factorizes the moderator)

    if(is.numeric(data[[varname]]) & is.numeric(data[[byvar]])) {

      varname_seq <- seq(min(data[[varname]]), max(data[[varname]]), length = nvals)

      bylevs <- seq(min(data[[byvar]]), max(data[[byvar]]), length = bynvals)
      bydata <- lapply(1:length(bylevs), function(i) i <- data)
      bydata <- lapply(1:length(bylevs), function(i) bydata[[i]] %>%
                         mutate(!!sym(byvar) := bylevs[[i]]))

      byres <- NULL
      byres <- lapply(1:length(bylevs), function(j) {
        res <- NULL
        for(i in seq_along(varname_seq)) {
          tmp_d <- bydata[[j]]
          tmp_d[[varname]] <- varname_seq[i]
          B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))
          X <- model.matrix(formula(obj), tmp_d)
          Xb <- X %*% t(B)
          p <- plogis(Xb)
          m <- apply(p, 2, weighted.mean, w = data[[weightvar]])
          res <- rbind(res, data.frame(
            z = bylevs[[j]],
            x = varname_seq[i],
            predicted = mean(m),
            conf.low = low(m),
            conf.high = high(m)))
          rownames(res) <- NULL
        }
        byres[[j]] <- as.data.frame(res)
      })

      preds <- do.call(rbind, byres)
      preds <- dplyr::rename(preds, !!sym(varname) := x)
      preds <- dplyr::rename(preds, !!sym(byvar) := z)
      preds[[byvar]] <- factor(preds[[byvar]], levels = bylevs)

    }

    ## OUTPUT ==================================================================

    output <- list(
      preds = dplyr::as_tibble(preds),
      seed = seed,
      sims = sims)
    class(output) <- "svyEffects"
    attributes(output)$predvar <- varname
    attributes(output)$byvar <- byvar
    attributes(output)$depvar <- colnames(obj$model)[1]
    return(output)

  }

}
