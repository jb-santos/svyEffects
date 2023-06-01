#' Average Marginal Effects For Ordered Logit Models Of Survey-Weighted Data
#'
#'
#' @description Calculates predicted probabilities and differences in probabilities
#' for a predictor holding all other observations at observed values
#' (i.e. a true average marginal effect). Uses simulations (the parametric
#' bootstrap) to derive 95% confidence intervals.
#'
#'
#' @param obj Model object of class \code{survey::svyolr}.
#' @param varname Character string denoting the name of the predictor variable
#' for which effects are to be calculated.
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
#' CONLDR <- svyolr(ftconldr ~ agegrp + gender + educ + region + marketlib,
#'   design = ces19_svy)
#' svyAME(CONLDR, varname = "region", seed = 2019)
#' svyAME(CONLDR, varname = "marketlib", seed = 2019)
#' }
#'
svyAME.svyolr <- function(obj,
                          varname,
                          nvals = 11,
                          diffchange = c("range", "unit", "sd"),
                          byvar = NULL,
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

    #========== NUMERIC / CONTINUOUS predictors ==================================

    if(is.numeric(data[[varname]])) {

      ## Predicted probabilities -----------------------------------------------

      varname_seq <- seq(min(data[[varname]]), max(data[[varname]]), length = nvals)

      preds <- NULL

      for(n in seq_along(varname_seq)){
        new <- dplyr::mutate(data, !!sym(varname) := varname_seq[n])
        Xmat <- model.matrix(formula(obj), data = new)[,-1]
        B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))
        tau_ind <- grep("\\|", names(coef(obj)))
        Tau <- B[,tau_ind]
        B <- B[,-tau_ind]

        CP <- NULL
        for (i in 1:ncol(Tau)) {
          CP[[i]] <- plogis(matrix(Tau[,i],
                                   ncol = nrow(Tau),
                                   nrow = nrow(new),
                                   byrow = TRUE) - Xmat %*% t(B))
        }

        P <- NULL
        for (i in 1:ncol(Tau)) {
          if(i == 1) {P[[i]] <- CP[[i]]}
          else{P[[i]] <- CP[[i]] - CP[[(i - 1)]]}}
        P[[(ncol(Tau)+1)]] <- 1 - CP[[ncol(Tau)]]

        WgtP <- lapply(1:length(P), function(i){
          apply(P[[i]], 2, function(x) weighted.mean(x, data$`(weights)`))
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

      ## Differences in predicted probabilities --------------------------------

      diffchange <- match.arg(diffchange)
      tmp0 <- switch(diffchange,
                     NULL = min(data[varname]),
                     range = min(data[varname]),
                     unit = survey::svymean(data[varname], svydata) - .5,
                     sd = survey::svymean(data[varname], svydata) -
                       (sqrt(survey::svyvar(data[varname], svydata)) / 2)
                     )
      tmp1 <- switch(diffchange,
                     NULL = max(data[varname]),
                     range = max(data[varname]),
                     unit = survey::svymean(data[varname], svydata) + .5,
                     sd = survey::svymean(data[varname], svydata) +
                       (sqrt(survey::svyvar(data[varname], svydata)) / 2)
                     )

      delta <- c(tmp0, tmp1)
      D_WgtP <- NULL

      for(n in seq_along(delta)) {
        new <- dplyr::mutate(data, !!sym(varname) := delta[n])
        Xmat <- model.matrix(formula(obj), data = new)[,-1]
        B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))
        tau_ind <- grep("\\|", names(coef(obj)))
        Tau <- B[,tau_ind]
        B <- B[,-tau_ind]

        CP <- NULL
        for (i in 1:ncol(Tau)) {
          CP[[i]] <- plogis(matrix(Tau[,i], ncol = nrow(Tau),
                                   nrow = nrow(new), byrow = TRUE) - Xmat %*% t(B))
        }
        P <- NULL
        for (i in 1:ncol(Tau)) {
          if(i == 1) {P[[i]] <- CP[[i]]}
          else{P[[i]] <- CP[[i]]- CP[[(i - 1)]]}}
        P[[(ncol(Tau) + 1)]] <- 1 - CP[[ncol(Tau)]]

        WgtP <- lapply(1:length(P), function(i){
          apply(P[[i]], 2, function(x)weighted.mean(x, data$`(weights)`))
        })

        D_WgtP[[n]] <- WgtP
      }

      D_WgtP <- mapply(x = D_WgtP[[2]], y = D_WgtP[[1]],
                       function(x, y) unlist(x) - unlist(y))

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

      ## Output ----------------------------------------------------------------

      preds <- rename(preds, !!sym(varname) := x)
      diffs <- rename(diffs, !!sym(varname) := x)
      output <- list(
        preds = tibble::as_tibble(preds),
        diffs = tibble::as_tibble(diffs),
        seed = seed,
        sims = sims,
        formula = formula(obj))
      class(output) <- "svyEffects"
      attributes(output)$predvar <- varname
      attributes(output)$depvar <- colnames(model.frame(obj))[1]
      attributes(output)$method <- "AME"
      return(output)

    } else {

      #========== FACTOR / CATEGORICAL predictors ==============================

      ## Predicted probabilities -----------------------------------------------

      nlev <- as.numeric(nlevels(data[[varname]]))
      levs <- levels(data[[varname]])
      nlevC2 <- (factorial(nlev) / (2 * factorial(nlev - 2)))

      # Create dataframes for simulations
      DFs <- lapply(1:nlev, function(x) x <- data)
      for (i in 1:nlev) {DFs[[i]] <- DFs[[i]] %>%
        dplyr::mutate(!!sym(varname) := factor(i, levels=1:nlev, labels = levs))}
      X <- NULL
      for (i in 1:nlev) {X[[i]] <- model.matrix(formula(obj), data = DFs[[i]])[,-1]}

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
          CP[[i]] <- plogis(matrix(Tau[,i], ncol = nrow(Tau),
                                   nrow = nrow(X[[k]]), byrow = TRUE) - X[[k]] %*% t(B))
        }
        for (j in 1:ncol(Tau)) {
          if(j == 1) {
            P[[j]] <- CP[[j]]
          }
          else{
            P[[j]] <- CP[[j]]- CP[[(j - 1)]]
          }
          P[[(ncol(Tau) + 1)]] <- 1 - CP[[(ncol(Tau))]]
        }
        Pr[[k]] <- P
      }

      # Apply weights
      Pr_ul <- unlist(Pr, recursive=FALSE)
      WgtPr <- lapply(1:length(Pr_ul), function(i)
        apply(Pr_ul[[i]], 2, function(x)weighted.mean(x, data$`(weights)`)))

      # Assemble table of predicted probabilities
      preds <- expand.grid(y = obj$lev, x = levs)
      preds$predicted <- sapply(WgtPr, mean)
      preds$conf.low <- sapply(WgtPr, low)
      preds$conf.high <- sapply(WgtPr, high)
      preds$type <- "Probability"

      ## Differences in predicted probabilities --------------------------------

      WgtPr <- as.data.frame(WgtPr)
      i <- 1
      j <- 1
      P <- NULL
      while (i <= nlev) {
        P[[i]] <- as.data.frame(WgtPr[,j:(j + length(obj$lev) - 1)])
        i <- i + 1
        j <- j + length(obj$lev)
      }

      P_combs <- combn(P, 2, simplify = FALSE)
      D <- lapply(1:length(P_combs), function(i) P_combs[[i]][[2]] - P_combs[[i]][[1]])
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

      ## Output ----------------------------------------------------------------

      preds <- rename(preds, !!sym(varname) := x)
      diffs <- rename(diffs, !!sym(varname) := x)
      output <- list(
        preds = tibble::as_tibble(preds),
        diffs = tibble::as_tibble(diffs),
        seed = seed,
        sims = sims,
        formula = formula(obj))
      class(output) <- "svyEffects"
      attributes(output)$predvar <- varname
      attributes(output)$depvar <- colnames(model.frame(obj))[1]
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

      # Setup
      levs <- levels(data[[varname]])
      bylevs <- levels(data[[byvar]])

      # Set up simulation dataframes
      df_list <- lapply(1:length(levs), function(i) i <- model.frame(obj))
      df_list <- lapply(1:length(levs), function(i)
        df_list[[i]] <- df_list[[i]] %>%
          dplyr::mutate(!!sym(varname) := factor(i, levels = 1:length(levs), labels = levs)))

      # Replicate original simulation dataframes across levels of byvar
      # and set byvar levels
      by_list <- lapply(1:length(bylevs), function(i) i <- df_list)
      i <- 1
      j <- 1
      while(j <= length(bylevs)) {
        if (i <= length(levs)) {
          by_list[[j]][[i]] <- by_list[[j]][[i]] %>%
            dplyr::mutate(!!sym(byvar) := factor(bylevs, levels = bylevs)[[j]])
          i <- i + 1
        } else {
          i <- 1
          j <- j + 1
        }
      }

      # Create model matrices
      X_list <- by_list
      i <- 1
      j <- 1
      while(j <= length(bylevs)) {
        if (i <= length(levs)) {
          X_list[[j]][[i]] <- model.matrix(formula(obj), data = by_list[[j]][[i]])[,-1]
          i <- i + 1
        } else {
          i <- 1
          j <- j + 1
        }
      }

      # Parametric bootstrap
      B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))
      tau_ind <- grep("\\|", names(coef(obj)))
      Tau <- B[,tau_ind]
      B <- B[,-tau_ind]

      # Generate predictions
      ByPr_list <- NULL
      ByPr_list <- lapply(1:length(bylevs), function(l) {
        CP <- NULL
        P <- NULL
        Pr_list <- NULL
        for(k in 1:length(levs)) {
          for(i in 1:ncol(Tau)) {
            CP[[i]] <- plogis(matrix(Tau[,i],
                                     ncol = nrow(Tau),
                                     nrow = nrow(X_list[[1]][[1]]),
                                     byrow = TRUE) - X_list[[l]][[k]] %*% t(B))
          }
          for (j in 1:ncol(Tau)) {
            if(j == 1) {
              P[[j]] <- CP[[j]]
            }
            else{
              P[[j]] <- CP[[j]]- CP[[(j - 1)]]
            }
            P[[(ncol(Tau) + 1)]] <- 1 - CP[[(ncol(Tau))]]
          }
          Pr_list[[k]] <- P
        }
        ByPr_list[[l]] <- Pr_list
      })

      # Weight predictions
      WgtByPr_list <- NULL
      for(j in 1:length(ByPr_list)) {
        Pr_list_ul <- NULL
        Pr_list_ul <- unlist(ByPr_list[[j]], recursive = FALSE)
        WgtPr_list <- NULL
        for (i in 1:length(Pr_list_ul)) {
          WgtPr_list[[i]] <- apply(Pr_list_ul[[i]], 2, function(x)
            weighted.mean(x, data$`(weights)`))
        }
        WgtByPr_list[[j]] <- WgtPr_list
      }


      # Assemble predictions
      preds <- NULL
      for(i in 1:length(WgtByPr_list)) {
        preds_tmp <- expand.grid(
          y = obj$lev,
          x = levs)
        preds_tmp$predicted <- sapply(WgtByPr_list[[i]], mean)
        preds_tmp$conf.low <- sapply(WgtByPr_list[[i]], low)
        preds_tmp$conf.high <- sapply(WgtByPr_list[[i]], high)
        preds_tmp$type <- "Probability"
        preds_tmp$z <- factor(bylevs[[i]], levels = bylevs)
        preds[[i]] <- preds_tmp
      }
      preds <- do.call("rbind", preds)
      preds <- preds %>%
        dplyr::relocate(z, .before = predicted) %>%
        dplyr::rename(!!sym(byvar) := z) %>%
        dplyr::rename(!!sym(varname) := x)

    }

    #========== CONTINUOUS * CATEGORICAL =======================================
    # (this is, basically, a loop of the univariate continuous function)

    if(is.numeric(data[[varname]]) & is.factor(data[[byvar]])) {

      # Define byvar
      bylevs <- levels(data[[byvar]])

      # Setup list of dataframes with levels of byvar set
      by_list <- lapply(1:length(bylevs), function(i) i <- model.frame(obj))
      by_list <- lapply(1:length(bylevs), function(i)
        by_list[[i]] <- by_list[[i]] %>%
          dplyr::mutate(!!sym(byvar) := factor(i, levels = 1:length(bylevs), labels = bylevs)))

      # Define varname_seq
      varname_seq <- seq(min(data[[varname]]), max(data[[varname]]), length = nvals)

      # Loop simulations over levels of varname and across bylist
      res_list <- NULL
      # Note: for loops and lapply appear to be equivalent in terms of speed
      # for(j in 1:length(by_list)) { ... }
      # res_list <- lapply(1:length(by_list), function(j) { ... })
      for(j in 1:length(by_list)) {
        res <- NULL
        for(k in seq_along(varname_seq)) {
          new <- by_list[[j]] %>% dplyr::mutate(!!sym(varname) := varname_seq[k])
          Xmat <- model.matrix(formula(obj), data = new)[,-1]
          B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))
          tau_ind <- grep("\\|", names(coef(obj)))
          Tau <- B[,tau_ind]
          B <- B[,-tau_ind]

          CP <- NULL
          for (i in 1:ncol(Tau)) {
            CP[[i]] <- plogis(matrix(Tau[,i],
                                     ncol = nrow(Tau),
                                     nrow = nrow(new),
                                     byrow = TRUE) - Xmat %*% t(B))
          }

          P <- NULL
          for (i in 1:ncol(Tau)) {
            if(i == 1) {P[[i]] <- CP[[i]]}
            else{P[[i]] <- CP[[i]] - CP[[(i - 1)]]}}
          P[[(ncol(Tau) + 1)]] <- 1 - CP[[ncol(Tau)]]

          WgtP <- lapply(1:length(P), function(i){
            apply(P[[i]], 2, function(x) weighted.mean(x, data$`(weights)`))
          })

          tmp <- data.frame(
            y = obj$lev,
            x = varname_seq[k],
            predicted = sapply(WgtP, mean),
            conf.low = sapply(WgtP, low),
            conf.high = sapply(WgtP, high),
            type = "Probability")

          res <- rbind(res, tmp)
        }
        res$z <- factor(bylevs, levels = bylevs)[[j]]
        res_list[[j]] <- res
      }
      preds <- do.call("rbind", res_list)
      preds$y <- factor(preds$y, levels = obj$lev)
      preds <- preds %>%
        dplyr::relocate(z, .before = predicted) %>%
        dplyr::rename(!!sym(byvar) := z) %>%
        dplyr::rename(!!sym(varname) := x)

    }

    # CATEGORICAL * CONTINUOUS =================================================
    # (loops the univariate categorical function over multiple levels of moderator)

    if(is.factor(data[[varname]]) & is.numeric(data[[byvar]])) {

      levs <- levels(data[[varname]])
      bylevs <- seq(min(data[[byvar]]), max(data[[byvar]]), length = bynvals)

      # Set up simulation dataframes
      df_list <- lapply(1:length(levs), function(i) i <- model.frame(obj))
      df_list <- lapply(1:length(levs), function(i)
        df_list[[i]] <- df_list[[i]] %>%
          dplyr::mutate(!!sym(varname) := factor(i, levels = 1:length(levs), labels = levs)))

      # Replicate original simulation dataframes across levels of byvar
      # and set byvar levels
      by_list <- lapply(1:length(bylevs), function(i) i <- df_list)
      i <- 1
      j <- 1
      while(j <= length(bylevs)) {
        if (i <= length(levs)) {
          by_list[[j]][[i]] <- by_list[[j]][[i]] %>%
            dplyr::mutate(!!sym(byvar) := bylevs[[j]])
          i <- i + 1
        } else {
          i <- 1
          j <- j + 1
        }
      }

      # Create model matrices
      X_list <- by_list
      i <- 1
      j <- 1
      while(j <= length(bylevs)) {
        if (i <= length(levs)) {
          X_list[[j]][[i]] <- model.matrix(formula(obj), data = by_list[[j]][[i]])[,-1]
          i <- i + 1
        } else {
          i <- 1
          j <- j + 1
        }
      }

      # Parametric bootstrap
      B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))
      tau_ind <- grep("\\|", names(coef(obj)))
      Tau <- B[,tau_ind]
      B <- B[,-tau_ind]

      # Generate predictions
      ByPr_list <- NULL
      ByPr_list <- lapply(1:length(bylevs), function(l) {
        CP <- NULL
        P <- NULL
        Pr_list <- NULL
        for(k in 1:length(levs)) {
          for(i in 1:ncol(Tau)) {
            CP[[i]] <- plogis(matrix(Tau[,i],
                                     ncol = nrow(Tau),
                                     nrow = nrow(X_list[[1]][[1]]),
                                     byrow = TRUE) - X_list[[l]][[k]] %*% t(B))
          }
          for (j in 1:ncol(Tau)) {
            if(j == 1) {
              P[[j]] <- CP[[j]]
            }
            else{
              P[[j]] <- CP[[j]]- CP[[(j - 1)]]
            }
            P[[(ncol(Tau) + 1)]] <- 1 - CP[[(ncol(Tau))]]
          }
          Pr_list[[k]] <- P
        }
        ByPr_list[[l]] <- Pr_list
      })

      # Weight predictions
      WgtByPr_list <- NULL
      for(j in 1:length(ByPr_list)) {
        Pr_list_ul <- NULL
        Pr_list_ul <- unlist(ByPr_list[[j]], recursive = FALSE)
        WgtPr_list <- NULL
        for (i in 1:length(Pr_list_ul)) {
          WgtPr_list[[i]] <- apply(Pr_list_ul[[i]], 2, function(x) weighted.mean(x, data$`(weights)`))
        }
        WgtByPr_list[[j]] <- WgtPr_list
      }

      # Assemble predictions
      preds <- NULL
      for(i in 1:length(WgtByPr_list)) {
        preds_tmp <- expand.grid(
          y = obj$lev,
          x = levs)
        preds_tmp$predicted <- sapply(WgtByPr_list[[i]], mean)
        preds_tmp$conf.low <- sapply(WgtByPr_list[[i]], low)
        preds_tmp$conf.high <- sapply(WgtByPr_list[[i]], high)
        preds_tmp$type <- "Probability"
        preds_tmp$z <- bylevs[[i]]
        preds[[i]] <- preds_tmp
      }
      preds <- do.call("rbind", preds)
      preds$z <- factor(round(preds$z, 3))
      preds <- preds %>%
        dplyr::relocate(z, .before = predicted) %>%
        dplyr::rename(!!sym(byvar) := z) %>%
        dplyr::rename(!!sym(varname) := x)

    }

    # CONTINUOUS * CONTINUOUS ==================================================
    # (same as continues * categorical, but first factorizes the moderator)

    if(is.numeric(data[[varname]]) & is.numeric(data[[byvar]])) {

      # Define byvar
      bylevs <- seq(min(data[[byvar]]), max(data[[byvar]]), length = bynvals)

      # Setup list of dataframes with levels of byvar set
      by_list <- lapply(1:length(bylevs), function(i) i <- model.frame(obj))
      by_list <- lapply(1:length(bylevs), function(i)
        by_list[[i]] <- by_list[[i]] %>%
          dplyr::mutate(!!sym(byvar) := bylevs[[i]]))

      # Define varname_seq
      varname_seq <- seq(min(data[[varname]]), max(data[[varname]]), length = nvals)

      # Loop simulations over levels of varname and across bylist
      res_list <- NULL
      # Note: for loops and lapply appear to be equivalent in terms of speed
      # for(j in 1:length(by_list)) { ... }
      # res_list <- lapply(1:length(by_list), function(j) { ... })
      for(j in 1:length(by_list)) {
        res <- NULL
        for(k in seq_along(varname_seq)) {
          new <- by_list[[j]] %>% dplyr::mutate(!!sym(varname) := varname_seq[k])
          Xmat <- model.matrix(formula(obj), data = new)[,-1]
          B <- MASS::mvrnorm(sims, coef(obj), vcov(obj))
          tau_ind <- grep("\\|", names(coef(obj)))
          Tau <- B[,tau_ind]
          B <- B[,-tau_ind]

          CP <- NULL
          for (i in 1:ncol(Tau)) {
            CP[[i]] <- plogis(matrix(Tau[,i],
                                     ncol = nrow(Tau),
                                     nrow = nrow(new),
                                     byrow = TRUE) - Xmat %*% t(B))
          }

          P <- NULL
          for (i in 1:ncol(Tau)) {
            if(i == 1) {P[[i]] <- CP[[i]]}
            else{P[[i]] <- CP[[i]] - CP[[(i - 1)]]}}
          P[[(ncol(Tau) + 1)]] <- 1 - CP[[ncol(Tau)]]

          WgtP <- lapply(1:length(P), function(i){
            apply(P[[i]], 2, function(x) weighted.mean(x, data$`(weights)`))
          })

          tmp <- data.frame(
            y = obj$lev,
            x = varname_seq[k],
            predicted = sapply(WgtP, mean),
            conf.low = sapply(WgtP, low),
            conf.high = sapply(WgtP, high),
            type = "Probability")

          res <- rbind(res, tmp)
        }
        res$z <- factor(round(bylevs, 3))[[j]]
        res_list[[j]] <- res
      }

      # Assemble table of predictions
      preds <- do.call("rbind", res_list)
      preds$y <- factor(preds$y, levels = obj$lev)
      preds <- preds %>%
        dplyr::relocate(z, .before = predicted) %>%
        dplyr::rename(!!sym(byvar) := z) %>%
        dplyr::rename(!!sym(varname) := x)

    }

    # OUTPUT ===================================================================

    output <- list(
      preds = tibble::as_tibble(preds),
      seed = seed,
      sims = sims,
      formula = formula(obj))
    class(output) <- "svyEffects"
    attributes(output)$predvar <- varname
    attributes(output)$byvar <- byvar
    attributes(output)$depvar <- colnames(obj$model)[1]
    attributes(output)$method <- "AME"
    return(output)

  }

}
