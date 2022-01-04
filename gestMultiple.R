#' G-Estimation for a Time-Varying Outcome
#'
#' Performs g-estimation of a structural nested mean model (SNMM), based on the outcome regression methods described
#' in Sjolander and Vansteelandt (2016) and Dukes and Vansteelandt (2018). We assume a dataset with a time-varying outcome that is either binary or continuous,
#' time-varying and/or baseline confounders, and a time-varying exposure that is either binary, continuous or categorical.
#'
#' @param data A data frame in long format containing the data to be analysed. See description for details.
#' @param idvar Character string specifying the name of the ID variable in data.
#' @param timevar Character string specifying the name of the time variable in the data.
#'  Note that timevar must specify time periods as integer values starting from 1 (must not begin at 0).
#' @param Yn Character string specifying the name of the time-varying outcome variable.
#' @param An Character string specifying the name of the time-varying exposure variable.
#' @param Cn Optional character string specifying the name of the censoring indicator variable. The variable specified in Cn should be a numeric vector taking values 0 or 1, with 1 indicating censored.
#' @param outcomemodels a list of formulas or formula objects specifying the outcome models for Yn prior to adjustment by propensity score.
#' The i'th entry of the list specifies the outcome model for the i step counterfactuals. See description for details.
#' @param propensitymodel A formula or formula object specifying the propensity score model for An.
#' @param censoringmodel A formula or formula object specifying the censoring model for Cn.
#' @param type Value from 1-4 specifying SNMM type to fit. See details.
#' @param EfmVar Character string specifying the name of the effect modifying variable for types 2 or 4.
#' @param cutoff An integer taking value from 1 up to T, where T is the maximum value of timevar.
#' Instructs the function to estimate causal effects based only on exposures up to \code{cutoff} time periods prior to
#' the outcome.
#' @param ... Additional arguments, currently not in use.
#'
#' @details Suppose a series of time periods \eqn{1,\ldots,T+1} whereby a time-varying exposure and confounder (\eqn{A_t} and \eqn{L_t}) are measured over times \eqn{t=1,\ldots,T} and
#' a time varying outcome \eqn{Y_s} is measured over times \eqn{s=2,\ldots,T+1}. Define \eqn{c=s-t} as the step length, that is the number of time periods separating an exposure measurement, and subsequent outcome measurement.
#' By using the transform \eqn{t=s-c}, \code{gestmult} estimates the causal parameters \eqn{\psi} of a SNMM of the form
#'  \deqn{E\{Y_s(\bar{a}_{s-c},0)-Y_s(\bar{a}_{s-c-1},0)|\bar{a}_{s-c-1},\bar{l}_{s-c}\}=\psi z_{sc}a_{s-c} \; \forall c=1,\ldots,T\; and\; \forall s>c}
#'  if Y is continuous or
#'  \deqn{\frac{E(Y_s(\bar{a}_{s-c},0)|\bar{a}_{s-c-1},\bar{l}_{s-c})}{E(Y_s(\bar{a}_{s-c-1},0)|\bar{a}_{s-c-1},\bar{l}_{s-c})}=exp(\psi z_{sc}a_{s-c}) \; \forall c=1,\ldots,T\; and \; \forall s>c }
#'  if Y is binary. The SNMMs form is defined by the parameter \eqn{z_{sc}}, which can be controlled by the input \code{type} as follows
#'  \itemize{
#'  \item{\code{type=1} }{sets \eqn{z_{sc}=1}. This implies that \eqn{\psi} is now the effect of exposure at any time t on all subsequent outcome periods.}
#'  \item{\code{type=2} }{sets \eqn{z_{sc}=c(1,l_{s-c})} and adds affect modification by the variable named in \code{EfmVar}, which we denote \eqn{l_t}. Now \eqn{\psi=c(\psi_0,\psi_1)} where \eqn{\psi_0} is the effect of exposure at any time t on all subsequent outcome periods, when \eqn{l_{s-c}=0} at all times t, modified by
#'  \eqn{\psi_1} for each unit increase in \eqn{l_{s-c}} at all times t. Note that effect modification
#'  is currently only supported for binary or continuous confounders.}
#'  \item{\code{type=3} }{can posit a time-varying causal effect for each value of \eqn{c}, that is the causal effect for the exposure on outcome \eqn{c} time periods later.
#'  We set \eqn{z_{sc}} to a vector of zeros of length T with a 1 in the \eqn{c=s-t}'th position. Now \eqn{\psi=c(\psi_{1},\ldots,\psi_{T})}
#'  where \eqn{\psi_(c)} is the effect of exposure on outcome \eqn{c} time periods later for all outcome periods \eqn{s>c} that is \eqn{A_{s-c}} on \eqn{Y_s}.}
#'  \item{\code{type=4} }{allows for a time-varying causal effect that can be modified by \code{EfmVar}, denoted \eqn{l_t}, that is it allows for both time-varying effects and effect modification. It sets \eqn{z_{sc}} to a vector of zeros of length T with \eqn{c(1,l_{s-c})} in the \eqn{c=s-t}'th position.
#'  Now \eqn{\psi=(\underline{\psi_1},\ldots,\underline{\psi_T})} where \eqn{\underline{\psi_c}=c(\psi_{0c},\psi_{1c})}. Here \eqn{\psi_{0c}} is the effect of exposure on outcome \eqn{c} time periods later, given \eqn{l_{s-c}=0} for all \eqn{s>c}, modified by
#'  \eqn{\psi_{1c}} for each unit increase in \eqn{l_{s-c}} for all \eqn{s>c}. Note that effect modification
#'  is currently only supported for binary or continuous confounders.}
#'  }
#'  The data must be in long format, where we assume the convention that each row with \code{time=t} contains \eqn{A_t,L_t} and \eqn{C_{t+1},Y_{t+1}}. That is the censoring indicator for each row
#'  should indicate that a user is censored AFTER time t and the outcome indicates the first outcome that occurs AFTER \eqn{A_t} and \eqn{L_t} are measured.
#'  For example, data at time 1, should contain \eqn{A_1}, \eqn{L_1}, \eqn{Y_{2}}, and optionally \eqn{C_2}. If either A or Y are binary, they must be written as numeric vectors taking values either 0 or 1.
#'  The same is true for any covariate that is used for effect modification.\cr
#'  The data must be rectangular with a row entry for every individual for each exposure time 1 up to T. Data rows after censoring should be empty apart from the ID and time variables. This can be done using the function \code{\link{FormatData}}.\cr
#'  The input outcomemodels should be a list with T elements (the number of exposure times), where element i describes the outcome model for up to the i step counterfactual outcomes, that is the model is fitted to all counterfactuals up to \code{Y_{s-i}}.
#'
#'
#' @return List of the fitted causal parameters of the posited SNMM. These are labeled as follows for each SNMM type, where An is
#' set to the name of the exposure variable, i is the current value of c, and EfmVar is the effect modifying variable.
#' \item{\code{type=1} }{An: The effect of exposure at any time t on outcome at all subsequent times.}
#' \item{\code{type=2} }{An: The effect of exposure on outcome at any time t, when EfmVar is set to zero, on all subsequent outcome times.\cr
#' An:EfmVar: The effect modification by EfmVar, the additional effect of A on all subsequent Y for each unit increase in EfmVar at all times t. }
#' \item{\code{type=3} }{c=i.An: The effect of exposure at any time t on outcome \code{c=i} time periods later.}
#' \item{\code{type=4} }{c=i.An: The effect of exposure at any time t on outcome \code{c=i} time periods later, when EfmVar is set to zero.\cr
#' c=i.An:EfmVar: The effect modification by EfmVar, the additional effect of exposure on outcome c=i time periods later for each unit increase in EfmVar.}
#' The function also returns a summary of the propensity scores and censoring scores via \code{PropensitySummary} and \code{CensoringSummary},
#' along with \code{Data}, holding the original dataset with the propensity and censoring scores as a tibble dataset.
#' @references Vansteelandt, S., & Sjolander, A. (2016). Revisiting g-estimation of the Effect of a Time-varying Exposure Subject to Time-varying Confounding, Epidemiologic Methods, 5(1), 37-56. <doi:10.1515/em-2015-0005>.
#' @references Dukes, O., & Vansteelandt, S. (2018). A Note on g-Estimation of Causal Risk Ratios, American Journal of Epidemiology, 187(5), 1079–1084. <doi:10.1093/aje/kwx347>.
#'
#' @examples
#' datas <- dataexamples(n = 1000, seed = 123, Censoring = FALSE)
#' data <- datas$datagestmult
#' data <- FormatData(
#'   data = data, idvar = "id", timevar = "time", An = "A",
#'   varying = c("Y", "A", "L"), GenerateHistory = TRUE, GenerateHistoryMax = 1
#' )
#' idvar <- "id"
#' timevar <- "time"
#' Yn <- "Y"
#' An <- "A"
#' Cn <- NA
#' outcomemodels <- list("Y~A+L+U+Lag1A", "Y~A+L+U+Lag1A", "Y~A+L+U")
#' propensitymodel <- c("A~L+U+as.factor(time)+Lag1A")
#' censoringmodel <- NULL
#' EfmVar <- NA
#' gestMultiple(data, idvar, timevar, Yn, An, Cn, outcomemodels, propensitymodel,
#'   censoringmodel = NULL, type = 1, EfmVar,
#'   cutoff = NA
#' )
#'
#' # Example with censoring
#' datas <- dataexamples(n = 1000, seed = 123, Censoring = TRUE)
#' data <- datas$datagestmult
#' data <- FormatData(
#'   data = data, idvar = "id", timevar = "time", An = "A", Cn = "C",
#'   varying = c("Y", "A", "L"), GenerateHistory = TRUE, GenerateHistoryMax = 1
#' )
#' Cn <- "C"
#' EfmVar <- "L"
#' outcomemodels <- list("Y~A+L+U+A:L+Lag1A", "Y~A+L+U+A:L+Lag1A", "Y~A+L+U+A:L")
#' censoringmodel <- c("C~L+U+as.factor(time)")
#' gestMultiple(data, idvar, timevar, Yn, An, Cn, outcomemodels, propensitymodel,
#'   censoringmodel = censoringmodel, type = 2, EfmVar,
#'   cutoff = 2
#' )
#' @importFrom DataCombine slide
#' @importFrom geeM geem
#' @importFrom nnet multinom
#' @importFrom tidyr unnest_legacy
#' @importFrom tidyr nest_legacy
#' @importFrom tibble as_tibble
#' @importFrom rsample bootstraps
#' @importFrom stats reshape rnorm rbinom gaussian glm predict reformulate
#' terms Gamma complete.cases quantile formula
#' @import testthat
#'
#' @export

gestMultiple <- function(data, idvar, timevar, Yn, An, Cn = NA, outcomemodels, propensitymodel, censoringmodel = NULL, type, EfmVar = NA, cutoff = NA, ...) {
  # Error messages
  if (!is.data.frame(data)) (stop("Either no data set has been given, or it is not in a data frame."))
  if (is.na(EfmVar) && type %in% c(2, 4)) (stop("Type 2 or 4 is specified but argument EfmVar not specified."))
  if (!is.na(EfmVar) && !is.numeric(data[, EfmVar]) && type %in% c(2, 4)) (stop("Effect modification is only supported for a continuous covariate, or binary covariate written as an as.numeric() 0,1 vector"))
  if (!is.na(Cn) == TRUE && !is.numeric(data[, Cn])) (stop("A censoring indicator must be written as an as.numeric() 0,1 vector, with 1 indicating censoring."))
  if (!is.null(censoringmodel)) (warning("Variables included in censoringmodel should ideally be included in propensitymodel else propensity scores may be invalid."))
  if (!is.factor(data[, Yn]) && !is.numeric(data[, Yn])) (stop("Outcome Yn must be an as.numeric() continuous variable, or if binary, an as.numeric() 0 1 variable."))
  if (!is.factor(data[, Yn]) && !is.numeric(data[, Yn])) (stop("Exposure An must be either an as.factor() categorical variable, or an as.numeric() variable. If Binary, it must be set either as a two category as.factor() variable or a numeric 0 1 variable."))

  # Check if Yn and An are either binary, or  categorical
  Ybin <- FALSE
  Abin <- FALSE
  Acat <- FALSE
  if (setequal(unique(data[, Yn][!is.na(data[, Yn])]), c(0, 1)) && is.numeric(data[, Yn])) (Ybin <- TRUE)
  if (setequal(unique(data[, An][!is.na(data[, An])]), c(0, 1)) && is.numeric(data[, An])) (Abin <- TRUE)
  if (is.factor(data[, An])) (Acat <- TRUE)

  # Check that timevar is correctly specified, and there are the correct number of data rows

  if (!is.numeric(data[, timevar])) (stop("timevar must be as as.numeric() variable starting at 1"))
  if (is.na(min(data[, idvar]))) (stop("idvar must not contain any missing values"))
  if (min(data[, timevar]) != 1) (stop("timevar must be as as.numeric() variable starting at 1. It must also not contain any missing values"))
  if (nrow(data) != (length(unique(data[, idvar])) * max(data[, timevar]))) (stop("There must a a row entry for each individual at each time period. For those with entries missing or censored at a time point, add
                                                                       rows of missing values except for the time and id variable. Consider using the function FormatData."))


  T <- max(data[, timevar])
  if (is.na(cutoff) == TRUE) {
    cutoff <- T
  }

  data$int <- 1
  # Define propensity score model and relevant glm family if given
  lmp <- formula(propensitymodel)
  # Fit the propensity Score Model
  if (Acat == TRUE) {
    modp <- multinom(lmp, data = data)
  } else if (Abin == TRUE) {
    modp <- glm(lmp, family = "binomial", data = data)
  } else {
    modp <- glm(lmp, family = "gaussian", data = data)
  }


  # Calculate propensity scores
  if (Acat == TRUE) {
    props <- predict(modp, type = "probs", newdata = data)
    if (nlevels(data[, An]) == 2) {
      data$prs <- props
    } else {
      data$prs <- props[, -1]
    }
  } else {
    props <- predict(modp, type = "response", newdata = data)
    data$prs <- props
  }

  cps <- NA

  # Calculate initial censoring weights if Cn and censoring model are supplied
  if (is.na(Cn)) {
    data$w <- 1
  } else {
    lmc <- formula(censoringmodel)
    modc <- glm(lmc, family = "binomial", data = data)
    # Censoring scores
    cps <- 1 - predict(modc, type = "response", newdata = data)
    data$cps <- cps
    # Indicator Function of whether C=0
    data[, paste(Cn, "0", sep = "")] <- as.integer(!data[, Cn])
    # Setting up the denominator product of censoring scores "cprod" and initial weights "w"
    data$cprod <- data$cps
    data[is.na(data$cprod) == TRUE, "cprod"] <- 1
    data$w <- data[, paste(Cn, "0", sep = "")] / data$cps
  }

  # Set up counterfactual outcomes for exposure at the preceding time period
  # that is the "1 step" counterfactuals H_{s(s-1)}
  data$H <- data[, Yn]
  dc <- data
  dc$cntstep <- 1

  # Create a copy of the data set with complete cases for initial estimate of psi
  dcom <- data[complete.cases(data), ]

  # Set up an augmented dataset 'dc', adding additonal rows to hold the
  # 2-step, up to c-step counterfactuals, with c=cutoff, that is the counterfactuals for time
  # s-c under outcome s H_{s(s-c)}. Generate variable 'cntstep' to indicate the step length
  # of the counterfactual the row holds.
  for (i in 2:cutoff) {
    d2 <- data[data[, timevar] %in% seq(1, T - (i - 1), by = 1), ]
    d2$cntstep <- i
    dc <- rbind(dc, d2)
  }
  dc <- dc[order(dc[, idvar], dc[, timevar]), ]

  # Define values of internal variables 'z' and 'timevarying' based on input 'type'
  if (type == 1) {
    z <- c("int")
    timevarying <- FALSE
  } else if (type == 2) {
    z <- c("int", EfmVar)
    timevarying <- FALSE
  } else if (type == 3) {
    z <- c("int")
    timevarying <- TRUE
  } else if (type == 4) {
    z <- c("int", EfmVar)
    timevarying <- TRUE
  }

  # Add the relevant elements of the propensity score, that is psi*z to the outcome models
  par1 <- paste(eval(An), eval(z), sep = ":")
  par1[par1 == paste(eval(An), "int", sep = ":")] <- paste(eval(An))
  par2 <- paste("prs", eval(z), sep = ":")
  par2[par2 == paste("prs", "int", sep = ":")] <- paste("prs")




  # Get correct family for outcome models
  if (Ybin == TRUE) {
    family <- Gamma(link = "log")
  } else {
    family <- gaussian
  }

  for (i in 1:length(outcomemodels)) {
    outcomemodels[[i]] <- formula(outcomemodels[[i]])
    termlabs <- attr(terms(outcomemodels[[i]]), which = "term.labels")
    if (identical(as.numeric(length(which(termlabs == An))), 0)) {
      stop("Every formula in outcomemodels must have an An term")
    }
    if (type %in% c(2, 4)) {
      # Ensure that the An:EfmVar effect is displayed in the correct way

      if (identical(as.numeric(length(which(termlabs == paste(An, EfmVar, sep = ":")))), 0)) {
        stop("For types 2 and 4. Every formula in outcomemodels must have an An term, Efm term, and
                                 an An:Efm term. The An term must appear before any EfmVar term in each formula in outcomemodels.
                                 Or there must be an An*EfmVar term")
      }
      if (!identical(as.numeric(length(which(termlabs == EfmVar))), 0) && which(termlabs == EfmVar) < which(termlabs == An)) (stop("For types 2 and 4. Every formula in outcomemodels must have an An term, Efm term, and
                                  an An:Efm term. The An term must appear before any EfmVar term in each formula in outcomemodels. Or there must be an An*EfmVar term"))
    }

    outcomemodels[[i]] <- reformulate(c(termlabs, par2), response = "H")
  }


  # Perform g-estimation for SNMM types 1 or 2 (non time-varying psi)
  if (timevarying == FALSE) {
    # Obtain the "1-step" causal effect as the initial estimate of psi
    # that is the effect of A_{s-1} on Y_s.
    lmy <- outcomemodels[[1]]
    out <- summary(geem(terms(lmy), family = family, id = dcom[, idvar], data = dcom, weights = dcom$w))
    if (Acat == TRUE) {
      nam1 <- paste(An, levels(data[, An])[-1], sep = "")
      nam2 <- apply(expand.grid(nam1, z[-1]), 1, paste, collapse = ":")
      Acoef <- c(nam1, nam2)
      psi0 <- out$beta[match(Acoef, out$coefnames)]
      names(psi0) <- Acoef
      # Place the estimate of psi for each non reference category of A (psi^l)
      # into a separate element in list psicat
      psicat <- as.list(NULL)
      for (l in 2:nlevels(data[, An])) {
        psicat[[l]] <- psi0[grep(levels(data[, An])[l], Acoef)]
      }
      # set the value of psi at the reference level to zero
      psicat[[1]] <- rep(0, length(psicat[[2]]))
    } else {
      psi <- out$beta[match(par1, out$coefnames)]
      names(psi) <- par1
    }

    if (Acat == TRUE) {
      # Need to obtain psi^l*z_{sc}*A^l for all relevant data points before calculating H
      # call this variable 'psiZA'.
      # i=Maximum value of 'cntstep' for which 'psiZA' is being calculated
      # j=Current value of 'cntstep' for which 'psiZA' is being calculated
      # l<-category level of exposure variable
      i <- 2
      while (i <= cutoff && i <= T) {
        j <- 1
        dc$psiZA <- 0
        while (j <= (i - 1)) {
          if (length(z) == 1) {
            for (l in 1:nlevels(dc[, An])) {
              dc[dc$cntstep == j & dc[, An] == levels(dc[, An])[l] & !is.na(dc[, An]), "psiZA"] <- psicat[[l]]
            }
          } else {
            for (l in 1:nlevels(dc[, An])) {
              dc[dc$cntstep == j & dc[, An] == levels(dc[, An])[l] & !is.na(dc[, An]), "psiZA"] <- rowSums(
                sweep(dc[dc$cntstep == j & dc[, An] == levels(dc[, An])[l] & !is.na(dc[, An]), z], 2, psicat[[l]], "*")
              )
            }
          }
          j <- j + 1
        }
        # Obtain the 2-step counterfactuals and associated censoring weights
        # Set i as the maximum step length for which counterfactuals may be calculated
        # Set j as the current step length for which counterfactuals are being calculated, that is H_{s(s-j)}
        # set k as the current time period for which the j-step counterfactual is being calculated
        j <- 2
        while (j <= i) {
          for (k in 1:(T - (j - 1))) {
            # Censoring Weights
            if (is.na(Cn) == FALSE) {
              dc[dc$cntstep == j & dc[, timevar] == k, "cprod"] <- dc[dc$cntstep == (j - 1) & dc[, timevar] == k, "cps"] *
                dc[dc$cntstep == (j - 1) & dc[, timevar] == (k + 1), "cprod"]
              dc[dc$cntstep == j & dc[, timevar] == k, "w"] <- data[data[, timevar] == (k + j - 1), paste(Cn, "0", sep = "")] / dc[dc$cntstep == j & dc[, timevar] == k, "cprod"]
            }

            # Counterfactuals
            if (Ybin == FALSE) {
              dc[dc$cntstep == j & dc[, timevar] == k, "H"] <- dc[dc$cntstep == (j - 1) & dc[, timevar] == (k + 1), "H"] -
                dc[dc$cntstep == (j - 1) & dc[, timevar] == (k + 1), "psiZA"]
            }
            if (Ybin == TRUE) {
              dc[dc$cntstep == j & dc[, timevar] == k, "H"] <- dc[dc$cntstep == (j - 1) & dc[, timevar] == (k + 1), "H"] *
                exp(-dc[dc$cntstep == (j - 1) & dc[, timevar] == (k + 1), "psiZA"])
            }
          }
          # Increase j and calculate the counterfactuals of the next step length, until the maximum step
          # length i is reached
          j <- j + 1
        }
        # Obtain relevant rows for which counterfactuals have been calculated, that is
        # initially the 1 and 2 step counterfactuals
        dt <- dc[dc$cntstep %in% seq(1, i, by = 1), ]
        dtcom <- dt[complete.cases(dt), ]
        lmH <- outcomemodels[[i]]
        # Update the estimates of psi^l including both the 1 step and 2 step counterfactuals
        out <- summary(geem(terms(lmH), family = family, id = dtcom[, idvar], data = dtcom, weights = dtcom$w))
        psi0 <- out$beta[match(Acoef, out$coefnames)]
        names(psi0) <- Acoef
        # Update psicat
        psicat <- as.list(NULL)
        for (l in 2:nlevels(data[, An])) {
          psicat[[l]] <- psi0[grep(levels(data[, An])[l], Acoef)]
        }
        psicat[[1]] <- rep(0, length(psicat[[2]]))
        # Increase maximum step length by and repeat calculation of counterfactuals
        # up to step length 3 (and relevant censoring weights), and update psi. This is repeated until
        # either step length T or step length equal to cutoff is reached and the final estimate of psi is obtained.
        i <- i + 1
      }
      # Output results
      results <- list(psi = unlist(psicat[-1]), Data = as_tibble(data[, !names(data) %in% c("H", "psiZA")]), PropensitySummary = summary(data$prs), CensoringSummary = summary(cps))
      class(results) <- "Results"
      return(results)
    } else {
      # Obtain the 2-step counterfactuals and associated censoring weights
      # Set i as the maximum step length for which counterfactuals may be calculated
      # Set j as the current step length for which counterfactuals are being calculated, that is H_{s(s-j)}
      # set k as the current time period for which the j-step counterfactual is being calculated
      i <- 2
      while (i <= cutoff && i <= T) {
        j <- 2
        while (j <= i) {
          for (k in 1:(T - (j - 1))) {
            if (is.na(Cn) == FALSE) {
              # Censoring weights
              dc[dc$cntstep == j & dc[, timevar] == k, "cprod"] <- dc[dc$cntstep == (j - 1) & dc[, timevar] == k, "cps"] *
                dc[dc$cntstep == (j - 1) & dc[, timevar] == (k + 1), "cprod"]
              dc[dc$cntstep == j & dc[, timevar] == k, "w"] <- data[data[, timevar] == (k + j - 1), paste(Cn, "0", sep = "")] / dc[dc$cntstep == j & dc[, timevar] == k, "cprod"]
            }

            # Counterfactuals
            if (Ybin == FALSE) {
              if (length(z) == 1) {
                dc[dc$cntstep == j & dc[, timevar] == k, "H"] <- dc[dc$cntstep == (j - 1) & dc[, timevar] == (k + 1), "H"] -
                  dc[dc$cntstep == (j - 1) & dc[, timevar] == (k + 1), An] * psi * data[data[, timevar] == (k + 1), z]
              } else {
                dc[dc$cntstep == j & dc[, timevar] == k, "H"] <- dc[dc$cntstep == (j - 1) & dc[, timevar] == (k + 1), "H"] -
                  rowSums(dc[dc$cntstep == (j - 1) & dc[, timevar] == (k + 1), An] * sweep(dc[dc$cntstep == (j - 1) & dc[, timevar] == (k + 1), z], 2, psi, "*"))
              }
            }
            if (Ybin == TRUE) {
              if (length(z) == 1) {
                dc[dc$cntstep == j & dc[, timevar] == k, "H"] <- dc[dc$cntstep == (j - 1) & dc[, timevar] == (k + 1), "H"] *
                  exp(dc[dc$cntstep == (j - 1) & dc[, timevar] == (k + 1), An] * -psi * data[data[, timevar] == (k + 1), z])
              } else {
                dc[dc$cntstep == j & dc[, timevar] == k, "H"] <- dc[dc$cntstep == (j - 1) & dc[, timevar] == (k + 1), "H"] *
                  exp(-rowSums(dc[dc$cntstep == (j - 1) & dc[, timevar] == (k + 1), An] * sweep(dc[dc$cntstep == (j - 1) & dc[, timevar] == (k + 1), z], 2, psi, "*")))
              }
            }
          }
          # Increase j and calculate the counterfactuals of the next step length, until the maximum step
          # length i is reached
          j <- j + 1
        }
        # Obtain relevant rows for which counterfactuals have been calculated, that is
        # initially the 1 and 2 step counterfactuals
        dt <- dc[dc$cntstep %in% seq(1, i, by = 1), ]
        dtcom <- dt[complete.cases(dt), ]
        # Update the estimate of psi including both the 1 step and 2 step counterfactuals
        lmH <- outcomemodels[[i]]
        out <- summary(geem(terms(lmH), family = family, id = dtcom[, idvar], data = dtcom, weights = dtcom$w))
        psi <- out$beta[match(par1, out$coefnames)]
        names(psi) <- par1
        # Increase maximum step length by and repeat calculation of counterfactuals
        # up to step length 3 (and relevant censoring weights), and update psi. This is repeated until
        # either step length T or step length equal to cutoff is reached and the final estimate of psi is obtained.
        i <- i + 1
      }
      results <- list(psi = psi, Data = as_tibble(data[, !names(data) %in% c("H")]), PropensitySummary = summary(data$prs), CensoringSummary = summary(cps))
      class(results) <- "Results"
      return(results)
    }
  }

  # Perform g-estimation for SNMM types 3 or 4 (time-varying psi)
  else if (timevarying == TRUE) {
    if (Acat == TRUE) {
      # Obtain initial estimates of psi^l_{s-1}, the causal effects of exposure
      # on the subsequent outcome.
      lmy <- outcomemodels[[1]]
      out <- summary(geem(terms(lmy), family = family, id = dcom[, idvar], data = dcom, weights = dcom$w))
      nam1 <- paste(An, levels(data[, An])[-1], sep = "")
      nam2 <- apply(expand.grid(nam1, z[-1]), 1, paste, collapse = ":")
      Acoef <- c(nam1, nam2)
      psi0 <- out$beta[match(Acoef, out$coefnames)]
      names(psi0) <- Acoef
      # Create psicatlist to hold the time-varying values of psicat for
      # each step length c and psicatresult to hold the same list
      # with the reference category estimates (set to 0) removed.
      psicatlist <- as.list(NULL)
      psicatresult <- as.list(NULL)
      psicat <- as.list(NULL)
      for (l in 2:nlevels(data[, An])) {
        psicat[[l]] <- psi0[grep(levels(data[, An])[l], Acoef)]
      }
      psicat[[1]] <- rep(0, length(psicat[[2]]))
      psicatlist[[1]] <- psicat
      psicatresult[[1]] <- psicat[-1]

      i <- 2
      # Need to obtain psi^l*z_{sc}*A^l for all relevant data points before calculating H
      # call this variable 'psiZA'.
      # i=Maximum value of 'cntstep' for which 'psiZA' is being calculated
      # j=Current value of 'cntstep' for which 'psiZA' is being calculated
      # l<-category level of exposure variable
      while (i <= cutoff && i <= T) {
        j <- 1
        dc$psiZA <- 0
        while (j <= (i - 1)) {
          if (length(z) == 1) {
            for (l in 1:nlevels(dc[, An])) {
              dc[dc$cntstep == j & dc[, An] == levels(dc[, An])[l] & !is.na(dc[, An]), "psiZA"] <- psicatlist[[j]][[l]]
            }
          } else {
            for (l in 1:nlevels(dc[, An])) {
              dc[dc$cntstep == j & dc[, An] == levels(dc[, An])[l] & !is.na(dc[, An]), "psiZA"] <- rowSums(
                sweep(dc[dc$cntstep == j & dc[, An] == levels(dc[, An])[l] & !is.na(dc[, An]), z], 2, psicatlist[[j]][[l]], "*")
              )
            }
          }
          j <- j + 1
        }

        # Obtain the 2-step counterfactuals and associated censoring weights
        # Set i as the maximum step length for which counterfactuals may be calculated
        # Set j as the current step length for which counterfactuals are being calculated, that is H_{s(s-j)}
        # set k as the current time period for which the j-step counterfactual is being calculated
        j <- 2
        while (j <= i) {
          for (k in 1:(T - (j - 1))) {
            # Censoring Weights
            if (is.na(Cn) == FALSE) {
              dc[dc$cntstep == j & dc[, timevar] == k, "cprod"] <- dc[dc$cntstep == (j - 1) & dc[, timevar] == k, "cps"] *
                dc[dc$cntstep == (j - 1) & dc[, timevar] == (k + 1), "cprod"]
              dc[dc$cntstep == j & dc[, timevar] == k, "w"] <- data[data[, timevar] == (k + j - 1), paste(Cn, "0", sep = "")] / dc[dc$cntstep == j & dc[, timevar] == k, "cprod"]
            }
            # Counterfactuals
            if (Ybin == FALSE) {
              dc[dc$cntstep == j & dc[, timevar] == k, "H"] <- dc[dc$cntstep == (j - 1) & dc[, timevar] == (k + 1), "H"] -
                dc[dc$cntstep == (j - 1) & dc[, timevar] == (k + 1), "psiZA"]
            }
            if (Ybin == TRUE) {
              dc[dc$cntstep == j & dc[, timevar] == k, "H"] <- dc[dc$cntstep == (j - 1) & dc[, timevar] == (k + 1), "H"] *
                exp(-dc[dc$cntstep == (j - 1) & dc[, timevar] == (k + 1), "psiZA"])
            }
          }
          j <- j + 1
        }
        # Obtain relevant data and obtain an estimates of psi^l_{s-2}
        dt <- dc[dc$cntstep %in% i, ]
        dtcom <- dt[complete.cases(dt), ]
        lmH <- outcomemodels[[i]]
        out <- summary(geem(terms(lmH), family = family, id = dtcom[, idvar], data = dtcom, weights = dtcom$w))
        psi0 <- out$beta[match(Acoef, out$coefnames)]

        names(psi0) <- Acoef
        psicat <- as.list(NULL)
        for (l in 2:nlevels(data[, An])) {
          psicat[[l]] <- psi0[grep(levels(data[, An])[l], Acoef)]
        }
        psicat[[1]] <- rep(0, length(psicat[[2]]))
        psicatlist[[i]] <- psicat
        psicatresult[[i]] <- psicat[-1]
        # Set maximum step length up by one and repeat algorithm until
        # psi_{s-T} or psi_{s-c} where c is set to the value of cutoff is calculated.
        i <- i + 1
      }
      # Obtain relevant names for output
      nam <- as.vector(NULL)
      for (p in 1:cutoff) {
        nam[p] <- paste("c=", p, sep = "")
      }
      names(psicatresult) <- nam
      results <- list(psi = unlist(psicatresult), Data = as_tibble(data[, !names(data) %in% c("H", "psiZA")]), PropensitySummary = summary(data$prs), CensoringSummary = summary(cps))
      class(results) <- "Results"
      return(results)
    } else {
      lmy <- outcomemodels[[1]]
      out <- summary(geem(terms(lmy), family = family, id = dcom[, idvar], data = dcom, weights = dcom$w))
      psi <- out$beta[match(par1, out$coefnames)]
      names(psi) <- par1
      psilist <- as.list(NULL)
      psilist[[1]] <- psi
      # Obtain counterfactuals and censoring weights as before
      # where each j step counterfactual is estimated using the estimate of
      # psi_{s-j+1}
      i <- 2
      while (i <= cutoff && i <= T) {
        j <- 2
        while (j <= i) {
          for (k in 1:(T - (j - 1))) {
            if (is.na(Cn) == FALSE) {
              dc[dc$cntstep == j & dc[, timevar] == k, "cprod"] <- dc[dc$cntstep == (j - 1) & dc[, timevar] == k, "cps"] *
                dc[dc$cntstep == (j - 1) & dc[, timevar] == (k + 1), "cprod"]
              dc[dc$cntstep == j & dc[, timevar] == k, "w"] <- data[data[, timevar] == (k + j - 1), paste(Cn, "0", sep = "")] / dc[dc$cntstep == j & dc[, timevar] == k, "cprod"]
            }
            if (Ybin == FALSE) {
              if (length(z) == 1) {
                dc[dc$cntstep == j & dc[, timevar] == k, "H"] <- dc[dc$cntstep == (j - 1) & dc[, timevar] == (k + 1), "H"] -
                  dc[dc$cntstep == (j - 1) & dc[, timevar] == (k + 1), An] * psilist[[j - 1]] * data[data[, timevar] == (k + 1), z]
              } else {
                dc[dc$cntstep == j & dc[, timevar] == k, "H"] <- dc[dc$cntstep == (j - 1) & dc[, timevar] == (k + 1), "H"] -
                  rowSums(dc[dc$cntstep == (j - 1) & dc[, timevar] == (k + 1), An] * sweep(dc[dc$cntstep == (j - 1) & dc[, timevar] == (k + 1), z], 2, psilist[[j - 1]], "*"))
              }
            }
            if (Ybin == TRUE) {
              if (length(z) == 1) {
                dc[dc$cntstep == j & dc[, timevar] == k, "H"] <- dc[dc$cntstep == (j - 1) & dc[, timevar] == (k + 1), "H"] *
                  exp(dc[dc$cntstep == (j - 1) & dc[, timevar] == (k + 1), An] * -psilist[[j - 1]] * data[data[, timevar] == (k + 1), z])
              } else {
                dc[dc$cntstep == j & dc[, timevar] == k, "H"] <- dc[dc$cntstep == (j - 1) & dc[, timevar] == (k + 1), "H"] *
                  exp(-rowSums(dc[dc$cntstep == (j - 1) & dc[, timevar] == (k + 1), An] * sweep(dc[dc$cntstep == (j - 1) & dc[, timevar] == (k + 1), z], 2, psilist[[j - 1]], "*")))
              }
            }
          }

          j <- j + 1
        }
        # Obtain relevant data and calculate psi_{s-2}
        dt <- dc[dc$cntstep %in% i, ]
        dtcom <- dt[complete.cases(dt), ]
        lmH <- outcomemodels[[i]]
        out <- summary(geem(terms(lmH, keep.order = T), family = family, id = dtcom[, idvar], data = dtcom, weights = dtcom$w))
        psi <- out$beta[match(par1, out$coefnames)]
        names(psi) <- par1
        psilist[[i]] <- psi

        # Set maximum step length up by one and repeat algorithm until
        # psi_{s-T} or psi_{s-c} where c is set to the value of cutoff is calculated.
        i <- i + 1
      }
      # Create relevant names for outputs
      nam <- as.vector(NULL)
      for (p in 1:cutoff) {
        nam[p] <- paste("c=", p, sep = "")
      }
      names(psilist) <- nam
      results <- list(psi = unlist(psilist), Data = as_tibble(data[, !names(data) %in% c("H")]), PropensitySummary = summary(data$prs), CensoringSummary = summary(cps))
      class(results) <- "Results"
      return(results)
    }
  }
}
