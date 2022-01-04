#' G-Estimation for an End of Study Outcome
#'
#' Performs g-estimation of a structural nested mean model (SNMM), based on the outcome regression methods described
#' in Sjolander and Vansteelandt (2016) and Dukes and Vansteelandt (2018). We expect a dataset that holds an end of study outcome that is either binary or continuous,
#' time-varying and/or baseline confounders, and a time-varying exposure that is either binary, continuous or categorical.
#'
#' @param data A data frame in long format containing the data to be analysed. See description for details.
#' @param idvar Character string specifying the name of the ID variable in the data.
#' @param timevar Character string specifying the name of the time variable in the data.
#'  Note that timevar must specify time periods as integer values starting from 1 (must not begin at 0).
#' @param Yn Character string specifying the name of the end of study outcome variable.
#' @param An Character string specifying the name of the time-varying exposure variable.
#' @param Cn Optional character string specifying the name of the censoring indicator variable. The variable specified in Cn should be a numeric vector taking values 0 or 1, with 1 indicating censored.
#' @param outcomemodels a list of formulas or formula objects specifying the outcome models for Yn prior to adjustment by propensity score.
#' The i'th entry of the list specifies the outcome model for the counterfactuals up to time i. See description for details.
#' @param propensitymodel A formula or formula object specifying the propensity score model for An.
#' @param censoringmodel A formula or formula object specifying the censoring model for Cn.
#' @param type Value from 1-4 specifying SNMM type to fit. See details.
#' @param EfmVar Character string specifying the name of the effect modifying variable for types 2 or 4.
#' @param ... Additional arguments, currently not in use.
#'
#' @details Given a time-varying exposure variable, \eqn{A_t} and time-varying confounders, \eqn{L_t} measured over time periods \eqn{t=1,\ldots,T}, and an end of study outcome \eqn{Y}
#' measured at time \eqn{T+1}, \code{gest} estimates the causal parameters \eqn{\psi} of a SNMM of the form
#'  \deqn{E(Y(\bar{a}_{t},0)-Y(\bar{a}_{t-1},0)|\bar{a}_{t-1},\bar{l}_{t})=\psi z_ta_t \;\forall\; t=1,\ldots,T}
#'  if Y is continuous or
#'  \deqn{\frac{E(Y(\bar{a}_{t},0)|\bar{a}_{t-1},\bar{l}_{t})}{E(Y(\bar{a}_{t-1},0)|\bar{a}_{t-1},\bar{l}_{t})}=exp(\psi z_ta_t)\;\forall\; t=1,\ldots,T }
#'  if Y is binary. The SNMMs form is defined by the parameter \eqn{z_t}, which can be controlled by the input \code{type} as follows
#'  \itemize{
#'  \item{\code{type=1} }{sets \eqn{z_t=1}. This implies that \eqn{\psi} is the effect of exposure at any time t on Y.}
#'  \item{\code{type=2} }{sets \eqn{z_t=c(1,l_t)}, and adds affect modification by \code{EfmVar}, which we denote \eqn{L_t}.
#'  Now \eqn{\psi=c(\psi_0,\psi_1)} where \eqn{\psi_0} is the effect of exposure at any time t on Y when \eqn{l_t=0} for all t, modified by
#'  \eqn{\psi_1} for each unit increase in \eqn{l_t} at all times t. Note that effect modification
#'  is currently only supported for binary (written as a numeric 0,1 vector) or continuous confounders.}
#'  \item {\code{type=3} }{allows for time-varying causal effects. It sets \eqn{z_t} to a vector of zeros of length T with a 1 in the t'th position. Now \eqn{\psi=c(\psi_1,\ldots,\psi_T)}
#'  where \eqn{\psi_t} is the effect of \eqn{A_t} on Y.}
#'  \item{\code{type=4} }{allows for a time-varying causal effect that can be modified by \code{EfmVar}, denoted \eqn{l_t}, that is it allows for both time-varying effects and effect modification. It sets \eqn{z_t} to a vector of zeros of length T with \eqn{c(1,l_t)} in the t'th position.
#'  Now \eqn{\psi=(\underline{\psi_1},\ldots,\underline{\psi_T})} where \eqn{\underline{\psi_t}=c(\psi_{0t},\psi_{1t})}. Here \eqn{\psi_{0t}} is the effect of exposure at time t on Y when \eqn{l_t=0} modified by
#'  \eqn{\psi_{1t}} for each unit increase in \eqn{l_t}. Note that effect modification
#'  is currently only supported for binary (written as a numeric 0,1 vector) or continuous confounders.}
#'  }
#'  The data must be in long format, where we assume the convention that each row with \code{time=t} contains \eqn{A_t,L_t} and \eqn{C_{t+1}} and \eqn{Y_{T+1}}. Thus the censoring indicator for each row
#'  should indicate that a user is censored AFTER time t. The end of study outcome \eqn{Y_{T+1}} should be repeated on each row. If either A or Y are binary, they must be written as numeric vectors taking values either 0 or 1.
#'  The same is true for any covariate that is used for effect modification.\cr
#'  The data must be rectangular with a row entry for every individual for each exposure time 1 up to T. Data rows after censoring should be empty apart from the ID and time variables. This can be done using the function \code{\link{FormatData}}.\cr
#'  The input outcomemodels should be a list with T elements (the number of exposure times), where element i describes the outcome model for the counterfactuals at time i.
#'
#' @return List of the fitted causal parameters of the posited SNMM. These are labeled as follows for each SNMM type, where \code{An} is
#' set to the name of the exposure variable, i is the current time period, and and EfmVar is the effect modifying variable.
#' \item{\code{type=1} }{\code{An}: The effect of exposure at any time t on outcome. }
#' \item{\code{type=2} }{\code{An}: The effect of exposure at any time t on outcome, when \code{EfmVar} is set to zero.\cr
#' \code{An:EfmVar}: The effect modification by \code{EfmVar}, the additional effect of A on Y for each unit increase in \code{EfmVar}}.
#' \item{\code{type=3} }{\code{t=i.An}: The effect of exposure at time t=i on outcome.}
#' \item{\code{type=4} }{\code{t=i.An}: The effect of exposure at time t=i on outcome, when \code{EfmVar} is set to zero.\cr
#' \code{t=i.An:EfmVar}: The effect modification by \code{EfmVar}, the additional effect of A on Y at time t=i for each unit increase in \code{EfmVar}.}
#' The function also returns a summary of the propensity scores and censoring scores via \code{PropensitySummary} and \code{CensoringSummary},
#' along with \code{Data}, holding the original dataset with the propensity and censoring scores as a tibble dataset.
#'
#'
#' @examples
#' datas <- dataexamples(n = 1000, seed = 123, Censoring = FALSE)
#' data <- datas$datagest
#' data <- FormatData(
#'   data = data, idvar = "id", timevar = "time", An = "A",
#'   varying = c("Y", "A", "L"), GenerateHistory = TRUE, GenerateHistoryMax = 1
#' )
#' idvar <- "id"
#' timevar <- "time"
#' Yn <- "Y"
#' An <- "A"
#' Cn <- NA
#' outcomemodels <- list("Y~A+L+U+Lag1A", "Y~A+L+U+Lag1A", "Y~A+L+U+Lag1A")
#' propensitymodel <- c("A~L+U+as.factor(time)+Lag1A")
#' censoringmodel <- NULL
#' EfmVar <- NA
#' gestSingle(data, idvar, timevar, Yn, An, Cn, outcomemodels, propensitymodel,
#' censoringmodel = NULL, type = 1, EfmVar)
#'
#' # Example with censoring
#' datas <- dataexamples(n = 1000, seed = 123, Censoring = TRUE)
#' data <- datas$datagest
#' data <- FormatData(
#'   data = data, idvar = "id", timevar = "time", An = "A", Cn = "C",
#'   varying = c("Y", "A", "L"), GenerateHistory = TRUE, GenerateHistoryMax = 1
#' )
#' Cn <- "C"
#' EfmVar <- "L"
#' outcomemodels <- list("Y~A+L+U+A:L+Lag1A", "Y~A+L+U+A:L+Lag1A", "Y~A+L+U+A:L")
#' censoringmodel <- c("C~L+U+as.factor(time)")
#' gestSingle(data, idvar, timevar, Yn, An, Cn, outcomemodels, propensitymodel,
#' censoringmodel = censoringmodel, type = 2, EfmVar)
#' @importFrom DataCombine slide
#' @importFrom geeM geem
#' @importFrom nnet multinom
#' @importFrom tidyr nest_legacy
#' @importFrom tibble as_tibble
#' @importFrom rsample bootstraps
#' @importFrom stats reshape rnorm rbinom gaussian glm predict reformulate
#' terms Gamma complete.cases quantile
#' @importFrom magrittr %>%
#' @importFrom tidyselect all_of
#'
#' @references Vansteelandt, S., & Sjolander, A. (2016). Revisiting g-estimation of the Effect of a Time-varying Exposure Subject to Time-varying Confounding, Epidemiologic Methods, 5(1), 37-56. <doi:10.1515/em-2015-0005>.
#' @references Dukes, O., & Vansteelandt, S. (2018). A Note on g-Estimation of Causal Risk Ratios, American Journal of Epidemiology, 187(5), 1079–1084. <doi:10.1093/aje/kwx347>.
#'
#' @export

gestSingle <- function(data, idvar, timevar, Yn, An, Cn = NA, outcomemodels, propensitymodel, censoringmodel = NULL, type, EfmVar = NA, ...) {
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


  data$int <- 1

  # Define propensity score model and relevant glm family
  lmp <- formula(propensitymodel)

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

  # Calculate initial censoring weights if Cn and LnC are supplied
  if (is.na(Cn) == TRUE) {
    data$w <- 1
  } else {
    lmc <- formula(censoringmodel)
    modc <- glm(lmc, family = "binomial", data = data)
    # Censoring scores
    cps <- 1 - predict(modc, type = "response", newdata = data)
    data$cps <- cps
    # Create Indicator Function of whether C=0
    data[, paste(Cn, "0", sep = "")] <- as.integer(!data[, Cn])
    # Setting up the initial denominator product of censoring scores "cprod" and initial censoring weights "w"
    data$cprod <- data$cps
    data[is.na(data$cprod) == TRUE, "cprod"] <- 1
    data$w <- data[, paste(Cn, "0", sep = "")] / data$cps
  }

  # Set up a variable "H" to hold the counterfactual outcome values
  # and set the value of H at the final exposure time T to the outcome
  data$H <- 0
  data[data[, timevar] == T, "H"] <- data[data[, timevar] == T, Yn]
  # Create a copy of the data set with complete cases for initial estimate of psi
  dcom <- data[complete.cases(data), ]

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

  # Get family for outcome model
  if (Ybin == TRUE) {
    family <- Gamma(link = "log")
  } else {
    family <- gaussian
  }

  # Obtain the correct formulas for the outcome model
  par1 <- paste(eval(An), eval(z), sep = ":")
  par1[par1 == paste(eval(An), "int", sep = ":")] <- paste(eval(An))
  par2 <- paste("prs", eval(z), sep = ":")
  par2[par2 == paste("prs", "int", sep = ":")] <- paste("prs")

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
    # Obtain the causal effect of A_T on Y as the initial estimate of psi
    lmy <- outcomemodels[[T]]
    dt <- dcom[dcom[, timevar] == T, ]
    out <- summary(geem(terms(lmy), family = family, id = dt[, idvar], data = dt, weights = dt$w))
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

    # Calculate Counterfactuals
    if (Acat == TRUE) {

      # Need to obtain psi^l*z_{sc}*A^l for all relevant data points before calculating H
      # call this variable 'psiZA'.
      # i=Miminum timeperiod for which 'psiZA' is being calculated
      # j=Current times period for which 'psiZA' is being calculated
      # l=category level of exposure variable
      i <- T - 1
      while (i >= 1) {
        j <- T
        data$psiZA <- 0
        while (j >= (i + 1)) {
          if (length(z) == 1) {
            for (l in 1:nlevels(data[, An])) {
              data[data[, timevar] == j & data[, An] == levels(data[, An])[l] & !is.na(data[, An]), "psiZA"] <- psicat[[l]]
            }
          } else {
            for (l in 1:nlevels(data[, An])) {
              data[data[, timevar] == j & data[, An] == levels(data[, An])[l] & !is.na(data[, An]), "psiZA"] <- rowSums(
                sweep(data[data[, timevar] == j & data[, An] == levels(data[, An])[l] & !is.na(data[, An]), z], 2, psicat[[l]], "*")
              )
            }
          }
          j <- j - 1
        }
        # Obtain counterfactual outcomes and censoring weights for time T-1
        # Set i as the minimum time period for which counterfactuals and censoring weights are calculated
        # Set k to the current time period being calculated.
        k <- (T - 1)
        while (k >= i) {
          if (is.na(Cn) == FALSE) {
            # Censoring weights
            data[data[, timevar] == k, "cprod"] <- data[data[, timevar] == k + 1, "cprod"] * data[data[, timevar] == k, "cps"]
            data[data[, timevar] == k, "w"] <- data[data[, timevar] == T, paste(Cn, "0", sep = "")] / data[data[, timevar] == k, "cprod"]
          }
          # Counterfactuals
          if (Ybin == FALSE) {
            data[data[, timevar] == k, "H"] <- data[data[, timevar] == (k + 1), "H"] -
              data[data[, timevar] == (k + 1), "psiZA"]
          }
          if (Ybin == TRUE) {
            data[data[, timevar] == k, "H"] <- data[data[, timevar] == (k + 1), "H"] *
              exp(-data[data[, timevar] == (k + 1), "psiZA"])
          }
          # Current time period is lowered and censoring weights and counterfactuals
          # calculated until current value of i is reached.
          k <- k - 1
        }
        # Obtain relevant time periods for which counterfactuals have been calculated
        dt <- data[data[, timevar] %in% seq(i, T, by = 1), ]
        dtcom <- dt[complete.cases(dt), ]
        lmH <- outcomemodels[[i]]
        out <- summary(geem(terms(lmH), family = family, id = dtcom[, idvar], data = dtcom, weights = dtcom$w))
        # Update the estimates of psi^l including data at time T-1.
        psi0 <- out$beta[match(Acoef, out$coefnames)]
        names(psi0) <- Acoef
        psicat <- as.list(NULL)
        # Update psicat
        for (l in 2:nlevels(data[, An])) {
          psicat[[l]] <- psi0[grep(levels(data[, An])[l], Acoef)]
        }
        psicat[[1]] <- rep(0, length(psicat[[2]]))
        # Set minimum time period to time T-2 and repeat calculation of counterfactuals
        # and censoring weights up to time T-2 and update psi. This is repeated until
        # time period 1 is reached and the final estimate of psi is obtained.
        i <- i - 1
      }
      results <- list(psi = unlist(psicat[-1]), Data = as_tibble(data[, !names(data) %in% c("H", "psiZA")]), PropensitySummary = summary(data$prs), CensoringSummary = summary(cps))
      class(results) <- "Results"

      return(results)
    } else {
      # Obtain counterfactual outcomes and censoring weights for time T-1
      # Set i as the minimum time period for which counterfactuals and censoring weights are calculated
      # Set k to the current time period being calculated.
      i <- T - 1
      while (i >= 1) {
        k <- T - 1
        while (k >= i) {
          if (is.na(Cn) == FALSE) {
            # Censoring weights
            data[data[, timevar] == k, "cprod"] <- data[data[, timevar] == k + 1, "cprod"] * data[data[, timevar] == k, "cps"]
            data[data[, timevar] == k, "w"] <- data[data[, timevar] == T, paste(Cn, "0", sep = "")] / data[data[, timevar] == k, "cprod"]
          }
          # Counterfactuals
          if (Ybin == FALSE) {
            if (length(z) == 1) {
              data[data[, timevar] == k, "H"] <- data[data[, timevar] == (k + 1), "H"] -
                data[data[, timevar] == (k + 1), An] * psi * data[data[, timevar] == (k + 1), z]
            } else {
              data[data[, timevar] == k, "H"] <- data[data[, timevar] == (k + 1), "H"] -
                rowSums(data[data[, timevar] == (k + 1), An] * sweep(data[data[, timevar] == (k + 1), z], 2, psi, "*"))
            }
          }
          if (Ybin == TRUE) {
            if (length(z) == 1) {
              data[data[, timevar] == k, "H"] <- data[data[, timevar] == (k + 1), "H"] *
                exp(-data[data[, timevar] == (k + 1), An] * psi * data[data[, timevar] == (k + 1), z])
            } else {
              data[data[, timevar] == k, "H"] <- data[data[, timevar] == (k + 1), "H"] *
                exp(-rowSums(data[data[, timevar] == (k + 1), An] * sweep(data[data[, timevar] == (k + 1), z], 2, psi, "*")))
            }
          }
          # Current time period is lowered and censoring weights and counterfactuals
          # calculated until current value of i is reached.
          k <- k - 1
        }
        # Obtain relevant time periods for which counterfactuals have been calculated
        dt <- data[data[, timevar] %in% seq(i, T, by = 1), ]
        dtcom <- dt[complete.cases(dt), ]
        lmH <- outcomemodels[[i]]
        # Update the estimate of psi including data at time T-1.
        out <- summary(geem(terms(lmH), family = family, id = dtcom[, idvar], data = dtcom, weights = dtcom$w))
        psi <- out$beta[match(par1, out$coefnames)]
        names(psi) <- par1
        # Set minimum time period to time T-2 and repeat calculation of counterfactuals
        # and censoring weights up to time T-2 and update psi. This is repeated until
        # time period 1 is reached and the final estimate of psi is obtained.
        i <- i - 1
      }
      results <- list(psi = psi, Data = as_tibble(data[, !names(data) %in% c("H")]), PropensitySummary = summary(data$prs), CensoringSummary = summary(cps))
      class(results) <- "Results"
      return(results)
    }
  }
  # Perform g-estimation for SNMM types 3 or 4 (time-varying psi)
  if (timevarying == TRUE) {
    if (Acat == TRUE) {
      dt <- dcom[dcom[, timevar] == T, ]
      lmy <- outcomemodels[[T]]
      out <- summary(geem(terms(lmy), family = family, id = dt[, idvar], data = dt, weights = dt$w))
      nam1 <- paste(An, levels(data[, An])[-1], sep = "")
      nam2 <- apply(expand.grid(nam1, z[-1]), 1, paste, collapse = ":")
      Acoef <- c(nam1, nam2)
      psi0 <- out$beta[match(Acoef, out$coefnames)]
      names(psi0) <- Acoef
      # Create psicatlist to hold the time-varying values of psicat for
      # each step length c and psicatresult to hold the same list
      # with the reference category estimates (set to 0) removed.
      psicat <- as.list(NULL)
      psicatlist <- as.list(NULL)
      psicatresult <- as.list(NULL)
      for (l in 2:nlevels(data[, An])) {
        psicat[[l]] <- psi0[grep(levels(data[, An])[l], Acoef)]
      }
      psicat[[1]] <- rep(0, length(psicat[[2]]))
      psicatlist[[T]] <- psicat
      psicatresult[[T]] <- psicat[-1]

      # Need to obtain psi^l*z_{sc}*A^l for all relevant data points before calculating H
      # call this variable 'psiZA'.
      # i=Miminum time period for which 'psiZA' is being calculated
      # j=Current time period for which 'psiZA' is being calculated
      # l=category level of exposure variable
      i <- (T - 1)
      while (i >= 1) {
        j <- T
        data$psiZA <- 0
        while (j >= (i + 1)) {
          if (length(z) == 1) {
            for (l in 1:nlevels(data[, An])) {
              data[data[, timevar] == j & data[, An] == levels(data[, An])[l] & !is.na(data[, An]), "psiZA"] <- psicatlist[[j]][[l]]
            }
          } else {
            for (l in 1:nlevels(data[, An])) {
              data[data[, timevar] == j & data[, An] == levels(data[, An])[l] & !is.na(data[, An]), "psiZA"] <- rowSums(
                sweep(data[data[, timevar] == j & data[, An] == levels(data[, An])[l] & !is.na(data[, An]), z], 2, psicatlist[[j]][[l]], "*")
              )
            }
          }
          j <- j - 1
        }
        # Obtain counterfactual outcomes and censoring weights for time T-1
        # Set i as the minimum time period for which counterfactuals and censoring weights are calculated
        # Set k to the current time period being calculated.
        k <- (T - 1)
        while (k >= i) {
          if (is.na(Cn) == FALSE) {
            data[data[, timevar] == k, "cprod"] <- data[data[, timevar] == k + 1, "cprod"] * data[data[, timevar] == k, "cps"]
            data[data[, timevar] == k, "w"] <- data[data[, timevar] == T, paste(Cn, "0", sep = "")] / data[data[, timevar] == k, "cprod"]
          }
          if (Ybin == FALSE) {
            data[data[, timevar] == k, "H"] <- data[data[, timevar] == (k + 1), "H"] -
              data[data[, timevar] == (k + 1), "psiZA"]
          }
          if (Ybin == TRUE) {
            data[data[, timevar] == k, "H"] <- data[data[, timevar] == (k + 1), "H"] *
              exp(-data[data[, timevar] == (k + 1), "psiZA"])
          }
          k <- k - 1
        }
        # Obtain data at time T-1 for calculation of psi
        dt <- data[data[, timevar] %in% i, ]
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
        psicatresult[[i]] <- psicat[-1]
        psicatlist[[i]] <- psicat
        # Set minimum time period to time T-2 and repeat calculation of counterfactuals
        # and censoring weights to time T-2 and update psi. This is repeated until
        # time period 1 is reached and the final estimate of psi is obtained.
        i <- i - 1
      }
      # Create relevant names for output and create results
      nam <- as.vector(NULL)
      for (p in 1:T) {
        nam[p] <- paste("t=", p, sep = "")
      }
      names(psicatresult) <- nam
      results <- list(psi = unlist(psicatresult), Data = as_tibble(data[, !names(data) %in% c("H", "psiZA")]), PropensitySummary = summary(data$prs), CensoringSummary = summary(cps))
      class(results) <- "Results"
      return(results)
    } else {
      dt <- dcom[dcom[, timevar] == T, ]
      lmy <- outcomemodels[[T]]
      out <- summary(geem(terms(lmy), family = family, id = dt[, idvar], data = dt, weights = dt$w))
      psi <- out$beta[match(par1, out$coefnames)]
      names(psi) <- par1
      psilist <- as.list(NULL)
      psilist[[T]] <- psi

      # Obtain counterfactuals and censoring weights for time T-1 as before
      # using the estimate of psi_T
      i <- (T - 1)
      while (i >= 1) {
        k <- (T - 1)
        while (k >= i) {
          if (is.na(Cn) == FALSE) {
            data[data[, timevar] == k, "cprod"] <- data[data[, timevar] == k + 1, "cprod"] * data[data[, timevar] == k, "cps"]
            data[data[, timevar] == k, "w"] <- data[data[, timevar] == T, paste(Cn, "0", sep = "")] / data[data[, timevar] == k, "cprod"]
          }
          if (Ybin == FALSE) {
            if (length(z) == 1) {
              data[data[, timevar] == k, "H"] <- data[data[, timevar] == (k + 1), "H"] -
                data[data[, timevar] == (k + 1), An] * psilist[[k + 1]] * data[data[, timevar] == (k + 1), z]
            } else {
              data[data[, timevar] == k, "H"] <- data[data[, timevar] == (k + 1), "H"] -
                rowSums(data[data[, timevar] == (k + 1), An] * sweep(data[data[, timevar] == (k + 1), z], 2, psilist[[k + 1]], "*"))
            }
          }
          if (Ybin == TRUE) {
            if (length(z) == 1) {
              data[data[, timevar] == k, "H"] <- data[data[, timevar] == (k + 1), "H"] *
                exp(-data[data[, timevar] == (k + 1), An] * psilist[[k + 1]] * data[data[, timevar] == (k + 1), z])
            } else {
              data[data[, timevar] == k, "H"] <- data[data[, timevar] == (k + 1), "H"] *
                exp(-rowSums(data[data[, timevar] == (k + 1), An] * sweep(data[data[, timevar] == (k + 1), z], 2, psilist[[k + 1]], "*")))
            }
          }
          k <- k - 1
        }
        # Obtain data at time T-1 for calculation of psi
        dt <- data[data[, timevar] %in% i, ]
        dtcom <- dt[complete.cases(dt), ]
        lmH <- outcomemodels[[i]]
        out <- summary(geem(terms(lmH), family = family, id = dtcom[, idvar], data = dtcom, weights = dtcom$w))
        psi <- out$beta[match(par1, out$coefnames)]
        names(psi) <- par1
        psilist[[i]] <- psi

        # Set minimum time period to time T-2 and repeat calculation of counterfactuals
        # and censoring weights to time T-2 and update psi. This is repeated until
        # time period 1 is reached and the final estimate of psi is obtained.
        i <- i - 1
      }
      # Create relevant names for output and create results
      nam <- as.vector(NULL)
      for (p in 1:T) {
        nam[p] <- paste("t=", p, sep = "")
      }
      names(psilist) <- nam
      results <- list(psi = unlist(psilist), Data = as_tibble(data[, !names(data) %in% c("H")]), PropensitySummary = summary(data$prs), CensoringSummary = summary(cps))
      class(results) <- "Results"
      return(results)
    }
  }
}
