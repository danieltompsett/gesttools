#' Percentile Based Bootstrap Confidence Intervals
#'
#' Generates percentile based confidence intervals for the causal parameters
#' of a fitted SNMM. Bonferroni corrected confidence
#' intervals are also reported for multiple comparisons.
#'
#' @param gestfunc Name (without quotations) of the g-estimation function to run.
#' One of \code{gestSingle} or \code{gestMultiple}.
#' @param data,idvar,timevar,Yn,An,Cn,outcomemodels,propensitymodel,censoringmodel,type,EfmVar,cutoff
#' Same arguments as in gest functions, to be input into gestfunc.
#' @param bn Number of bootstrapped datasets.
#' @param alpha Confidence level of confidence intervals.
#' @param onesided Controls the type of confidence interval generated. Takes one of three inputs, \code{"upper"} for upper one-sided confidence intervals,
#' \code{"lower"} for lower one-sided confidence intervals, and \code{"twosided"} for two-sided confidence intervals. Defaults to \code{"twosided"}.
#' @param seed Integer specifying the random seed for generation of bootstrap samples.
#' @param ... additional arguments.
#'
#' @return Returns a list of the following four elements.
#' \item{original }{The value of the causal parameters estimated on the original data \code{data}.}
#' \item{mean.boot }{The average values of the causal parameters estimated on the bootstrapped datasets.}
#' \item{conf }{The upper and/or lower bounds of \eqn{1-\alpha} confidence intervals for each element of \eqn{\psi}.
#' For example, if \code{type=2}, and \eqn{\psi=(\psi_0,\psi_1)}, a separate confidence interval is fitted for \eqn{\psi_0} and \eqn{\psi_1}.}
#' \item{conf.Bonferroni }{The upper and/or lower bounds of Bonferroni corrected confidence
#' intervals for \eqn{\psi}, used for multiple comparisons.}
#' \item{boot.results}{A tibble containing the result for each bootstrapped dataset}
#'
#' @examples
#' datas <- dataexamples(n = 1000, seed = 123, Censoring = FALSE)
#' data <- datas$datagest
#' data <- FormatData(
#'   data = data, idvar = "id", timevar = "time", An = "A",
#'   varying = c("A", "L"), GenerateHistory = TRUE, GenerateHistoryMax = 1
#' )
#' idvar <- "id"
#' timevar <- "time"
#' Yn <- "Y"
#' An <- "A"
#' Cn <- NA
#' outcomemodels <- list("Y~A+L+U+Lag1A", "Y~A+L+U+Lag1A", "Y~A+L+U+Lag1A")
#' propensitymodel <- c("A~L+U+as.factor(time)+Lag1A")
#' censoringmodel <- NULL
#' type <- 1
#' EfmVar <- NA
#' bn <- 5
#' alpha <- 0.05
#' gestfunc <- gestSingle
#' gestboot(gestfunc, data, idvar, timevar, Yn, An, Cn, outcomemodels, propensitymodel,
#'   censoringmodel = NULL, type = 1, EfmVar,
#'   bn = bn, alpha = alpha, onesided = "twosided", seed = 123
#' )
#' @export
gestboot <- function(gestfunc, data, idvar, timevar, Yn, An, Cn, outcomemodels, propensitymodel, censoringmodel = NULL, type, EfmVar = NA,
                     cutoff, bn, alpha = 0.05, onesided = "twosided", seed = NULL, ...) {
  t0 <- gestfunc(data = data, idvar = idvar, timevar = timevar, Yn = Yn, An = An, Cn = Cn, outcomemodels = outcomemodels, propensitymodel = propensitymodel, censoringmodel = censoringmodel, type = type, EfmVar = EfmVar, cutoff = cutoff)$psi
  nams <- names(t0)
  # Create tibble data based on ID
  data <- data[order(data[, idvar], data[, timevar]), ]
  Data <- data %>% nest_legacy(-all_of(idvar))
  if (!is.null(seed)) set.seed(seed)
  bs <- bootstraps(Data, times = bn)


  results1 <- as.list(NULL)
  results <- as.list(NULL)
  for (j in 1:bn) {
    tryCatch(
      {
        b <- as.data.frame(as_tibble(bs$splits[[j]]) %>% unnest_legacy())
        b[, idvar] <- rep(1:length(unique(data[, idvar])), each = max(b[, timevar]))
        results1[[j]] <- gestfunc(data = b, idvar = idvar, timevar = timevar, Yn = Yn, An = An, Cn = Cn, outcomemodels = outcomemodels, propensitymodel = propensitymodel, censoringmodel = censoringmodel, type = type, EfmVar = EfmVar, cutoff = cutoff)$psi
        results[[j]] <- unlist(results1[[j]])
      },
      error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
      }
    )
  }

  if (length(unlist(results)) < bn) (warning("One or more bootstrapped datasets failed to obtain a fitted causal parameter. Consider removing terms from Lny to avoid collinearity, or assess the sparseness of the data."))

  mean <- colMeans(do.call(rbind, results))

  resultsmat <- do.call(rbind, results)

  ci.quant <- function(x = NA) {
    return(quantile(x, probs = c(alpha / 2, 1 - alpha / 2)))
  }

  ci.quant.bonf <- function(x = NA) {
    return(quantile(x, probs = c(alpha / (2 * length(unlist(t0))), 1 - alpha / (2 * length(unlist(t0))))))
  }

  ci.quant.upper <- function(x = NA) {
    return(quantile(x, probs = c(1 - alpha)))
  }

  ci.quant.bonf.upper <- function(x = NA) {
    return(quantile(x, probs = c(1 - alpha / (length(unlist(t0))))))
  }

  ci.quant.lower <- function(x = NA) {
    return(quantile(x, probs = c(alpha)))
  }

  ci.quant.bonf.lower <- function(x = NA) {
    return(quantile(x, probs = c(alpha / (length(unlist(t0))))))
  }


  results.sort <- apply(resultsmat, 2, sort)

  if (onesided == "twosided") {
    conf.quant <- t(apply(results.sort, 2, ci.quant))
    conf.quant.bonf <- t(apply(results.sort, 2, ci.quant.bonf))
  } else if (onesided == "upper") {
    conf.quant <- t(apply(results.sort, 2, ci.quant.upper))
    conf.quant.bonf <- t(apply(results.sort, 2, ci.quant.bonf.upper))
  } else if (onesided == "lower") {
    conf.quant <- t(apply(results.sort, 2, ci.quant.lower))
    conf.quant.bonf <- t(apply(results.sort, 2, ci.quant.bonf.lower))
  }

  results <- list(
    original = t0, mean.boot = mean, conf = conf.quant,
    conf.Bonferroni = conf.quant.bonf, boot.results = as_tibble(resultsmat)
  )

  class(results) <- "Results"

  return(results)
}
