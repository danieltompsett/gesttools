#' Formats Data Into Correct Form
#'
#' Takes a dataset in long format and puts it into the required format for use
#' with the g-estimation functions. Specifically it ensures there exists a data
#' entry for each individual at each time period, by adding empty rows, and orders the dataset by
#' time and identifier. It can also create variables for the exposure histories of all time-varying
#' variables in the data.
#'
#' @param data A data frame in long format containing the data to be analysed.
#' @param idvar A character string specifying the name of of the variable specifying
#' an individuals identifier.
#' @param timevar A character string specifying the name of the time variable.
#' Note that time periods must be labeled as integers starting from 1
#' (\eqn{1,2,\ldots}).
#' @param An A character string specifying the name of the exposure variable
#' @param varying A vector of character strings specifying the names of the variables
#' to be included in the analysis which are time-varying. Specifically
#' the exposure, time-varying confounders and (if applicable) the time-varying outcome.
#' If \code{Cn} is specified, it is added to \code{varying} automatically.
#' @param Cn Optional character string specifying the name of the censoring indicator if present.
#' @param GenerateHistory A TRUE or FALSE indicator. If set to TRUE, variables are generated
#' corresponding to the lagged histories of all variables included in \code{varying}.
#' These will be labeled as \code{LagVari} where \code{Var} is the variable name and \code{i}
#' indicates how much the variable is lagged by. For example \code{LagAn2} is the value of \code{An}, 2
#' time periods prior.
#' @param GenerateHistoryMax An optional positive integer specifying \code{GenerateHistory} to generate exposure histories
#' up to \code{GenerateHistoryMax} time periods prior.
#'
#' @return A data frame in long format with additional rows added as necessary. If
#' \code{data} is already in the correct format then no additional rows will be added.
#'
#' @details Note that any variable in \code{varying} that is strictly categorical MUST be declared as
#' an \code{as.factor()} variable. Binary or continuous variables should be declared as an
#' \code{as.numeric()} variable.
#'
#' @examples
#' data <- dataexamples(n = 1000, seed = 3456, Censoring = TRUE)$datagest
#' # To demonstrate the function we
#' # Delete the third row, corresponding to the entry for ID 1 at time 3
#' data <- data[-3, ]
#' datanew <- FormatData(
#'   data = data, idvar = "id", timevar = "time", An = "A",
#'   Cn = "C", varying = c("A", "L"), GenerateHistory = TRUE, GenerateHistoryMax = 1
#' )
#' head(datanew)
#' # Note that the missing entry has been re-added,
#' # with missing values for A and L in the third row
#' # An example with lagged history of time varying variables created.
#' data <- dataexamples(n = 1000, seed = 3456, Censoring = TRUE)$datagestmultcat
#' datanew <- FormatData(
#'   data = data, idvar = "id", timevar = "time", An = "A",
#'   Cn = "C", varying = c("Y","A", "L"), GenerateHistory = TRUE, GenerateHistoryMax = NA
#' )
#' head(datanew)
#' @export


FormatData <- function(data, idvar, timevar, An, varying, Cn = NA, GenerateHistory = FALSE,
                       GenerateHistoryMax = NA) {
  if (!is.data.frame(data)) (stop("Either no data set has been given, or it is not in a data frame."))
  if (!is.na(Cn)) {
    varying <- c(varying, Cn)
  }

  datwide <- reshape(data, direction = "wide", timevar = timevar, idvar = idvar, v.names = varying)
  datrec <- reshape(datwide, direction = "long", timevar = timevar, idvar = idvar)
  datrec <- datrec[order(datrec[, idvar], datrec[, timevar]), ]

  if (GenerateHistory == TRUE) {
    T <- max(datrec[, timevar])
    if (T <= 1) (stop("Lagged variables cannot be created with only 1 time period"))
    if (is.na(GenerateHistoryMax) == TRUE) {
      GenerateHistoryMax <- T - 1
    }
    varying <- varying[!(varying %in% Cn)]
    histmax <- min(T - 1, GenerateHistoryMax)

    # Function to generate lagged variables
    lagged <- function(name) {
      for (i in 1:histmax) {
        datrec <- slide(data = datrec, Var = name, GroupVar = "id", slideBy = -i, NewVar = paste("Lag", i, name, sep = ""), reminder = F)
        # Set Value of lagged variable that do not exists to either 0, or the reference category.
        if (is.factor(datrec[, name]) == TRUE) {
          datrec[, paste("Lag", i, name, sep = "")] <- as.factor(datrec[, paste("Lag", i, name, sep = "")])
          levels(datrec[, paste("Lag", i, name, sep = "")]) <- levels(datrec[, name])
          datrec[datrec[, timevar] %in% seq(1, i, by = 1), paste("Lag", i, name, sep = "")] <- levels(datrec[, name])[1]
        } else {
          datrec[datrec[, timevar] %in% seq(1, i, by = 1), paste("Lag", i, name, sep = "")] <- 0
        }
      }
      return(datrec)
    }

    datrec <- Reduce(merge, lapply(varying, lagged))
  }


  return(datrec)
}
