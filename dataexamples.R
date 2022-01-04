#' Generate Simulated Example Datasets
#'
#' The code simulates four datasets designed to demonstrate the
#' g-estimation functions of the package. These are used in the examples
#' in the user manual. Each dataset comprises of an outcome Y (time-varying or end of study), time-varying exposure A, time-varying confounder L,
#' a baseline confounder U, and optionally a censoring indicator C over 3 time periods.
#'
#' @param n Number of individuals in the dataset.
#' @param seed Random seed used for data generation.
#' @param Censoring TRUE or FALSE indicator of whether to include a censoring indicator C.
#' If \code{Censoring=TRUE}, data entries for A, Y, L and U are set to missing after censoring.
#'
#' @return Returns a list of four datasets labeled \code{datagest}, \code{datagestmult},
#' \code{datagestcat}, and \code{datagestmultcat}, designed to demonstrate an end of study outcome with a binary exposure (\code{datagest}), a time varying outcome study with a binary exposure (\code{datagestmult}),
#' or an end of study or time varying outcome with a categorical exposure (\code{datagestcat} or \code{datagestmultcat}).
#'
#' @examples
#' datas <- dataexamples(n = 1000, seed = 34567, Censoring = FALSE)
#' data <- datas$datagest
#' head(data, n = 20)
#' # Multiple outcome data with censoring
#' datas <- dataexamples(n = 100, seed = 34567, Censoring = TRUE)
#' data <- datas$datagestmultcat
#' head(data, n = 20)
#' @importFrom stats plogis
#'
#' @export

dataexamples <- function(n = 1000, seed = NULL, Censoring = FALSE) {
  if (!is.null(seed)) set.seed(seed)

  n <- n
  id <- seq(1, n, by = 1)
  U <- rnorm(n, 0, 1)
  L.1 <- rnorm(n, 1 + U, 1)
  a.1 <- plogis(1 + 0.1 * L.1)
  A.1 <- rbinom(n, 1, a.1)
  Y.1 <- rnorm(n, 1 + A.1 + L.1 + U, 1)

  L.2 <- rnorm(n, 1 + (A.1 / 2) + L.1 + U, 1)
  a.2 <- plogis(1 + 0.1 * L.2 + 0.1 * A.1)
  A.2 <- rbinom(n, 1, a.2)
  Y.2 <- rnorm(n, 1 + (A.1 / 2) + A.2 + L.1 + L.2 + U, 1)

  L.3 <- rnorm(n, 1 + (A.2 / 2) + L.2 + U, 1)
  a.3 <- plogis(1 + 0.1 * L.3 + 0.1 * A.2)
  A.3 <- rbinom(n, 1, a.3)
  Y.3 <- rnorm(n, 1 + (A.2 / 2) + A.3 + L.1 + L.2 + L.3 + U, 1)
  # End of study Outcome
  Y <- rnorm(n, 1 + (A.2 / 2) + A.3 + L.1 + L.2 + L.3 + U, 1)

  if (Censoring == TRUE) {
    C.1 <- rbinom(n, 1, plogis(-1 + 0.001 * A.1 + 0.001 * L.1))
    C.2 <- rbinom(n, 1, plogis(-1 + 0.001 * A.2 + 0.001 * L.2))
    C.3 <- rbinom(n, 1, plogis(-1 + 0.001 * A.3 + 0.001 * L.3))
    C.2[C.1 == 1] <- 1
    C.3[C.2 == 1] <- 1
    Y[C.3 == 1] <- NA
    Y.3[C.3 == 1] <- NA
    A.3[C.2 == 1] <- NA
    L.3[C.2 == 1] <- NA
    Y.2[C.2 == 1] <- NA
    A.2[C.1 == 1] <- NA
    L.2[C.1 == 1] <- NA
    Y.1[C.1 == 1] <- NA
    # Create end of study outcome data and set as long format
    dw <- as.data.frame(cbind(id, Y, A.1, A.2, A.3, L.1, L.2, L.3, C.1, C.2, C.3, U))
    dl <- reshape(dw, direction = "long", varying = c(
      "A.1", "A.2", "A.3",
      "L.1", "L.2", "L.3",
      "C.1", "C.2", "C.3"
    ))
    # Order data by time and ID
    dl <- dl[order(dl$id, dl$time), ]
    datagest <- dl

    dw <- as.data.frame(cbind(id, Y.1, Y.2, Y.3, A.1, A.2, A.3, L.1, L.2, L.3, C.1, C.2, C.3, U))
    dl <- reshape(dw, direction = "long", varying = c(
      "Y.1", "Y.2", "Y.3",
      "A.1", "A.2", "A.3",
      "L.1", "L.2", "L.3",
      "C.1", "C.2", "C.3"
    ))
    # Order data by time and ID
    dl <- dl[order(dl$id, dl$time), ]
    datagestmult <- dl
  } else {
    # Create end of study outcome data and set as long format
    dw <- as.data.frame(cbind(id, Y, A.1, A.2, A.3, L.1, L.2, L.3, U))
    dl <- reshape(dw, direction = "long", varying = c(
      "A.1", "A.2", "A.3",
      "L.1", "L.2", "L.3"
    ))
    # Order data by time and ID
    dl <- dl[order(dl$id, dl$time), ]
    datagest <- dl

    dw <- as.data.frame(cbind(id, Y.1, Y.2, Y.3, A.1, A.2, A.3, L.1, L.2, L.3, U))
    dl <- reshape(dw, direction = "long", varying = c(
      "Y.1", "Y.2", "Y.3",
      "A.1", "A.2", "A.3",
      "L.1", "L.2", "L.3"
    ))
    # Order data by time and ID
    dl <- dl[order(dl$id, dl$time), ]
    datagestmult <- dl
  }


  ### Categorical A###########

  n <- n
  set.seed(seed)
  id <- seq(1, n, by = 1)
  U <- rnorm(n, 0, 1)
  L.1 <- rnorm(n, 1 + U, 1)
  a.1 <- plogis(1 + 0.1 * L.1)

  A.1 <- as.vector(NULL)
  for (i in 1:n) {
    A.1[i] <- sample(letters[1:3], 1,
      replace = TRUE,
      prob = c(1 - (3 * a.1[i]) / 5, a.1[i] / 5, 2 * (a.1[i] / 5))
    )
  }
  A.1 <- as.factor(A.1)
  A.1par <- as.vector(NULL)
  A.1par[A.1 == letters[1]] <- 1
  A.1par[A.1 == letters[2]] <- 2
  A.1par[A.1 == letters[3]] <- 3

  Y.1 <- rnorm(n, 1 + A.1par + L.1 + U, 1)

  L.2 <- rnorm(n, 1 + A.1par / 2 + L.1 + U, 1)
  a.2 <- plogis(1 + 0.1 * L.2 + A.1par)
  A.2 <- as.vector(NULL)
  for (i in 1:n) {
    A.2[i] <- sample(letters[1:3], 1,
      replace = TRUE,
      prob = c(1 - (3 * a.2[i] / 5), a.2[i] / 5, 2 * (a.2[i] / 5))
    )
  }
  A.2 <- as.factor(A.2)
  A.2par <- as.vector(NULL)
  A.2par[A.2 == letters[1]] <- 1
  A.2par[A.2 == letters[2]] <- 2
  A.2par[A.2 == letters[3]] <- 3

  Y.2 <- rnorm(n, 1 + A.1par / 2 + A.2par + L.1 + L.2 + U, 1)

  L.3 <- rnorm(n, 1 + A.2par / 2 + L.2 + U, 1)
  a.3 <- plogis(1 + 0.1 * L.3 + A.2par)
  A.3 <- as.vector(NULL)
  for (i in 1:n) {
    A.3[i] <- sample(letters[1:3], 1,
      replace = TRUE,
      prob = c(1 - (3 * a.3[i] / 5), a.3[i] / 5, 2 * (a.3[i] / 5))
    )
  }
  A.3 <- as.factor(A.3)
  A.3par <- as.vector(NULL)
  A.3par[A.3 == letters[1]] <- 1
  A.3par[A.3 == letters[2]] <- 2
  A.3par[A.3 == letters[3]] <- 3

  Y.3 <- rnorm(n, 1 + A.2par / 2 + A.3par + L.1 + L.2 + L.3 + U, 1)

  Y <- rnorm(n, 1 + A.2par / 2 + A.3par + L.1 + L.2 + L.3 + U, 1)

  if (Censoring == TRUE) {
    C.1 <- rbinom(n, 1, plogis(-1 + 0.001 * A.1par + 0.001 * L.1))
    C.2 <- rbinom(n, 1, plogis(-1 + 0.001 * A.2par + 0.001 * L.2))
    C.3 <- rbinom(n, 1, plogis(-1 + 0.001 * A.3par + 0.001 * L.3))
    C.2[C.1 == 1] <- 1
    C.3[C.2 == 1] <- 1
    Y[C.3 == 1] <- NA
    Y.3[C.3 == 1] <- NA
    A.3[C.2 == 1] <- NA
    L.3[C.2 == 1] <- NA
    Y.2[C.2 == 1] <- NA
    A.2[C.1 == 1] <- NA
    L.2[C.1 == 1] <- NA
    Y.1[C.1 == 1] <- NA


    # Create end of study outcome data and set as long format
    dw <- as.data.frame(cbind(id, Y, A.1, A.2, A.3, L.1, L.2, L.3, C.1, C.2, C.3, U))
    dl <- reshape(dw, direction = "long", varying = c(
      "A.1", "A.2", "A.3",
      "L.1", "L.2", "L.3",
      "C.1", "C.2", "C.3"
    ))
    # Order data by time and ID
    dl <- dl[order(dl$id, dl$time), ]
    datagestcat <- dl
    datagestcat$A <- as.factor(datagestcat$A)
    levels(datagestcat$A) <- letters[1:3]

    dw <- as.data.frame(cbind(id, Y.1, Y.2, Y.3, A.1, A.2, A.3, L.1, L.2, L.3, C.1, C.2, C.3, U))
    dl <- reshape(dw, direction = "long", varying = c(
      "Y.1", "Y.2", "Y.3",
      "A.1", "A.2", "A.3",
      "L.1", "L.2", "L.3",
      "C.1", "C.2", "C.3"
    ))
    # Order data by time and ID
    dl <- dl[order(dl$id, dl$time), ]
    datagestmultcat <- dl
    datagestmultcat$A <- as.factor(datagestmultcat$A)
    levels(datagestmultcat$A) <- letters[1:3]
  } else {
    # Create end of study outcome data and set as long format
    dw <- as.data.frame(cbind(id, Y, A.1, A.2, A.3, L.1, L.2, L.3, U))
    dl <- reshape(dw, direction = "long", varying = c(
      "A.1", "A.2", "A.3",
      "L.1", "L.2", "L.3"
    ))
    # Order data by time and ID
    dl <- dl[order(dl$id, dl$time), ]
    datagestcat <- dl
    datagestcat$A <- as.factor(datagestcat$A)
    levels(datagestcat$A) <- letters[1:3]


    dw <- as.data.frame(cbind(id, Y.1, Y.2, Y.3, A.1, A.2, A.3, L.1, L.2, L.3, U))
    dl <- reshape(dw, direction = "long", varying = c(
      "Y.1", "Y.2", "Y.3",
      "A.1", "A.2", "A.3",
      "L.1", "L.2", "L.3"
    ))
    # Order data by time and ID
    dl <- dl[order(dl$id, dl$time), ]
    datagestmultcat <- dl
    datagestmultcat$A <- as.factor(datagestmultcat$A)
    levels(datagestmultcat$A) <- letters[1:3]
  }


  datas <- list(
    datagest = datagest, datagestmult = datagestmult, datagestcat = datagestcat,
    datagestmultcat = datagestmultcat
  )


  return(datas)
}
