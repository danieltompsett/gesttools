% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gestmultcat.R
\name{gestmultcat}
\alias{gestmultcat}
\title{G-Estimation for a Time-Varying Outcome and Categorical Time-Varying Exposure}
\usage{
gestmultcat(
  data,
  idvar,
  timevar,
  Yn,
  An,
  Ybin,
  Lny,
  Lnp,
  type = 1,
  Cn = NA,
  LnC = NA,
  cutoff = NA,
  ...
)
}
\arguments{
\item{data}{A data frame in long format containing the data to be analysed. See description for details.}

\item{idvar}{Character string specifying the name of the ID variable in data.}

\item{timevar}{Character string specifying the name of the time variable in the data.
Note that timevar must specify time periods as integer values starting from 1 (must not begin at 0).}

\item{Yn}{Character string specifying the name of the time-varying outcome variable.}

\item{An}{Character string specifying the name of the time-varying exposure variable.}

\item{Ybin}{TRUE or FALSE indicator of whether the outcome is binary.}

\item{Lny}{Vector of character strings specifying the names of the
confounders to be included in the outcome model in quotations.}

\item{Lnp}{Vector of character strings specifying the names of the confounders
to be included in the model calculating the propensity scores.}

\item{type}{Value from 1-4 specifying SNMM type to fit. See details.}

\item{Cn}{Optional character string specifying the name of the censoring indicator variable.
Cn should be a numeric vector taking values 0 or 1, with 1 indicating censored.}

\item{LnC}{Vector of character strings specifying the names of the covariates to be used in the censoring score model to calculate
the censoring weights. Note that any variable in \code{LnC} should also be in \code{Lnp} for the validity of the censoring and propensity weights.}

\item{cutoff}{An integer taking value from 1 up to T, where T is the maximum value of timevar.
Instructs the function to estimate causal effects based only on exposures up to \code{cutoff} time periods prior to
outcomes. See details.}

\item{...}{Additional arguments, currently not in use.}
}
\value{
List of the fitted causal parameters of the posited SNMM. These are labeled as follows for each SNMM type, where An is
set to the name of the exposure variable, i is the current value of c, j is the category level, and Lny[1] is set to the name of the first confounder in \code{Lny}.
\item{\code{type=1} }{Anj: The effect of exposure at category j, at any time t on all subsequent outcome times.}
\item{\code{type=2} }{Anj: The effect of exposure at category j on outcome at any time t, when Ln[1] is set to zero, on all subsequent outcome times.\cr
Anj:Ln[1]: The effect modification by Lny[1], the additional effect of A at category j on all subsequent Y for each unit increase in Lny[1].}
\item{\code{type=3} }{c=i.Anj: The effect of exposure at any time t on outcome \code{c=i} time periods later.}
\item{\code{type=4} }{c=i.Anj: The effect of exposure at any time t on outcome \code{c=i} time periods later, when Ln[1] is set to zero.\cr
c=i.Anj:Ln[1]: The effect modification by Lny[1], the additional effect of exposure at category j, on outcome c=i time periods later for each unit increase in Lny[1]. }
}
\description{
Performs g-estimation of a structural nested mean model (SNMM), based on the outcome regression methods described
in Sjolander and Vansteelandt (2016) and Dukes and Vansteelandt (2018). We assume a dataset with a time-varying outcome that is either binary or continuous,
time-varying and/or baseline confounding, and a categorical time-varying exposure of three or more categories.
}
\details{
Suppose a series of time periods \eqn{1,\ldots,T+1} whereby a time-varying exposure and confounder (\eqn{A_t=a_t} and \eqn{L_t=l_t}) are measured over times \eqn{t=1,\ldots,T} and
a time varying outcome \eqn{Y_s=y_s} is measured over times \eqn{s=2,\ldots,T+1}. Define \eqn{c=s-t} as the step length, that is number of time periods separating an exposure measurement, and subsequent outcome measurement. Also suppose that \eqn{A_t=a_t} is a categorical variable consisting of \eqn{k>2} categories.
These categories may take any arbitrary list of names, but we assume for theory purposes they are labeled as \eqn{j=0,1\ldots,k}
where \eqn{j=0} denotes no exposure, or some reference category. Define binary variables \eqn{A_t^j=a_t^j} \eqn{j=0,1,\ldots,k} where \eqn{A_t^j=1} if \eqn{A_t=j} and 0 otherwise.
By using the transform \eqn{t=s-c}, \code{gestmultcat} estimates the causal parameters \eqn{\psi} of a SNMM of the form
 \deqn{E\{Y_s(\bar{a}_{s-c},a^0)-Y_s(\bar{a}_{s-c-1},a^0)|\bar{a}_{s-c-1},\bar{l}_{s-c}\}=\sum_{j=1}^{k}\psi^j z_{sc}a^j_{s-c} \; \forall c=1,\ldots,T\; and\; \forall s>c}
 if Y is continuous or
 \deqn{\frac{E(Y_s(\bar{a}_{s-c},a^0)|\bar{a}_{s-c-1},\bar{l}_{s-c-1})}{E(Y_s(\bar{a}_{s-c-1},a^0)|\bar{a}_{s-c-1},\bar{l}_{s-c-1})}=exp(\sum_{j=1}^{k}\psi^j z_{sc}a^j_{s-c}) \; \forall c=1,\ldots,T\; and \; \forall s>c }
 if Y is binary. The SNNM fits a separate set of causal parameters \eqn{\psi^j}, for the effect of exposure at category \eqn{j} on outcome,
 compared to exposure at the reference category 0, for each non-reference category. The models form is defined by the parameter \eqn{z_{sc}}, which can be controlled by the input \code{type} as follows
 \itemize{
 \item{\code{type=1} }{sets \eqn{z_{sc}=1}. This implies that \eqn{\psi^j} is now the effect of exposure when set to category \eqn{j}, compared to when set to the reference category, at any time t on all subsequent outcome periods.}
 \item{\code{type=2} }{sets \eqn{z_{sc}=c(1,l_{s-c})} and adds affect modification by the first named variable in \code{Lny}, which we denote \eqn{L_t}. Now \eqn{\psi^j=c(\psi^j_0,\psi^j_1)} where \eqn{\psi^j_0} is the effect of exposure when set to category \eqn{j}, compared to when set to the reference category, at any time t on all subsequent outcome periods when \eqn{l_{s-c}=0} for all t, modified by
 \eqn{\psi^j_1} for each unit increase in \eqn{l_{s-c}} at all times t. Note that effect modification
 is currently only supported for binary or continuous confounders.}
 \item {\code{type=3} }{can posit a time-varying causal effect for each value of \eqn{c}, that is the causal effect for the exposure on outcome \eqn{c} time periods later.
 We set \eqn{z_{sc}} to a vector of zeros of length T with a 1 in the \eqn{c=s-t}'th position. Now \eqn{\psi^j=c(\psi^j_{1},\ldots,\psi^j_{T})}
 where \eqn{\psi^j_c} is the effect of exposure, when set to category \eqn{j}, on outcome \eqn{c} time periods later for all \eqn{s>c} that is \eqn{A^j_{s-c}} on \eqn{Y_s} for all \eqn{s>c}.}
 \item{\code{type=4} }{allows for a time-varying causal effect that can be modified by the first named variable in \code{Lny}, that is it allows for both time-varying effects and effect modification. It sets \eqn{z_{sc}} to a vector of zeros of length T with \eqn{c(1,l_{s-c})} in the \eqn{c=s-t}'th position.
 Now \eqn{\psi^j=(\underline{\psi^j_1},\ldots,\underline{\psi^j_T})} where \eqn{\underline{\psi^j_c}=c(\psi^j_{0c},\psi^j_{1c})}. Here \eqn{\psi^j_{0c}} is the effect of exposure when set to category \eqn{j} on outcome \eqn{c} time periods later, given \eqn{l_{s-c}=0}, for all \eqn{s>c}, modified by
 \eqn{\psi^j_{1c}} for each unit increase in \eqn{l_{s-c}} for all \eqn{s>c}. Note that effect modification
 is currently only supported for binary or continuous confounders.}
 }
The data must be in long format, where we assume the convention that each row with \code{time=t} contains \eqn{A_t,L_t} and \eqn{C_{t+1},Y_{t+1}}. That is the censoring indicator for each row
should indicate that a user is censored AFTER time t, and the outcome the first outcome that occurs AFTER \eqn{a_t} and \eqn{l_t} are measured.
For example, data at time 1, should contain \eqn{A_1}, \eqn{L_1}, \eqn{Y_{2}}, and optionally \eqn{C_2}. If Y is binary, it must be written as a numeric vector taking values either 0 or 1.
The same is true for any covariate that is used for effect modification.\cr
The data must be rectangular with a row entry for every individual for each exposure time 1 up to T. Data rows after censoring should be empty apart from the ID and time variables. This can be done using the function \code{\link{FormatData}}.\cr
By default the censoring, propensity and outcome models include the exposure history at the previous time as a variable. One may consider also including all previous exposure and confounder history as variables in \code{Lny},\code{Lnp}, and \code{LnC} if necessary.\cr
Censoring weights are handled as described in Sjolander and Vansteelandt (2016). Note that it is necessary that any variable included in \code{LnC} must also be in \code{Lnp}. Missing data not due to censoring is automatically handled by removing rows with missing data prior to fitting the model. If outcome models fail to fit, consider removing covariates from \code{Lny} but keeping
them in \code{Lnp} to reduce collinearity issues, or consider the sparseness of the data.
}
\examples{
datas<-dataexamples(n=10000,seed=123,Censoring=FALSE)
data=datas$datagestmultcat
#A is a categorical variable with categories labeled 1,2 and 3, with 1 the
#reference category
idvar="id"
timevar="time"
Yn="Y"
An="A"
Ybin=FALSE
#Remove U from Y to avoid collinearity
Lny=c("L","U")
Lnp=c("L","U")
Cn<-NA
LnC<-NA
type=NA

gestmultcat(data,idvar,timevar,Yn,An,Ybin,Lny,Lnp,type=1)
gestmultcat(data,idvar,timevar,Yn,An,Ybin,Lny,Lnp,type=3)

#Example with censoring
datas<-dataexamples(n=10000,seed=123,Censoring=TRUE)
data=datas$datagestmultcat
Cn="C"
LnC=c("L","U")
gestmultcat(data,idvar,timevar,Yn,An,Ybin,Lny,Lnp,type=1,Cn,LnC)

}
\references{
Vansteelandt, S., & Sjolander, A. (2016). Revisiting g-estimation of the Effect of a Time-varying Exposure Subject to Time-varying Confounding, Epidemiologic Methods, 5(1), 37-56. <doi:10.1515/em-2015-0005>.

Dukes, O., & Vansteelandt, S. (2018). A Note on g-Estimation of Causal Risk Ratios, American Journal of Epidemiology, 187(5), 1079–1084. <doi:10.1093/aje/kwx347>.
}
