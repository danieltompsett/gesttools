#' G-Estimation for an End of Study Outcome
#'
#' Performs g-estimation of a structural nested mean model (SNMM), based on the outcome regression methods described
#' in Sjolander and Vansteelandt (2016) and Dukes and Vansteelandt (2018). We expect a dataset that holds an end of study outcome that is either binary or continuous,
#' time-varying and/or baseline confounding, and a time-varying exposure that is either binary or continuous.
#'
#'@param data A data frame in long format containing the data to be analysed. See description for details.
#'@param idvar Character string specifying the name of the ID variable in the data.
#'@param timevar Character string specifying the name of the time variable in the data.
#'  Note that timevar must specify time periods as integer values starting from 1 (must not begin at 0).
#'@param Yn Character string specifying the name of the end of study outcome variable.
#'@param An Character string specifying the name of the time-varying exposure variable.
#'@param Ybin TRUE or FALSE indicator of whether the outcome is binary.
#'@param Abin TRUE or FALSE indicator of whether the exposure is binary.
#'Note that if \code{Abin==TRUE} then \code{An} MUST be written as a numeric variable
#'taking values 0 or 1. If not use \code{\link{gestcat}}
#'@param Lny Vector of character strings specifying the names of the
#'  time-varying and/or baseline confounders to be included in the outcome model in quotations.
#'@param Lnp Vector of character strings specifying the names of the time-varying and/or baseline confounders
#'  to be included in the model calculating the propensity scores.
#'@param type Value from 1-4 specifying SNMM type to fit. See details
#'@param Cn Optional character string specifying the name of the censoring indicator variable.
#'  Cn should be a numeric variable taking values 0 or 1, with 1 indicating censored.
#'@param LnC Vector of character strings specifying the names of the time-varying and/or baseline covariates to be used in the censoring score model to calculate
#'  the censoring weights. Note that any variable in \code{LnC} should also be in \code{Lnp} for the validity of the censoring and propensity weights.
#'@param ... Additional arguments, currently not in use.
#'
#' @details Given a time-varying exposure variable, \eqn{A_t=a_t} and time-varying confounders, \eqn{L_t=l_t} measured over time periods \eqn{t=1,\ldots,T}, and an end of study outcome \eqn{Y}
#' measured at time \eqn{T+1}, \code{gest} estimates the causal parameters \eqn{\psi} of a SNMM of the form
#'  \deqn{E(Y(\bar{a}_{t},0)-Y(\bar{a}_{t-1},0)|\bar{a}_{t-1},\bar{l}_{t-1})=\psi z_ta_t \;\forall\; t=1,\ldots,T}
#'  if Y is continuous or
#'  \deqn{\frac{E(Y(\bar{a}_{t},0)|\bar{a}_{t-1},\bar{l}_{t-1})}{E(Y(\bar{a}_{t-1},0)|\bar{a}_{t-1},\bar{l}_{t-1})}=exp(\psi z_ta_t)\;\forall\; t=1,\ldots,T }
#'  if Y is binary. The SNMMs form is defined by the parameter \eqn{z_t}, which can be controlled by the input \code{type} as follows
#'  \itemize{
#'  \item{\code{type=1} }{sets \eqn{z_t=1}. This implies that \eqn{\psi} is the effect of exposure at any time t on Y.}
#'  \item{\code{type=2} }{sets \eqn{z_t=c(1,l_t)}, and adds affect modification by the first named variable in \code{Lny}, which we denote \eqn{L_t}.
#'  Now \eqn{\psi=c(\psi_0,\psi_1)} where \eqn{\psi_0} is the effect of exposure at any time t on Y when \eqn{l_t=0} for all t, modified by
#'  \eqn{\psi_1} for each unit increase in \eqn{l_t} at all times t. Note that effect modification
#'  is currently only supported for binary (written as a numeric 0,1 vector) or continuous confounders.}
#'  \item {\code{type=3} }{allows for time-varying causal effects. It sets \eqn{z_t} to a vector of zeros of length T with a 1 in the t'th position. Now \eqn{\psi=c(\psi_1,\ldots,\psi_T)}
#'  where \eqn{\psi_t} is the effect of \eqn{A_t} on Y.}
#'  \item{\code{type=4} }{allows for a time-varying causal effect that can be modified by the first named variable in \code{Lny}, that is it allows for both time-varying effects and effect modification. It sets \eqn{z_t} to a vector of zeros of length T with \eqn{c(1,l_t)} in the t'th position.
#'  Now \eqn{\psi=(\underline{\psi_1},\ldots,\underline{\psi_T})} where \eqn{\underline{\psi_t}=c(\psi_{0t},\psi_{1t})}. Here \eqn{\psi_{0t}} is the effect of exposure at time t on Y when \eqn{l_t=0} modified by
#'  \eqn{\psi_{1t}} for each unit increase in \eqn{l_t}. Note that effect modification
#'  is currently only supported for binary (written as a numeric 0,1 vector) or continuous confounders.}
#'  }
#'  The data must be in long format, where we assume the convention that each row with \code{time=t} contains \eqn{A_t,L_t} and \eqn{C_{t+1}}. Thus the censoring indicator for each row
#'  should indicate that a user is censored AFTER time t. The end of study outcome should be repeated on each row. If either A or Y are binary, they must be written as numeric vectors taking values either 0 or 1.
#'  The same is true for any covariate that is used for effect modification.\cr
#'  The data must be rectangular with a row entry for every individual for each exposure time 1 up to T. Data rows after censoring should be empty apart from the ID and time variables. This can be done using the function \code{\link{FormatData}}.\cr
#'  By default the censoring, propensity and outcome models include the exposure history at the previous time as a variable. One may consider also including all previous exposure and confounder history as variables in \code{Lny},\code{Lnp}, and \code{LnC} using \code{\link{FormatData}} if necessary.\cr
#'  Censoring weights are handled as described in Sjolander and Vansteelandt (2016). Note that it is necessary that any variable included in \code{LnC} must also be in \code{Lnp}. Missing data not due to censoring is handled automatically by removing rows with missing data prior to fitting the model. If outcome models fail to fit, consider removing covariates from \code{Lny} but keeping
#'  them in \code{Lnp} to reduce collinearity issues.
#'
#'
#' @return List of the fitted causal parameters of the posited SNMM. These are labeled as follows for each SNMM type, where \code{An} is
#' set to the name of the exposure variable, i is the current time period, and \code{Lny[1]} is set to the name of the first confounder in \code{Lny}.
#' \item{\code{type=1} }{\code{An}: The effect of exposure at any time t on outcome. }
#' \item{\code{type=2} }{\code{An}: The effect of exposure at any time t on outcome, when \code{Ln[1]} is set to zero.\cr
#' \code{An:Ln[1]}: The effect modification by \code{Lny[1]}, the additional effect of A on Y for each unit increase in \code{Lny[1]}. }
#' \item{\code{type=3} }{\code{t=i.An}: The effect of exposure at time t=i on outcome.}
#' \item{\code{type=4} }{\code{t=i.An}: The effect of exposure at time t=i on outcome, when \code{Ln[1]} is set to zero.\cr
#' \code{t=i.An:Ln[1]}: The effect modification by \code{Lny[1]}, the additional effect of A on Y at time t=i for each unit increase in \code{Lny[1]}.}
#'
#'
#'
#' @examples
#' datas<-dataexamples(n=10000,seed=123,Censoring=FALSE)
#' data=datas$datagest
#' idvar="id"
#' timevar="time"
#' Yn="Y"
#' An="A"
#' Ybin=FALSE
#' Abin=TRUE
#' Lny=c("L","U")
#' Lnp=c("L","U")
#' type=1
#' Cn=NA
#' LnC=NA
#' gest(data,idvar=idvar,timevar,Yn,An,Ybin,Abin,Lny,Lnp,type=1)
#' gest(data,idvar=idvar,timevar,Yn,An,Ybin,Abin,Lny,Lnp,type=4)
#'
#' #Example with censoring
#' datas<-dataexamples(n=10000,seed=123,Censoring=TRUE)
#' data=datas$datagest
#' Cn="C"
#' LnC=c("L","U")
#' gest(data,idvar,timevar,Yn,An,Ybin,Abin,Lny,Lnp,Cn,LnC,type=1)
#'
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

gest<-function(data,idvar,timevar,Yn,An,Ybin,Abin,Lny,Lnp,type=1,Cn=NA,LnC=NA,...){
# Error messages
if(!is.data.frame(data))(stop("Either no data set has been given, or it is not in a data frame."))
if(any(!sapply(unlist(c("idvar","timevar","Yn","An","Ybin","Abin","Lny","Lnp","type")),exists)))(stop('Missing input'))
if(!is.numeric(data[,Lny[1]]) && type %in% c(2,4))(stop("Effect modification is only supported for a continuous covariate, or binary covariate written as an as.numeric() 0,1 vector"))
if(!is.numeric(data[,An]) && Abin==TRUE)(stop("A binary exposure A must be written as an as.numeric() 0,1 vector, with 1 indicating exposure. If not, consider using gestcat()"))
if(!is.numeric(data[,Yn]) && Ybin==TRUE)(stop("A binary outcome Y must be written as an as.numeric() 0,1 vector, with 1 indicating an event."))
if(!is.na(Cn)==TRUE && !is.numeric(data[,Cn]))(stop("A censoring indicator must be written as an as.numeric() 0,1 vector, with 1 indicating censoring."))
if(!is.na(LnC[1]))(warning("Be sure that any variable in LnC is also in Lnp, else propensity scores may be invalid"))
# Define the ID and time variable
data$id<-data[,idvar]
data$time<-data[,timevar]
T<-max(data$time)

# Obtain the necessary variables of the data and create lagged exposure variable (A_{t-1})
if(is.na(Cn)==TRUE){
  d<-data[,unique(c(An,Yn,Lnp,Lny,idvar,timevar))]
  }else{
  d<-data[,unique(c(An,Yn,Lnp,Lny,idvar,timevar,LnC,Cn))]}

d<-slide(data=d,Var=An,GroupVar="id",slideBy=-1,NewVar=paste("L",1,An,sep=""),reminder=F)
Alagn<-paste("L",1,An,sep="")

# Set lagged exposure at time 1 to zero
d[d$time %in% 1,Alagn]<-0
d$int<-1

# Define propensity score model and relevant glm family
if (T==1){
  lmp<-reformulate(termlabels=c(Lnp),response=An)
  } else {
  lmp<-reformulate(termlabels=c(Alagn,"as.factor(time)",Lnp),response=An)
  }
if (Abin==TRUE){
  modp<-glm(lmp,family="binomial",data=d)
  } else {
  modp<-glm(lmp,family="gaussian",data=d)
  }

# Calculate propensity scores
props<-predict(modp,type="response",newdata=d)
d$prs<-props

# Calculate initial censoring weights if Cn and LnC are supplied
if(is.na(Cn)==TRUE){
  d$w<-1} else{
  lmc<-reformulate(termlabels=c(LnC,"as.factor(time)"),response=Cn)
  modc<-glm(lmc,family="binomial",data=d)
  # Censoring scores
  cps<-1-predict(modc,type="response",newdata=d)
  d$cps<-cps
  # Create Indicator Function of whether C=0
  d[,paste(Cn,"0",sep="")]<-as.integer(!d[,Cn])
  # Setting up the initial denominator product of censoring scores "cprod" and initial censoring weights "w"
  d$cprod<-d$cps
  d[is.na(d$cprod)==TRUE,"cprod"]<-1
  d$w<-d[,paste(Cn,"0",sep="")]/d$cps
  }

# Set up a variable "H" to hold the counterfactual outcome values
# and set the value of H at the final exposure time T to the outcome
d$H<-0
d[d$time==T,"H"]<-d[d$time==T,Yn]
# Create a copy of the data set with complete cases for initial estimate of psi
dcom<-d[complete.cases(d),]

# Define values of internal variables 'z' and 'timevarying' based on input 'type'
if (type==1){
  z<-c("int")
  timevarying<-FALSE
  }else if (type==2){
  z<-c("int",Lny[1])
  timevarying<-FALSE
  }else if (type==3){
  z<-c("int")
  timevarying<-TRUE
  }else if (type==4){
  z<-c("int",Lny[1])
  timevarying<-TRUE
  }

# Get family for outcome model
if(Ybin==TRUE){
  family<-Gamma(link="log")
  }else{
  family<-gaussian
  }

# Obtain the correct formulas for the outcome model
par1<-paste(eval(An),eval(z),sep=":")
par1[par1==paste(eval(An),"int",sep=":")]<-paste(eval(An))
par2<-paste("prs",eval(z),sep=":")
par2[par2==paste("prs","int",sep=":")]<-paste("prs")

if (T==1){
  lmy<-reformulate(termlabels=c(par1,par2,Lny),response=Yn)
  } else {
  lmy<-reformulate(termlabels=c(par1,par2,Alagn,Lny),response=Yn)
  lmH<-reformulate(termlabels=c(par1,par2,Alagn,Lny),response="H")
  lmH1<-reformulate(termlabels=c(par1,par2,Lny),response="H")
  }

# Perform g-estimation for SNMM types 1 or 2 (non time-varying psi)
if (timevarying==FALSE){
  # Obtain the causal effect of A_T on Y as the initial estimate of psi
  dt<-dcom[dcom$time==T,]
  out<-summary(geem(terms(lmy),family=family,id=dt$id,data=dt,weights=dt$w))
  psi<-out$beta[match(par1,out$coefnames)]
  names(psi)<-par1

  # Obtain counterfactual outcomes and censoring weights for time T-1
  # Set i as the minimum time period for which counterfactuals and censoring weights are calculated
  # Set k to the current time period being calculated.
  i<-T-1
  while(i>=1){
    k<-T-1
	  while(k>=i){
	    if(is.na(Cn)==FALSE){
	    # Censoring weights
	    d[d$time==k,"cprod"]<-d[d$time==k+1,"cprod"]*d[d$time==k,"cps"]
	    d[d$time==k,"w"]<-d[d$time==T,paste(Cn,"0",sep="")]/d[d$time==k,"cprod"]
	    }
	    # Counterfactuals
	    if(Ybin==FALSE){
	      if (length(z)==1){
	        d[d$time==k,"H"]<-d[d$time==(k+1),"H"]-
	        d[d$time==(k+1),An]*psi*d[d$time==(k+1),z]
          }else{
          d[d$time==k,"H"]<-d[d$time==(k+1),"H"]-
	        rowSums(d[d$time==(k+1),An]*sweep(d[d$time==(k+1),z],2,psi,"*"))
	        }
	      }
      if(Ybin==TRUE){
	      if (length(z)==1){
	        d[d$time==k,"H"]<-d[d$time==(k+1),"H"]*
	        exp(-d[d$time==(k+1),An]*psi*d[d$time==(k+1),z])
          }else{
          d[d$time==k,"H"]<-d[d$time==(k+1),"H"]*
	        exp(-rowSums(d[d$time==(k+1),An]*sweep(d[d$time==(k+1),z],2,psi,"*")))
	        }
        }
	  # Current time period is lowered and censoring weights and counterfactuals
	  # calculated until current value of i is reached.
    k<-k-1
	  }
    # Obtain relevant time periods for which counterfactuals have been calculated
    dt<-d[d$time %in% seq(i,T,by=1),]
    dtcom<-dt[complete.cases(dt),]
    # Update the estimate of psi including data at time T-1.
    out<-summary(geem(terms(lmH),family=family,id=dtcom$id,data=dtcom,weights=dtcom$w))
    psi<-out$beta[match(par1,out$coefnames)]
    names(psi)<-par1
    # Set minimum time period to time T-2 and repeat calculation of counterfactuals
    # and censoring weights up to time T-2 and update psi. This is repeated until
    # time period 1 is reached and the final estimate of psi is obtained.
    i<-i-1
    }

  return(list(psi=psi))
  }

# Perform g-estimation for SNMM types 3 or 4 (time-varying psi)
if(timevarying==TRUE){
  # Obtain estimate of psi_T
  dt<-dcom[dcom$time==T,]
  out<-summary(geem(terms(lmy),family=family,id=dt$id,data=dt,weights=dt$w))
  psi<-out$beta[match(par1,out$coefnames)]
  names(psi)<-par1
  psilist<-as.list(NULL)
  psilist[[T]]<-psi

  # Obtain counterfactuals and censoring weights for time T-1 as before
  # using the estimate of psi_T
  i<-(T-1)
  while(i>=1){
    k<-(T-1)
	  while(k>=i){
      if(is.na(Cn)==FALSE){
	    d[d$time==k,"cprod"]<-d[d$time==k+1,"cprod"]*d[d$time==k,"cps"]
      d[d$time==k,"w"]<-d[d$time==T,paste(Cn,"0",sep="")]/d[d$time==k,"cprod"]
      }
	    if (Ybin==FALSE){
	      if (length(z)==1){
	        d[d$time==k,"H"]<-d[d$time==(k+1),"H"]-
	        d[d$time==(k+1),An]*psilist[[k+1]]*d[d$time==(k+1),z]
          }else{
          d[d$time==k,"H"]<-d[d$time==(k+1),"H"]-
	        rowSums(d[d$time==(k+1),An]*sweep(d[d$time==(k+1),z],2,psilist[[k+1]],"*"))
	       }
      }
	    if (Ybin==TRUE){
	      if (length(z)==1){
	        d[d$time==k,"H"]<-d[d$time==(k+1),"H"]*
	        exp(-d[d$time==(k+1),An]*psilist[[k+1]]*d[d$time==(k+1),z])
          }else{
          d[d$time==k,"H"]<-d[d$time==(k+1),"H"]*
	        exp(-rowSums(d[d$time==(k+1),An]*sweep(d[d$time==(k+1),z],2,psilist[[k+1]],"*")))
	        }
	       }
    k<-k-1
	  }
    # Obtain data at time T-1 for calculation of psi
	  dt<-d[d$time %in% i,]
	  dtcom<-dt[complete.cases(dt),]
	  if (i==1){
	    out<-summary(geem(terms(lmH1),family=family,id=dtcom$id,data=dtcom,weights=dtcom$w))
	    psi<-out$beta[match(par1,out$coefnames)]
	    names(psi)<-par1
	    psilist[[i]]<-psi
      }else{
	    out<-summary(geem(terms(lmH),family=family,id=dtcom$id,data=dtcom,weights=dtcom$w))
	    psi<-out$beta[match(par1,out$coefnames)]
	    names(psi)<-par1
	    psilist[[i]]<-psi
      }
	  # Set minimum time period to time T-2 and repeat calculation of counterfactuals
	  # and censoring weights to time T-2 and update psi. This is repeated until
	  # time period 1 is reached and the final estimate of psi is obtained.
	  i<-i-1
  }
  # Create relevant names for output and create results
  nam<-as.vector(NULL)
  for (p in 1:T){
  nam[p]<-paste("t=",p,sep="")
  }
  names(psilist)<-nam
  results<-unlist(psilist)
  return(list(psi=results))
  }
}




