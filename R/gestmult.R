#' G-estimation for a time-varying outcome.
#'
#' Performs G-estimation for a structural nested mean model (SNMM), based on the outcome regression methods described
#' in Sjolander and Vansteelendt (2016) and Dukes and Vansteelendt (2018). We assume a dataset with a time-varying outcome that is either binary or continuous,
#' time-varying counfounding, and a time-varying exposure that is either binary or continuous.
#'
#'@param data A data table in long format containing the data to be analysed. See description for details.
#'@param idvar Character string specifying the name of the ID variable in data.
#'@param timevar Character string specifying the name of the time variable in the data.
#'  Note that timevar must specify time periods as integer values starting from 1 (must not begin at 0).
#'@param Yn Character string specifying the name of the outcome variable.
#'@param An Character string specifying the name of the exposure variable.
#'@param Ybin True or False indicator of whether the outcome is binary.
#'@param Abin True or false indicator of whether the exposure is binary.
#'@param Lny Vector of character strings specifying the names of the
#'  confounders to be included in the outcome model in quotations.
#'@param Lnp Vector of character strings specifying the names of the confounders
#'  to be included in the model calculating the propensity scores.
#'@param type Value from 1-4 specifying SNMM type to fit. See details
#'@param Cn Optional character string specifying the name of the censoring indicator variable.
#'  Cn should be a vector taking values 0 or 1, with 1 indicating censored.
#'@param LnC Vector of character strings specifying the names of the covariates to be used in the censoring score model to calculate
#'  the censoring weights.
#'@param cutoff An integer c taking values from 1 up to T, where T is the maximum value of timevar.
#'Instructs the function to estimate causal effects based only on exposures up to c time periods prior to
#'outcomes. See details
#'
#'@param ... Additional arguments, currently not in use.
#'
#' @details Suppose a series of time periods \eqn{1,\ldots,T+1} whereby a time-varying exposure and confounder (\eqn{A_t=a_t} and \eqn{L_t=l_t}) are measured over times \eqn{t=1,\ldots,T} and
#' a time varying outcome \eqn{Y_s=y_s} is measured over times \eqn{s=2,\ldots,T+1}. Define \eqn{c=s-t} as the number of time periods separating an exposure measurement, and subsequent outcome measurement.
#' By using the transform \eqn{t=s-c} gestmult estimates the causal parameters \eqn{\psi} of a SNMM of the form
#'  \deqn{E\{Y_s(\bar{a}_{s-c},0)-Y_s(\bar{a}_{s-c-1},0)|\bar{a}_{s-c-1},\bar{l}_{s-c}\}=\psi z_{sc}a_{s-c} \; \forall c=1,\ldots,T\; and\; \forall s>c}
#'  if Y is continuous or
#'  \deqn{\frac{E(Y_s(\bar{a}_{s-c},0)|\bar{a}_{s-c-1},\bar{l}_{s-c-1})}{E(Y_s(\bar{a}_{s-c-1},0)|\bar{a}_{s-c-1},\bar{l}_{s-c-1})}=exp(\psi z_{sc}a_{s-c}) \; \forall c=1,\ldots,T\; and \; \forall s>c }
#'  if Y is binary. The SNMMs form is defined by the parameter \eqn{z_{sc}}, which can be controlled by the input type as follows
#'  \itemize{
#'  \item{\code{type=1} }{Sets \eqn{z_{sc}=1}. Now \eqn{\psi} is now the effect of exposure (or unit increase in exposure) at any time t on all subsequent outcome periods.}
#'  \item{\code{type=2} }{Sets \eqn{z_{sc}=c(1,l_{t=s-c})}. Now \eqn{\psi=c(\psi_0,\psi_1)} where \eqn{\psi_0} is the effect of exposure (or unit increase in exposure) at any time t on all subsequent outcome periods, when \eqn{l_{t=s-c}=0} for all t, modified by
#'  \eqn{\psi_1} for each unit increase in \eqn{l_{t=s-c}} at all times t. If there is more than one confounder, the first confounder named in Lny is taken as the effect modifier. Note that effect modification
#'  is currently only supported for binary or continuous confounders.  }
#'  \item{\code{type=3} }{A time-varying causal effect can be posited for each value of \eqn{c}, that is the causal effect for the exposure on outcome \eqn{c} time periods later.
#'  We set \eqn{z_{sc}} to a vector of zeros of length T with a 1 in the \eqn{c=s-t}'th position. Now \eqn{\psi=c(\psi_{1},\ldots,\psi_{T})}
#'  where \eqn{\psi_(c)} is the effect of exposure on outcome \eqn{c} time periods later for all outcome periods \eqn{s>c} that is \eqn{A_{s-c}} on \eqn{Y_s}.}
#'  \item{\code{type=4} }{Sets \eqn{z_{sc}} to a vector of zeros of length T with \eqn{c(1,l_{s-c})} in the \eqn{c=s-t}'th position.
#'  Now \eqn{\psi=(\underline{\psi_1},\ldots,\underline{\psi_T})} where \eqn{\underline{\psi_c}=c(\psi_{0c},\psi_{1c})}. Here \eqn{\psi_{0c}} is the effect of exposure on outcome \eqn{c} time periods later, given \eqn{l_{s-c}=0} for all outcome periods \eqn{s>c}, modified by
#'  \eqn{\psi_{1c}} for each unit increase in \eqn{l_{s-c}} for all \eqn{s>c}. If there is more than one confounder, the first confounder named in Lny is taken as the effect modifier. Note that effect modification
#'  is currently only supported for binary or continuous confounders.}
#'  }
#'  The data must be in long format, where we assume the convention that each row with \code{time=t} contains \eqn{A_t,L_t} and \eqn{C_{t+1},Y_{t+1}}. That is the censoring indicator for each row
#'  should indicate that a user is censored AFTER time t and the outcome indicates the first outcome that occurs AFTER \eqn{a_t} and \eqn{l_t} are measured.
#'  For example, data at time 1, should contain \eqn{a_1}, \eqn{l_1}, \eqn{y_{s=2}}, and optionally \eqn{C_2}. If either A or Y are binary, they must be written as numeric vectors taking values either 0 or 1.
#'  The same is true for any covariate that is used for effect modification.\cr
#'  The data must be rectangular with a row entry for every individual for each exposure time 1 up to T. Data rows after censoring should be empty apart from the ID and time variables. This can be done using the function \code{\link{FormatData}}.\cr
#'  By default the censoring, propensity and outcome models include the exposure history at the previous time as a variable. One may consider also including all previous exposure and confounder history as variables in \code{Lny},\code{Lnp}, and \code{LnC} if necessary.\cr
#'  Censoring weights are handled as described in Sjolander and Vansteelendt (2016). Note that it is necessary that any variable included in \code{LnC} must also be in \code{Lnp}. Missing data not due to censoring is handled by removing rows with missing data prior to fitting the model. If outcome models fail to fit, consider removing covariates from \code{Lny} but keeping
#'  them in \code{Lnp} to reduce collinearity issues, or consider the sparseness of the data.
#'
#'
#' @return List of the fitted causal parameters of the posited SNMM. These are labeled as follows for each SNMM type, where An is
#' set to the name of the exposure variable, i is the current value of c, and Lny[1] is set to the name of the first confounder in \code{Lny}.
#' \item{\code{type=1} }{An: The effect of exposure at any time t on outcome at all subsequent times.}
#' \item{\code{type=2} }{An: The effect of exposure on outcome at any time t, when Ln[1] is set to zero, on all subsequent outcome times.\cr
#' An:Ln[1]: The effect modification by Lny[1], the additional effect of A on all subsequent Y for each unit increase in Lny[1] at all times t. }
#' \item{\code{type=3} }{c=i.An: The effect of exposure at any time t on outcome \code{c=i} time periods later.}
#' \item{\code{type=4} }{c=i.An: The effect of exposure at any time t on outcome \code{c=i} time periods later, when Ln[1] is set to zero.\cr
#' c=i.An:Ln[1]: The effect modification by Lny[1], the additional effect of exposure on outcome c=i time periods later for each unit increase in Lny[1]. }
#'
#' @references Vansteelandt, S., & Sjolander, A. (2016). Revisiting g-estimation of the Effect of a Time-varying Exposure Subject to Time-varying Confounding, Epidemiologic Methods, 5(1), 37-56. doi: \url{https://doi.org/10.1515/em-2015-0005}.
#' @references Dukes, O., & Vansteelandt, S. (2018). A Note on G-Estimation of Causal Risk Ratios, American Journal of Epidemiology, 187(5), 1079–1084. doi: \url{https://doi.org/10.1093/aje/kwx347}.
#'
#' @examples
#' datas<-dataexamples(n=10000,seed=123)
#' data=datas$datagestmult
#' idvar="id"
#' timevar="time"
#' Yn="Y"
#' An="A"
#' Ybin=FALSE
#' Abin=TRUE
#' Lny=c("L","U")
#' Lnp=c("L","U")
#' type=NA
#' cutoff=NA
#'
#' gestmult(data,idvar,timevar,Yn,An,Ybin,Abin,Lny,Lnp,type=1)
#' gestmult(data,idvar,timevar,Yn,An,Ybin,Abin,Lny,Lnp,type=2)
#' gestmult(data,idvar,timevar,Yn,An,Ybin,Abin,Lny,Lnp,type=3)
#' gestmult(data,idvar,timevar,Yn,An,Ybin,Abin,Lny,Lnp,type=4)
#' gestmult(data,idvar,timevar,Yn,An,Ybin,Abin,Lny,Lnp,type=3,cutoff=2)
#'
#' @importFrom DataCombine slide
#' @importFrom geeM geem
#' @importFrom nnet multinom
#' @importFrom tidyr unnest_legacy
#' @importFrom tidyr nest_legacy
#' @importFrom tibble as_tibble
#' @importFrom rsample bootstraps
#' @importFrom stats reshape rnorm rbinom gaussian glm predict reformulate
#' terms Gamma complete.cases quantile
#'
#' @export


gestmult<-function(data,idvar,timevar,Yn,An,Ybin,Abin,Lny,Lnp,type=1,Cn=NA,LnC=NA,cutoff=NA,...){
# Error messages
if(!is.data.frame(data))(stop("Either no data set has been given, or it is not in a data frame."))
if(any(!sapply(unlist(c("idvar","timevar","Yn","An","Ybin","Abin","Lny","Lnp","type")),exists)))(stop('Missing input'))
if(!is.numeric(data[,Lny[1]]) && type %in% c(2,4))(stop("Effect modification is only supported for a continuous covariate, or binary covariate written as an as.numeric() 0,1 vector"))
if(!is.numeric(data[,An]) && Abin==TRUE)(stop("A binary exposure A must be written as an as.numeric() 0,1 vector, with 1 indicating exposure. If not, consider using gestcat()"))
if(!is.numeric(data[,Yn]) && Ybin==TRUE)(stop("A binary outcome Y must be written as an as.numeric() 0,1 vector, with 1 indicating an event."))

# Define the ID and time variable
data$id<-data[,idvar]
data$time<-data[,timevar]
T<-max(data$time)
if(is.na(cutoff)==TRUE){
  cutoff<-T}

# Obtain the necessary variables of the data and create lagged exposure variable (A_{t-1})
if(is.na(Cn)==TRUE){
  d<-data[,unique(c(An,Yn,Lnp,Lny,idvar,timevar))]
  }else{
  d<-data[,unique(c(An,Yn,Lnp,Lny,idvar,timevar,LnC))]}

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
  d$w<-1
  } else{
  lmc<-reformulate(termlabels=c(LnC,"as.factor(time)"),response=Cn)
  modc<-glm(lmc,family="binomial",data=d)
  # Censoring scores
  cps<-1-predict(modc,type="response",newdata=d)
  d$cps<-cps
  #Indicator Function of whether C=0
  d[,paste(Cn,"0",sep="")]<-as.integer(!d[,Cn])
  #Setting up the denominator product of censoring scores "cprod" and initial weights "w"
  d$cprod<-d$cps
  d[is.na(d$cprod)==TRUE,"cprod"]<-1
  d$w<-d[,paste(Cn,"0",sep="")]/d$cps
  }

# Set up counterfactual outcomes for exposure at the preceding time period
# that is the "1 step" counterfactuals H_{s(s-1)}
d$H<-d[,Yn]
dc<-d
dc$cntstep<-1

# Create a copy of the data set with complete cases for initial estimate of psi
dcom<-d[complete.cases(d),]

# Set up an augmented dataset 'dc', adding additonal rows to hold the
# 2-step, up to c-step counterfactuals, with c=cutoff, that is the counterfactuals for time
# s-c under outcome s H_{s(s-c)}. Generate variable 'cntstep' to indicate the step length
# of the counterfactual the row holds.
for (i in 2:cutoff){
d2<-d[d$time %in% seq(1,T-(i-1),by=1),]
d2$cntstep<-i
dc<-rbind(dc,d2)
}
dc<-dc[order(dc$id,dc$time),]

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

# Get correct family for outcome models
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
  # Obtain the "1-step" causal effect as the initial estimate of psi
  # that is the effect of A_{s-1} on Y_s.
  out<-summary(geem(terms(lmy),family=family,id=dcom$id,data=dcom,weights=dcom$w))
  psi<-out$beta[match(par1,out$coefnames)]
  names(psi)<-par1

  # Obtain the 2-step counterfactuals and associated censoring weights
  # Set i as the maximum step length for which counterfactuals may be calculated
  # Set j as the current step length for which counterfactuals are being calculated, that is H_{s(s-j)}
  # set k as the current time period for which the j-step counterfactual is being calculated
  i<-2
	while(i<=cutoff && i<=T){
	  j<-2
  	while(j<=i) {
    	for(k in 1:(T-(j-1))){
    	  if(is.na(Cn)==FALSE){
          #Censoring weights
        	dc[dc$cntstep==j & dc$time==k,"cprod"]<-dc[dc$cntstep==(j-1) & dc$time==k,"cps"]*
        	dc[dc$cntstep==(j-1) & dc$time==(k+1),"cprod"]
        	dc[dc$cntstep==j & dc$time==k,"w"]<-d[d$time==(k+j-1),paste(Cn,"0",sep="")]/dc[dc$cntstep==j & dc$time==k,"cprod"]
        	}

    	  # Counterfactuals
        if (Ybin==FALSE){
      	  if (length(z)==1){
      	    dc[dc$cntstep==j & dc$time==k,"H"]<-dc[dc$cntstep==(j-1) & dc$time==(k+1),"H"]-
      	    dc[dc$cntstep==(j-1) & dc$time==(k+1),An]*psi*d[d$time==(k+1),z]
      	    }else{
      	    dc[dc$cntstep==j & dc$time==k,"H"]<-dc[dc$cntstep==(j-1) & dc$time==(k+1),"H"]-
      	    rowSums(dc[dc$cntstep==(j-1) & dc$time==(k+1),An]*sweep(dc[dc$cntstep==(j-1) & dc$time==(k+1),z],2,psi,"*"))
      	    }
      	  }
    	  if (Ybin==TRUE){
    	    if (length(z)==1){
    	      dc[dc$cntstep==j & dc$time==k,"H"]<-dc[dc$cntstep==(j-1) & dc$time==(k+1),"H"]*
    	      exp(dc[dc$cntstep==(j-1) & dc$time==(k+1),An]*-psi*d[d$time==(k+1),z])
            }else{
    	      dc[dc$cntstep==j & dc$time==k,"H"]<-dc[dc$cntstep==(j-1) & dc$time==(k+1),"H"]*
    	      exp(-rowSums(dc[dc$cntstep==(j-1) & dc$time==(k+1),An]*sweep(dc[dc$cntstep==(j-1) & dc$time==(k+1),z],2,psi,"*")))
    	      }
    	    }
    	}
  	# Increase j and calculate the counterfactuals of the next step length, until the maximum step
  	# length i is reached
  	j<-j+1
  	}
	# Obtain relevant rows for which counterfactuals have been calculated, that is
	# initially the 1 and 2 step counterfactuals
	dt<-dc[dc$cntstep %in% seq(1,i,by=1),]
	dtcom<-dt[complete.cases(dt),]
	# Update the estimate of psi including both the 1 step and 2 step counterfactuals
	out<-summary(geem(terms(lmH),family=family,id=dtcom$id,data=dtcom,weights=dtcom$w))
	psi<-out$beta[match(par1,out$coefnames)]
	names(psi)<-par1
	# Increase maximum step length by and repeat calculation of counterfactuals
	# up to step length 3 (and relevant censoring weights), and update psi. This is repeated until
	# either step length T or step length equal to cutoff is reached and the final estimate of psi is obtained.
	i<-i+1
	}

  return(list(psi=psi))
}

# Perform g-estimation for SNMM types 3 or 4 (time-varying psi)
else if (timevarying==TRUE){
  # Obtain estimate of psi_{s-1}
  out<-summary(geem(terms(lmy),family=family,id=dcom$id,data=dcom,weights=dcom$w))
  psi<-out$beta[match(par1,out$coefnames)]
  names(psi)<-par1
  psilist<-as.list(NULL)
  psilist[[1]]<-psi
  # Obtain counterfactuals and censoring weights as before
  # where each j step counterfactual is estimated using the estimate of
  # psi_{s-j+1}
  i<-2
	while(i<=cutoff && i<=T){
	  j<-2
	  while(j<=i) {
	    for(k in 1:(T-(j-1))){
	      if(is.na(Cn)==FALSE){
          dc[dc$cntstep==j & dc$time==k,"cprod"]<-dc[dc$cntstep==(j-1) & dc$time==k,"cps"]*
	        dc[dc$cntstep==(j-1) & dc$time==(k+1),"cprod"]
          dc[dc$cntstep==j & dc$time==k,"w"]<-d[d$time==(k+j-1),paste(Cn,"0",sep="")]/dc[dc$cntstep==j & dc$time==k,"cprod"]
	        }
		    if (Ybin==FALSE){
	        if (length(z)==1){
	          dc[dc$cntstep==j & dc$time==k,"H"]<-dc[dc$cntstep==(j-1) & dc$time==(k+1),"H"]-
	          dc[dc$cntstep==(j-1) & dc$time==(k+1),An]*psilist[[j-1]]*d[d$time==(k+1),z]
            }else{
	          dc[dc$cntstep==j & dc$time==k,"H"]<-dc[dc$cntstep==(j-1) & dc$time==(k+1),"H"]-
	          rowSums(dc[dc$cntstep==(j-1) & dc$time==(k+1),An]*sweep(dc[dc$cntstep==(j-1) & dc$time==(k+1),z],2,psilist[[j-1]],"*"))
	          }
	        }
	      if (Ybin==TRUE){
	        if (length(z)==1){
	          dc[dc$cntstep==j & dc$time==k,"H"]<-dc[dc$cntstep==(j-1) & dc$time==(k+1),"H"]*
	          exp(dc[dc$cntstep==(j-1) & dc$time==(k+1),An]*-psilist[[j-1]]*d[d$time==(k+1),z])
            }else{
	          dc[dc$cntstep==j & dc$time==k,"H"]<-dc[dc$cntstep==(j-1) & dc$time==(k+1),"H"]*
	          exp(-rowSums(dc[dc$cntstep==(j-1) & dc$time==(k+1),An]*sweep(dc[dc$cntstep==(j-1) & dc$time==(k+1),z],2,psilist[[j-1]],"*")))
	          }
	        }
	      }

	  j<-j+1
	  }
	# Obtain relevant data and calculate psi_{s-2}
	dt<-dc[dc$cntstep %in% i,]
	dtcom<-dt[complete.cases(dt),]
	if (i==T){
	  out<-summary(geem(terms(lmH1),family=family,id=dtcom$id,data=dtcom,weights=dtcom$w))
	  psi<-out$beta[match(par1,out$coefnames)]
	  names(psi)<-par1
	  psilist[[i]]<-psi
	  }else{
	  out<-summary(geem(terms(lmH,keep.order=T),family=family,id=dtcom$id,data=dtcom,weights=dtcom$w))
	  psi<-out$beta[match(par1,out$coefnames)]
	  names(psi)<-par1
	  psilist[[i]]<-psi
	  }
	# Set maximum step length up by one and repeat algorithm until
	# psi_{s-T} or psi_{s-c} where c is set to the value of cutoff is calculated.
	i<-i+1
	}
  # Create relevant names for outputs
  nam<-as.vector(NULL)
  for (p in 1:cutoff){
    nam[p]<-paste("c=",p,sep="")
  }
  names(psilist)<-nam
  results<-unlist(psilist)
  return(list(psi=results))
  }
}

