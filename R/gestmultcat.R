#' G-estimation for a time-varying outcome and categorical time-varying exposure.
#'
#' Performs G-estimation for a structural nested mean model (SNMM), based on the outcome regression methods described
#' in Sjolander and Vansteelendt (2016) and Dukes and Vansteelendt (2018). We assume a dataset with a time-varying outcome that is either binary or continuous,
#' time-varying counfounding, and a categorical time-varying exposure of three or more categories.
#'
#'@param data A data table in long format containing the data to be analysed. See description for details.
#'@param idvar Character string specifying the name of the ID variable in data.
#'@param timevar Character string specifying the name of the time variable in the data.
#'  Note that timevar must specify time periods as integer values starting from 1 (must not begin at 0).
#'@param Yn Character string specifying the name of the outcome variable.
#'@param An Character string specifying the name of the exposure variable.
#'@param Ybin True or False indicator of whether the outcome is binary.
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
#' a time varying outcome \eqn{Y_s=y_s} is measured over times \eqn{s=2,\ldots,T+1}. Define \eqn{c=s-t} as the number of time periods separating an exposure measurement, and subsequent outcome measurement. Also suppose that \eqn{A_t=a_t} is a categorical variable consisting of \eqn{k>2} categories.
#' These categories may take any arbitrary list of names, but we assume for theory purposes they are labeled as \eqn{j=0,1\ldots,k}
#' where \eqn{j=0} denotes no exposure, or some reference category. Define binary variables \eqn{A_t^j=a_t^j} \eqn{j=0,1,\ldots,k} where \eqn{A_t^j=1} if \eqn{A_t=j} and 0 otherwise.
#' By using the transform \eqn{t=s-c} gestmult estimates the causal parameters \eqn{\psi} of a SNMM of the form
#'  \deqn{E\{Y_s(\bar{a}_{s-c},a^0)-Y_s(\bar{a}_{s-c-1},a^0)|\bar{a}_{s-c-1},\bar{l}_{s-c}\}=\sum_{j=1}^{k}\psi^j z_{sc}a^j_{s-c} \; \forall c=1,\ldots,T\; and\; \forall s>c}
#'  if Y is continuous or
#'  \deqn{\frac{E(Y_s(\bar{a}_{s-c},a^0)|\bar{a}_{s-c-1},\bar{l}_{s-c-1})}{E(Y_s(\bar{a}_{s-c-1},a^0)|\bar{a}_{s-c-1},\bar{l}_{s-c-1})}=exp(\sum_{j=1}^{k}\psi^j z_{sc}a^j_{s-c}) \; \forall c=1,\ldots,T\; and \; \forall s>c }
#'  if Y is binary. The SNNM fits a separate set of causal parameters \eqn{\psi^j}, for the effect of exposure at category \eqn{j} on outcome,
#'  compared to exposure at the reference category 0, for each non-reference category. The models form is defined by the parameter \eqn{z_{sc}}, which can be controlled by the input type as follows
#'  \itemize{
#'  \item{\code{type=1:} }{Sets \eqn{z_{sc}=1}. Now \eqn{\psi^j} is now the effect of exposure when set to category \eqn{j}, compared to when set to the reference category, at any time t on all subsequent outcome periods.}
#'  \item{\code{type=2:} }{Sets \eqn{z_{sc}=c(1,l_{t=s-c})}. Now \eqn{\psi^j=c(\psi^j_0,\psi^j_1)} where \eqn{\psi^j_0} is the effect of exposure when set to category \eqn{j}, compared to when set to the reference category, at any time t on all subsequent outcome periods when \eqn{l_{t=s-c}=0} for all t, modified by
#'  \eqn{\psi^j_1} for each unit increase in \eqn{l_{t=s-c}} at all times t. If there is more than one confounder, the first confounder named in Lny is taken as the effect modifier. Note that effect modification
#'  is currently only supported for binary or continuous confounders.}
#'  \item {\code{type=3:} }{A time-varying causal effect can be posited for each value of \eqn{c}, that is the causal effect for the exposure on outcome \eqn{c} time periods later.
#'  We set \eqn{z_{sc}} to a vector of zeros of length T with a 1 in the \eqn{c=s-t}'th position. Now \eqn{\psi^j=c(\psi^j_{1},\ldots,\psi^j_{T})}
#'  where \eqn{\psi^j_(c)} is the effect of exposure, when set to category \eqn{j}, on outcome \eqn{c} time periods later for all \eqn{s>c} that is \eqn{A^j_{s-c}} on \eqn{Y_s} for all \eqn{s>c}.}
#'  \item{\code{type=4:} }{Sets \eqn{z_{sc}} to a vector of zeros of length T with \eqn{c(1,l_{s-c})} in the \eqn{c=s-t}'th position.
#'  Now \eqn{\psi^j=(\underline{\psi^j_1},\ldots,\underline{\psi^j_T})} where \eqn{\underline{\psi^j_c}=c(\psi^j_{0c},\psi^j_{1c})}. Here \eqn{\psi^j_{0c}} is the effect of exposure when set to category \eqn{j} on outcome \eqn{c} time periods later for all \eqn{s>c}, given \eqn{l_{s-c}=0}, modified by
#'  \eqn{\psi^j_{1c}} for each unit increase in \eqn{l_{s-c}} for all \eqn{s>c}. If there is more than one confounder, the first confounder named in Lny is taken as the effect modifier. Note that effect modification
#'  is currently only supported for binary or continuous confounders.}
#'  }
#' The data must be in long format, where we assume the convention that each row with \code{time=t} contains \eqn{A_t,L_t} and \eqn{C_{t+1},Y_{t+1}}. That is the censoring indicator for each row
#' should indicate that a user is censored AFTER time t, and the outcome the first outcome that occurs AFTER \eqn{a_t} and \eqn{l_t} are measured.
#' For example, data at time 1, should contain \eqn{a_1}, \eqn{l_1}, \eqn{y_{2}}, and optionally \eqn{C_2}. If Y is binary, it must be written as a numeric vector taking values either 0 or 1.
#' The same is true for any covariate that is used for effect modification.\cr
#' The data must be rectangular with a row entry for every individual for each exposure time 1 up to T. Data rows after censoring should be empty apart from the ID and time variables. This can be done using the function \code{\link{FormatData}}.\cr
#' By default the censoring, propensity and outcome models include the exposure history at the previous time as a variable. One may consider also including all previous exposure and confounder history as variables in \code{Lny},\code{Lnp}, and \code{LnC} if necessary.\cr
#' Censoring weights are handled as described in Sjolander and Vansteelendt (2016). Note that it is necessary that any variable included in \code{LnC} must also be in \code{Lnp}. Missing data not due to censoring is handled by removing rows with missing data prior to fitting the model. If outcome models fail to fit, consider removing covariates from \code{Lny} but keeping
#' them in \code{Lnp} to reduce collinearity issues, or consider the sparseness of the data.
#'
#' @return List of the fitted causal parameters of the posited SNMM. These are labeled as follows for each SNMM type, where An is
#' set to the name of the exposure variable, i is the current value of c, j is the category level, and Lny[1] is set to the name of the first confounder in \code{Lny}.
#' \item{\code{type=1} }{Anj: The effect of exposure at category j, at any time t on all subsequent outcome times.}
#' \item{\code{type=2} }{Anj: The effect of exposure at category j on outcome at any time t, when Ln[1] is set to zero, on all subsequent outcome times.\cr
#' Anj:Ln[1]: The effect modification by Lny[1], the additional effect of A at category j on all subsequent Y for each unit increase in Lny[1].}
#' \item{\code{type=3} }{c=i.Anj: The effect of exposure at any time t on outcome \code{c=i} time periods later.}
#' \item{\code{type=4} }{c=i.Anj: The effect of exposure at any time t on outcome \code{c=i} time periods later, when Ln[1] is set to zero.\cr
#' c=i.Anj:Ln[1]: The effect modification by Lny[1], the additional effect of exposure at category j, on outcome c=i time periods later for each unit increase in Lny[1]. }
#'
#' @references Vansteelandt, S., & Sjolander, A. (2016). Revisiting g-estimation of the Effect of a Time-varying Exposure Subject to Time-varying Confounding, Epidemiologic Methods, 5(1), 37-56. doi: \url{https://doi.org/10.1515/em-2015-0005}.
#' @references Dukes, O., & Vansteelandt, S. (2018). A Note on G-Estimation of Causal Risk Ratios, American Journal of Epidemiology, 187(5), 1079–1084. doi: \url{https://doi.org/10.1093/aje/kwx347}.
#'
#' @examples
#' datas<-dataexamples(n=10000,seed=123)
#' data=datas$datagestmultcat
#' #A is a categorical variables with categories labeled 1,2 and 3, with 1 the
#' #reference category
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
#' gestmultcat(data,idvar,timevar,Yn,An,Ybin,Lny,Lnp,type=1)
#' gestmultcat(data,idvar,timevar,Yn,An,Ybin,Lny,Lnp,type=2)
#' gestmultcat(data,idvar,timevar,Yn,An,Ybin,Lny,Lnp,type=3)
#' gestmultcat(data,idvar,timevar,Yn,An,Ybin,Lny,Lnp,type=4)
#' gestmultcat(data,idvar,timevar,Yn,An,Ybin,Lny,Lnp,type=3,cutoff=2)
#'
#' @importFrom DataCombine slide
#' @importFrom geeM geem
#' @importFrom nnet multinom
#' @importFrom tidyr nest_legacy
#' @importFrom tibble as_tibble
#' @importFrom rsample bootstraps
#' @importFrom stats reshape rnorm rbinom gaussian glm predict reformulate
#' terms Gamma complete.cases quantile
#'
#' @export

gestmultcat<-function(data,idvar,timevar,Yn,An,Ybin,Lny,Lnp,type=1,Cn=NA,LnC=NA,cutoff=NA,...){
# Error messages
if(!is.data.frame(data))(stop("Either no data set has been given, or it is not in a data frame."))
if(any(!sapply(unlist(c("idvar","timevar","Yn","An","Ybin","Lny","Lnp","type")),exists)))(stop('Missing input'))
if(!is.numeric(data[,Lny[1]]) && type %in% c(2,4))(stop("Effect modification is only supported for a continuous covariate, or binary covariate written as an as.numeric() 0,1 vector"))
if(!is.factor(data[,An]))(stop("Exposure A must be an as.factor() categorical variable of 2 or more categories."))
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

# Set lagged exposure at time 1 to the reference category
d[d$time %in% 1,Alagn]<-levels(d[,An])[1]
d[,Alagn]<-as.factor(d[,Alagn])
d$int<-1

#Propensity Scores
if (T==1){
  lmp<-reformulate(termlabels=c(Lnp),response=An)
} else {
  lmp<-reformulate(termlabels=c(Alagn,"as.factor(time)",Lnp),response=An)
}
modp<-multinom(lmp,data=d)
props<-predict(modp,type="probs",newdata=d)
d$prs<-props[,-1]


# Calculate initial censoring weights if Cn and LnC are supplied
if(is.na(Cn)==TRUE){
  d$w<-1
  }else{
  lmc<-reformulate(termlabels=c(LnC,"as.factor(time)"),response=Cn)
  modc<-glm(lmc,family="binomial",data=d)
  # Censoring scores
  cps<-1-predict(modc,type="response",newdata=d)
  d$cps<-cps
  # Indicator Function of whether C=0
  d[,paste(Cn,"0",sep="")]<-as.integer(!d[,Cn])
  # Setting up the denominator product of censoring weight "cprod" and weights "w"
  d$cprod<-d$cps
  # Set cprod to 1 where it is missing to allow for censoring weights
  d[is.na(d$cprod)==TRUE,"cprod"]<-1
  d$w<-d[,paste(Cn,"0",sep="")]/d$cps
  }

# Create a copy of the data set with complete cases for initial estimate of psi
dcom<-d[complete.cases(d),]

##Set up 1 step counterfactuals
d$H<-d[,Yn]
dc<-d
dc$cntstep<-1

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
#Create types
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

#Get correct family for outcome models
if(Ybin==TRUE){
  family<-Gamma(link="log")
  } else {
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
  # Obtain initial estimates of psi^l from the one step counterfactuals
  # and draw the correct parameters corresponding to psi from the outcome model
  out<-summary(geem(terms(lmy),family=family,id=dcom$id,data=dcom,weights=dcom$w))
  nam1<-paste(An,levels(d[,An])[-1],sep="")
  nam2<-apply(expand.grid(nam1,z[-1]), 1, paste, collapse=":")
  Acoef<-c(nam1,nam2)
  psi0<-out$beta[match(Acoef,out$coefnames)]
  names(psi0)<-Acoef
  #Place the estimate of psi for each non reference category of A (psi^l)
  #into a separate element in list psicat
  psicat<-as.list(NULL)
  for (l in 2:nlevels(d[,An])){
    psicat[[l]]<-psi0[grep(levels(d[,An])[l],Acoef)]
    }
  #set the value of psi at the reference level to zero
  psicat[[1]]<-rep(0,length(psicat[[2]]))

  # Need to obtain psi^l*z_{sc}*A^l for all relevant data points before calculating H
  # call this variable 'psiZA'.
  # i=Maximum value of 'cntstep' for which 'psiZA' is being calculated
  # j=Current value of 'cntstep' for which 'psiZA' is being calculated
  # l<-category level of exposure variable
  i<-2
  while(i<=cutoff && i<=T){
    j<-1
    dc$psiZA<-0
    while(j<=(i-1)) {
      if (length(z)==1){
	      for (l in 1:nlevels(dc[,An])){
	      dc[dc$cntstep==j & dc[,An]==levels(dc[,An])[l] & !is.na(dc[,An]),"psiZA"]<-psicat[[l]]}
	      }else{
	      for (l in 1:nlevels(dc[,An])){
	      dc[dc$cntstep==j & dc[,An]==levels(dc[,An])[l] & !is.na(dc[,An]),"psiZA"]<-rowSums(
	      sweep(dc[dc$cntstep==j & dc[,An]==levels(dc[,An])[l] & !is.na(dc[,An]),z],2,psicat[[l]],"*"))}
	      }
    j<-j+1
    }
    # Obtain the 2-step counterfactuals and associated censoring weights
    # Set i as the maximum step length for which counterfactuals may be calculated
    # Set j as the current step length for which counterfactuals are being calculated, that is H_{s(s-j)}
    # set k as the current time period for which the j-step counterfactual is being calculated
	  j<-2
    while(j<=i) {
	    for(k in 1:(T-(j-1))){
	      # Censoring Weights
        if(is.na(Cn)==FALSE){
	        dc[dc$cntstep==j & dc$time==k,"cprod"]<-dc[dc$cntstep==(j-1) & dc$time==k,"cps"]*
	        dc[dc$cntstep==(j-1) & dc$time==(k+1),"cprod"]
	        dc[dc$cntstep==j & dc$time==k,"w"]<-d[d$time==(k+j-1),paste(Cn,"0",sep="")]/dc[dc$cntstep==j & dc$time==k,"cprod"]
      	}

	      # Counterfactuals
	      if (Ybin==FALSE){
	        dc[dc$cntstep==j & dc$time==k,"H"]<-dc[dc$cntstep==(j-1) & dc$time==(k+1),"H"]-
	        dc[dc$cntstep==(j-1) & dc$time==(k+1),"psiZA"]
	      }
    	  if (Ybin==TRUE){
	        dc[dc$cntstep==j & dc$time==k,"H"]<-dc[dc$cntstep==(j-1) & dc$time==(k+1),"H"]*
	        exp(-dc[dc$cntstep==(j-1) & dc$time==(k+1),"psiZA"])
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
	  # Update the estimates of psi^l including both the 1 step and 2 step counterfactuals
	  out<-summary(geem(terms(lmH),family=family,id=dtcom$id,data=dtcom,weights=dtcom$w))
	  psi0<-out$beta[match(Acoef,out$coefnames)]
	  names(psi0)<-Acoef
	  # Update psicat
	  psicat<-as.list(NULL)
	  for (l in 2:nlevels(d[,An])){
	    psicat[[l]]<-psi0[grep(levels(d[,An])[l],Acoef)]
	  }
	  psicat[[1]]<-rep(0,length(psicat[[2]]))
	  # Increase maximum step length by and repeat calculation of counterfactuals
	  # up to step length 3 (and relevant censoring weights), and update psi. This is repeated until
	  # either step length T or step length equal to cutoff is reached and the final estimate of psi is obtained.
	  i<-i+1
	}
  # Output results
  results<-unlist(psicat[-1])
  return(list(psi=results))
}
# Perform g-estimation for SNMM types 3 or 4 (time-varying psi)
else if (timevarying==TRUE){
  # Obtain initial estimates of psi^l_{s-1}, the causal effects of exposure
  # on the subsequent outcome.
  out<-summary(geem(terms(lmy),family=family,id=dcom$id,data=dcom,weights=dcom$w))
  nam1<-paste(An,levels(d[,An])[-1],sep="")
  nam2<-apply(expand.grid(nam1,z[-1]), 1, paste, collapse=":")
  Acoef<-c(nam1,nam2)
  psi0<-out$beta[match(Acoef,out$coefnames)]
  names(psi0)<-Acoef
  # Create psicatlist to hold the time-varying values of psicat for
  # each step length c and psicatresult to hold the same list
  # with the reference category estimates (set to 0) removed.
  psicatlist<-as.list(NULL)
  psicatresult<-as.list(NULL)
  psicat<-as.list(NULL)
  for (l in 2:nlevels(d[,An])){
    psicat[[l]]<-psi0[grep(levels(d[,An])[l],Acoef)]
  }
  psicat[[1]]<-rep(0,length(psicat[[2]]))
  psicatlist[[1]]<-psicat
  psicatresult[[1]]<-psicat[-1]

  i<-2
  # Need to obtain psi^l*z_{sc}*A^l for all relevant data points before calculating H
  # call this variable 'psiZA'.
  # i=Maximum value of 'cntstep' for which 'psiZA' is being calculated
  # j=Current value of 'cntstep' for which 'psiZA' is being calculated
  # l<-category level of exposure variable
	while(i<=cutoff && i<=T){
    j<-1
    dc$psiZA<-0
    while(j<=(i-1)) {
      if (length(z)==1){
	      for (l in 1:nlevels(dc[,An])){
	      dc[dc$cntstep==j & dc[,An]==levels(dc[,An])[l] & !is.na(dc[,An]),"psiZA"]<-psicatlist[[j]][[l]]}
	      }else{
	      for (l in 1:nlevels(dc[,An])){
	      dc[dc$cntstep==j & dc[,An]==levels(dc[,An])[l] & !is.na(dc[,An]),"psiZA"]<-rowSums(
	      sweep(dc[dc$cntstep==j & dc[,An]==levels(dc[,An])[l] & !is.na(dc[,An]),z],2,psicatlist[[j]][[l]],"*"))}
	      }
    j<-j+1
    }

  # Obtain the 2-step counterfactuals and associated censoring weights
  # Set i as the maximum step length for which counterfactuals may be calculated
  # Set j as the current step length for which counterfactuals are being calculated, that is H_{s(s-j)}
  # set k as the current time period for which the j-step counterfactual is being calculated
	j<-2
	while(j<=i) {
	  for(k in 1:(T-(j-1))){
	    # Censoring Weights
	    if(is.na(Cn)==FALSE){
	      dc[dc$cntstep==j & dc$time==k,"cprod"]<-dc[dc$cntstep==(j-1) & dc$time==k,"cps"]*
	      dc[dc$cntstep==(j-1) & dc$time==(k+1),"cprod"]
	      dc[dc$cntstep==j & dc$time==k,"w"]<-d[d$time==(k+j-1),paste(Cn,"0",sep="")]/dc[dc$cntstep==j & dc$time==k,"cprod"]
	    }
	    # Counterfactuals
	    if (Ybin==FALSE){
	      dc[dc$cntstep==j & dc$time==k,"H"]<-dc[dc$cntstep==(j-1) & dc$time==(k+1),"H"]-
	      dc[dc$cntstep==(j-1) & dc$time==(k+1),"psiZA"]
	    }
	    if (Ybin==TRUE){
	      dc[dc$cntstep==j & dc$time==k,"H"]<-dc[dc$cntstep==(j-1) & dc$time==(k+1),"H"]*
	      exp(-dc[dc$cntstep==(j-1) & dc$time==(k+1),"psiZA"])
	    }
	  }
  j<-j+1
	}
	# Obtain relevant data and obtain an estimates of psi^l_{s-2}
	dt<-dc[dc$cntstep %in% i,]
	dtcom<-dt[complete.cases(dt),]

	if (i==T){
	  out<-summary(geem(terms(lmH1),family=family,id=dtcom$id,data=dtcom,weights=dtcom$w))
	  psi0<-out$beta[match(Acoef,out$coefnames)]
	  psimar<-out$beta[match(z[-1],out$coefnames)]
	  }else{
	  out<-summary(geem(terms(lmH),family=family,id=dtcom$id,data=dtcom,weights=dtcom$w))
	  psi0<-out$beta[match(Acoef,out$coefnames)]
	}
	names(psi0)<-Acoef
	psicat<-as.list(NULL)
	for (l in 2:nlevels(d[,An])){
	  psicat[[l]]<-psi0[grep(levels(d[,An])[l],Acoef)]
	}
	psicat[[1]]<-rep(0,length(psicat[[2]]))
	psicatlist[[i]]<-psicat
	psicatresult[[i]]<-psicat[-1]
	# Set maximum step length up by one and repeat algorithm until
	# psi_{s-T} or psi_{s-c} where c is set to the value of cutoff is calculated.
	i<-i+1
	}
  # Obtain relevant names for output
  nam<-as.vector(NULL)
  for (p in 1:cutoff){
    nam[p]<-paste("s-",p,sep="")
  }
  names(psicatresult)<-nam
  results<-unlist(psicatresult)
  return(list(psi=results))
	}
}
