#' G-Estimation for an End of Study Outcome and Categorical Exposure Variable
#'
#' Performs g-estimation of a structural nested mean model (SNMM), based on the outcome regression methods described
#' in Sjolander and Vansteelandt (2016) and Dukes and Vansteelandt (2018). We expect a dataset with an end of study outcome that is either binary or continuous,
#' time-varying and/or baseline confounding, and a categorical time-varying exposure of two of more categories.
#'
#'@param data A data frame in long format containing the data to be analysed. See description for details.
#'@param idvar Character string specifying the name of the ID variable in data.
#'@param timevar Character string specifying the name of the time variable in the data.
#'  Note that timevar must specify time periods as integer values starting from 1 (must not begin at 0).
#'@param Yn Character string specifying the name of the end of study outcome variable.
#'@param An Character string specifying the name of the time-varying exposure variable.
#'@param Ybin TRUE or FALSE indicator of whether the outcome is binary.
#'@param Lny Vector of character strings specifying the names of the
#'  variables to be included in the outcome model in quotations.
#'@param Lnp Vector of character strings specifying the names of the variables
#'  to be included in the model calculating the propensity scores.
#'@param type Value from 1-4 specifying SNMM type to fit. See details.
#'@param Cn Optional character string specifying the name of the censoring indicator variable.
#'  Cn should be a numeric vector taking values 0 or 1, with 1 indicating censored.
#'@param LnC Vector of character strings specifying the names of the variables to be used in the censoring score model to calculate
#'  the censoring weights. Note that any variable in \code{LnC} should also be in \code{Lnp} for the validity of the censoring and propensity weights.
#'@param ... Additional arguments, currently not in use.
#'
#' @details Suppose a set of time varying confounders \eqn{L_t=l_t}, and a time-varying categorical exposure variable \eqn{A_t=a_t}, measured over time periods \eqn{t=1,\ldots,T}, with an end of study outcome \eqn{Y} measured at time \eqn{T+1}.
#' Also suppose that \eqn{A_t} holds data consisting of \eqn{k\geq 2} categories. These categories may take any arbitrary list of names, but we assume for theory purposes they are labeled as \eqn{j=0,1\ldots,k},
#' where \eqn{j=0} denotes no exposure, or some reference category. Define binary variables \eqn{A_t^j=a_t^j} \eqn{j=0,1,\ldots,k} where \eqn{A_t^j=1} if \eqn{A_t=j} and 0 otherwise.
#' Then \code{gestcat} fits a SNMM of the form
#'  \deqn{E(Y(\bar{a}_{t},a^0)-Y(\bar{a}_{t-1},a^0)|\bar{a}_{t-1},\bar{l}_{t-1})=\sum_{j=1}^{k}\psi^j z_ta^j_{t} \;\forall\; t=1,\ldots,T}
#'  if Y is continuous or
#'  \deqn{\frac{E(Y(\bar{a}_{t},a^0)|\bar{a}_{t-1},\bar{l}_{t-1})}{E(Y(\bar{a}_{t-1},a^0)|\bar{a}_{t-1},\bar{l}_{t-1})}=exp(\sum_{j=1}^{k}\psi^j z_ta^j_{t})\;\forall\; t=1,\ldots,T }
#'  if Y is binary. The SNNM fits a separate set of causal parameters \eqn{\psi^j}, for the effect of exposure at category \eqn{j} on outcome,
#'  compared to exposure at the reference category 0, for each non-reference category. The models form is defined by the parameter \eqn{z_t}, which can be controlled by the input \code{type} as follows
#'  \itemize{
#'  \item{\code{type=1} }{sets \eqn{z_t=1}. This implies that \eqn{\psi^j} is now the effect of exposure when set to category \eqn{j}, (compared the reference category) at any time t on Y.}
#'  \item{\code{type=2} }{sets \eqn{z_t=c(1,l_t)}, and adds affect modification by the first named variable in \code{Lny}, which we denote \eqn{L_t}.
#'  Now \eqn{\psi^j=c(\psi^j_0,\psi^j_1)} where \eqn{\psi^j_0} is the effect of exposure when set to category \eqn{j}, (compared the reference category) at any time t on Y when \eqn{l_t=0} for all t, modified by
#'  \eqn{\psi^j_1} for each unit increase in \eqn{l_t} for all t. Note that effect modification
#'  is currently only supported for binary or continuous confounders.}
#'  \item {\code{type=3} }{sets \eqn{z_t} to a vector of zeros of length T with a 1 in the t'th position. Now \eqn{\psi^j=c(\psi^j_1,\ldots,\psi^j_T)}
#'  where \eqn{\psi^j_t} is the effect of exposure, when set to category \eqn{j}, at time t on Y.}
#'  \item{\code{type=4} }{allows for a time-varying causal effect that can be modified by the first named variable in \code{Lny}, that is it allows for both time-varying effects and effect modification. It sets \eqn{z_t} to a vector of zeros of length T with \eqn{c(1,l_t)} in the t'th position.
#'  Now \eqn{\psi^j=(\underline{\psi^j_1},\ldots,\underline{\psi^j_T})} where \eqn{\underline{\psi^j_t}=c(\psi^j_{0t},\psi^j_{1t})}. Here \eqn{\psi^j_{0t}} is the effect of exposure when set to category \eqn{j} at time t on Y at \eqn{l_t=0}, modified by
#'  \eqn{\psi^j_{1t}} for each unit increase in \eqn{l_t}. Note that effect modification
#'  is currently only supported for binary or continuous confounders. }
#'  }
#'  The data must be in long format, where we assume the convention that each row with \code{time=t} contains \eqn{A_t,L_t} and \eqn{C_{t+1}}. That is the censoring indicator for each row
#'  should indicate that a user is censored AFTER time t. The end of study outcome should be repeated on each row. If Y is binary, it must be written as a numeric vector taking values either 0 or 1.
#'  The same is true for any covariate that is used for effect modification.\cr
#'  The data must be rectangular with a row entry for every individual for each exposure time 1 up to T. Data rows after censoring should be empty apart from the ID and time variables. This can be done using the function \code{\link{FormatData}}.\cr
#'  By default the censoring, propensity and outcome models include the exposure history at the previous time as a variable. One may consider also including all previous exposure and confounder history as variables in \code{Lny},\code{Lnp}, and \code{LnC} if necessary.\cr
#'  Censoring weights are handled as described in Sjolander and Vansteelandt (2016). Note that it is necessary that any variable included in \code{LnC} must also be in \code{Lnp}. Missing data not due to censoring is handled automatically by removing rows with missing data prior to fitting the model. If outcome models fail to fit, consider removing covariates from \code{Lny} but keeping
#'  them in \code{Lnp} to reduce collinearity issues, or consider the sparseness of the data.
#'
#' @return List of the fitted causal parameters of the posited SNMM. These are labeled as follows for each SNMM type, where An is
#' set to the name of the exposure variable, i is the current time period, j is the category level, and Lny[1] is set to the name of the first confounder in \code{Lny}.
#' \item{\code{type=1} }{Anj: The effect of exposure at category j at any time t on outcome. }
#' \item{\code{type=2} }{Anj: The effect of exposure at category j, at any time t on outcome, when Ln[1] is set to zero.\cr
#' Anj:Ln[1]: The effect modification by Lny[1], the additional effect of A at category j on Y for each unit increase in Lny[1]. }
#' \item{\code{type=3} }{t=i.Anj: The effect of exposure at category j, at time t=i on outcome.}
#' \item{\code{type=4} }{t=i.Anj: The effect of exposure at category j, at time t=i on outcome when Ln[1] is set to zero.\cr
#' t=i.Anj:Ln[1]: The effect modification by Lny[1], the additional effect of A at category j at time t=i on Y, for each unit increase in Lny[1].}
#'
#' @references Vansteelandt, S., & Sjolander, A. (2016). Revisiting g-estimation of the Effect of a Time-varying Exposure Subject to Time-varying Confounding, Epidemiologic Methods, 5(1), 37-56. <doi:10.1515/em-2015-0005>.
#' @references Dukes, O., & Vansteelandt, S. (2018). A Note on g-Estimation of Causal Risk Ratios, American Journal of Epidemiology, 187(5), 1079–1084. <doi:10.1093/aje/kwx347>.
#'
#' @examples
#' datas<-dataexamples(n=10000,seed=123,Censoring=FALSE)
#' data=datas$datagestcat
#' #A is a categorical variable with categories labeled 1,2 and 3, with 1 the
#' #reference category
#' idvar="id"
#' timevar="time"
#' Yn="Y"
#' An="A"
#' Ybin=FALSE
#' Lny=c("L","U")
#' Lnp=c("L","U")
#' Cn<-NA
#' LnC<-NA
#' type=NA
#'
#' gestcat(data,idvar,timevar,Yn,An,Ybin,Lny,Lnp,type=1)
#' gestcat(data,idvar,timevar,Yn,An,Ybin,Lny,Lnp,type=3)
#'
#' #Example with censoring
#' datas<-dataexamples(n=10000,seed=123,Censoring=TRUE)
#' data=datas$datagestcat
#' Cn="C"
#' LnC=c("L","U")
#' gestcat(data,idvar,timevar,Yn,An,Ybin,Lny,Lnp,type=3,Cn,LnC)
#'
#' @export

gestcat<-function(data,idvar,timevar,Yn,An,Ybin,Lny,Lnp,type=1,Cn=NA,LnC=NA,...){
# Error messages
if(!is.data.frame(data))(stop("Either no data set has been given, or it is not in a data frame."))
if(any(!sapply(unlist(c("idvar","timevar","Yn","An","Ybin","Lny","Lnp","type")),exists)))(stop('Missing input'))
if(!is.numeric(data[,Lny[1]]) && type %in% c(2,4))(stop("Effect modification is only supported for a continuous covariate, or binary covariate written as an as.numeric() 0,1 vector"))
if(!is.factor(data[,An]))(stop("Exposure A must be an as.factor() categorical variable of 2 or more categories."))
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

# Set lagged exposure at time 1 to the reference category

d[,Alagn]<-as.factor(d[,Alagn])
levels(d[,paste("L",1,An,sep="")])<-levels(d[,An])
d[d$time %in% 1,Alagn]<-levels(d[,An])[1]

d$int<-1

# Propensity Scores
if (T==1){
  lmp<-reformulate(termlabels=c(Lnp),response=An)
} else {
  lmp<-reformulate(termlabels=c(Alagn,"as.factor(time)",Lnp),response=An)
}
modp<-multinom(lmp,data=d)
props<-predict(modp,type="probs",newdata=d)
if(nlevels(d[,An])==2){
  d$prs<-props
  }else{d$prs<-props[,-1]}

# Calculate initial censoring weights if Cn and LnC are supplied
if(is.na(Cn)==TRUE){
  d$w<-1} else{
  lmc<-reformulate(termlabels=c(LnC,"as.factor(time)"),response=Cn)
  modc<-glm(lmc,family="binomial",data=d)
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
  # Obtain the causal effects of A^l_T on Y as the initial estimates of psi^l
  # and draw the correct parameters corresponding to psi^l from the outcome model
  dt<-dcom[dcom$time==T,]
  out<-summary(geem(terms(lmy),family=family,id=dt$id,data=dt,weights=dt$w))
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
  # i=Miminum time period for which 'psiZA' is being calculated
  # j=Current time period for which 'psiZA' is being calculated
  # l=category level of exposure variable
  i<-T-1
  while(i>=1){
    j<-T
    d$psiZA<-0
    while(j>=(i+1)){
      if(length(z)==1){
        for (l in 1:nlevels(d[,An])){
          d[d$time==j & d[,An]==levels(d[,An])[l] & !is.na(d[,An]),"psiZA"]<-psicat[[l]]}
          }else{
          for (l in 1:nlevels(d[,An])){
          d[d$time==j & d[,An]==levels(d[,An])[l] & !is.na(d[,An]),"psiZA"]<-rowSums(
          sweep(d[d$time==j & d[,An]==levels(d[,An])[l] & !is.na(d[,An]),z],2,psicat[[l]],"*"))}
          }
      j<-j-1
    }
    # Obtain counterfactual outcomes and censoring weights for time T-1
    # Set i as the minimum time period for which counterfactuals and censoring weights are calculated
    # Set k to the current time period being calculated.
	  k<-(T-1)
	    while(k>=i){
	      if(is.na(Cn)==FALSE){
	        # Censoring weights
          d[d$time==k,"cprod"]<-d[d$time==k+1,"cprod"]*d[d$time==k,"cps"]
	        d[d$time==k,"w"]<-d[d$time==T,paste(Cn,"0",sep="")]/d[d$time==k,"cprod"]
	        }
	      # Counterfactuals
	      if (Ybin==FALSE){
	        d[d$time==k,"H"]<-d[d$time==(k+1),"H"]-
	        d[d$time==(k+1),"psiZA"]
	        }
	      if (Ybin==TRUE){
	        d[d$time==k,"H"]<-d[d$time==(k+1),"H"]*
        	exp(-d[d$time==(k+1),"psiZA"])
	        }
	      # Current time period is lowered and censoring weights and counterfactuals
	      # calculated until current value of i is reached.
	      k<-k-1
	      }
	  # Obtain relevant time periods for which counterfactuals have been calculated
    dt<-d[d$time %in% seq(i,T,by=1),]
    dtcom<-dt[complete.cases(dt),]
    out<-summary(geem(terms(lmH),family=family,id=dtcom$id,data=dtcom,weights=dtcom$w))
    # Update the estimates of psi^l including data at time T-1.
    psi0<-out$beta[match(Acoef,out$coefnames)]
	  names(psi0)<-Acoef
	  psicat<-as.list(NULL)
    # Update psicat
	  for (l in 2:nlevels(d[,An])){
	  psicat[[l]]<-psi0[grep(levels(d[,An])[l],Acoef)]
	  }
	  psicat[[1]]<-rep(0,length(psicat[[2]]))
	  # Set minimum time period to time T-2 and repeat calculation of counterfactuals
	  # and censoring weights up to time T-2 and update psi. This is repeated until
	  # time period 1 is reached and the final estimate of psi is obtained.
    i<-i-1
    }
  results<-unlist(psicat[-1])
  return(list(psi=results))
  }
# Perform g-estimation for SNMM types 3 or 4 (time-varying psi)
else if (timevarying==TRUE){
  dt<-dcom[dcom$time==T,]
  out<-summary(geem(terms(lmy),family=family,id=dt$id,data=dt,weights=dt$w))
  nam1<-paste(An,levels(d[,An])[-1],sep="")
  nam2<-apply(expand.grid(nam1,z[-1]), 1, paste, collapse=":")
  Acoef<-c(nam1,nam2)
  psi0<-out$beta[match(Acoef,out$coefnames)]
  names(psi0)<-Acoef
  # Create psicatlist to hold the time-varying values of psicat for
  # each step length c and psicatresult to hold the same list
  # with the reference category estimates (set to 0) removed.
  psicat<-as.list(NULL)
  psicatlist<-as.list(NULL)
  psicatresult<-as.list(NULL)
  for (l in 2:nlevels(d[,An])){
  psicat[[l]]<-psi0[grep(levels(d[,An])[l],Acoef)]
  }
  psicat[[1]]<-rep(0,length(psicat[[2]]))
  psicatlist[[T]]<-psicat
  psicatresult[[T]]<-psicat[-1]

  # Need to obtain psi^l*z_{sc}*A^l for all relevant data points before calculating H
  # call this variable 'psiZA'.
  # i=Miminum time period for which 'psiZA' is being calculated
  # j=Current time period for which 'psiZA' is being calculated
  # l=category level of exposure variable
  i<-(T-1)
  while(i>=1){
  j<-T
  d$psiZA<-0
  while(j>=(i+1)){
    if(length(z)==1){
      for (l in 1:nlevels(d[,An])){
        d[d$time==j & d[,An]==levels(d[,An])[l] & !is.na(d[,An]),"psiZA"]<-psicatlist[[j]][[l]]}
        }else{
        for (l in 1:nlevels(d[,An])){
        d[d$time==j & d[,An]==levels(d[,An])[l] & !is.na(d[,An]),"psiZA"]<-rowSums(
        sweep(d[d$time==j & d[,An]==levels(d[,An])[l] & !is.na(d[,An]),z],2,psicatlist[[j]][[l]],"*"))}
        }
    j<-j-1
    }
    # Obtain counterfactual outcomes and censoring weights for time T-1
    # Set i as the minimum time period for which counterfactuals and censoring weights are calculated
    # Set k to the current time period being calculated.
	  k<-(T-1)
	  while(k>=i){
	    if(is.na(Cn)==FALSE){
	      d[d$time==k,"cprod"]<-d[d$time==k+1,"cprod"]*d[d$time==k,"cps"]
	      d[d$time==k,"w"]<-d[d$time==T,paste(Cn,"0",sep="")]/d[d$time==k,"cprod"]
        }
	    if (Ybin==FALSE){
	      d[d$time==k,"H"]<-d[d$time==(k+1),"H"]-
	      d[d$time==(k+1),"psiZA"]
	      }
	    if (Ybin==TRUE){
	      d[d$time==k,"H"]<-d[d$time==(k+1),"H"]*
	      exp(-d[d$time==(k+1),"psiZA"])
	      }
	    k<-k-1
	    }
	  # Obtain data at time T-1 for calculation of psi
    dt<-d[d$time %in% i,]
    dtcom<-dt[complete.cases(dt),]
    if (i==1){
      out<-summary(geem(terms(lmH1),family=family,id=dtcom$id,data=dtcom,weights=dtcom$w))
    }else{
      out<-summary(geem(terms(lmH),family=family,id=dtcom$id,data=dtcom,weights=dtcom$w))
      }
    psi0<-out$beta[match(Acoef,out$coefnames)]
	  names(psi0)<-Acoef
	  psicat<-as.list(NULL)
	  for (l in 2:nlevels(d[,An])){
	  psicat[[l]]<-psi0[grep(levels(d[,An])[l],Acoef)]
	  }
	  psicat[[1]]<-rep(0,length(psicat[[2]]))
    psicatresult[[i]]<-psicat[-1]
    psicatlist[[i]]<-psicat
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
  names(psicatresult)<-nam
  results<-unlist(psicatresult)
  return(list(psi=results))
  }
}


