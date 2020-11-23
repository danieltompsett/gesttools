#' Percentile based bootstrap confidence intervals
#'
#' Generates percentile based confidence intervals for the causal parameters
#' of a fitted SNMM by gest functions. Bonferroni corrected confidence
#' intervals are also reported for multiple comparisons.
#'
#' @param gestfunc Name (without quotations) of the gest function to run.
#' One of gest,gestmult,gestcat or gestmultcat.
#' @param data,idvar,timevar,Yn,An,Ybin,Abin,Lny,Lnp,type,Cn,LnC,cutoff
#' Same arguments as in gest functions, to be input into gestfunc.
#' @param bn Number of bootstrapped datasets
#' @param alpha Confidence level of confidence intervals
#' @param onesided Controls the type of confidence interval generated. Takes one of three inputs, \code{"upper"} for upper one-sided confidence intervals,
#' \code{"lower"} for lower one-sided confidence intervals, and \code{"twosided"} for two-sided confidence intervals. Defaults to \code{"twosided"}.
#' @param seed Integer specifying the random seed for generation of bootstrap samples.
#' @param ... additional arguments
#'
#' @return Returns a list of the following four elements.
#' \item{t0 }{The fitted causal parameters of the original data.}
#' \item{mean.boot }{The average fitted values of the causal parameters of the bootstrapped datasets.}
#' \item{conf }{The upper and/or lower bounds of a \eqn{1-\alpha} confidence interval for each causal parameter
#' comprising \eqn{\psi}.}
#' \item{conf.Bonferroni }{The upper and/or lower bounds of a Bonferroni corrected confidence
#' interval for each causal parameter comprising \eqn{psi}, used for multiple comparisons.}
#'
#' @examples
#' datas<-dataexamples(n=1000,seed=1234567)
#' data=datas$datagest
#' idvar="id"
#' timevar="time"
#' Yn="Y"
#' An="A"
#' Ybin=FALSE
#' Abin=TRUE
#' Lny=c("L")
#' Lnp=c("L")
#' gestfunc<-gest
#' type=NA
#'
#' gest.boot(gest,data,idvar,timevar,Yn,An,Ybin,Abin,Lny,Lnp,type=1,bn=10,alpha=0.05,onesided="lower",seed=123)
#'
#'
#' @export
gest.boot<-function(gestfunc,data,idvar,timevar,Yn,An,Ybin,Abin,Lny,Lnp,
type=1,Cn=NA,LnC=NA,cutoff=NA,bn,alpha,onesided="twosided",seed=123,...){


t0<-gestfunc(data=data,idvar=idvar,timevar=timevar,Yn=Yn,An=An,Ybin=Ybin,Abin=Abin,Lny=Lny,Lnp=Lnp,type=type,Cn=Cn,LnC=LnC,cutoff=cutoff)$psi
nams<-names(t0)
#Create tibble data based on ID
Data<-data %>% nest_legacy(-all_of(idvar))
set.seed(seed)
bs <- bootstraps(Data, times = bn)


results1<-as.list(NULL)
results<-as.list(NULL)
for (j in 1:bn){
tryCatch({
b<-as.data.frame(as_tibble(bs$splits[[j]]) %>% unnest_legacy())

b<-b[order(b[,idvar],b[,timevar]),]

results1[[j]]<-gestfunc(data=b,idvar=idvar,timevar=timevar,Yn=Yn,An=An,Ybin=Ybin,Abin=Abin,Lny=Lny,Lnp=Lnp,type=type,Cn=Cn,LnC=LnC,cutoff=cutoff)$psi
results[[j]]<-unlist(results1[[j]])
},error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

mean<-colMeans(do.call(rbind,results))

resultsmat<-do.call(rbind,results)

ci.quant<-function(x=NA){
  return(quantile(x,probs=c(alpha/2,1-alpha/2)))
}

ci.quant.bonf<-function(x=NA){
  return(quantile(x,probs=c(alpha/(2*length(unlist(t0))),1-alpha/(2*length(unlist(t0))))))
}

ci.quant.upper<-function(x=NA){
  return(quantile(x,probs=c(1-alpha)))
}

ci.quant.bonf.upper<-function(x=NA){
  return(quantile(x,probs=c(1-alpha/(length(unlist(t0))))))
}

ci.quant.lower<-function(x=NA){
  return(quantile(x,probs=c(alpha)))
}

ci.quant.bonf.lower<-function(x=NA){
  return(quantile(x,probs=c(alpha/(length(unlist(t0))))))
}


results.sort<-apply(resultsmat,2,sort)

if(onesided=="twosided"){
  conf.quant<-t(apply(results.sort,2,ci.quant))
  conf.quant.bonf<-t(apply(results.sort,2,ci.quant.bonf))

}else if(onesided=="upper"){
  conf.quant<-t(apply(results.sort,2,ci.quant.upper))
  conf.quant.bonf<-t(apply(results.sort,2,ci.quant.bonf.upper))

}else if(onesided=="lower"){
  conf.quant<-t(apply(results.sort,2,ci.quant.lower))
  conf.quant.bonf<-t(apply(results.sort,2,ci.quant.bonf.lower))
}




return(list(original=t0,mean.boot=mean,conf=conf.quant,
conf.Bonferroni=conf.quant.bonf))

}
