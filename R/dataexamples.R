#' Generates simulated example datasets.
#'
#' The code simulates four datasets designed to demonstrate each of the four
#' G-estimation functions of the package. Used in the examples section for each function
#' in the user manual. Each dataset comprises of an outcome Y (time-varying or end of study), time-varying exposure A, time-varying confounder L,
#' a baseline unmeasured confounder U, and optionally a censoring indicator C over 3 time periods.
#'
#' @param n Number of individuals in the dataset.
#' @param seed Random seed used for data generation.
#' @param Censoring TRUE or FALSE indicator of whether to include a censoring indicator \code{C}.
#' If \code{Censoring=TRUE}, data entries for A, Y, L and U are set to missing after censoring.
#'
#' @return Returns a list of four datasets labeled datagest,datagestmult,
#' datagestcat, and datagestmultcat, designed to demonstrate the respective
#' functions
#'
#' @examples
#' datas<-dataexamples(n=100,seed=34567,Censoring=FALSE)
#' data<-datas$datagest
#' #Multiple outcome data with censoring
#' datas<-dataexamples(n=100,seed=34567,Censoring=TRUE)
#' data<-datas$datagestmultcat
#'
#' @export

dataexamples<-function(n=1000,seed=1234,Censoring=FALSE){

expit<-function(x){
exp(x)/(1+exp(x))
}

n<-n
set.seed(seed)
id<-seq(1,n,by=1)
U<-rnorm(n,0,1)
L.1<-rnorm(n,1+U,1)
a.1<-expit(1+0.1*L.1)
A.1<-rbinom(n,1,a.1)
Y.1<-rnorm(n,1+A.1+L.1+U,1)

L.2<-rnorm(n,1+(A.1/2)+L.1+U,1)
a.2<-expit(1+0.1*L.2+0.1*A.1)
A.2<-rbinom(n,1,a.2)
Y.2<-rnorm(n,1+(A.1/2)+A.2+L.1+L.2+U,1)

L.3<-rnorm(n,1+(A.2/2)+L.2+U,1)
a.3<-expit(1+0.1*L.3+0.1*A.2)
A.3<-rbinom(n,1,a.3)
Y.3<-rnorm(n,1+(A.2/2)+A.3+L.1+L.2+L.3+U,1)
#End of study Outcome
Y<-rnorm(n,1+(A.2/2)+A.3+L.1+L.2+L.3+U,1)

if(Censoring==TRUE){
  C.1<-rbinom(n,1,expit(-1+0.001*A.1+0.001*L.1))
  C.2<-rbinom(n,1,expit(-1+0.001*A.2+0.001*L.2))
  C.3<-rbinom(n,1,expit(-1+0.001*A.3+0.001*L.3))
  C.2[C.1==1]<-1
  C.3[C.2==1]<-1
  Y[C.3==1]<-NA
  Y.3[C.3==1]<-NA
  A.3[C.2==1]<-NA
  L.3[C.2==1]<-NA
  Y.2[C.2==1]<-NA
  A.2[C.1==1]<-NA
  L.2[C.1==1]<-NA
  Y.1[C.1==1]<-NA
  #Create end of study outcome data and set as long format
  dw<-as.data.frame(cbind(id,Y,A.1,A.2,A.3,L.1,L.2,L.3,C.1,C.2,C.3,U))
  dl<-reshape(dw,direction="long",varying=c("A.1","A.2","A.3",
                                            "L.1","L.2","L.3",
                                            "C.1","C.2","C.3"))
  #Order data by time and ID
  dl<-dl[order(dl$id,dl$time),]
  datagest<-dl

  dw<-as.data.frame(cbind(id,Y.1,Y.2,Y.3,A.1,A.2,A.3,L.1,L.2,L.3,C.1,C.2,C.3,U))
  dl<-reshape(dw,direction="long",varying=c("Y.1","Y.2","Y.3",
                                            "A.1","A.2","A.3",
                                            "L.1","L.2","L.3",
                                            "C.1","C.2","C.3"))
#Order data by time and ID
  dl<-dl[order(dl$id,dl$time),]
  datagestmult<-dl
  }else{
  #Create end of study outcome data and set as long format
  dw<-as.data.frame(cbind(id,Y,A.1,A.2,A.3,L.1,L.2,L.3,U))
  dl<-reshape(dw,direction="long",varying=c("A.1","A.2","A.3",
  "L.1","L.2","L.3"))
  #Order data by time and ID
  dl<-dl[order(dl$id,dl$time),]
  datagest<-dl

  dw<-as.data.frame(cbind(id,Y.1,Y.2,Y.3,A.1,A.2,A.3,L.1,L.2,L.3,U))
  dl<-reshape(dw,direction="long",varying=c("Y.1","Y.2","Y.3",
                                            "A.1","A.2","A.3",
                                            "L.1","L.2","L.3"))
  #Order data by time and ID
  dl<-dl[order(dl$id,dl$time),]
  datagestmult<-dl
  }


###Categorical A###########

n<-n
set.seed(seed)
id<-seq(1,n,by=1)
U<-rnorm(n,0,1)
L.1<-rnorm(n,1+U,1)
a.1<-expit(1+0.1*L.1)

A.1<-as.vector(NULL)
for (i in 1:n){
  A.1[i]<-sample(letters[1:3],1, replace=TRUE,
        prob=c(1-(3*a.1[i])/5,a.1[i]/5,2*(a.1[i]/5)))}
A.1<-as.factor(A.1)
A.1par<-as.vector(NULL)
A.1par[A.1==letters[1]]<-1
A.1par[A.1==letters[2]]<-2
A.1par[A.1==letters[3]]<-3

Y.1<-rnorm(n,1+A.1par+L.1+U,1)

L.2<-rnorm(n,1+A.1par+L.1+U,1)
a.2<-expit(1+0.1*L.2+A.1par)
A.2<-as.vector(NULL)
for (i in 1:n){
  A.2[i]<-sample(letters[1:3],1, replace=TRUE,
          prob=c(1-(3*a.2[i]/5),a.2[i]/5,2*(a.2[i]/5)))}
A.2<-as.factor(A.2)
A.2par<-as.vector(NULL)
A.2par[A.2==letters[1]]<-1
A.2par[A.2==letters[2]]<-2
A.2par[A.2==letters[3]]<-3

Y.2<-rnorm(n,1+A.1par+A.2par+L.1+L.2+U,1)

L.3<-rnorm(n,1+A.2par+L.2+U,1)
a.3<-expit(1+0.1*L.3+A.2par)
A.3<-as.vector(NULL)
for (i in 1:n){
  A.3[i]<-sample(letters[1:3],1, replace=TRUE,
          prob=c(1-(3*a.3[i]/5),a.3[i]/5,2*(a.3[i]/5)))}
A.3<-as.factor(A.3)
A.3par<-as.vector(NULL)
A.3par[A.3==letters[1]]<-1
A.3par[A.3==letters[2]]<-2
A.3par[A.3==letters[3]]<-3

Y.3<-rnorm(n,1+A.1par+A.2par+A.3par+L.1+L.2+L.3+U,1)

Y<-rnorm(n,1+A.1par+A.2par+A.3par+L.1+L.2+L.3+U,1)

if(Censoring==TRUE){
  C.1<-rbinom(n,1,expit(-1+0.05*A.1par+0.05*L.2))
  C.2<-rbinom(n,1,expit(-1+0.05*A.2par+0.05*L.2))
  C.3<-rbinom(n,1,expit(-1+0.05*A.3par+0.05*L.3))
  C.2[C.1==1]<-1
  C.3[C.2==1]<-1
  Y[C.3==1]<-NA
  Y.3[C.2==1]<-NA
  A.3[C.2==1]<-NA
  L.3[C.2==1]<-NA
  Y.2[C.2==1]<-NA
  A.2[C.1==1]<-NA
  L.2[C.1==1]<-NA
  Y.1[C.1==1]<-NA


  #Create end of study outcome data and set as long format
  dw<-as.data.frame(cbind(id,Y,A.1,A.2,A.3,L.1,L.2,L.3,C.1,C.2,C.3,U))
  dl<-reshape(dw,direction="long",varying=c("A.1","A.2","A.3",
                                            "L.1","L.2","L.3",
                                            "C.1","C.2","C.3"))
  #Order data by time and ID
  dl<-dl[order(dl$id,dl$time),]
  datagestcat<-dl
  datagestcat$A<-as.factor(datagestcat$A)

  dw<-as.data.frame(cbind(id,Y.1,Y.2,Y.3,A.1,A.2,A.3,L.1,L.2,L.3,C.1,C.2,C.3,U))
  dl<-reshape(dw,direction="long",varying=c("Y.1","Y.2","Y.3",
                                            "A.1","A.2","A.3",
                                            "L.1","L.2","L.3",
                                            "C.1","C.2","C.3"))
  #Order data by time and ID
  dl<-dl[order(dl$id,dl$time),]
  datagestmultcat<-dl
  datagestmultcat$A<-as.factor(datagestmultcat$A)
}else{
  #Create end of study outcome data and set as long format
  dw<-as.data.frame(cbind(id,Y,A.1,A.2,A.3,L.1,L.2,L.3,U))
  dl<-reshape(dw,direction="long",varying=c("A.1","A.2","A.3",
                                            "L.1","L.2","L.3"))
  #Order data by time and ID
  dl<-dl[order(dl$id,dl$time),]
  datagestcat<-dl
  datagestcat$A<-as.factor(datagestcat$A)

  dw<-as.data.frame(cbind(id,Y.1,Y.2,Y.3,A.1,A.2,A.3,L.1,L.2,L.3,U))
  dl<-reshape(dw,direction="long",varying=c("Y.1","Y.2","Y.3",
                                            "A.1","A.2","A.3",
                                            "L.1","L.2","L.3"))
  #Order data by time and ID
  dl<-dl[order(dl$id,dl$time),]
  datagestmultcat<-dl
  datagestmultcat$A<-as.factor(datagestmultcat$A)
}


datas<-list(datagest=datagest,datagestmult=datagestmult,datagestcat=datagestcat,
         datagestmultcat=datagestmultcat)


return(datas)

}

