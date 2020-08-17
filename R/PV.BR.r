#' @name PV.BR
#' @aliases PV.BR
#' @title Function to compute PPV and NPV with specified base rates
#' @description This function computes positive predictive values (PPV) and negative predictive values (NPV)
#'  with provided base rates (or known prevalence).
#' @param outcome The outcome variable indicating the status in the form of a data frame or matrix. 
#' This variable is typically coded as 0 (positive) and 1 (negative). 
#' @param predictor  A numerical vector of scores used to predict the status of the outcome. This variable 
#' should be of the same length as the outcome variable (i.e., two variables are 
#' from the same data set and also of the same number of data rows).
#' @param cut.off Specification of the criterion used to select the optimal cut score. 
#' Three options available: (1) 'max.Youden' returns the cut score that maximizes the Youden Index (the default);
#' (2) 'max.sen' returns the cut score that maximizes the sensitivity; and (3) 'max.spe' returns
#' the cut score that maximizes the specificity.
#' @param BR  Base rates or known prevalence. Multiple values can be specified simultaneously. 
#' By default BR=1.
#' @return An object that contains results of classification statistics.
#' \item{Result}{
#'  * Cut.off, the optimal cut score.\cr
#'  * Sensitivity, also true positive rate, the y-axis of the ROC.\cr
#'  * Specificity, also true negative rate. \cr
#'  * Youden.Index. \cr
#'  * PPV or positive predictive value for each specified base rate.\cr
#'  * NPV or negative predictive value for each specified base rate.\cr
#'  * PPV for the sample.\cr
#'  * NPV for the sample.}
#' @examples  
#' #read the example data
#' data(ROC.data.ex)
#' #run the function
#' PV.BR(ROC.data.ex$outcome, ROC.data.ex$predictor,
#'       cut.off='max.Youden', BR=1)
#' @export
#' @usage  PV.BR(outcome, predictor,cut.off='max.Youden', BR=1)
#' @references {
#'  McCaffrey R.J., Palav A.A., O’Bryant S.E., Labarge A.S. (2003). 
#' "A Brief Overview of Base Rates. 
#' In: McCaffrey R.J., Palav A.A., O’Bryant S.E., Labarge A.S. (eds) 
#' Practitioner’s Guide to Symptom Base Rates in Clinical Neuropsychology. Critical Issues in Neuropsychology. ."
#' Springer, Boston, MA. doi:10.1007/978-1-4615-0079-7_1. 
#' }

PV.BR<-function (outcome,predictor,cut.off='max.Youden', BR=1) {
  
  PPV.NPV.BR <- function(SEN,SPE,BR=1) {
    PPV <- (SEN*BR)/((SEN*BR)+((1-SPE)*(1-BR)))
    NPV <-(SPE*(1-BR))/((SPE*(1-BR))+((1-SEN)*BR))
    result<-c(PPV,NPV)
    return(result)
  }
  
  temp<-cutscores(outcome,predictor)
  if (cut.off=="max.sen") {
  cutoff<-temp$Summary[1,1] 
  } else if (cut.off=="max.spe") {
      cutoff<-temp$Summary[2,1] 
      } else  {
        cutoff<-temp$Summary[3,1] 
        }
  
  Result<-matrix(NA,1,(4+length(BR)*2)) 

  test<-rep(0,length(outcome))
  test[which(predictor>=cutoff)]<-1
  TP<-length(outcome[which(test==1 & outcome==1)])
  FP<-length(outcome[which(test==1 & outcome==0)])
  FN<-length(outcome[which(test==0 & outcome==1)])
  TN<-length(outcome[which(test==0 & outcome==0)])
  SEN<-TP/(TP+FN)
  SPE<-TN/(FP+TN)
  Youden<-SEN+SPE-1
  ppv<-matrix(NA,1,length(BR)) 
  npv<-matrix(NA,1,length(BR))
  for (b in 1:length(BR)) {
    ppv[b]<-PPV.NPV.BR (SEN,SPE,BR[b])[1]
    npv[b]<-PPV.NPV.BR (SEN,SPE,BR[b])[2]
  }
  Result[1,]<-round(c(cutoff,SEN,SPE,Youden,ppv,npv),3)
  ppv.BR.name<-paste('PPV(BR=',BR,')',sep='')
  npv.BR.name<-paste('NPV(BR=',BR,')',sep='')
  colnames(Result)<-c('Cut.off','Sensitivity','Specificity', 'Youden.Index',
                      ppv.BR.name,npv.BR.name)
  return(as.data.frame(Result))
  }
  