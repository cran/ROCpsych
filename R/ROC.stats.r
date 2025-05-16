#' @name ROC.stats
#' @aliases ROC.stats
#' @title Function to compute statistics from a confusion matrix
#' @description This function computes all diagnostic statistics from a confusion matrix.
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
#' @return An object that contains the results.
#' \item{ROC.stats}{Summary and classification statistics for all participants and
#' all the consecutive groups. Specifically. \cr
#'  * N, sample size for each category. \cr
#'  * TP, true positives. \cr
#'  * FP, false positives. \cr 
#'  * FN, false negatives.\cr
#'  * TN, true negatives.\cr
#'  * Cut.off, the optimal cut score.\cr
#'  * AUC, Area under the ROC curve.\cr
#'  * AUC.SE, Standard error of AUC.\cr
#'  * AUC.low & AUC.up, '95%' confidence interval of AUC.\cr
#'  * Sensitivity, also true positive rate, the y-axis of the ROC.\cr
#'  * Specificity, also true negative rate. \cr
#'  * Youden.Index. \cr
#'  * PPV or positive predictive value for each specified base rate.\cr
#'  * NPV or negative predictive value for each specified base rate.\cr
#'  * PPV for the sample.\cr
#'  * NPV for the sample.\cr
#'  * FNR, false negative rate, or miss rate.\cr
#'  * FPR, false positive rate, or fall-out rate.\cr
#'  * FOR, false omission rate.\cr
#'  * FDR, false discovery rate.\cr
#'  * Prevalence.\cr
#'  * Accuracy.\cr
#'  * PLR, positive likelihood ratio.\cr
#'  * NLR, negative likelihood ratio.\cr
#'  * DOR, Diagnostic odds ratio.}
#' @import stats
#' @importFrom reportROC reportROC
#' @examples  
#' #read the example data
#' data(ROC.data.ex)
#' #run the function
#' ROC.stats(ROC.data.ex$outcome, ROC.data.ex$predictor,
#'           cut.off='max.Youden',BR=1)
#' @export
#' @usage ROC.stats(outcome, predictor,cut.off='max.Youden',BR=1)

ROC.stats<-function (outcome, predictor,cut.off='max.Youden',BR=1) {
  
  temp.1<-reportROC::reportROC(gold=outcome,predictor=predictor,important="se",plot=F)
  temp.1.to.use<-as.matrix(temp.1[1,c(2:5)])
  temp.2<-PV.BR(outcome, predictor,cut.off,BR)
  cutoff<-temp.2$Cut.off
  test<-rep(0,length(outcome))
  test[which(predictor>=cutoff)]<-1
  TP<-length(outcome[which(test==1 & outcome==1)])
  FP<-length(outcome[which(test==1 & outcome==0)])
  FN<-length(outcome[which(test==0 & outcome==1)])
  TN<-length(outcome[which(test==0 & outcome==0)])
  SEN<-TP/(TP+FN)
  SPE<-TN/(FP+TN)
  Youden<-SEN+SPE-1
  SE.SP<-cbind(TP/(TP+FN),TN/(FP+TN))
  P.N.T<-cbind(TP/(TP+FP),TN/(FN+TN))
  FNR<-FN/(TP+FN)
  FPR<-FP/(FP+TN)
  FOR<-FN/(FN+TN)
  FDR<-FP/(TP+FP)
  Prevalance<-(TP+FN)/(TP+FP+FN+TN)
  Accuracy<-(TP+TN)/(TP+FP+FN+TN)
  PIC<-(FN+FP)/(TP+FP+FN+TN)
  PLR<-SEN/FPR
  NLR<-FNR/SPE
  DOR<-PLR/NLR
  
result.matrix<-cbind(TP,FP,FN,TN,round(cutoff,3),
                     temp.1.to.use,temp.2[,-1],
                      round(P.N.T,3),round(FNR,3),round(FPR,3),round(FOR,3),round(FDR,3),
                      round(Prevalance,3),round(Accuracy,3),round(PLR,3),round(NLR,3),round(DOR,3))
names(result.matrix)<-c('TP','FP','FN','TN',"Cut.off",
                           "AUC","AUC.SE","AUC.low", "AUC.up",
                           names(temp.2[,-1]),
                           'PPV.sample','NPV.sample','FNR/MissRate','FPR/FallOut',
                           'FOR','FDR','Prevalance','Accuracy','PLR','NLR','DOR')
ROC.stats<-as.data.frame(result.matrix)
  return(ROC.stats)
}
