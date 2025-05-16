#' @name group.auc.test
#' @aliases group.auc.test
#' @title Function to compare AUC across all consecutive categories of an ordinal scale
#' @description This function computes commonly used classification statistics of a confusion matrix
#' and compares the area under the curve (AUC) across all consecutive categories of an ordinal variable. 
#' The function of roc.test () from the pROC package 
#' (\url{https://cran.r-project.org/package=pROC}) is used for AUC comparison. 
#' @param outcome The outcome variable indicating the status in the form of a data frame or matrix. 
#' This variable is typically coded as 0 (positive) and 1 (negative). 
#' @param predictor  A numerical vector of scores used to predict the status of the outcome. This variable 
#' should be of the same length as the outcome variable (i.e., two variables are 
#' from the same data set and also of the same number of data rows).
#' @param groups  A data frame that contains all created indicator variables using the function
#' group.to.vars () in this package.
#' @param cut.off Specification of the criterion used to select the optimal cut score. 
#' Three options available: (1) 'max.Youden' returns the cut score that maximizes the Youden Index (the default);
#' (2) 'max.sen' returns the cut score that maximizes the sensitivity; and (3) 'max.spe' returns
#' the cut score that maximizes the specificity.
#' @param BR  Base rates or known prevalence. Multiple values can be specified simultaneously. 
#' By default BR=1.
#' @return A list of two objects: (1) descriptive and classification statistics, and
#' (2) results of the AUC comparison for each pair of the consecutive categories.
#' \item{Summary.Stats}{Summary and classification statistics for all participants and
#' all the consecutive groups. The first row is the results of the entire sample and has a row name of "All", 
#' followed by results for each pair of the groups specified by group.to.vars (). For example, 
#' if the first indicator of age is age.40, then the second row of results will have the row name of "age.40" and
#' includes results for participants with age at or below 40, the third row will have the row name of
#' "age.40.1" and includes results for those with age beyond 40. \cr
#' The results include the following statistics: \cr
#'  * N, the sample size for each category. \cr
#'  * TP, true positives. \cr
#'  * FP, false positives. \cr 
#'  * FN, false negatives.\cr
#'  * TN, true negatives.\cr
#'  * Cut.off, the optimal cut score.\cr
#'  * AUC, Area under the ROC curve.\cr
#'  * AUC.SE, Standard error of AUC.\cr
#'  * AUC.low & AUC.up, '95%' confidence interval of AUC. \cr
#'  * Sensitivity, also true positive rate, y-axis of the ROC.\cr
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
#' \item{AUC.test}{Results of the AUC comparison for each pair of the consecutive categories.}
#' @import stats
#' @importFrom pROC roc roc.test
#' @importFrom reportROC reportROC
#' @examples  
#' #read the example data
#' data(ROC.data.ex)
#' #create new binary variables for the ordinal variable
#' data.new.age<-group.to.vars(ROC.data.ex,
#'                             ROC.data.ex$age,
#'                            root.name='age')
#' #run the function
#' result.age<-group.auc.test(ROC.data.ex$outcome,ROC.data.ex$predictor, 
#'                            groups=data.new.age[,5:ncol(data.new.age)],
#'                            cut.off='max.Youden', BR=1)
#' #obtain results
#' result.age$Summary.Stats
#' result.age$AUC.test
#' @export
#' @usage group.auc.test(outcome,predictor, 
#'                       groups, cut.off='max.Youden',BR=1)

group.auc.test<-function (outcome,predictor, groups, cut.off='max.Youden', BR=1) {
  if (is.null(dim(groups))){
  groups<-as.data.frame(groups)
  names(groups)<-'Group'
  }
  ROC.sig.matrix<-matrix(NA,ncol(groups),10)
  #overall
  result.matrix<-cbind('All',length(predictor),
                       ROC.stats(outcome, predictor,cut.off='max.Youden',BR=BR))

  #subsetting data
  for (d in 1:ncol(groups)) {
    
    #ROC for first data
    outcome.c<-outcome[which(groups[,d]==0)]
    predictor.c<-predictor[which(groups[,d]==0)]
    result.matrix.sub.1<-matrix(NA,1,ncol(result.matrix))
    
    if (!(sum(outcome.c==0)>1 & sum(outcome.c==1)>1)) {
        result.matrix.sub.1[,1:2]<-cbind(colnames(groups)[d],as.matrix(length(predictor.c)))
      }  else { result.matrix.sub.1<-data.frame(colnames(groups)[d],length(predictor.c),
                                   ROC.stats(outcome.c, predictor.c,cut.off='max.Youden',BR=BR))
                          } 
  
    #ROC for second data
    outcome.r<-outcome[which(groups[,d]==1)]
    predictor.r<-predictor[which(groups[,d]==1)]
    result.matrix.sub.2<-matrix(NA,1,ncol(result.matrix))
    
    if (!(sum(outcome.r==0)>1 & sum(outcome.r==1)>1)) {
      result.matrix.sub.2[,1:2]<-cbind(colnames(groups)[d],as.matrix(length(predictor.r)))
    } else {
      result.matrix.sub.2<-data.frame(colnames(groups)[d],length(predictor.r),
                                    ROC.stats(outcome.r, predictor.r,cut.off='max.Youden',BR=BR))
      } 

    if (is.null(colnames(result.matrix.sub.1))) {
      colnames(result.matrix.sub.2)[1:2]<-c('Group','N')
      colnames(result.matrix.sub.1)<-colnames(result.matrix.sub.2)
    } else {
      colnames(result.matrix.sub.1)[1:2]<-c('Group','N')
      colnames(result.matrix.sub.2)<-colnames(result.matrix.sub.1)
    }
    colnames(result.matrix)[1:2]<-c('Group','N')
    
    result.matrix<-as.matrix(result.matrix)
    result.matrix.sub.1<-as.matrix(result.matrix.sub.1)
    result.matrix.sub.2<-as.matrix(result.matrix.sub.2)
    
    result.matrix<-rbind(result.matrix,result.matrix.sub.1,
                             result.matrix.sub.2)
    
    #compare two ROCs
    if (!(sum(outcome.c==0)>1 & sum(outcome.c==1)>1)) {
      roc1<-NULL
    } else {
      roc1 <- pROC::roc(outcome.c~predictor.c,data.frame(outcome.c,predictor.c))
    }
    
    if (!(sum(outcome.r==0)>1 & sum(outcome.r==1)>1)) {
      roc2<-NULL
    } else {
      roc2 <- pROC::roc(outcome.r~predictor.r,data.frame(outcome.r,predictor.r))
    }
    
    if (!is.null(roc1) & !is.null(roc2)){
      result.compare<-pROC::roc.test(roc1, roc2, method=c("delong"),
                               paired=NULL, alternative = c("two.sided"))
      ROC.sig.matrix[d,]<-cbind(colnames(groups)[d],
                                t(result.matrix.sub.1[,c(2,8,9)]),
                                t(result.matrix.sub.2[,c(2,8,9)]),
                                diff(result.compare$estimate),
                                result.compare$statistic,result.compare$p.value)
    } else {  ROC.sig.matrix[d,]<-cbind(colnames(groups)[d],
                                        t(result.matrix.sub.1[,c(2,8,9)]),
                                        t(result.matrix.sub.2[,c(2,8,9)]),
                                        NA,NA,NA)}
    
  } #end d
  rownames(result.matrix)<-result.matrix[,1]
  Summary.Stats<-as.data.frame(result.matrix[,-1])
  
  colnames(ROC.sig.matrix)<-c('Group','N.0',"AUC.0","AUC.SE.0",'N.1',"AUC.1","AUC.SE.1",
                              'AUC.diff','Z','p.value')
  ROC.sig.matrix[,8:10] = round(as.numeric(ROC.sig.matrix[,8:10]),3)
  AUC.test<-as.data.frame(ROC.sig.matrix)
  
  return(list(Summary.Stats=Summary.Stats,AUC.test=AUC.test))
  }
  