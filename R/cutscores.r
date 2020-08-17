#' @name cutscores
#' @aliases cutscores
#' @title Function to compute optimal cut-off scores
#' @description This function computes the optimal cut-off scores based on sensitivity, specificity, and the
#' Youden Index (Youden, 1950) <doi:10.1002/1097-0142(1950)3:1<32::AID-CNCR2820030106>3.0.CO;2-3>.
#' @param outcome The outcome variable indicating the status in the form of a data frame or matrix. 
#' This variable is typically coded as 0 (positive) and 1 (negative). 
#' @param predictor  A numerical vector of scores used to predict the status of the outcome. This variable 
#' should be of the same length as the outcome variable (i.e., two variables are 
#' from the same data set and also of the same number of data rows).
#' @return A list of two objects: (1) summary statistics of selected cut scores, and (2) detailed
#' information of each used cut score and corresponding classification statistics.
#' \item{Summary}{Summary statistics of selected cut scores. Specifically, \cr
#'  * Cut.off, the select cut-off scores according to different criteria\cr
#'  * SEN, Sensitivity, also true positive rate, the y-axis of the ROC.\cr
#'  * SPE, Specificity, also true negative rate. \cr
#'  * 1-SPE, the x-axis of the ROC. \cr
#'  * Youden.Index. \cr
#'  * TP, true positives. \cr
#'  * FP, false positives. \cr 
#'  * FN, false negatives.\cr
#'  * TN, true negatives.}
#' \item{Details}{Detailed information of each used cut score and corresponding classification statistics.}
#' @import stats
#' @examples  
#' #read the example data
#' data(ROC.data.ex)
#' #run the function
#' result<-cutscores(ROC.data.ex$outcome, ROC.data.ex$predictor)
#' #obtain results
#' result$Summary
#' result$Details 
#' @export
#' @references {
#'  Youden, W.J. (1950). 
#' "Index for rating diagnostic tests."
#'  Cancer,3, 32-35. doi:10.1002/1097-0142(1950)3:1<32::AID-CNCR2820030106>3.0.CO;2-3. 
#' }
#' @usage cutscores(outcome, predictor)

cutscores<-function (outcome, predictor) {
  
  cut.off.init<-sort(unique(predictor))
  cut.off.low<-cut.off.init[1]-1
  cut.off.high<-cut.off.init[length(cut.off.init)]+1
  cut.off.others<-(cut.off.init[1:(length(cut.off.init)-1)]+cut.off.init[2:length(cut.off.init)])/2
  cut.off.temp<-c(cut.off.low,cut.off.others,cut.off.high)
  
  
  Result<-matrix(NA,length(cut.off.temp),9)   
  for (i in 1:length(cut.off.temp)) {
  cut.off.value<-cut.off.temp[i]
  test<-rep(0,length(outcome))
  test[which(predictor>=cut.off.value)]<-1
  
  TP<-length(outcome[which(test==1 & outcome==1)])
  FP<-length(outcome[which(test==1 & outcome==0)])
  FN<-length(outcome[which(test==0 & outcome==1)])
  TN<-length(outcome[which(test==0 & outcome==0)])
  SEN<-TP/(TP+FN)
  SPE<-TN/(FP+TN)
  SPE.rv<-1-SPE
  Youden<-SEN+SPE-1
  Result[i,]<-c(cut.off.value,SEN,SPE,SPE.rv,Youden,TP,FP,FN,TN)
  }
  colnames(Result)<-c('Cut.off','SEN','SPE', '1-SPE', 'Youden.Index',
                      'TP','FP','FN','TN')
    Result<-as.data.frame(Result)
    
    cut.max.sen<-max(Result$SEN[which(Result$SEN<1)])
    if (sum(Result$SEN==cut.max.sen)==1) {
      cut.max.sen.row<-Result[which(Result$SEN==cut.max.sen),]
    } else {
      result.sub.sen<-subset(Result, SEN==cut.max.sen )
      cut.max.yd<-max(result.sub.sen$Youden.Index)
      cut.max.sen.row<-result.sub.sen[which(result.sub.sen$Youden.Index==cut.max.yd),]
    }
    
    cut.max.spe<-max(Result$SPE[which(Result$SPE<1)])
    if (sum(Result$SPE==cut.max.spe)==1) {
      cut.max.spe.row<-Result[which(Result$SPE==cut.max.spe),]
    } else {
      result.sub.spe<-subset(Result, SPE==cut.max.spe )
      cut.max.yd<-max(result.sub.spe$Youden.Index)
      cut.max.spe.row<-result.sub.spe[which(result.sub.spe$Youden.Index==cut.max.yd),]
    }
    
    cut.max.yd<-max(abs(Result$Youden.Index[which(Result$Youden.Index<1)]))
    cut.max.yd.row<-Result[which(abs(Result$Youden.Index)==cut.max.yd),]
    if (nrow(cut.max.yd.row)>1) {
      cut.max.yd.row<-cut.max.yd.row[which(cut.max.yd.row$SEN==max(cut.max.yd.row$SEN)),]
      }
    
    cut.off.summary<-rbind(cut.max.sen.row,cut.max.spe.row,cut.max.yd.row)
    row.names(cut.off.summary)<-c('max.sen','max.spe','max.Youden')
    return (list(Summary=cut.off.summary,
                 Details=Result))
}