#' @name group.to.vars
#' @aliases group.to.vars
#' @title Function to create new variables from the ordinal variable for further analysis 
#' @description This function collapses group memberships or categories of the ordinal variable
#'  into binary variables (or indicators) for each category and appends the new variables to the end of the original data.
#'  For each new variable, 0 repsrents participants at or below the selected category and 1 reprents 
#'  participants above the selected category. For example, age.40 = 0 means participants with age at or below 40, whereas
#'  age.40 = 1 indicates participants with age beyond 40.
#' @param data A data frame or matrix that contains the ordinal variable.
#' @param group The ordinal variable in the 'data' object.
#' @param root.name Indicate whether a root name is used to name the new variables. 
#' If not specified (by default, root.name=NULL), the variable name will be used as the root.
#' @return A data frame with the original data and newly created variables.
#' @import stats  
#' @examples  
#' #read the example data
#' data(ROC.data.ex)
#' #create new binary variables for the ordinal variable
#' data.new.age<-group.to.vars(ROC.data.ex,
#'                             ROC.data.ex$age,
#'                            root.name='age')
#' @export
#' @usage group.to.vars(data, group, root.name=NULL)

group.to.vars<-function (data, group, root.name=NULL) {
  
  if (is.null(root.name)){
    dt<-subset(data,select=group)
    root.name<-names(dt)
  } 
  unique.values<- sort(unique(group))
  print(paste('A total number of ',length(unique.values),' unique categories are identified',sep=''))
  print(unique.values)

  dummies<-matrix(NA,nrow(data),length(unique.values)-1)
  for (c in 1:(length(unique.values)-1)) {
    dummies[which(group > unique.values[c]),c]<-1
    dummies[which(group <= unique.values[c]),c]<-0
  }
  dummies.names<-paste(root.name,".",unique.values[1:(length(unique.values)-1)],sep='')
  colnames(dummies)<-dummies.names
  data.new<-data.frame(data,dummies)
  
  return(data.new)
}
