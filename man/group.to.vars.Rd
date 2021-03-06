% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/group.to.vars.r
\name{group.to.vars}
\alias{group.to.vars}
\title{Function to create new variables from the ordinal variable for further analysis}
\usage{
group.to.vars(data, group, root.name=NULL)
}
\arguments{
\item{data}{A data frame or matrix that contains the ordinal variable.}

\item{group}{The ordinal variable in the 'data' object.}

\item{root.name}{Indicate whether a root name is used to name the new variables. 
If not specified (by default, root.name=NULL), the variable name will be used as the root.}
}
\value{
A data frame with the original data and newly created variables.
}
\description{
This function collapses group memberships or categories of the ordinal variable
 into binary variables (or indicators) for each category and appends the new variables to the end of the original data.
 For each new variable, 0 repsrents participants at or below the selected category and 1 reprents 
 participants above the selected category. For example, age.40 = 0 means participants with age at or below 40, whereas
 age.40 = 1 indicates participants with age beyond 40.
}
\examples{
 
#read the example data
data(ROC.data.ex)
#create new binary variables for the ordinal variable
data.new.age<-group.to.vars(ROC.data.ex,
                            ROC.data.ex$age,
                           root.name='age')
}
