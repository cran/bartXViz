\name{census06}
\alias{census06}
\docType{data}
\title{
Census-based Population Proportions for Covariate Bins (2006)
}
\description{
The data frame has the following components:

\itemize{
\item \emph{Y} This dataset provides population counts in covariate bins based on the 2006 U.S. Census,
Each row represents a unique combination of demographic covariates within a state. 
A data frame with 2940 rows and 9 variables:
  \describe{
  \item{stateid}{Numeric identifier for the state}
  \item{region}{Region code}
  \item{age}{Age group (1 = 18-30, 2=31-50, 3= 51-65, 4 =65+)}
  \item{gXr}{Gender and race interaction}
  \item{educ}{Education level (1 = LTHS,2 = HS,3 = Some Coll,4 = Coll+)}
  \item{pvote}{Republican presidential vote share in the previous election}
  \item{religcon}{Proportion of population identifying as religious conservatives}
  \item{libcon}{State-level ideology score (liberal to conservative)}
  \item{n}{Population count for the given covariate bin within the state}
}

}
}

\references{
 Bisbee, James. "Barp: Improving mister p using bayesian additive regression trees." American Political Science Review 113.4 (2019): 1060-1065.
}

\keyword{datasets}
