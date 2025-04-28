\name{svy}
\alias{svy}
\docType{data}
\title{
Survey Data on Support for Gay Marriage (2006)
}
\description{
A dataset used for modeling support for gay marriage in the United States, 
combining individual- and state-level covariates from a 2006 survey.

\itemize{
\item \emph{Y} A data frame with 5000 rows and 11 variables: 
  \describe{
  \item{id}{Unique observation identifier}
  \item{state}{Two-letter abbreviation for U.S. state}
  \item{stateid}{Numeric identifier for the state}
  \item{region}{Region code}
  \item{age}{Age group (1 = 18-30, 2 = 31-50, 3 = 51-65, 4 = 65+)}
  \item{gXr}{Gender and race interaction}
  \item{educ}{Education level (1 = LTHS,2 = HS,3 = Some Coll,4 = Coll+)}
  \item{supp_gaymar}{Support for gay marriage (0 = oppose, 1 = support)}
  \item{pvote}{Republican presidential vote share in the previous election}
  \item{religcon}{Proportion of population identifying as religious conservatives}
  \item{libcon}{State-level ideology score (liberal to conservative)}
}

}
}

\references{
 Bisbee, James. "Barp: Improving mister p using bayesian additive regression trees." American Political Science Review 113.4 (2019): 1060-1065.
}

\keyword{datasets}
