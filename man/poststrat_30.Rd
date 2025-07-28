\name{poststrat_30}
\alias{poststrat_30}
\docType{data}
\title{
Post-Stratification Table of 2014-2018 American Community Survey (ACS)
}
\description{
A dataset used for poststratification modeling based on the 2014–2018 American Community Survey (ACS), 
including 12,000 combinations of demographic and geographic strata. The table includes more combinations than the number of observed units, so some strata are not represented in the CCES sample. This dataset includes only the 30 states that are present in the \code{cces_30_df}.


A data frame with 7200 rows and 6 variables:
\describe{
  \item{state}{A factor variable representing the 30 selected U.S. states. The same states as those included in cces_30_df.}
  \item{eth}{An ethnicity variable stored as a factor (White, Black, Hispanic, Other).}
  \item{gender}{A factor variable representing gender (Female, Male).}
  \item{age}{A factor variable dividing respondents into age groups (18–29, 30–39, 40–49, 50–59, 60–69, 70+).}
  \item{educ}{A factor variable representing educational attainment (No HS, HS, Some college, 4-Year College, Post-grad).}
  \item{n}{Number of individuals in each demographic-geographic combination}
}
}

\references{
 Juan Lopez-Martin, Justin H. Phillips, and Andrew Gelman, Multilevel Regression and Poststratification Case Studies
 \url{https://bookdown.org/jl5522/MRP-case-studies/introduction-to-mrp.html#data}
}

\keyword{datasets}
