\name{cces_30_df}
\alias{cces_30_df}
\docType{data}
\title{
Survey Data on Public Opinion about Abortion Coverage in Insurance Plans(2018)
}
\description{
A dataset used for modeling support for abortion coverage in insurance plans in the United States, 
combining individual- and district-level covariates from a 2018 survey.
Out of all U.S. states, only data from 30 states are included in this dataset.

A data frame with 49,095 rows and 6 variables:
\describe{
  \item{abortion}{Allow employers to decline coverage of abortions in insurance plans (0: Oppose, 1: Support).}
  \item{state}{A factor variable representing the 30 selected U.S. states. The raw data includes states with state_fips values ranging from 1 to 30.}
  \item{eth}{An ethnicity variable stored as a factor (White, Black, Hispanic, Other).}
  \item{gender}{A factor variable representing gender (Female, Male).}
  \item{age}{A factor variable dividing respondents into age groups (18–29, 30–39, 40–49, 50–59, 60–69, 70+).}
  \item{educ}{A factor variable representing educational attainment (No HS, HS, Some college, 4-Year College, Post-grad).}
}
}

\references{
 Juan Lopez-Martin, Justin H. Phillips, and Andrew Gelman, Multilevel Regression and Poststratification Case Studies
 \url{https://bookdown.org/jl5522/MRP-case-studies/introduction-to-mrp.html#data}
}

\keyword{datasets}
