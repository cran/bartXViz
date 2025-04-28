#' Approximate Shapley values computed from the BARP model
#'
#' This function is implemented to calculate the contribution of each variable
#' in the BARP (Bayesian Additive Regression Tree with post-stratification) model
#' using the permutation method.
#'
#' @param object A BARP model (Bayesian Additive Regression Tree) estimated
#' using the `barps` function, a modified version of the `barp` function from the BARP library with a fixed seed.
#' @param feature_names The name of the variable for which you want to check the contribution.
#' The default value is set to `NULL`, which means the contribution of all variables in `X` will be calculated.
#' @param X The dataset containing all independent variables used as input when estimating the BART model.T he explanatory variables X included in the model must be converted to factors prior to input.
#' @param nsim The number of Monte Carlo sampling iterations, which is fixed at `1` by default in the case of the BARP model.
#' @param pred_wrapper A function used to estimate the predicted values of the model.
#' @param census  Census data containing the names of the `X` columns. It should also have the same format as `X` and include a variable named `proportion`, which indicates the number of individuals corresponding to each combination.
#' @param geo.unit  Enter the name of the stratification variable used in post stratification.
#' @param parallel The default value is set to `FALSE`, but it can be changed to `TRUE` for parallel computation.
#' @param ... Additional arguments to be passed
#' @return Returns of class `Explainbarp` with consisting of a list with the following components:
#' \item{phis}{A list containing the Shapley values for each variable.}
#' \item{newdata}{The data used to check the contribution of variables. If a variable has two categories, it is dummy-coded, and if it has three or more categories, categorical variables are one-hot encoded.}
#' \item{fnull}{The expected value of the model's predictions.}
#' \item{fx}{The prediction value for each observation.}
#' \item{factor_names}{The name of the categorical variable. If the data contains only continuous or dummy variables, it is set to NULL.}
#' @export

Explain.barp <- function(object, feature_names = NULL, X = NULL, nsim = 1, pred_wrapper = NULL,
                         census = NULL, geo.unit = NULL,   parallel = FALSE,   ...) {
  
  i <- 0;
  if (is.null(census)) {
    stop("Please enter your census data.", call. = FALSE)
  }
  
  # Only nsim = 1
  if (nsim > 1) stop ("It stops because nsim > 1.", call. = FALSE)
  
  # prediction function check
  if (is.null(pred_wrapper)) {
    stop("Prediction function required for approximate Shapley values.", call. = FALSE)
  }
  
  if (is.null(X)) {
    X <- object $ trees $ X
  }
  # Convert to data frame
  if (inherits(X, what = "tbl_df")) {
    X <- as.data.frame(X)
  }

    # expacted value
    if( sum(names (object$barp.dat) == "SL.bartMachine_1_All")==1  ){
      fx <-  object$barp.dat$SL.bartMachine_1_All
    } else if(  sum(names (object$barp.dat)== "SL.bartMachine_All")==1 ){
      fx <-  object$barp.dat$SL.bartMachine_All
    }
  

    # predicted value
    fnull <-  object $ pred.opn [,1:2]
    names (  fnull ) [2] <- "f_null"

    # variable names
    feature_names <- names(object $ trees $ X)

    # Variable sorting
    censusdata <- as.data.frame(census[,  feature_names])
    X <- X [,  feature_names]

    # Categorical variable processing
    newdata <- as.data.frame( pre_process_new_data(censusdata,object $ trees))

  # Set up the 'foreach' "do" operator
  `%.do%` <- if (isTRUE(parallel)) `%dopar%` else `%do%`

  set.seed(object $ setSeed)
  # Compute approximate Shapley values
  phis_temp <-  foreach(i = feature_names  ) %.do% {
    Explain_column ( object$trees, X = X, column = i,
                              pred_wrapper = pred_wrapper,newdata = censusdata)
  }
  names( phis_temp) <- feature_names


  # Deal with other NULL arguments
  geo_check <- stringr::str_detect(colnames(newdata),paste0("^", geo.unit))
  featurenames <- colnames(newdata) [-which(geo_check)]


  phis  <-  foreach(i = 1:length(featurenames) ) %.do% {
    ind <- featurenames[i]
    phi_idx <- which( str_detect( ind , feature_names ))

    temp <- matrix (0, nrow = dim (newdata)[1],ncol= object$trees$num_iterations_after_burn_in )
    temp [which(newdata[,ind]==1),] <- phis_temp[[phi_idx]] [which(newdata[,ind]==1),]
    temp
  }
  names( phis) <-  featurenames
  
 
  split_last_underscore <- function(x) {
    idx <- regexpr("_[^_]*$", x) 
    if (idx > 0) {
      left <- substr(x, 1, idx - 1)
      right <- substr(x, idx + 1, nchar(x))
      return(c(left, right))
    } else {
      return(c(x, NA)) 
    }
  }
  
 
  split_names <- lapply(featurenames, split_last_underscore)
  name_temp <- as.data.frame(do.call(rbind, split_names), stringsAsFactors = FALSE)
 
  check <- name_temp[(name_temp[,1] == name_temp [,2]),1] 
  
  if (length(check) >=1) {
    factor_names <- featurenames [-which(featurenames %in% check)]
  } else if (  identical(check, character(0)) ){
    factor_names <- featurenames
  }
 
    out <- list (phis = phis, newdata = newdata,   fnull =  fnull, fx = fx ,factor_names = factor_names )

  class( out ) <- "Explainbarp"

  return( out )
}


