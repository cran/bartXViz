#' Approximate Shapley values computed from a BART model fitted using \code{bartMachine}
#'
#' This function is used to calculate the contribution of each variable
#' in the Bayesian Additive Regression Trees (BART) model using permutation.
#' It is used to compute the Shapley values of models estimated 
#' using the \code{bartMachine} function from the \pkg{bartMachine}.
#'
#' @param object A BART model (Bayesian Additive Regression Tree) estimated
#' using the \code{bartMachine} function from the \pkg{bartMachine}.
#' @param feature_names The name of the variable for which you want to check the contribution.
#' The default value is set to \code{NULL}, which means the contribution of all variables in \code{X} will be calculated.
#' @param X The dataset containing all independent variables used as input when estimating the BART model. Categorical or character variables must not contain an underscore ("_") in their values or labels.
#' @param nsim The number of Monte Carlo repetitions used for estimating each Shapley value is set to \code{1} by default for the BART model.
#' @param pred_wrapper A function used to estimate the predicted values of the model.
#' @param newdata New data containing the variables included in the model.
#' This is used when checking the contribution of newly input data using the model.
#' The default value is set to \code{NULL}, meaning that the input \code{X} data,
#' i.e., the data used for model estimation, will be used by default.
#' @param parallel The default value is set to \code{FALSE},
#' but it can be changed to \code{TRUE} for parallel computation.
#' @param ... Additional arguments to be passed
#' @return An object of class \code{ExplainbartMachine} with the following components :
#' \item{phis}{A list containing the Shapley values for each variable.}
#' \item{newdata}{The data used to check the contribution of variables. If a variable has categories, categorical variables are one-hot encoded.}
#' \item{fnull}{The expected value of the model's predictions.}
#' \item{fx}{The prediction value for each observation.}
#' \item{factor_names}{The name of the categorical variable. If the data contains only continuous or dummy variables, it is set to \code{NULL}.}
#' @export
#' @examples 
#' \donttest{
#'## Friedman data
#'set.seed(2025)
#'n = 200
#'p = 5
#'X = data.frame(matrix(runif(n * p), ncol = p))
#'y = 10 * sin(pi* X[ ,1] * X[,2]) +20 * (X[,3] -.5)^2 + 10 * X[ ,4] + 5 * X[,5] + rnorm(n)
#'
#'##  Using the bartMachine library 
#'model = bartMachine::bartMachine(X,y, seed = 2025, num_iterations_after_burn_in =200 )
#'
#'## prediction wrapper function
#'pfun <- function (object, newdata) {
#'   bartMachine::bart_machine_get_posterior(object,newdata) $ y_hat_posterior_samples
#'   }
#'   
#'## Calculate shapley values
#'model_exp =  Explain  ( model, X = X,  pred_wrapper =  pfun )
#'}

Explain.bartMachine <- function(object, feature_names = NULL,  X = NULL,
                                nsim = 1, pred_wrapper = NULL,
                                newdata = NULL,   parallel = FALSE, ...) {
  
  i<- 0; n <- 0;
  # Only nsim = 1
  if (nsim > 1) stop ("It stops because nsim > 1.",
                      "Because the BART model uses posterior samples,",
                      "it is used by setting nsim=1.", call. = FALSE)


  # Compute baseline/average training prediction (fnull) and predictions
  # associated with each explanation (fx); if `adjust = FALSE`, then the
  # baseline is not needed and defaults to zero.
  if (is.null(X)) {
    stop("Training features required for approximate Shapley values.", call. = FALSE)
  }

  if (inherits(X, what = "tbl_df")) {
    X <- as.data.frame(X)
  }


  if (is.null(newdata)) {  # Explain all rows of background data set
    newdata <- X
  }

  if (is.null(pred_wrapper)) {
    stop("Prediction function required for approximate Shapley values. Please ",
         "supply one via the `pred_wrapper` argument", call. = FALSE)
  }


  fx <-   predict(object, new_data = newdata)
  # baseline value (i.e., avg training prediction)
  fnull <- mean(as.numeric(object $ y_hat_train))

  # Deal with other NULL arguments
  if (is.null(feature_names)) {
    feature_names <- colnames(X)
  }

  # Set up the 'foreach' "do" operator
  `%.do%` <- if (isTRUE(parallel)) `%dopar%` else `%do%`
  # Compute approximate Shapley values # ,  ...
  phis_temp <- foreach(i = feature_names ) %.do% {
    Explain_column(object, X = X, column = i, pred_wrapper = pred_wrapper,
                              newdata = newdata)
    # number of post by obs = n matrix , list = number of variable
  }
  names(phis_temp) <- feature_names

  temp_new <-  as.data.frame( pre_process_new_data(newdata,object))
  featurenames <- names(  temp_new )


  count_matches <- sapply(feature_names, function(fn) {
    sum(startsWith(featurenames, fn))
  })
  
 
  tmp_var <- data.frame(
    var = names(count_matches),
    n = as.integer(count_matches),
    row.names = NULL
  )
  
  if ( (sum(sapply(newdata, is.factor)) > 0  | sum(sapply(newdata, is.character)) > 0) & (
    sum( tmp_var $ n  >=2) >= 1)) {
 
    phis  <-  foreach(i = 1:length(  featurenames ) ) %.do% {

      ind <- names(temp_new)[i]
      phi_idx <- which( stringr::str_detect( ind , feature_names ))

      temp_var <- NULL
      value_factor <- NULL
      
      if( tmp_var$n[tmp_var$ var ==   feature_names[ phi_idx]] >=  2 ){
        value_factor <-  tail(unlist(strsplit( ind , split = "_")),1)
        temp_var <- stringr::str_replace (  ind , paste0("_",value_factor),"")
      }


      temp <- matrix (0, nrow = dim (newdata)[1],
                     ncol= object$num_iterations_after_burn_in )

      if(is.null(temp_var) ==FALSE & length(which(stringr::str_detect( featurenames, paste0("^",temp_var)  ))) >  2) {

        temp [which(newdata[,phi_idx ]==value_factor),] <-
          phis_temp[[phi_idx]] [which(newdata[,phi_idx] == value_factor),]
        temp

      } else  if(is.null(temp_var) ==FALSE & length(which(stringr::str_detect( featurenames, paste0("^",temp_var)  ))) == 2) {

        temp [which(newdata[,phi_idx ]==1),]  <-
          phis_temp[[phi_idx]] [which(newdata[,phi_idx]==  1),]
        temp

      } else {
        temp  <- phis_temp[[phi_idx]]
        temp
      }
    }
    names( phis) <-  featurenames

    factor_names <- NULL
    
    if ( sum( tmp_var$n >=2) >= 1  ){
      
      names( phis) <- names(temp_new)
      factor_names <-  names(temp_new ) [which((names(temp_new )%in% names(X) ) ==FALSE)]

    } 
    }else if ( sum(sapply(newdata, is.factor)) == 0 ){
    phis <- phis_temp
    factor_names <- NULL
  }

  out <- list (phis = phis, newdata = temp_new,   fnull =  fnull,  fx = fx, factor_names = factor_names)

  class(out) <- "ExplainbartMachine"
  return (out)


}
