#' Approximate Shapley values computed from a BART model fitted using \code{bart}
#'
#' `Explain.bart` function is used to calculate the contribution of each variable
#' in the Bayesian Additive Regression Trees (BART) model using permutation.
#' It is used to compute the Shapley values of models estimated using the \code{bart} function from the \code{dbarts}.
#'
#' @param object A BART model (Bayesian Additive Regression Tree) estimated
#' using the \code{bart} function from the \pkg{dbarts}.
#' @param feature_names The name of the variable for which you want to check the contribution.
#' The default value is set to \code{NULL}, which means the contribution of all variables in \code{X} will be calculated.
#' @param X The dataset containing all independent variables used as input when estimating the BART model.
#' @param nsim The number of Monte Carlo sampling iterations, which is fixed at \code{1} by default
#' in the case of the BART model.
#' @param pred_wrapper A function used to estimate the predicted values of the model.
#' @param newdata New data containing the variables included in the model.
#' This is used when checking the contribution of newly input data using the model.
#' The default value is set to \code{NULL}, meaning that the input \code{X} data,
#' i.e., the data used for model estimation, will be used by default.
#' @param parallel The default value is set to \code{FALSE}, but it can be changed to \code{TRUE} for parallel computation.
#' @param ... Additional arguments to be passed
#' @return Returns of class \code{ExplainBART} with consisting of a list with the following components:
#' \item{phis}{A list containing the Shapley values for each variable.}
#' \item{newdata}{The data used to check the contribution of variables. If a variable has categories, categorical variables are one-hot encoded.}
#' \item{fnull}{The expected value of the model's predictions.}
#' \item{fx}{The prediction value for each observation.}
#' \item{factor_names}{The name of the categorical variable. If the data contains only continuous or dummy variables, it is set to \code{NULL}.}
#' @export
#' @examples
#' \donttest{
#' ## Friedman data
#' set.seed(2025)
#' n = 200
#' p = 5
#' X = data.frame(matrix(runif(n * p), ncol = p))
#' y = 10 * sin(pi* X[ ,1] * X[,2]) +20 * (X[,3] -.5)^2 + 10 * X[ ,4] + 5 * X[,5] + rnorm(n)
#' 
#' ## Using the dbarts library
#' model = dbarts::bart(X,y,keeptrees = TRUE ,  ndpost = 200)
#' 
#' ## prediction wrapper function
#' pfun <- function(object, newdata) {
#'        predict(object , newdata)
#'        }
#'        
#'## Calculate shapley values
#'model_exp =  Explain  ( model, X = X,  pred_wrapper =  pfun )
#'}
Explain.bart <- function(object, feature_names = NULL, X = NULL, nsim = 1,  pred_wrapper = NULL,
                         newdata = NULL,   parallel = FALSE, ...) {

  i<-0;
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

  if (is.null(pred_wrapper)) {
    stop("Prediction function required for approximate Shapley values. Please ",
         "supply one via the `pred_wrapper` argument", call. = FALSE)
  }

  if (is.null(newdata)) {
       newdata <-  X
  }

  fx  <-   colMeans ( predict(object, newdata =  newdata ))
  fnull <-  mean( colMeans (object$yhat.train)) # baseline value (i.e., avg training prediction)

  # Deal with other NULL arguments
  if (is.null(feature_names)) {
    feature_names = colnames(X)
  }

  # Set up the 'foreach' "do" operator
  `%.do%` <- if (isTRUE(parallel)) `%dopar%` else `%do%`

  # Compute approximate Shapley values # ,  ...
  phis_temp <- foreach(i = feature_names ) %.do% {
     t(Explain_column(object, X = X, column = i, pred_wrapper = pred_wrapper,
                                 newdata = newdata))
    # number of post by obs = n matrix , list = number of variable
  }
  names( phis_temp ) <- feature_names


  featurenames <- names(as.data.frame(makeModelMatrixFromDataFrame (newdata)))
  
  count_matches <- sapply(feature_names, function(fn) {
    sum(startsWith(featurenames, fn))
  })
  
  
  tmp_var <- data.frame(
    var = names(count_matches),
    n = as.integer(count_matches),
    row.names = NULL
  )


  if (( sum(sapply(newdata, is.factor)) > 0 | sum(sapply(newdata, is.character)) > 0) &
      sum( tmp_var $ n  >=2) >= 1 ) {

    factor_names <- names(X) [which( ( names(X) %in% featurenames  ) ==FALSE)]

    phis  <-  foreach(i = 1:length(featurenames) ) %.do% {
      
      ind <- featurenames[i]
      phi_idx <- which( stringr::str_detect( ind , feature_names ))
      
      temp_var <- NULL
      value_factor <- NULL
      
      if( tmp_var$n[tmp_var$ var ==   feature_names[ phi_idx]] >=  2 ){
        value_factor <-  tail(unlist(strsplit( ind , split = "\\.")),1)
        temp_var <- stringr::str_replace (  ind , paste0("\\.",value_factor),"")
      }
      
        temp <- matrix (0, nrow = dim (newdata)[1],ncol= dim(object $ yhat.train) [1] )
       
       if( is.null(temp_var) ==FALSE & length(which(stringr::str_detect( featurenames, paste0("^",temp_var)  ))) >  2 ) {
 
         temp [which(newdata[,phi_idx ]==value_factor),] <- phis_temp[[phi_idx]] [which(newdata[,phi_idx]==  value_factor),]
         temp
       }  else  if(is.null(temp_var) ==FALSE & length(which(stringr::str_detect( featurenames, paste0("^",temp_var)  ))) == 2) {
         
         temp [which(newdata[,phi_idx ]==1),]  <-
           phis_temp[[phi_idx]] [which(newdata[,phi_idx]==  1),]
         temp
         
       } else {
         temp <- phis_temp[[ phi_idx]]
         temp
       }

    }
    names( phis) <-  featurenames
    newdata <- as.data.frame(makeModelMatrixFromDataFrame (newdata))

  } else if ( sum(sapply(newdata, is.factor)) == 0 ){
    phis <- phis_temp

    if( sum(featurenames%in%feature_names) != length(featurenames) ){
      factor_names <- names(X) [which( ( names(X) %in% featurenames ) ==FALSE)]
    }else{
      factor_names <- NULL
    }

  }

  out <- list (phis = phis, newdata = newdata,   fnull =  fnull,  fx = fx, factor_names = factor_names)

  class(out) <- "ExplainBART"
  return (out)

}
