#' Approximate Shapley values
#' 
#' Compute fast (approximate) Shapley values for a set of features using the 
#' Monte Carlo algorithm described in Strumbelj and Igor (2014).
#' An efficient algorithm for tree-based models, commonly referred to as Tree SHAP, is also 
#' supported for [lightgbm](https://cran.r-project.org/package=lightgbm) and
#' [xgboost](https://cran.r-project.org/package=xgboost) models; see Lundberg 
#' et. al. (2020) for details.
#' 
#' @param object A fitted model object (e.g., a
#' `ranger::ranger()`, or` xgboost::xgboost()`,object, to name a few).
#'
#' @param feature_names Character string giving the names of the predictor
#' variables (i.e., features) of interest. If `NULL` (default) they will be
#' taken from the column names of `X`. 
#'
#' @param X A matrix-like R object (e.g., a data frame or matrix) containing 
#' ONLY the feature columns from the training data (or suitable background data 
#' set). If the input includes categorical variables that need to be one-hot encoded, 
#' please input data that has been processed using '`data.table::one_hot()`.
#' In XGBoost, the input should be the raw dataset containing only the explanatory variables, 
#' not the data created using `xgb.DMatrix.`
#'  **NOTE:** This argument is required whenever `exact = FALSE`.
#'
#' @param pred_wrapper Prediction function that requires two arguments,
#' `object` and `newdata`. **NOTE:** This argument is required 
#' whenever `exact = FALSE`. The output of this function should be 
#' determined according to:
#'
#' \describe{
#'   \item{Regression}{A numeric vector of predicted outcomes.}
#'   \item{Binary classification}{A vector of predicted class probabilities
#'   for the reference class.}
#'   \item{Multiclass classification}{A vector of predicted class probabilities
#'   for the reference class.}
#' }
#' 
#' @param nsim The number of Monte Carlo repetitions to use for estimating each 
#' Shapley value (only used when `exact = FALSE`). Default is 1. 
#' **NOTE:** To obtain the most accurate results, `nsim` should be set 
#' as large as feasibly possible.
#' 
#' @param newdata A matrix-like R object (e.g., a data frame or matrix) 
#' containing ONLY the feature columns for the observation(s) of interest; that 
#' is, the observation(s) you want to compute explanations for. Default is 
#' `NULL` which will produce approximate Shapley values for all the rows in 
#' `X` (i.e., the training data). 
#' If the input includes categorical variables that need to be one-hot encoded, 
#' please input data that has been processed using `data.table::one_hot()`.
#' 
#' 
#' @param exact Logical indicating whether to compute exact Shapley values. 
#' Currently only available for [stats::lm()](https://cran.r-project.org/web/packages/STAT/STAT.pdf), 
#' [xgboost::xgboost()](https://cran.r-project.org/web/packages/xgboost/xgboost.pdf), and [lightgbm::lightgbm()](https://cran.rstudio.com/web/packages/lightgbm/lightgbm.pdf) objects. 
#' Default is `FALSE`. Note that setting `exact = TRUE` will return 
#' explanations for each of the `stats::terms()` in an 
#'`stats::lm()` object. Default is `FALSE`.
#' 
#' 
#' @param parallel Logical indicating whether or not to compute the approximate
#' Shapley values in parallel across features; default is `FALSE`. **NOTE:**
#' setting `parallel = TRUE` requires setting up an appropriate (i.e., 
#' system-specific) *parallel backend* as described in the 
#' [foreach](https://cran.r-project.org/package=foreach); for details, see
#' `vignette("foreach", package = "foreach")` in R.
#' 
#' @param ... Additional arguments to be passed 
#' 
#' @return An object of class `Explain` with the following components :
#' \item{newdata}{The data frame formatted dataset employed for the estimation of Shapley values.
#'              If a variable has categories, categorical variables are one-hot encoded.} 
#' \item{phis}{A list format containing Shapley values for individual variables.}
#' \item{fnull}{The expected value of the model's predictions.}
#' \item{fx}{The prediction value for each observation.}
#' \item{factor_names}{The name of the categorical variable.
#' If the data contains only continuous or dummy variables, it is set to `NULL`.}
#' 
#' 
#' @note 
#' Setting `exact = TRUE` with a linear model (i.e., an 
#' `stats::lm()` or `stats::glm()` object) assumes that the
#' input features are independent.
#' 
#' @references 
#' Strumbelj, E., and Igor K. (2014). Explaining prediction models and 
#' individual predictions with feature contributions. Knowledge and information 
#' systems, 41(3), 647-665.
#' 
#' Lundberg, S. M., Erion, G., Chen, H., DeGrave, A., Prutkin, J. M., Nair, B.,
#' Katz, R., Himmelfarb, J., Bansal, N., and Lee, Su-In (2020). From local 
#' explanations to global understanding with Explainable AI for trees. 
#' Nature Machine Intelligence, 2(1), 2522–5839.
#' 
#' @rdname Explain.default
#' 
#' @export
#' 
#' @examples 
#' \donttest{
#' #
#' # A projection pursuit regression (PPR) example
#' #
#' 
#' # Load the sample data; see ?datasets::mtcars for details
#' data(mtcars)
#' 
#' # Fit a projection pursuit regression model
#' fit <- ppr(mpg ~ ., data = mtcars, nterms = 5)
#' 
#' # Prediction wrapper
#' pfun <- function(object, newdata) {  # needs to return a numeric vector
#'   predict(object, newdata = newdata)  
#' }
#' 
#' # Compute approximate Shapley values using 10 Monte Carlo simulations
#' set.seed(101)  # for reproducibility
#' shap <- Explain(fit, X = subset(mtcars, select = -mpg), nsim = 10, 
#'                 pred_wrapper = pfun)
#'}
#' 
Explain <- function(object, ...) {
  UseMethod("Explain")
} 


#' @rdname Explain.default
#' 
#' @export
#' 
Explain.default  <- function(object, feature_names = NULL, X = NULL,
                             nsim = 1, pred_wrapper = NULL,
                             newdata = NULL, parallel = FALSE,...) {

  i <- 0;
  # Compute baseline/average training prediction (fnull) and predictions
  # associated with each explanation (fx); if `adjust = FALSE`, then the
  # baseline is not needed and defaults to zero.
  if (is.null(X)) {
    stop("Training features required for approximate Shapley values. Please ",
         "supply them via the `X` argument; see `?fastshap::Explain` for ",
         "details.", call. = FALSE)
  }
  if (inherits(X, what = "tbl_df")) {
    X <- as.data.frame(X)
  }
  
  
  if (is.null(newdata)) {
      newdata <- as.data.frame( X )  # Explain all rows of background data set
  }
  
  if( inherits(object, "xgb.Booster")) {
    X <- as.matrix( X) 
    newdata <- as.matrix(newdata) 
  }  

  fx <- pred_wrapper(object, newdata = newdata)
  fnull <-  mean(pred_wrapper(object, newdata = X))


  # Deal with other NULL arguments
  if (is.null(feature_names)) {
    feature_names = colnames(X)
  }
  if (is.null(pred_wrapper)) {
    stop("Prediction function required for approximate Shapley values. Please ",
         "supply one via the `pred_wrapper` argument; see ",
         "`?fastshap::Explain` for details.", call. = FALSE)
  }

  # Set up the 'foreach' "do" operator
  `%.do%` <- if (isTRUE(parallel)) `%dopar%` else `%do%`

   # decoded data ----------------- # mltools: one_hot
   temp_new <- label_data  ( newdata) $  decoded_data
   temp_X <- label_data  ( X )$  decoded_data

   temp_feature <- names( temp_X )
   
   if( inherits(object,"xgb.Booster")) {
     temp_new <- as.matrix(temp_new) 
     temp_X <- as.matrix(temp_X) 
   }  

  # Compute approximate Shapley values # ,  ...
   phis_temp <- foreach(i = temp_feature ) %.do% {
    replicate(nsim, {  # replace later with vapply()
      Explain_column_default (object, X = temp_X , column = i,
                              pred_wrapper = pred_wrapper, newdata =  temp_new)
    })
    # number of post by obs = n matrix , list = number of variable
  }
  names( phis_temp) <- temp_feature


  tmp_var <- label_data  ( newdata) $ factor_check

  idx <- NULL
  if ( sum (tmp_var $ n >= 2) > 1 ){

  for (j in  unique(  tmp_var$var_tmp [  tmp_var $ uniq_len ==2 &   tmp_var $ n >=2]) ){
    idx <- c(idx ,which(stringr::str_detect( feature_names, paste0("^", j) )))
  }
 }

  phis  <-  foreach(i = 1:length( feature_names ) ) %.do% {

    temp <- matrix (0, nrow = dim (newdata)[1], ncol = nsim )

    if ( i %in% idx ) {
      temp [which(newdata[, i]==1),]  <-
        phis_temp[[tmp_var$var_tmp[i]]] [which(newdata[,i]==  1),]
      temp

    } else {
      temp <-  phis_temp[[feature_names[i]]]

    }

  }
  names( phis) <-  feature_names

  factor_names <- NULL
  factor_names <- names(X) [which((names(X) %in% names( temp_new)) ==FALSE)]

  out <- list (phis = phis, newdata = newdata,  fnull =  fnull, fx = fx, factor_names = factor_names )

  class(out) <- "Explain"
  return (out)

}

#' @rdname Explain.default
#'
#' @export
#' 
Explain.lm <- function(object, feature_names = NULL, X, nsim = 1, pred_wrapper,
                       newdata = NULL, exact = FALSE,  
                       parallel = FALSE, ...) {
  if (isTRUE(exact)) {  # use Linear SHAP
    phis <- if (is.null(newdata)) {
      predict(object, type = "terms" )
    } else {
      predict(object, newdata = newdata, type = "terms" )
    }
    
    df_phis <- mapply(function(col, name) {
      setNames(data.frame(col, stringsAsFactors = FALSE), name)
    },   data.frame(phis), names( data.frame(phis)), SIMPLIFY = FALSE)
    
    baseline <- attr(phis, which = "constant")  # mean response for all training data
    
    return(structure(list(
      "phis" = df_phis,
      "newdata" = newdata[, feature_names, drop = FALSE],
      "fnull" =  baseline,
      "fx" = predict(object,newdata = newdata, type = "response") ,
      "factor_names" = NULL
    ), class = "Explain"))
    
   }else{
    Explain.default (object, feature_names = feature_names, X = X, nsim = nsim,
              pred_wrapper = pred_wrapper, newdata = newdata, parallel = parallel )
  }
}

#' @rdname Explain.default
#'
#' @export

Explain.xgb.Booster <- function(object, feature_names = NULL, X = NULL, nsim = 1,
                                pred_wrapper, newdata = NULL, exact = FALSE,   parallel = FALSE, ...) {
  if (isTRUE(exact)) {  # use Tree SHAP
    if (is.null(X) && is.null(newdata)) {
      stop("Must supply `X` or `newdata` argument (but not both).",
           call. = FALSE)
    }
    X <- if (is.null(X)) newdata else X
    phis <- predict(object, newdata = X, predcontrib = TRUE,
                           approxcontrib = FALSE )

    pred_y <- predict(object, X)
    
    df_phis <- mapply(function(col, name) {
      setNames(data.frame(col, stringsAsFactors = FALSE), name)
    },   data.frame(phis), names( data.frame(phis)), SIMPLIFY = FALSE)
    

      return(structure(list(
        "phis" = df_phis,
        "newdata" = newdata[, feature_names, drop = FALSE],
        "fnull" =  unique(phis[, "BIAS"]),
        "fx" =  pred_y,
        "factor_names" =NULL
      ), class = "Explain"))

  } else {
    Explain.default (object, feature_names = feature_names, X = X, nsim = nsim,
                    pred_wrapper = pred_wrapper, newdata = newdata,
                    parallel = parallel, ...)
  }
}

#' @rdname Explain.default
#'
#' @export


Explain.lgb.Booster <- function(object, feature_names = NULL, X = NULL, nsim = 1,
                                pred_wrapper, newdata = NULL,
                                exact = FALSE,  parallel = FALSE, ...) {
  if (isTRUE(exact)) {  # use Tree SHAP
    if (is.null(X) && is.null(newdata)) {
      stop("Must supply `X` or `newdata` argument (but not both).",
           call. = FALSE)
    }
    X <- if (is.null(X)) newdata else X

    # Adapt LightGBM predict() interface
    if (packageVersion("lightgbm") > package_version("3.3.2")) {
      phis <- predict(object, X, type = "contrib", ...)
    } else {
      phis <- predict(object, X, predcontrib = TRUE, ...)
    }

    colnames(phis) <- c(colnames(X), "BIAS")
    
    df_phis <- mapply(function(col, name) {
      setNames(data.frame(col, stringsAsFactors = FALSE), name)
    },   data.frame(phis), names( data.frame(phis)), SIMPLIFY = FALSE)
    

    pred_y = predict(object, X)

      return(structure(list(
        "phis" =  df_phis ,
        "newdata" = newdata[, feature_names, drop = FALSE],
        "fnull" =  unique(phis[, "BIAS"]),
        "fx" =  pred_y,
        "factor_names" = NULL
      ), class = "Explain"))

  } else {
    Explain.default (object, feature_names = feature_names, X = X, nsim = nsim,
                    pred_wrapper = pred_wrapper, newdata = newdata,
                    parallel = parallel, ...)
  }
}
