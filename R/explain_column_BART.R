Explain_column_BART <- function(object, X, column, pred_wrapper, newdata ) {

  # Find column position if feature name given
  if (is.character(column)) {
    column <- which(column == colnames(X))
  }


  # Extract dimensions of X (before possible subsetting)
  n <- nrow(X)  # number of training instances
  p <- ncol(X)  # number of features

  # Generate original and sampled feature instances
  if (is.null(newdata)) {  # FIXME: Should sampling be done with replacement?
    W <-   X[sample(n, replace = TRUE), ]
    O <-  genOMat(n, p)
  } else {
    W <- X[sample(n, size = nrow(newdata), replace = TRUE), , drop = FALSE]  # randomly sample rows from full X
    O <-  genOMat(n, p)[sample(n, size = nrow(newdata), replace = TRUE), , drop = FALSE]
    X <- newdata  # observations of interest
  }

  O <- if (column == 1) {  # case 1
    cbind(TRUE, O)
  } else if (column == p) {  # case 2
    cbind(O, TRUE)
  } else {  # case 3
    cbind(
      O[, 1:(column - 1), drop = FALSE],
      TRUE,
      O[, column:(p - 1), drop = FALSE]
    )
  }

  if (is.matrix(X)) {

    # Use RcppArmadillo for slight performance gain
    B <-  genFrankensteinMatrices(X, W, O, feature = column)
    colnames(B[[1L]]) <- colnames(B[[2L]]) <- colnames(X)

  } else {

    # Use base R's logical subsetting
    B <- list(X, X)
    B[[1L]][O] <- X[O]
    B[[1L]][!O] <- W[!O]
    O[, column] <- FALSE
    B[[2L]][O] <- X[O]
    B[[2L]][!O] <- W[!O]

  }

  # Make sure each component of `B` has the same column classes as `X`
  B[[1L]] <-  copy_classes(B[[1L]], y =  X)
  B[[2L]] <-  copy_classes(B[[2L]], y =  X)


  B[[1L]] <-  bartModelMatrix(B[[1L]] )
  B[[2L]] <-  bartModelMatrix(B[[2L]] )

  # Return differences in predictions
  pred_wrapper( object, newdata = B[[1L]]) -
    pred_wrapper( object, newdata = B[[2L]])


}
