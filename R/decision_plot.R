#' Decision plot
#'
#' The decision_plot function is a graph that visualizes how individual features
#' contribute to a model's prediction for a specific observation using Shapley values.
#'  It can be used to visualize one or multiple observations.
#'
#' @param object Enter the name of the object that contains the model's contributions
#' and results obtained using the Explain function.
#' @param obs_num single or multiple observation numbers
#' @param title plot title
#' @param geo.unit The name of the stratum variable in the BARP model as a character.
#' @param geo.id Enter a single value of the stratum variable as a character.
#' @param bar_default bar_default is an option for adjusting the legend's color scale to fit the window length, and its default value is set to TRUE. 
#' If plots fail to render in LaTeX documents, it is recommended to set this option to FALSE.
#' @return 
#' \item{plot_out}{The decision plot for one or multiple observations specified in `obs_num`.}
#' @examples
#' \donttest{
#' ## Friedman data
#'set.seed(2025)
#'n = 200
#'p = 5
#'X = data.frame(matrix(runif(n * p), ncol = p))
#'y = 10 * sin(pi* X[ ,1] * X[,2]) +20 * (X[,3] -.5)^2 + 10 * X[ ,4] + 5 * X[,5] + rnorm(n)
#'
#'## BART model
#' model = dbarts::bart (X,y, keeptrees = TRUE,ndpost = 200 )
#' 
#'# prediction wrapper function
#'pfun <- function (object, newdata) {
#'predict(object, newdata)
#'}
#'
#'# Calculate shapley values
#'model_exp =  Explain  ( model, X = X,  pred_wrapper =  pfun )
#'
#'# Single observation 
#'decision_plot(model_exp, obs_num=1 )
#'
#'#Multiple observation 
#'decision_plot(model_exp, obs_num=10:40 )
#'}
#' @export

decision_plot <- function( object, obs_num, title=NULL,
                          geo.unit=NULL, geo.id=NULL, bar_default = TRUE ){
  ID <- value <- shap_cumsum <- yhat <- 0;
  
  if (is.null (geo.unit) == FALSE & is.null (geo.unit) == FALSE){
    geo_check <-  str_split(colnames(object $ newdata), paste0( geo.unit,"_"),simplify = T)[,2]
    index <- which (object $ newdata [,  which(geo_check==geo.id)] == 1 )
    baseline  <-  object $fnull[ object $fnull[ ,geo.unit]==geo.id,2]

    feature_names <- names(object $newdata) [which(geo_check=="")]

      # Compute mean and variance of posterior samples
      acomb <- function(...) abind(..., along = 3)
      phis.stats <- foreach(i = feature_names, .combine = "acomb") %do% {
        cbind(rowMeans( object$phis[[i]][ index, ]  ),
              apply( object$phis[[i]][ index, ] , MARGIN = 1, FUN = var ))
      }

      for (i in seq_len(dim(phis.stats)[1L])) {  # loop through each observation
        err <- object$ fx [ index ]- sum(phis.stats[i, 1L, ]) - object$fnull[object$fnull[,geo.unit] ==geo.id ,2]
        if( sum(phis.stats[i, 2L, ])==0) { v <- 0
        }else { v <- (phis.stats[i, 2L, ] / max(phis.stats[i, 2L, ])) * 1e6 }
        adj <- err[i] * (v - (v * sum(v)) / (1 + sum(v)))
        phis.stats[i, 1L, ] <- phis.stats[i, 1L, ] + adj  # adjust Shapley estimate
      }

      phis_adj <- phis.stats[, 1L, ]


  } else {
    baseline <-  object $fnull
    feature_names <- names(object $newdata)

    if(inherits(object,"Explain")& dim (object$phis[[1]])[2]==1){

      phis_data <- as.matrix(do.call("cbind",object$phis))
      names (phis_data) <- names(object$phis)
      phis_adj <- phis_data

    } else {

    # Compute mean and variance of posterior samples
    acomb <- function(...) abind(..., along = 3)
    phis.stats <- foreach(i = feature_names, .combine = "acomb") %do% {
      cbind(rowMeans( object$phis[[i]] ), apply( object$phis[[i]] , MARGIN = 1, FUN = var))
    }
    for (i in seq_len(dim(phis.stats)[1L])) {  # loop through each observation
      err <- object$ fx - sum(phis.stats[i, 1L, ]) - object$fnull
      if( sum(phis.stats[i, 2L, ])==0) {v <- 0
      }else {v <- (phis.stats[i, 2L, ] / max(phis.stats[i, 2L, ])) * 1e6}
      adj <- err[i] * (v - (v * sum(v)) / (1 + sum(v)))

      phis.stats[i, 1L, ] <- phis.stats[i, 1L, ] + adj  # adjust Shapley estimate
    }
    phis_adj <- phis.stats[, 1L, ]
    }
  }


  phis_adj <- as.data.frame(phis_adj)
  names ( phis_adj) <- feature_names

  mean_data <- as.data.frame(colMeans(abs(phis_adj))[order(colMeans(abs(phis_adj)) )])
  names ( mean_data ) <- "mean"
  mean_data $ variable <-  rownames ( mean_data)
  mean_data$ rank <- seq_len(nrow( mean_data))


  dt  <- setDT(phis_adj)[, names(phis_adj)[1:dim(phis_adj)[2]], with = FALSE]
  dt [, ID:= .I]
  dt_long <- melt.data.table(dt, measure.vars = colnames(phis_adj))


  dt_long <- merge(dt_long,  mean_data[ c("variable", "rank")], all.x=TRUE)

  local_t <- dt_long %>%
    arrange(ID, rank)%>% group_by(ID) %>%
    mutate(shap_cumsum = cumsum(value) + baseline  ) %>%
    mutate(ID = factor(ID)) %>%   mutate(yhat = last(shap_cumsum))


  leng_ID <-  length(unique (local_t  $ID))
  leng_var <- length(unique (local_t  $variable))
  leng_yhat <-  local_t  %>% filter ( rank ==1 )  %>% select(ID,yhat)

  basevalue <- data.frame(variable = rep("",leng_ID),ID = factor(1:leng_ID),value=rep(0,leng_ID),
                         rank = rep(0,leng_ID), shap_cumsum = baseline  ,yhat = leng_yhat$yhat )

  lastvalue <-  data.frame(variable = rep("",leng_ID),ID = factor(1:leng_ID),value= leng_yhat$yhat ,
                           rank = rep(leng_var+1,leng_ID), shap_cumsum = leng_yhat$yhat,yhat = leng_yhat$yhat )

  local_t <- rbind (  local_t ,  basevalue, lastvalue)


  local_select <- local_t  %>% filter (ID %in% obs_num)
  local_select <- local_select  %>% arrange ( ID,  rank )
  local_select <- as.data.frame( local_select )

  Yhat_var <- var( local_t $ shap_cumsum)
  Yhat_range <- range( local_t$ shap_cumsum )
  Yhat_range <- round(c(Yhat_range[1]-0.1*Yhat_var,Yhat_range[2]+0.1*Yhat_var),4 )

  plot_out <- local_select %>% ggplot(aes(x= rank , y = shap_cumsum, group = ID, color=yhat)) +
    geom_hline(yintercept =  baseline, color = "gray") +           coord_flip ()+
    geom_line() +   theme_bw(base_size = 12) + labs(x ="", y="Model output value") +
    theme(  panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
            legend.position = "top",legend.title = element_blank() ) +
    scale_y_continuous (limits = Yhat_range, sec.axis = sec_axis(~. )) +
    scale_x_discrete(limits = mean_data$variable, labels = mean_data$variable )

  if(isTRUE(bar_default)) {
    plot_out <- plot_out +
      scale_colour_stepsn( limits = Yhat_range, colours = c("blue", "red") ,
                           guide= guide_colorbar(barheight = 0.8,
                                                 barwidth =  unit(1, "npc") - sum(ggplotGrob(plot_out)[["widths"]][-7]) - unit(1,"line") ))
  }else if(isFALSE(bar_default)){
    plot_out <- plot_out +
      scale_colour_stepsn( limits = Yhat_range, colours = c("blue", "red"),
                           guide= guide_colorbar( barheight = 0.8, barwidth = 12))
  }

  if(!is.null(title)) {
    plot_out <- plot_out + ggtitle(title)
  }


  print( plot_out)
}

