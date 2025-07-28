#' Waterfall Plot
#'
#' The \code{waterfall_plot} function is a bar chart that displays the positive and
#' negative contributions across sequential data points, visualizing how each
#' variable's contributions change for a single observation.
#'
#' @param object Enter the name of the object that contains the model's contributions
#' and results obtained using the Explain function.
#' @param obs_num  observation number (only one)
#' @param title plot title
#' @param geo.unit The name of the stratum variable in the BARP model as a character.
#' @param geo.id Enter a single value of the stratum variable as a character.
#' @param obs_name Enter the name of the vector containing observation IDs or names.
#' @return The function returns a waterfall plot.
#' \item{plot_out}{The waterfall plot of the observation at index \code{obs_num}.}
#' @examples
#' \donttest{
#' ## Friedman data
#' set.seed(2025)
#' n <- 200
#' p <- 5
#' X <- data.frame(matrix(runif(n * p), ncol = p))
#' y <- 10 * sin(pi* X[ ,1] * X[,2]) +20 * (X[,3] -.5)^2 + 10 * X[ ,4] + 5 * X[,5] + rnorm(n)
#' 
#'## Using dbarts 
#'model <- dbarts::bart (X, y, keeptrees = TRUE, ndpost = 200)
#'
#'# prediction wrapper function
#'pfun <- function (object, newdata) {
#'  predict(object, newdata)
#'  }
#'  
#'# Calculate shapley values
#'model_exp <-  Explain(model, X = X, pred_wrapper=pfun)
#'
#'# Waterfall plot of 100th observation
#'waterfall_plot(model_exp, obs_num=100) 
#'}
#'
#' @export
waterfall_plot <- function (object, obs_num = NULL, title = NULL,
                            geo.unit=NULL, geo.id = NULL, obs_name = NULL){
  
  if (missing(object)) {
    message("Please provide an object containing the results from the 'Explain' function.")
    return(invisible(NULL))
  }
  if (is.null(obs_num)) {
    stop("Input for obs_num is required and must be either a single number.")
  }

  from <- to <- fill <- values <- 0; 
  
  if (is.null(geo.unit) == FALSE & is.null (geo.unit) ==FALSE){
    geo_check <- stringr::str_split(colnames(object$newdata), paste0( geo.unit,"_"),simplify = T)[,2]
    index <- which(object$newdata[,which(geo_check==geo.id)]== 1 )
    baseline <- object$fnull[object$fnull[,geo.unit]==geo.id,2]

    feature_names <- names(object$newdata)[which(geo_check=="")]
    geo_data <- object$newdata[index,]


      # Compute mean and variance of posterior samples
      acomb <- function(...) abind(..., along = 3)
      phis.stats <- foreach(i = feature_names, .combine = "acomb") %do% {
        cbind(rowMeans(object$phis[[i]][index,]),
              apply(object$phis[[i]][index,], MARGIN = 1, FUN = var))
      }

      for (i in seq_len(dim(phis.stats)[1L])) {  # loop through each observation
        err <- object$fx [ index ]- sum(phis.stats[i, 1L, ]) - object$fnull[object$fnull[,geo.unit] ==geo.id ,2]
        if( sum(phis.stats[i, 2L, ])==0) { v <- 0
        }else { v <- (phis.stats[i, 2L, ] / max(phis.stats[i, 2L, ])) * 1e6 }
        adj <- err[i] * (v - (v * sum(v)) / (1 + sum(v)))
        phis.stats[i, 1L, ] <- phis.stats[i, 1L, ] + adj  # adjust Shapley estimate
      }
      phis_adj <- phis.stats[, 1L, ]


  }else {
    baseline <- object$fnull
    feature_names <- names(object$newdata)
    
    if(inherits(object,"Explain") & dim (object$phis[[1]])[2]==1){

      phis_data <- as.matrix(do.call("cbind",object$phis))
      names (phis_data) <- names(object$phis)
      phis_adj <- phis_data

    } else {

    # Compute mean and variance of posterior samples
    acomb <- function(...) abind(..., along = 3)
    phis.stats <- foreach(i = feature_names, .combine = "acomb") %do% {
      cbind(rowMeans( object$phis[[i]] ), apply( object$phis[[i]], MARGIN = 1, FUN = var))
    }
    for (i in seq_len(dim(phis.stats)[1L])) {  # loop through each observation
      err <- object$fx - sum(phis.stats[i, 1L, ]) - object$fnull
      if( max(phis.stats[i, 2L, ])==0) { v <- 0 }
      else { v <- (phis.stats[i, 2L, ] / max(phis.stats[i, 2L, ])) * 1e6 }
      adj <- err[i] * (v - (v * sum(v)) / (1 + sum(v)))
      phis.stats[i, 1L, ] <- phis.stats[i, 1L, ] + adj  # adjust Shapley estimate
    }
    phis_adj <- phis.stats[, 1L, ]
    }
  }


  phis_adj <- as.data.frame(phis_adj)
  names (phis_adj) <- feature_names
  
   obj <- as.data.frame(t(phis_adj[as.numeric(obs_num),]))
   names(obj) <- "values"
   obj$names <- rownames(obj)

   if (is.null(geo.unit) == FALSE & is.null(geo.unit)==FALSE){
     temp <- as.data.frame(t(geo_data[obs_num,]))
     var_idx <- rownames(temp)[temp[,1]==1]
     obj <- obj[which(obj$names %in% var_idx),]
     obj <- obj[order(abs(obj$values)),]

     obj$to <- cumsum(obj$values) + baseline
     obj$from <- dplyr::lag(obj$to, default = baseline)
     obj$fill <- obj$to <  obj$from
     obj$i <- seq_len(nrow(obj))

   } else if(is.null(geo.unit) & is.null(geo.unit)) {

     if(isFALSE(is.null(object$factor_names))){

       temp <- as.data.frame(t(object$newdata[obs_num,]))
         name_temp <- NULL
         num_factor <- NULL
         for (i in object$factor_names){
           name_temp <- c(name_temp, str_which(rownames(temp),i))
           num_factor <- c(num_factor, sum(str_count(rownames(temp), i)))
         }

         if(is.null(name_temp)==F) {
           num_ind <- rownames(temp)[-name_temp]
         } else if(is.null(name_temp)) {
           num_ind <- rownames(temp)
         }
         
         dumm_ind <- NULL
         if(!inherits(object, "ExplainbartMachine")) {
           dumm_tmp <- object$factor_names[num_factor==1]
           for (i in dumm_tmp) {
             dumm_ind <- rownames(temp)[str_which(rownames(temp),i )]
           }
         }
        
       var_idx <- c(rownames(temp) [temp[,1]==1],num_ind ,dumm_ind)
       obj <- obj[which(obj$names %in% var_idx),]

     }

     obj <- obj[order(abs(obj$values)),]
     obj$to  <- cumsum(obj$values) + baseline
     obj$from  <- dplyr::lag(obj$to, default = baseline)
     obj$fill <- obj$to <  obj$from
     obj$i <- seq_len(nrow(obj))
   }

   plot_out <- ggplot(obj,aes(xmin = from,xmax = to,y = reorder(names, i),
                              fill = factor(fill, levels = c(FALSE, TRUE)))) +
     gggenes::geom_gene_arrow( show.legend = FALSE,
                               arrowhead_width = unit(2, "mm"),
                               arrowhead_height = unit(1 / (1 + 2 * nrow(obj)),"npc"),
                               arrow_body_height = unit(1 / (1 + 2 * nrow(obj)),"npc"),alpha= 0.8) +
     geom_fit_text(aes(label = paste0(round(values,3))), show.legend = FALSE) +
     scale_fill_manual(values =c("#FF0051","#008BFB"), drop = FALSE) +
     theme_bw() +
     theme( panel.border = element_blank(), panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(), axis.line.x = element_line(),
            axis.ticks.y =  element_blank()) +
     labs(y = element_blank(), x ="Prediction")  +
     annotate("text", x = obj$to[1], y = obj$i[1] - 0.5, label = paste0("E[f(x)]=",round(obj$from [1],3)),
              vjust = "inward", hjust = "inward", size = 3)+
     annotate("text",x = obj$to[dim(obj)[1]],y = obj$i[dim(obj)[1]]-0.5 ,
              label = paste0("f(x)=", round(obj$to[dim(obj)[1]],3)),
              vjust = "inward", hjust = "inward" , size = 3 )

   if(!is.null(title)) {
     plot_out <- plot_out + ggtitle(title)
   }

  if(is.null(obs_name) ) {
    plot_out <- annotate_figure(plot_out,
                                 top = text_grob(paste0("Observation number = ", obs_num),
                                                  hjust = 1, x = 0.95,vjust = 1.05, size = 10))
  } else if(is.null(obs_name)==FALSE) {
    plot_out <- annotate_figure(plot_out,
                                 top = text_grob(paste0("Observation number = ", obs_num,"\n",obs_name[obs_num]),
                                                  hjust = 1, x = 0.95,vjust = 1.05, size = 10))
  }

  print(plot_out)
}
