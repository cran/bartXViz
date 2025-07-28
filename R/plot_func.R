#' @importFrom bartMachine bartMachine calc_credible_intervals 
#' @importFrom bartMachine bart_machine_get_posterior 
#' @importFrom SuperLearner create.Learner 
#' @importFrom SuperLearner CVFolds SuperLearner.control 
#' @importFrom SuperLearner SuperLearner.CV.control
#' @importFrom abind abind
#' @importFrom reshape2 melt
#' @importFrom data.table setDT melt.data.table as.data.table .I ':=' set dcast setcolorder
#' @importFrom BART  bartModelMatrix
#' @importFrom dbarts makeModelMatrixFromDataFrame
#' @importFrom utils packageVersion setTxtProgressBar tail txtProgressBar
#' @importFrom ggpubr annotate_figure text_grob
#' @importFrom forcats fct_rev
#' @importFrom gridExtra grid.arrange
#' @importFrom grid textGrob unit
#' @importFrom ggforce geom_sina
#' @importFrom gggenes geom_gene_arrow
#' @importFrom ggfittext  geom_fit_text
#' @importFrom tidyr pivot_longer 
#' @importFrom stringr str_detect str_replace str_split str_which str_count
#' @importFrom Rcpp sourceCpp 
#' @importFrom foreach foreach %do% %dopar%
#' @importFrom stats var binomial  gaussian median model.frame model.matrix 
#' @importFrom stats quantile predict setNames reorder lag
#' @importFrom dplyr `%>%`  arrange distinct filter left_join group_by 
#' @importFrom dplyr group_by_  select starts_with  summarise  mutate last
#' @importFrom ggplot2 aes annotate coord_flip element_blank element_line  
#' @importFrom ggplot2 geom_bar geom_boxplot geom_errorbar geom_hline geom_line
#' @importFrom ggplot2 geom_point geom_text ggplot ggplotGrob ggtitle 
#' @importFrom ggplot2 guide_colorbar  position_dodge  scale_fill_manual 
#' @importFrom ggplot2 scale_x_discrete scale_y_continuous sec_axis 
#' @importFrom ggplot2 scale_color_gradient scale_colour_stepsn  xlab  
#' @importFrom ggplot2 ylab ylim  theme_bw labs theme




# barplot -----
bar_plot <- function (object, probs, plot.flag,title, ...) {
  lower<-upper<-group<-0;
  if (is.null(probs)) {
    object_summary <- data.frame(mean = apply(object,2,mean, na.rm=TRUE)  )
    object_summary$names <- rownames(object_summary)
  } else if (is.null(probs) == FALSE){
    object_summary <- data.frame(mean = apply(object,2,mean, na.rm=TRUE)   ,
                                lower = apply(object, 2, quantile,probs=(1-probs)/2, na.rm=TRUE) ,
                                upper=apply(object, 2, quantile,probs=1-(1-probs)/2, na.rm=TRUE) ,
                                median = apply(object,2,median, na.rm=TRUE)  )

    object_summary <- object_summary [order(-object_summary$mean),]
    object_summary$names <- rownames(object_summary)
    p <- nrow(object_summary)
    object_summary$flag <- c(FALSE,object_summary$upper[-1] - object_summary$median[-p]<0)
    object_summary$group <- factor(cumsum(object_summary$flag)+1)
  }



  if(plot.flag){

    plot_bar <- ggplot( object_summary, aes(x = reorder(names, +mean), y = mean) ) +
      geom_bar( stat = "identity" ,fill="grey80")  +
      ylab("Mean(|SHAP|)") + xlab("") +
      geom_errorbar(aes(ymin = lower,ymax= upper,color=group), width = 0.2,linewidth=1)+
      geom_point(aes(y=median),color="blue",size=3) +
      geom_text(aes(label = sprintf("%.3f",round(mean,3))), vjust = -0.3, hjust = -0.3,  size=3.5)+
      theme_bw()  + annotate("text", label = paste0(probs*100,"% quantile") ,x= 1,
                             y = max(object_summary$upper) , vjust = "inward", hjust = "inward") +
      coord_flip(ylim = c(0 , max(object_summary$upper) ))+
      theme(legend.position = "none")

  } else {
    ymin <- (max (object_summary$mean) * (-0.1))

    plot_bar <- ggplot( object_summary, aes(x = reorder(names, +mean), y = mean) ) +
      geom_bar( stat = "identity" ,fill="grey")  +
      ylab("Mean(|SHAP|)") + xlab("")  +
      geom_text(aes(label = sprintf("%.3f",round(mean,3))), ,  hjust = "inward",  size=3.5)+
      coord_flip(ylim = c(ymin , max(object_summary$mean) ))  + theme_bw()
  }

  if(!is.null(title)) {
    plot_bar <- plot_bar + ggtitle(title)
  }
  plot_bar

}


# long format

long_data <- function(data, normalize , absolute ) {
  ID <- variable <- value <- 0;
  if (isTRUE(absolute)) {
    data_list <- list(raw_data = data.frame(data),
                      mean_data = colMeans(abs(data), na.rm = TRUE)[order(colMeans(abs(data),na.rm = TRUE), decreasing = TRUE)])
  } else if (isFALSE(absolute)) {

    data_list <- list(raw_data = data.frame(data),
                      mean_data = as.data.frame(colMeans(abs(data),na.rm = TRUE )[order(colMeans(abs(data), na.rm = TRUE))]))

    names ( data_list$mean_data ) <- "mean"
    data_list$mean_data$variable <- rownames (data_list$mean_data)
    data_list$mean_data$rank <- seq_len(nrow(data_list$mean_data))

  }

  names ( data_list$raw_data) <- colnames (data)
  dt  <- setDT(data_list$raw_data)[, names(data_list$raw_data)[1:dim(data_list$raw_data)[2]], with = FALSE]
  dt [, ID:= .I]


  dt_long <- melt.data.table(dt, measure.vars = colnames(data_list$raw_data))

  if(isTRUE(normalize)) {
    dt_long  <- dt_long  %>% group_by(variable) %>% mutate (normalize = min_max_normalization(value))
  }

  out <- list (long = dt_long, mean = data_list$mean_data)
  return(out)
}


min_max_normalization <- function(x) {
  return ((x - min(x,na.rm = TRUE)) / (max(x,na.rm = TRUE) - min(x,na.rm = TRUE)))
}




# summary plot ------------------

summary_plot <- function (object,title,...) {
  value <- variable <- normalize <- 0;
  plot_summary <- object$long %>% filter(is.na(value)==FALSE) %>%
     ggplot(aes(variable, value)) +
    geom_sina(aes(color= normalize*12),
                       method = "counts", maxwidth = 0.7  ) +
    scale_x_discrete(limits = rev(names (object$mean)),
                     labels = rev(names (object$mean))) +
    scale_color_gradient(low="blue", high= "red", breaks=c(0,12), labels=c("Low", "High"),
                         guide = guide_colorbar(barwidth = 0.3, barheight = 12))+
    theme_bw(base_size = 12) + coord_flip()  +
    theme(axis.line.y = element_blank(),axis.ticks.y = element_blank(),
          legend.position="right",legend.title = element_blank() ) +
    labs(y = "SHAP (Impact on model output)", x = "" )

  if(!is.null(title)) {
    plot_summary  <- plot_summary  + ggtitle(title)
  }
  plot_summary
}


# one-hot encoding -> decoded
label_data <- function(data) {
  n <- value <- 0;
  data <- as.data.frame(data)

  var_tmp <- uniq_len <- factor_value <- NULL
  for (i in  names(data) ) {
    temp <- unlist  (strsplit(i, split =  "_" ) )
    var_tmp <-  c( var_tmp, temp[1] )
    uniq_len <- c( uniq_len, length(unique(data[,i])))
    factor_value <- c(factor_value,temp[length(temp)] )
  }

  factor_check  <-  data.frame(varname =  names(data),
                              var_tmp = var_tmp, uniq_len=uniq_len,
                              factor_value = factor_value)
  factor_check  <-  factor_check %>% group_by(var_tmp) %>% dplyr::add_count(var_tmp, name = "n") 

  fac <- unique(factor_check$var_tmp [factor_check$uniq_len ==2 & factor_check$n >=2])

  decoded_data <- data
  for (i in fac){
    decoded_data <-  decoded_data  %>% 
      pivot_longer(cols=starts_with(i), names_to=i, 
                   names_prefix= paste0(i,"."))  %>%
      filter(value==1) %>% select (-value)
  }

  decoded_data <- as.data.frame(decoded_data)

  out <- list (decoded_data = decoded_data, factor_check = factor_check)
  return (out)
}

