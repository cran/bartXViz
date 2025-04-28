#' Visualization of Shapley values from the BARP model
#'
#'
#' This function is implemented to visualize the computed Shapley values in 
#' various ways for objects of the `Explainbarp` class. The type of plot 
#' generated depends on the input parameters.
#' Since the BARP model is designed to be visualized for a single stratum, 
#' the user must specify both the stratum variable and the value of the stratum to be visualized.
#'
#' @param x An `ExplainBARP` class object containing the Shapley values of the BARP model.
#' @param average  Input the reference value for calculating the mean of the object's phi list.
#' `"obs"` represents the average based on observations (#post by #variable),
#' while `"post"` represents the average based on posterior samples (#obs by #variable).
#' If `"both"` is entered, calculations are performed based on both observation and posterior sample criteria.
#' @param title The title of the plot, with a default value of `NULL`.
#' @param type `"bar"` represents a bar chart that includes the average contribution of each variable,
#'  while `"bee"` represents a summary plot, allowing you to determine the graph's format.
#' @param num_post To check the contribution of variables for a single posterior sample,
#'  enter a value within the number of posterior samples.
#' @param adjust The default value is `FALSE`.
#' Enter `TRUE` to check the Shapley values adjusted based on the model's average contribution.
#' @param plot.flag If `average = "obs"`, the quantile interval of each variable's is provided by default.
#' @param probs Enter the probability for the quantile interval. The default value is `0.95`.
#' @param geo.unit  Enter the name of the stratification variable used in post stratification.
#' @param geo.id Enter one value of interest among the values of the stratification variable.
#' @param ... Additional arguments to be passed
#' @return The plot is returned based on the specified option.: 
#' \item{out}{If average is `"obs"` or `"post"`, a bar plot or summary plot is generated based on
#'the selected averaging criterion. If average is `"both"`, a boxplot is displayed to
#'show the distribution of Shapley values computed using both criteria.
#' If adjust is `TRUE`, the adjusted Shapley values are displayed.
#' If `num_post` is specified, a bar plot or summary plot for the selected posterior sample is generated.}
#' @importFrom dplyr n
#' @importFrom dplyr mutate
#' @export

plot.Explainbarp <- function(x, average = NULL, type = NULL,  num_post= NULL,
                             plot.flag=TRUE, adjust= FALSE,  probs=0.95, title=NULL,geo.unit=NULL,geo.id=NULL,...){
  object <- x;
  value <- variable <- n <- Average_re <- 0;
  geo_check <- str_split(colnames(object $ newdata), paste0( geo.unit,"_"),simplify = T)[,2]
  index <- which (object $ newdata [,  which(geo_check==geo.id)] == 1 )
  
  if(isTRUE (adjust)) {
    
    baseline  <-  object $fnull[ object $fnull[ ,geo.unit]==geo.id,2]
    feature_names <- names(object $newdata) [which(geo_check=="")]
    new_data <- object $ newdata [index , geo_check =="" ]
    
    # Compute mean and variance of posterior samples
    acomb <- function(...) abind(..., along = 3)
    phis.stats <- foreach(i = feature_names, .combine = "acomb") %do% {
      cbind(rowMeans( object$phis[[i]][ index, ]  ),
            apply( object$phis[[i]][ index, ] , MARGIN = 1, FUN = var ))
    }
    
    for (i in seq_len(dim(phis.stats)[1L])) {  # loop through each observation
      err <- object$ fx [ index ]- sum(phis.stats[i, 1L, ]) - object$fnull[object$fnull[,geo.unit] ==geo.id ,2]
      if( max(phis.stats[i, 2L, ])==0) {  v <- 0  } else {
        v <- (phis.stats[i, 2L, ] / max(phis.stats[i, 2L, ])) * 1e6 }
      adj <- err[i] * (v - (v * sum(v)) / (1 + sum(v)))
      phis.stats[i, 1L, ] <- phis.stats[i, 1L, ] + adj  # adjust Shapley estimate
      
    }
    phis_adj <- phis.stats[, 1L, ]
    
    phis_adj <- as.data.frame(phis_adj)
    names ( phis_adj ) <- feature_names
    
    for ( i in  names (phis_adj) ){
      temp_i <- which ( new_data [,i] !=1)
      phis_adj[temp_i,i] <- NA
    }
    
    # Exclude cases where all values are not estimated.
    phis_adj <- phis_adj[, colSums(is.na(phis_adj)) != nrow(phis_adj)]
    
    local_obs <- long_data ( phis_adj, normalize = TRUE, absolute = TRUE)
    
    if (isTRUE(type =="bar")| is.null (type)){
      # bar plot
      Local_mean <- as.data.frame(abs( phis_adj ))
      
      
      sample_summary <-  local_obs $  long %>% filter (is.na (value)==F  )%>%
        select(c("variable", "value"))%>% group_by(variable) %>% mutate (n=n() ) %>%
        select (c("variable","n"))%>%distinct(.keep_all = T)
      
      sample_zero <- local_obs $  long %>% filter (is.na (value)==F  )%>%
        select(c("variable", "value"))%>% group_by(variable) %>% mutate (n=n() )%>%
        filter(value==0) %>% mutate (zero =n()) %>% distinct(.keep_all = T) %>% select(c("variable", "zero"))
      
      sample_summary <- left_join(sample_summary, sample_zero, by = "variable" )
      sample_summary [which( is.na(sample_summary$zero)), "zero"] <- 0
      
      sample_summary$ percent <- paste0(sprintf("%.1f", round(sample_summary$zero/sample_summary$n*100, 1)) ,"%")
      
      out <- bar_plot (Local_mean,probs = NULL, plot.flag = FALSE ,title  )
      ymin <-  (- diff(range(out$coordinates$limits$y)) * 0.02)
      
      out <- out +
        annotate("text",x=  sample_summary$variable , y =  ymin,
                 hjust = 1, label = sample_summary$percent , size = 3.5)
      
      print(out)
      
    } else if (isTRUE( type == "bees")){
      #  summary plot
      out <- summary_plot (local_obs,title )
      print(out)
    }
    
    
  }
  else {
    
    
    if (isFALSE (adjust) & is.null(average)){
      
      stop("If you want to calculate the average based on observations, enter 'obs'.",
           "If you want to calculate the average based on the posterior sample, enter 'post'.",
           "To compare the results based on observation and posterior sample, enter 'both'.",call. = FALSE)
      
    } else if ( isTRUE(average =="obs") ){
      
      
      new_data <- object $ newdata [index , geo_check =="" ]
      
      temp_phis  <- foreach(i = 1 : length(object $ phis) ) %do% {
        object $ phis[[i]][index,]
      }
      names (temp_phis) <- names (object $ phis)
      
      
      for ( i in  names (temp_phis) ){
        temp_i <- which ( new_data [,i] !=1)
        temp_phis[[i]][temp_i,] <- NA
        
      }
      
      
      # Average based on observation= #post by #variable
      obs_mean_shap <- foreach(i = 1 : length(object $ phis), .combine = 'cbind' ) %do% {
        apply( temp_phis[[i]] , 2, function(x) {mean(abs(x),na.rm = T)})
      }
      colnames(obs_mean_shap) <- names(object $ phis)
      
      
      # Exclude cases where all values are not estimated.
      obs_mean_shap <- obs_mean_shap[, colSums(is.nan(obs_mean_shap)) != nrow(obs_mean_shap)]
      
      
      out <- bar_plot (obs_mean_shap, probs, plot.flag= T,title )
      print ( out )
      
      if ( isTRUE(type =="bees") ) {
        message("In Average based on observation,\n a barplot is also provided by entering type='bees'.")
      }
      
    } else if ( isTRUE(average =="post") ){
      if (is.null(num_post)){
        
        new_data <- object $ newdata [index , geo_check =="" ]
        
        # Average based on posterior sample = #obs by #variable
        post_mean_shap <- foreach(i = 1: length(object $ phis), .combine = 'cbind' ) %do% {
          apply(object $ phis[[i]][index,], 1, function(x) {mean( x) })
        }
        colnames(post_mean_shap) <- names(object $ phis)
        
        for ( i in colnames(post_mean_shap) ){
          temp_i <- which ( new_data [,i] !=1)
          post_mean_shap[temp_i,i] <- NA
        }
        
        # Exclude cases where all values are not estimated.
        post_mean_shap <- post_mean_shap[, colSums(is.na (post_mean_shap)) != nrow(post_mean_shap)]
        
        
        # Change data format & mean values => local_total is List
        local_total <- long_data (post_mean_shap, normalize = TRUE, absolute = TRUE)
        
        if (isTRUE(type =="bar")| is.null (type)){
          
          # bar plot
          Local_mean <- as.data.frame(abs(post_mean_shap))
          
          sample_summary <-  local_total $  long %>% filter (is.na (value)==F  )%>%
            select(c("variable", "value"))%>% group_by(variable) %>% mutate (n=n() ) %>%
            select (c("variable","n"))%>%distinct(.keep_all = T)
          
          sample_zero <- local_total $  long %>% filter (is.na (value)==F  )%>%
            select(c("variable", "value"))%>% group_by(variable) %>% mutate (n=n() )%>%
            filter(value==0) %>% mutate (zero =n()) %>% distinct(.keep_all = T) %>% select(c("variable", "zero"))
          
          sample_summary <- left_join(sample_summary, sample_zero, by = "variable" )
          sample_summary [which( is.na(sample_summary$zero)), "zero"] <- 0
          
          sample_summary$ percent <- paste0(sprintf("%.1f", round(sample_summary$zero/sample_summary$n*100, 1)) ,"%")
          
          
          out <- bar_plot (Local_mean, probs = NULL, plot.flag = FALSE,title  )
          ymin <- (- diff(range(out$coordinates$limits$y)) * 0.02)
          
          out <- out +
            annotate("text",x= sample_summary$variable, y =  ymin,
                     hjust = 1, label = sample_summary$percent , size = 3.5)
          print ( out )
          
        } else if (isTRUE( type == "bees")){
          
          #  summary plot
          out <- summary_plot (local_total,title )
          print ( out )
        }
        
        
        
      } else if (isFALSE (is.null(num_post)) ){
        
        
        # num_post th posterior sample = #obs by #variable
        obs <- foreach(i = 1: length(object $ phis), .combine = 'cbind' ) %do% {
          object $ phis[[i]][index, num_post]
        }
        colnames(obs) <- names(object $ phis)
        
        new_data <- object $ newdata [index , geo_check =="" ]
        
        for ( i in colnames(obs) ){
          temp_i <- which ( new_data [,i] !=1)
          obs[temp_i,i] <- NA
        }
        
        # Exclude cases where all values are not estimated.
        obs <- obs[, colSums(is.na (obs)) != nrow(obs)]
        
        
        local_obs <- long_data (obs, normalize = TRUE, absolute = TRUE)
        
        if (isTRUE(type == "bees") | is.null (type)){
          # Local SHAP summary plot - i th posterior
          out <- summary_plot (local_obs ,title)
          annotate_figure(out, top =  text_grob(paste0("Sample number = ", num_post),
                                                hjust = 1, x = 0.95,vjust = 1.05, size = 10 ))
        } else if (isTRUE(type == "bar")) {
          
          out <- bar_plot (abs(obs),probs =NULL, plot.flag = FALSE,title  )
          
          sample_summary <-  local_obs $  long %>% filter (is.na (value)==F  )%>%
            select(c("variable", "value"))%>% group_by(variable) %>% mutate (n=n() ) %>%
            select (c("variable","n"))%>%distinct(.keep_all = T)
          
          sample_zero <- local_obs $  long %>% filter (is.na (value)==F  )%>%
            select(c("variable", "value"))%>% group_by(variable) %>% mutate (n=n() )%>%
            filter(value==0) %>% mutate (zero =n()) %>% distinct(.keep_all = T) %>% select(c("variable", "zero"))
          
          sample_summary <- left_join(sample_summary, sample_zero, by = "variable" )
          sample_summary [which( is.na(sample_summary$zero)), "zero"] <- 0
          
          sample_summary$ percent <- paste0(sprintf("%.1f", round(sample_summary$zero/sample_summary$n*100, 1)) ,"%")
          
          ymin <-  (- diff(range(out$coordinates$limits$y)) * 0.02)
          
          
          out <- out +
            annotate("text",x= sample_summary$variable, y =  ymin,
                     hjust = 1, label = sample_summary$percent , size = 3.5)
          annotate_figure(out, top =  text_grob(paste0("Sample number = ", num_post),
                                                hjust = 1, x = 0.95,vjust = 1.05, size = 10 ))
        }
      }
      
    }else if ( isTRUE(average =="both") ){
      
      new_data <- object $ newdata [index , geo_check =="" ]
      
      temp_phis <-  foreach(i = 1 : length(object $ phis) ) %do% {
        object $ phis[[i]][index,]
      }
      names (temp_phis) <- names (object $ phis)
      
      
      for ( i in  names (temp_phis) ){
        temp_i <- which ( new_data [,i] !=1)
        temp_phis[[i]][temp_i,] <- NA
        
      }
      
      # Average based on observation= #post by #variable
      obs_mean_shap <- as.data.frame( foreach(i = 1 : length(object $ phis), .combine = 'cbind' ) %do% {
        apply( temp_phis[[i]] , 2, function(x) {mean(abs(x),na.rm = T)})
      })
      colnames(obs_mean_shap) <- names(object $ phis)
      
      # Exclude cases where all values are not estimated.
      obs_mean_shap <- obs_mean_shap[, colSums(is.na (obs_mean_shap)) != nrow(obs_mean_shap)]
      
      mean_obs <- colMeans(obs_mean_shap) [order( colMeans(obs_mean_shap), decreasing = TRUE) ]
      mean_obs <- as.data.frame(mean_obs)
      names (mean_obs) <- "mean"
      mean_obs$ Average <- "Observation"
      mean_obs$ variable <- obs_variable <- rownames (mean_obs)
      
      obs_mean_shap $ Average <- "Observation"
      obs_mean_shap $ ID <-  rownames(obs_mean_shap)
      
      
      
      # Average based on posterior sample = #obs by #variable
      post_mean_shap <- as.data.frame(foreach(i = 1: length(object $ phis), .combine = 'cbind' ) %do% {
        apply( object $ phis[[i]][index,] , 1, function(x) {abs(mean( x ) )})
      })
      colnames(post_mean_shap) <- names(object $ phis)
      
      for ( i in colnames(post_mean_shap) ){
        temp_i <- which ( new_data [,i] !=1)
        post_mean_shap[temp_i,i] <- NA
      }
      
      # Exclude cases where all values are not estimated.
      post_mean_shap <- post_mean_shap[, colSums(is.na (post_mean_shap)) != nrow(post_mean_shap)]
      
      mean_post <- colMeans(post_mean_shap,na.rm = T) [order( colMeans(post_mean_shap,na.rm = T ), decreasing = TRUE) ]
      mean_post <- as.data.frame(mean_post)
      names (mean_post) <- "mean"
      mean_post $ Average <- "Samples"
      mean_post$ variable <- rownames (mean_post)
      
      post_mean_shap$ Average <- "Samples"
      post_mean_shap $ ID <- rownames(post_mean_shap)
      
      
      total_mean_shap <- rbind (obs_mean_shap, post_mean_shap )
      
      total_long <- melt(total_mean_shap, id =c('ID',"Average"))
      total_long $ Average_re <- fct_rev( total_long $ Average )
      
      mean_data <- rbind (mean_obs, mean_post)
      mean_data  $ Average_re <- fct_rev( mean_data  $ Average)
      
      data_list <- list(total_long = total_long, mean_data = mean_data  )
      
      
      if ( is.null (type) | isTRUE(type == "bees")) {
        
        ymax <- max (data_list $ total_long$value, na.rm = T) * (1+0.001)
        ymin <-  (- diff(range(data_list $ total_long$value, na.rm = T))*0.05)
        
        range_p <- (ymax-ymin)*0.02
        
        
        left_p <-  data_list $ total_long %>% filter (is.na (value)==F & Average_re=="Observation")%>%
          ggplot(  aes(x= variable, y=value, fill= Average_re )) +
          geom_boxplot() +coord_flip() + theme_bw(base_size = 11) +
          labs(fill = "", x="",y="Mean(|SHAP|)") +
          scale_x_discrete(limits = rev(obs_variable),
                           labels = rev(obs_variable))+
          scale_fill_manual(values=c(  "#F8766D"))+
          theme(legend.position =  c(0.8, 0.2), legend.background=element_blank()) + ylim(ymin, ymax )
        
        
        sample_summary <-  data_list $ total_long %>% filter (is.na (value)==F  & Average_re=="Samples")%>%
          select(c("variable", "value"))%>% group_by(variable) %>% mutate (n=n() ) %>%
          select (c("variable","n"))%>%distinct(.keep_all = T)
        
        sample_zero <- data_list $ total_long %>% filter (is.na (value)==F & Average_re=="Samples" )%>%
          select(c("variable", "value"))%>% group_by(variable) %>% mutate (n=n() )%>%
          filter(value==0) %>% mutate (zero =n()) %>% distinct(.keep_all = T) %>% select(c("variable", "zero"))
        
        sample_summary <- left_join(sample_summary, sample_zero, by = "variable" )
        sample_summary [which( is.na(sample_summary$zero)), "zero"] <- 0
        
        sample_summary $ percent <- paste0(sprintf("%.1f", round(sample_summary$zero/sample_summary$n*100, 1)) ,"%")
        
        right_p <-  data_list $ total_long %>% filter (is.na (value)==F & Average_re=="Samples")%>%
          ggplot(  aes(x= variable, y=value, fill= Average_re )) +
          geom_boxplot() +coord_flip() + theme_bw(base_size = 11) +
          labs(fill = "", x="",y="Mean(|SHAP|)") +
          scale_x_discrete(limits = rev(obs_variable),
                           labels = rev(obs_variable))+
          scale_fill_manual(values=c(  "#00BFC4"))+
          theme(legend.position =  c(0.8, 0.2), legend.background=element_blank()) + ylim(ymin, ymax )
        
        right_p <- right_p + annotate("text",x= sample_summary$variable,
                                      y =  - range_p,
                                      hjust = 1, label = sample_summary$percent , size = 3 )


       grid.arrange(left_p, right_p, ncol= 1, top=textGrob(title, x = 0.2, hjust = 0 ) )

      if ( isTRUE(type =="bees") ) {
        message("In Average = 'both',\n a baxplot is also provided by entering type='bees'.")
      }

    } else if (isTRUE(type == "bar")  ){
       out <- ggplot( data_list $ mean_data, aes(x= variable, y= mean, fill =  Average_re )) +
        geom_bar( stat = "identity" , position="dodge") +coord_flip() + theme_bw(base_size = 11) +
        scale_fill_manual(values=c( "#00BFC4","#F8766D"))+
        labs(fill = "Criteria for average ", x="",y="Mean(|SHAP|)") +
        scale_x_discrete(limits = rev(obs_variable), labels = rev(obs_variable)) +
        ylim(c(0,max( data_list $ mean_data$mean)  ))+
        geom_text(aes( label = sprintf("%.3f",round(mean,3)) ),
                  hjust ="inward",  size= 3.5, position = position_dodge(width=0.9))+
        theme(legend.position = "bottom")
       if(!is.null(title)) {
         out <- out + ggtitle(title)
       }

       print(out)
    }

  }
}
}

