#' A Function for Visualizing the Shapley Values of BART Models
#'
#' The \code{plot.ExplainbartMachine} function provides various visualization methods for Shapley values. 
#' It is designed to visualize \code{ExplainbartMachine} class objects, which contain Shapley values computed from models estimated using the \code{bartMachine} function from the \pkg{bartMachine}.
#' The values and format used in the graph are determined based on the input parameters.
#'
#' @param x An \code{ExplainbartMachine} class object containing the Shapley values of the BART model.
#' @param average  Input the reference value for calculating the mean of the object's \code{phi} list.
#' \code{"obs"} represents abind the average based on observations (#post by #variable),
#' while \code{"post"} represents the average based on posterior samples (#obs by #variable).
#' If \code{"both"} is entered, calculations are performed based on both observation and posterior sample criteria.
#' If no value is specified, "both" is used as the default.
#' @param title The title of the plot, with a default value of \code{NULL}.
#' @param type \code{"bar"} represents a bar chart that includes the average contribution of each variable,
#'  while \code{"bee"} represents a summary plot, allowing you to determine the graph's format.
#' @param num_post To check the contribution of variables for a single posterior sample,
#'  enter a value within the number of posterior samples.
#' @param adjust The default value is \code{FALSE}.
#' Enter \code{TRUE} to check the Shapley values adjusted based on the model's average contribution.
#' @param plot.flag If \code{average = "obs"}, the quantile interval of each variable's is provided by default.
#' @param probs Enter the probability for the quantile interval. The default value is \code{0.95}.
#' @param ... Additional arguments to be passed
#' @return The plot is returned based on the specified option.: 
#' \item{out}{If average is \code{"obs"} or \code{"post"}, a bar plot or summary plot is generated based on the selected averaging criterion. If average is \code{"both"}, a boxplot is displayed to
#'show the distribution of Shapley values computed using both criteria. If adjust is \code{TRUE}, the adjusted Shapley values are displayed.
#' If \code{num_post} is specified, a bar plot or summary plot for the selected posterior sample is generated.}
#' @export
#' 
plot.ExplainbartMachine <- function(x, average = NULL, type = NULL,  num_post= NULL,
plot.flag=TRUE, adjust= FALSE,  probs=0.95, title=NULL,...){

  if (!inherits(x, "ExplainbartMachine") ) {
    message("The input object must be of class 'ExplainbartMachine'.")
    return(invisible(NULL))
  } else if (missing(x) ) {
    stop("The object containing the Shapley values of models obtained from the 'Explain' function is missing.")
    return(invisible(NULL))
  }
  
    plot.ExplainBART (x= x, average=average , type = type,
                      adjust =  adjust, num_post =  num_post,
                      plot.flag = plot.flag , probs= as.numeric(probs),
                      title = title)

}
