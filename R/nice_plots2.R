#' @title Funciton to plot nice plots
#'
#' @description
#' This function takes a vector of xvalues and a vector or matrix of yvalues
#' then plot
#'
#' @param x_vec A numeric vector of x-axis values
#' @param y_matrix A numeric matrix of y-axis values
#' @param x_name Character, name of the x variable
#' @param y_names Character vector, name of the y variables
#' @param xlab Character, lebel of x-axis
#' @param ylab Character, label of y-axis
#' @param title Character, title of the plot
#'
#' @details
#' This function is desigend to plot nice curves
#'
#' @author Mohamad S. Hasan, \email{shakilmohamad7@gmail.com}
#' @export
#'
#'
#' @return \code{Data}
#' A plot of multiple curves
#'
#' @references Hasan and Schliekelman (2017)
#'
#' @examples
#' x_vec <- 1:10
#' y_matrix = cbind(sort(runif(10)), sort(rnorm(10)))
#' plt <- nice_plots2(x_vec = x_vec, y_matrix = y_matrix, x_name = "groups",
#'              y_names = c("uniform", "normal"),
#'              xlab = "groups", ylab = "Distribution", title = "uniform vs. normal")
#===============================================================================
# function to generate nice plots------------

# internal parameters:-----
# dat = data frame
# dat_melt = melted data
# g_plot = ggplot
#
#===============================================================================
nice_plots2 <- function(x_vec, y_matrix, x_name = NULL, y_names = NULL,
                        xlab = NULL, ylab = NULL, title = NULL)
{

    dat <- data.frame(x_vec, y_matrix)
    names(dat) <- c(x_name, y_names)
    dat_melt <- melt(dat, id.var = x_name, variable.name = "y_comn")
    g_plot <- ggplot(dat_melt, aes_string(x = names(dat_melt)[[1]], y = "value",
                                    group = "y_comn", colour = "y_comn")) +
        geom_line(size = 1.5) +
        labs(x = xlab, y = ylab, title = title) +
        theme(legend.position = "none")
    return(g_plot)
}





