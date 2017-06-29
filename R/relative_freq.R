#' @title compute relative frequency
#'
#' @description A function to comoute relative frequency efficiently
#'
#' @param bin_idx Integer, bin number of the histogram
#' @param h_breaks number of breaks for the histogram
#' @param obs A numeric vector of observations
#'
#' @details This function compute relative frequency effectively using tabulate
#'
#' @author Mohamad S. Hasan, mshasan@uga.edu
#'
#' @export
#'
#' @import stats
#'
#' @return \code{rf} A numerical vector of the relative frequencies
#'  \code{rf_one} A numerical value of the relative frequency
#' of a particular bin
#'
#' @examples
#'
#' # generating data (known in practice)
#' X = runif(1000, min = 0, max = 2)         # covariate
#' r_freq <- relative_freq(bin_idx = 1, h_breaks = 20, obs = X)
#'
#===============================================================================
# function to compute p(rank=k|filterEffect=ey) emperically

# internal parameters:-----
# bin = bin points
# bin.counts = number of point per bin
# rel_freq_all = relative frequency of all bins
# rel_freq_one = relative frequency of one bin

#===============================================================================

relative_freq <- function(bin_idx = 1L, h_breaks = 20L, obs)
{
    bin <- c(0, (1:h_breaks)/h_breaks)
    bin.counts <- tabulate(cut(obs, bin))
    rel_freq_all = bin.counts/sum(bin.counts)
    rel_freq_one = rel_freq_all[bin_idx]
    return(list(rf = rel_freq_all, rf_one = rel_freq_one))
}


