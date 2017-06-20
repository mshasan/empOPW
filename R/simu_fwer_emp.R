#' @title Simulate Family Wise Error Rate (FWER) empirically
#'
#' @description This function simulate family wise error rate or test type I error rate
#'
#' @param s number of replication in a simulation
#' @param m total number of hypothesis test
#' @param alphaVec a vector of significance levels
#' @param max.group maximum number of p-value groups to be used, minimum is 5
#'
#' @details
#' This function generate pvalues times form \code{uniform(0, 1)} then applying
#' OPWeight method to obtain the Familly Wise Error Rate (FWER)
#'
#' @author Mohamad S. Hasan, \email{mshasan@uga.edu}
#' @export
#'
#' @seealso \code{\link{qvalue}}
#' \code{\link{prob_rank_givenEffect}}
#' \code{\link{weight_binary}}
#' \code{\link{weight_continuous}}

#'
#' @return a matrix of fwer for different methods
#'
#' @references Hasan and Schliekelman (2017)
#'
#' @examples
#' alphaVec = seq(.01, .1, .02)
#' simVal = 1:3  # in actual case use at least simVal = 1000
#' typeIerror_mat = sapply(simVal, simu_fwer, m = 100, alphaVec = alphaVec)
#'
#===============================================================================
#
# internal parameters:-----
# pval = pvalues from null tests
# pval_filter = filter pvalues from null tests
# test = test statistics
# filter = filter test statistics
#
#===============================================================================

simu_fwer_emp <- function(s, m, alphaVec, max.group = 10L)
    {
    fwer_per_rep <- function(alpha)
        {
            pval <- runif(m)
            pval_filter <- runif(m)
            test = qnorm(pval, lower.tail = FALSE)
            filter = qnorm(pval_filter, lower.tail = FALSE)

            pro_bin <- empOPW(pvalue = pval, filter = filter,
                              alpha = alpha, max.group = max.group,
                              effectType = "binary", method = "BON")$rejections

            pro_cont<- empOPW(pvalue = pval, filter = filter,
                              alpha = alpha, max.group = max.group,
                           effectType = "continuous", method = "BON")$rejections

            ihw_fwer <- ihw(pval, filter, alpha = alpha,
                                            adjustment_type = "bonferroni")

            bon = sum(pval <= alpha/m, na.rm = TRUE)
            IHW <- rejections(ihw_fwer)

            return(c(bon, pro_bin, pro_cont, IHW))
        }

        fwer_per_rep_mat = sapply(alphaVec, fwer_per_rep)
        return(fwer_per_rep_mat)
    }







