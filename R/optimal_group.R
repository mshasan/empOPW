#' @title Obtain optimal number of groups and degrees of freedom
#'
#' @description A function to obtain the optimal number of groups and the spline
#' degrees of freedom from the p-values and the filters by applying the smooting
#' spline regression
#'
#' @param group number of groups to be used to split the p-values, default is five
#' @param pvalue a vector of pvalues of the test statistics
#' @param filter a vector of filter statistics
#' @param h_breaks number of breaks to be used for the histogram, default is 71
#' @param m total number of tests
#' @param m1 number of true alternatve tests
#' @param alpha significance level of the hypothesis test
#' @param mean_testEffect mean test effect of the true alterantives
#' @param effectType type of effect sizes; c("continuous", "binary")
#' @param method type of methods is used to obtain the results; c("BH", "BON"),
#' Benjemini-Hochberg or Bonferroni
#'
#' @details Optimal group and degrees of freedom of the spline regression is a
#' vital parameter to maximize the number of rejections. This function uses
#' smooting spline regresion to obtain that.
#'
#' @author Mohamad S. Hasan
#'
#' @export
#'
#' @return the number of rejected tests and the corresponding Optimal number of
#' groups and the degrees of freedom of the spline smo0thing
#'
#' @examples
#' # generate pvalues and filter statistics
#' m = 10000
#' set.seed(3)
#' filters = runif(m, min = 0, max = 2.5)          # filter statistics
#' H = rbinom(m, size = 1, prob = 0.1)             # hypothesis true or false
#' tests = rnorm(m, mean = H * filters)            # Z-score
#' pvals = 1 - pnorm(tests)                        # pvalue
#'
#' results <- optimal_group(group = 10, pvalue = pvals, filter = filters,
#'              h_breaks = 71, m = m, m1 = 8000, alpha = .05,
#'              mean_testEffect = 2.5, effectType = "continuous", method = "BH")
#'
#===============================================================================
# function to find the optimal number of groups and degrees of freedom---------
optimal_group <- function(group = 5L, pvalue, filter, h_breaks = 71L, m, m1,
           alpha = .05, mean_testEffect, effectType = c("continuous", "binary"),
           method = c("BH", "BON"))
{

    # ranks probability--------------
    ranksProb <- prob_rank_givenEffect_emp(pvalue = pvalue, filter = filter,
                group = group, h_breaks = h_breaks, effectType = effectType)

    # weights-----------
    if(effectType == "continuous"){
        wgt = weight_continuous(alpha = alpha, et = mean_testEffect,
                                m = group, ranksProb = ranksProb)
    } else {
        wgt = weight_binary(alpha = alpha, et = mean_testEffect, m = group,
                            m1 = m1/m*group, ranksProb = ranksProb)
    }

    # weight for all test----------
    grpSize <- ceiling(m/group)
    wgt_all = rep(wgt, each = grpSize)[1:m]

    # count number of rejections-----------
    if(method == "BH"){
        padj <- p.adjust(pvalue/wgt_all, method = method)
        n_rejections = sum(padj <= alpha, na.rm = TRUE)
    } else {
        n_rejections = sum(pvalue <= alpha*wgt_all/m, na.rm = TRUE)
    }

   return(c(group = group, n_rej = n_rejections))
}



