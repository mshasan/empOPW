#' @title Perform Empirical Optimal P-value Weighting
#'
#' @description A function to perform weighted p-value multiple hypothesis test.
#' This function compute the ranks probability of the test statistics by the filter
#' statistics given the effect sizes, and consequently the weights if neighter
#' the weights nor the probabilities are given. Then provides the number of rejected null
#' hypothesis and the list of the rejected pvalues as well as the corresponing
#' filter statistics.
#'
#' @param pvalue a vector of pvalues of the test statistics
#' @param filter a vector of filter statistics
#' @param weight optional weight vector not required
#' @param ranksProb ranks probabilities of the test-groups by the filters given
#' the mean effect. Note that, for each group of tests ranks probbaility would
#' be the same.
#' @param mean_testEffect mean test effect of the true alterantives
#' @param alpha significance level of the hypothesis test
#' @param tail right-tailed or two-tailed hypothesis test. default is right-tailed test.
#' @param delInterval interval between the \code{delta} values of a sequence.
#' Note that, \code{delta} is a LaGrange multiplier, necessary to normalize the weight
#' @param group number of groups to be used to split the p-values, default is five
#' @param  h_breaks number of breaks to be used for the histogram, default is 101
#' @param effectType type of effect sizes; c("continuous", "binary")
#' @param method type of methods is used to obtain the results; c("BH", "BON"),
#' Benjemini-Hochberg or Bonferroni
#' @param ... Arguments passed to internal functions
#'
#' @details If one wants to test \deqn{H_0: epsilon_i = 0 vs. H_a: epsilon_i > 0,}
#' then the \code{mean_testEffect}  and \code{mean_filterEffect} should be mean
#' of the test and filter effect sizes, respectively. This is called hypothesis
#' testing for the continuous effect sizes.\cr
#'
#' If one wants to test \deqn{H_0: epsilon_i = 0 vs. H_a: epsilon_i = epsilon,}
#' then \code{mean_testEffect} and \code{mean_filterEffect} should be median or
#' any discrete value of the test and filter effect sizes. This is called hypothesis
#' testing for the Binary effect sizes, where \code{epsilon} refers to a fixed value.\cr
#'
#' The main goal of the function is to compute the probabilities of the ranks from
#' the pvalues ranked by the filter statistics, consequently the weights.
#' Although \code{weights} \code{ranksProb} are optional, \code{empOPW} has the
#' options so that one can compute the probabilities and the weights externally
#' if necessary (see the examples).\cr
#'
#' Internally, \code{empOPW} function compute the \code{ranksProb} and consequently
#' the weights, then uses the p-values to make conclusions about hypotheses.
#' Although \code{ranksProb} is not required to the function,
#' One can compute \code{ranksProb} empirically by using the function
#' \code{\link{prob_rank_givenEffect_emp}}.\cr
#'
#' The function internally compute \code{mean_testEffect} from the test statistics,
#' which is obtainde from the p-values.
#'
#' @author Mohamad S. Hasan
#'
#' @export
#'
#' @import tibble tibble
#'
#' @seealso \code{\link{prob_rank_givenEffect_emp}} \code{\link{weight_binary}}
#' \code{\link{weight_continuous}}
#'
#'
#' @return \code{totalTests} total number of hypothesis tests evaluated
#' @return \code{nullProp} estimated propotion of the true null hypothesis
#' @return \code{ranksProb} probability of the ranks given the mean filter effect,
#' p(rank | ey = mean_filterEffect)
#' @return \code{weight} normalized weight
#' @return \code{rejections} total number of rejections
#' @return \code{rejections_list} list of rejected pvalues and the corresponding
#' filter statistics
#'
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
#' # general use
#' results <- empOPW(pvalue = pvals, filter = filters, effectType = "continuous",
#'                                               method = "BH")
#'
#' # supply the mean test effect externally
#' library(qvalue)
#' nullProp = qvalue(p = pvals, pi0.method = "bootstrap")$pi0
#' m0 = ceiling(nullProp*m)
#' m1 = m - m0
#'
#' et = mean(sort(tests, decreasing = TRUE)[1:m1])
#' results2 <- empOPW(pvalue = pvals, filter = filters, mean_testEffect = et,
#'                tail = 2, effectType = "continuous", method = "BH")
#'
#' # supply the ranks probability externally
#' grp = 5
#' probs = prob_rank_givenEffect_emp(pvalue = pvals, filter = filters, group = grp,
#'                                h_breaks = 101, effectType = "continuous")
#' results3 <- empOPW(pvalue = pvals, filter = filters, ranksProb = probs,
#'                  effectType = "continuous", tail = 2, method = "BH")
#'
#' # supply weight externally
#' wgt <- weight_continuous(alpha = .05, et = et, m = grp, ranksProb = probs)
#' results4 <- empOPW(pvalue = pvals, filter = filters, weight = wgt,
#'                         effectType = "continuous", alpha = .05, method = "BH")
#'
#===============================================================================
# # function to apply empOPW methods on data
#---------------------------------------------------
# Input:
#----------------------------
# pvalue = vector of pvalues
# filter = vector of filter statistics
# ranksProb = the probabilities of the filters or filters' ranks given
# the mean of the filter effects
# mean_filterEffect = filter effect size
# mean_testEffect = test effect size
# effectType = type of effect size c("binary","continuous")
# alpha = significance level of the hypotheis test
# nrep = number of replications for importance sampling
# tail = right-tailed or two-tailed hypothesis test. default is two-tailed test
# delInterval = interval between the delta values of a sequence
# method = Benjamini_HOchberg (BH) or Bonferroni (BON)

# internal parameters:-----
# m = number of hypothesis test
# nullProp = proportion of true null hypothesis
# m0 =  number of the true null tests
# m1 = number of the true alternative tests
# test =  compute test statistics from the pvalues if not given
# test_effect_vec = estiamted number of the true alternaitve test statistics
# mean_testEffect = mean test effect sizes of the true alternaive hypotheis
# mean_filterEffect = mean filter effect sizes of the true alternaive hypotheis
# ranksProb = probailities of the ranks given the mean effect size
# wgt = weights
# Data = create a data set
# OD = odered by covariate
# odered.pvalues = odered pvalues for all tests
# padj = adjusted pvalues for FDR uses

# Output:
#-------------------------
# totalTests = total number of hypothesis tests
# nullProp = estimated propotion of the true null hypothesis
# ranksProb = probability of the ranks given the mean filter effect,
#                                           p(rank | ey = mean_filterEffet)
# weight = normalized weight
# rejections = total number of rejections
# rejections_list = list of rejected pvalues and corresponding filter statistics
#-------------------------------------------------------------------------------

empOPW <- function(pvalue, filter, weight = NULL, ranksProb = NULL, mean_testEffect = NULL,
                alpha = .05, tail = 1L, delInterval = .0001, group = 5L, h_breaks = 101L,
                effectType = c("continuous", "binary"), method = c("BH", "BON"), ... )
{
    # compute the number of tests------------
    m = length(pvalue)
    nullProp = qvalue(p = pvalue, pi0.method = "bootstrap")$pi0
    m0 = ceiling(nullProp*m)
    m1 = m - m0


    # formulate a data set-------------
    Data = tibble(pvalue, filter)
    OD <- Data[order(Data$filter, decreasing = TRUE), ]
    OD_pvalue <- OD$pvalue


    #check whether weight is provided------------
    if(!is.null(weight)){
        wgt <- weight
    } else {

        # estimate the mean test effect size-------------
        if(!is.null(mean_testEffect)){
            mean_testEffect <- mean_testEffect
        } else {

            # compute test statistics from the pvalues---------
            test <- qnorm(pvalue/tail, lower.tail = FALSE)
            test[which(!is.finite(test))] <- NA

            # estimate the true alterantive test effect sizes----------------
            if(m1 == 0){
                test_effect_vec <- 0
            } else {
                test_effect_vec <-  sort(test, decreasing = TRUE)[1:m1]
            }

            # compute mean test effect----------
            if(effectType == "continuous"){
                mean_testEffect <- mean(test_effect_vec, na.rm = TRUE)
            } else {
                mean_testEffect <- median(test_effect_vec, na.rm = TRUE)
            }
        }

        #check whether filter ranks probability is provided------------
        if(!is.null(ranksProb)){
            ranksProb <- ranksProb
        } else {

            message("computing ranks probabilities")
            # compute the ranks probability of the tests given the mean effect
            ranksProb <- prob_rank_givenEffect_emp(pvalue = pvalue, filter = filter,
                            group = group, h_breaks = h_breaks, effectType = effectType)
            message("finished computing the ranks probabilities")
        }

        # compute the weights (always right-tailed)------------
        message("computing weights")
        if(effectType == "continuous"){
            wgt = weight_continuous(alpha = alpha, et = mean_testEffect, m = group,
                                    ranksProb = ranksProb)
        } else {
            wgt = weight_binary(alpha = alpha, et = mean_testEffect, m = group,
                                m1 = m1/m*group, ranksProb = ranksProb)
        }
        message("finished computing the weights")
    }

    wgt_all = rep(wgt, each = m/group)

    message("comparing pvalues with thresholds")
    if(method == "BH"){
        padj <- p.adjust(OD_pvalue/wgt_all, method = "BH")
        rejections_list = OD[which((padj <= alpha) == TRUE), ]
    } else {
        rejections_list = OD[which((OD_pvalue <= alpha*wgt_all/m) == TRUE), ]
    }


    # outputs--------------
    n_rejections = dim(rejections_list)[1]

    return(list(totalTests = m, nullProp = nullProp,
                ranksProb = ranksProb, weight = wgt_all,
                rejections = n_rejections, rejections_list = rejections_list))
}




