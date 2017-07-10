#' @title Emperical ranks probbaility of the test given the effect size
#'
#' @description Emperical comnputation of the ranks probability of a test being
#' higher than any other test given the effect size from the external information.
#'
#' @param alpha Nmeric, significance level of the hypothesis test
#' @param pvalue a vector of pvalues of the test statistics
#' @param filter a vector of filter statistics
#' @param N_current Integer vector, number of observations per test in the
#' current data
#' @param N_prior Integer vector, number of observations per test in the
#' prior data
#' @param tail right-tailed or two-tailed hypothesis test.
#' default is right-tailed test.
#' @param max.group maximum number of groups to be used to split the p-values,
#' default is five. Note that, it is better to keep approximately 1000 p-values
#' per group.
#' @param standarized Character of c("TRUE" or "FALSE") determine whether
#' standarization is required. Default is FALSE means filter is a vector of
#' zscore. Thus, standarization will not used.
#' @param effectType Character of type ("binary" or"continuous") of effect sizes
#'
#' @details Perform data analysis for the different methods such as proposed,
#' bonferroni, Benjamini and Hoghburgh, IHW, and Dorbibian methods.
#'
#' @author Mohamad S. Hasan, shakilmohamad7@gmail.com
#'
#' @export
#'
#'
#' @return \code{rejections} A numeric vector of the number of rejected test of
#' the different methods.
#'
#' @examples
#'
#' # generating data (known in practice)
#' set.seed(123)
#' m = 10000
#' X = runif(m, min = 0, max = 2.5)               # covariate
#' H = rbinom(length(X), size = 1, prob = 0.1)   # hypothesis true or false
#' Z = rnorm(length(X), mean = H * X)            # Z-score
#' p = 1 - pnorm(Z)
#' rejections <- data_analysis(alpha = .1, pvalue = p, filter = X, N_current = m,
#'          N_prior = m, tail = 2, max.group = 10, effectType = "continuous")
#'
#===============================================================================
# sigma = a vector of standard deviaitons
# z_prior = prior information for the bayes_weight. In the bayes_weight
# it is treated as the zscore with left-tailed test. Therefore, standarization
# is necessary if filter is not the zcore and also neagaitive(-) sign is used
# to indacate left-tailed test, because all other methods are based on
# right-tailed test
#===============================================================================
data_analysis <- function(alpha, pvalue, filter, N_current, N_prior, tail,
                          max.group, standarized = FALSE,
                          effectType = c("continuous", "binary"))
{
    m = length(pvalue)
    sigma <- sqrt(N_current/N_prior)

    if(standarized == TRUE){
        mn = mean(filter, na.rm = TRUE)
        sdv = sd(filter, na.rm = TRUE)
        z_prior = -(filter - mn)/sdv
    } else {
        z_prior = -filter
    }

    mu <- sqrt(N_current/N_prior)*z_prior
    dbn_wgt <- bayes_weights(mu = mu, sigma = sigma, q = alpha/m)$w

    # FWER cotrols-------
    pro_bon <- empOPW(pvalue = pvalue, filter = filter, alpha = alpha,
                      tail = tail, max.group = max.group,
                      effectType = effectType, method = "BON")$rejections
    bon <- sum(pvalue <= alpha/length(pvalue), na.rm = TRUE)
    dbn_bon <- sum(pvalue <= alpha*dbn_wgt/m, na.rm = TRUE)
    ihw_bon <- rejections(ihw(pvalue, filter, alpha = alpha,
                              adjustment_type = "bonferroni"))


    # FDR controls-------
    pro_bh <- empOPW(pvalue = pvalue, filter = filter, alpha = alpha,
                     tail = tail, max.group = max.group,
                     effectType = effectType, method = "BH")$rejections
    bh <- sum(p.adjust(pvalue, method = "BH") <= alpha, na.rm = TRUE)
    dbn_bh <- sum(p.adjust(pvalue/dbn_wgt, method = "BH") <= alpha, na.rm = TRUE)
    ihw_bh <- rejections(ihw(pvalue, filter, alpha = alpha))

    # results-----
    rejections <- c(pro_bon, bon, dbn_bon, ihw_bon, pro_bh, bh, dbn_bh, ihw_bh)

    return(rejections)
}










