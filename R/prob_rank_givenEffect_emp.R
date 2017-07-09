#' @title Emperical ranks probbaility of the test given the effect size
#'
#' @description Emperical comnputation of the ranks probability of a test being
#' higher than any other test given the effect size from the external information.
#' @param pvalue vector of test pvalues
#' @param filter vector of filter statistics
#' @param group number of groups, should be at least four
#' @param h_breaks number of breaks for the histogram
#' @param bin_idx Integer, bin number of the histogram. Almost alwyas it is one,
#' becasue we are interested to obtain the ranks probability of the alternative
#' tests
#' @param smooth Character of ("TRUE" or "FALSE") to apply spline smoothing,
#' default is TRUE
#' @param effectType type of effect size c("binary","continuous")
#'
#' @details If one wants to test \deqn{H_0: \epsilion_i=0 vs. H_a: \epsilion_i > 0,}
#' then \code{et}  and \code{ey} should effect type continuous effect size and
#' if one wants to test \deqn{H_0: \epsilion_i=0 vs. H_a: \epsilion_i = \epsilion,}
#' one should the binary effect size
#'
#' @author Mohamad S. Hasan, mshasan@uga.edu
#'
#' @export
#'
#' @import stats
#' @import limma propTrueNull
#'
#' @return \code{ranksProb} emperical probability of the rank of the test
#'
#' @examples
#'
#' # generating data (known in practice)
#' set.seed(123)
#' X = runif(10000, min = 0, max = 2.5)         # covariate
#' H = rbinom(length(X), size = 1, prob = 0.1)   # hypothesis true or false
#' Z = rnorm(length(X), mean = H * X)            # Z-score
#' p = 1 - pnorm(Z)
#'
#' # apply the function to compute the rank proabbility
#' grp = 10
#' ranksProb = prob_rank_givenEffect_emp(pvalue = p, filter = X, group = grp,
#'                                h_breaks = 71, effectType = "continuous")
#'
#' # plot the probability
#' plot(1:grp, ranksProb, type="l", xlab = "ranks", ylab = "P(rank | effect)")
#'
#===============================================================================
# function to compute p(rank=k|filterEffect=ey) emperically

# internal parameters:-----
# grpSize = number of pvalues per group
# Data = a data frame of pvalue and filter statistics
# OD = ordered data by the fitler statsistics
# OD_pvalue = ordered pvaluse by the fitler statistics
# pval_perGrp = vector of pvalues per group
# fun_prob = funtion to compute probability for each group
# prob = proability for each group
# hist_dens = compute density from histogram
# probAll = normalized density for all points but we need only the first

#===============================================================================
prob_rank_givenEffect_emp <- function(pvalue, filter, group = 5L, h_breaks = 100L,
                                      bin_idx = 1L, smooth = TRUE,
                                      effectType = c("continuous", "binary"))
    {
        Data = tibble(pvalue, filter)
        data_omit_na <- Data[which(!is.na(Data$pvalue)),]
        OD <- data_omit_na[order(data_omit_na$filter, decreasing = TRUE), ]
        OD_pvalue <- OD$pvalue

        grpSize <- ceiling(length(OD_pvalue)/group)

        # function to compute ranks probbaility per group--------------
        fun_prob <- function(grp)
            {
                pval_perGrp <- OD_pvalue[(grp*grpSize - grpSize + 1):(grp*grpSize)]

                if(effectType == "continuous"){

                    prob <- relative_freq(bin_idx = bin_idx, h_breaks = h_breaks,
                                  obs = pval_perGrp)$rf_one

                } else {

                    prob <- 1 - propTrueNull(p = pval_perGrp, method = "hist")
                }

                return(prob)
            }

        probVec <- sapply(1:group, fun_prob)

        # smooting and nomalizing the ranks probability-------------
        if(smooth == TRUE){

            probVec_smooth <- smooth.spline(x = 1:group, y = probVec, cv = FALSE)$y
            if(any(probVec_smooth < 0)){
                neg_val <- min(probVec_smooth[probVec_smooth < 0])
                probVec_smooth <- probVec_smooth - neg_val
            }

        } else {
            probVec_smooth <- probVec
        }

        probVec_smooth_norm <- probVec_smooth/sum(probVec_smooth, na.rm = TRUE)

        return(probVec_smooth_norm)
    }
