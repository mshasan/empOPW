#' @title Emperical probbaility of rank of test given effect size
#'
#' @description Emperical comnputation of the probbaility of rank of a test being
#' higher than any other test given the effect size from external information.
#' @param group number of groups
#' @param pvalue vector of test pvalues
#' @param filterStat vector of filter statistics
#' @param effectType type of effect size c("binary","continuous")
#'
#' @details If one wants to test \deqn{H_0: \epsilion_i=0 vs. H_a: \epsilion_i > 0,}
#' then \code{et}  and \code{ey} should be mean of the test and filter effect sizes,
#' respectively. This is called hypothesis testing for the continuous effect sizes.
#' If one wants to test \deqn{H_0: \epsilion_i=0 vs. H_a: \epsilion_i = \epsilion,}
#' then \code{et} and \code{ey} should be median or any discrete value of the
#' test and filter effect sizes. This is called hypothesis testing for the Binary
#' effect sizes
#' @author Mohamad S. Hasan, mshasan@uga.edu
#' @export
#' @import stats
#' @seealso \code{\link{qvalue}} \code{\link{runif}} \code{\link{rnorm}}
#' @return \code{prob} emperical probability of the rank of the test
#' @examples
#'
#' # generating data (known in practice)
#' pvalue <- runif(100000)
#' filterStat <- rnorm(100000)
#'
#' # apply the function to compute the rank proabbility
#' group=10
#' ranksProb=sapply(group,prob_rank_givenEffect_emp, pvalue, filterStat,
#'                        effectType="continuous")
#'
#' # plot the probability
#' plot(1:group,ranksProb,type="l",lwd=2,xlab="ranks",ylab="P(rank|effect)")
#'
#===============================================================================
# function to compute p(rank=k|filterEffect=ey) emperically

# Input:-----
# group = number of groups
# pvalue = vector of test pvalues
# filterStat = vector of filter statistics
# effectType = type of effect size c("binary","continuous")

# internal parameters:-----
# groupSize = number of pvalues per group
# Data = a data frame of pvalue and filter statistics
# OD = ordered data by the fitler statsistics
# OD_pvalue = ordered pvaluse by the fitler statistics
# pvalue_perGroup = vector of pvalues per group
# fun.prob = funtion to compute probability for each group
# prob = proability for each group
# h = compute density from histogram
# probAll = normalized density for all points but we need only the first

# output:-----
# probVec = normalized probability of rank given effect size, p(rank=k|effect=ey)
#===============================================================================
prob_rank_givenEffect_emp <- function(group = 5L, pvalue, filterStat,
                                      effectType = c("binary", "continuous"))
{
    groupSize <- ceiling(length(pvalue)/group)
    Data <- data.frame(pvalue, filterStat)
    OD <- Data[order(Data$filterStat, decreasing=T), ]
    OD_pvalue <- OD$pvalue

    fun.prob <- function(group)
    {
        pvalue_perGroup <- OD_pvalue[(group*groupSize-groupSize+1):(group*groupSize)]

        if(effectType=="binary"){
            prob <- 1-qvalue(p = pvalue_perGroup, pfdr = TRUE, pi0.method="bootstrap",
                             lambda = max(pvalue_perGroup))$pi0
        } else {
            h <- hist(pvalue_perGroup, freq = FALSE, breaks=seq(0,1,length=11))$density
            probAll = h/sum(h)
            prob = probAll[1]
        }
    return(prob)
    }
    probVec <- sapply(1:group, fun.prob)
    return(probVec/sum(probVec))
}
