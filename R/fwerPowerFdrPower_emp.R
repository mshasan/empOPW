#' @title Simulate FWER, POWER, FDR, and POWER
#'
#' @description This function simulate Family Wise Error Rate (FWER) and
#' corresponding Power, and False Discovery Rate (FDR) and the corresponding
#' Power for different effect sizes
#'
#' @param i i-th filter effect
#' @param simu number of replications
#' @param m number of hypothesis test
#' @param null proportion of the true null hypothesis
#' @param corr correlation between test statistics
#' @param cv determine whether the test mean effect and the filter mean effects
#' are the same
#' @param alpha significance threshold
#' @param groupSize number of test statistics per group
#' @param max.group maximum number of p-value groups to be used, minimum is 5
#' @param filterEffectVec a vector of different effect size
#' @param effectType type of effect sizes, c("continuous", "binary")
#'
#' @details
#' This function empirically simulate Family Wise Error Rate (FWER) and
#' corresponding Power, and False Discovery Rate (FDR) and the corresponding
#' Power for the different effect sizes
#'
#' @author Mohamad S. Hasan, \email{mshasan@uga.edu}
#'
#' @export
#'
#' @import pweight bayes_weights
#'
#' @seealso \code{\link{fwerPowerFdrPower}}
#'
#' @return a matrix of 16 rows containing information about FWER, POWER, FDR,
#' and POWER (4 rows for each)
#'
#' @references Hasan and Schliekelman (2017)
#'
#' @examples
#' # vector of effect sizes
#' effectVec <- c(1, 1.5, 2)
#' simuVal = 3
#' FwerPowerFdrPower <- sapply(1:length(effectVec), fwerPowerFdrPower_emp,
#'              simu = simuVal, m = 10000, null = .5, corr = 0,
#'              cv = 0, alpha = .05, groupSize = 100, max.group = 10,
#'              filterEffectVec = effectVec, effectType = "continuous")
#'
#===============================================================================
#----------------------fwerPowerFdrPower----------------------------
# function to compute Simulated FWER, POWER, and FDR by effect size
#
# internal parameters:-----
# ey = filter effect size
# m0 = no. of null hypothesis
# m1 = no. of alt hyp.
# xf = only alt. filter effect vector
# xt = alt. test effect vector
# Sigma = test correlation matrix
# H = alternative hypothesis true or false
# ef = filter effect vector (mixture of null and alt)
# et = test effect vector (mixture of null and alt)
# mGrp = subgroup of tests.
# test = filter test stat
# filter = actual test stat
# pval = filter test pvalues
# pro = proposed, bon = bonferroni, rdw = roeder and wasserman,
# IHW =  independent Hyp. Weight
#
#===============================================================================

fwerPowerFdrPower_emp <- function(i, simu, m, null, corr = 0, cv = 0, alpha = .05,
                        groupSize = 100L, max.group = 10L, filterEffectVec,
                        effectType = c("continuous", "binary"))
{
    # initial parameters---------
    ey <- filterEffectVec[i]
    m0 <- ceiling(m*null)
    m1 <- m - m0

    if(effectType == "continuous"){
        xf <- as.vector(runif_by_mean(n = m, mean = ey))
    } else {
        xf <- rep(ey, m)
    }

    xt <- if(cv == 0){xf} else {rnorm(m, ey, cv*ey)}
    Sigma <- matrix(corr, groupSize, groupSize) + diag(groupSize)*(1 - corr)


    fwerPowerFdrPower_simu <- function(s)
    {
        # generate synthetic data-------------
        H <- rbinom(m, 1, 1 - null)
        ef <- H*xf
        et <- H*xt
        mGrp = m/groupSize

        test <- if(corr == 0) {rnorm(m, et, 1)
        } else {as.vector(sapply(1:mGrp, test_by_block, eVec = et,
                                 groupSize = groupSize, Sigma = Sigma))}

        filter <- if(corr == 0) {rnorm(m, ef, 1)
        } else {as.vector(sapply(1:mGrp, test_by_block, eVec = ef,
                                 groupSize = groupSize, Sigma = Sigma))}

        pval <- pnorm(test, lower.tail = FALSE)

        dat = tibble(test, pval, et, filter)
        OD = dat[order(dat$filter, decreasing = TRUE), ]
        # end of synthetic data generation------------------


        # applying different methods----------------------
        pro_fwer<- empOPW(pvalue = dat$pval, filter = dat$filter, max.group = max.group,
                          alpha = alpha, effectType = effectType, method = "BON")
        pro_fdr <- empOPW(pvalue = dat$pval, filter = dat$filter, max.group = max.group,
                          alpha = alpha, effectType = effectType, method = "BH")

        weight_rdw <- roeder_wasserman_weight(pvalue = OD$pval, alpha = alpha)
        dbn_wgt <- bayes_weights(mu = OD$filter, sigma = rep(1, m), q = alpha/m)$w

        ihw_fwer <- ihw(OD$pval, OD$filter, alpha = alpha, adjustment_type = "bonferroni")
        ihw_fdr <-  ihw(OD$pval, OD$filter, alpha = alpha, adjustment_type = "BH")

        rej_pro <- OD$pval%in%pro_fwer$rejections_list$pvalue
        rej_bon <- OD$pval <= alpha/m
        rej_rdw <- OD$pval <= alpha*weight_rdw/m
        rej_dbn <- OD$pval <= alpha*dbn_wgt/m
        rej_ihwFwer <- adj_pvalues(ihw_fwer) <= alpha

        n_null <- max(1, sum(OD$et == 0, na.rm = TRUE))
        n_alt <-  max(1, sum(OD$et != 0, na.rm = TRUE))

        FWER_pro <- sum(rej_pro[OD$et == 0])
        FWER_bon <- sum(rej_bon[OD$et == 0])
        FWER_rdw <- sum(rej_rdw[OD$et == 0])
        FWER_dbn <- sum(rej_dbn[OD$et == 0])
        FWER_ihw <- sum(rej_ihwFwer[OD$et == 0])

        POWER_pro <- sum(rej_pro[OD$et != 0])/n_alt
        POWER_bon <- sum(rej_bon[OD$et != 0])/n_alt
        POWER_rdw <- sum(rej_rdw[OD$et != 0])/n_alt
        POWER_dbn <- sum(rej_dbn[OD$et != 0])/n_alt
        POWER_ihw <- sum(rej_ihwFwer[OD$et != 0])/n_alt


        # fdr based information----------------
        rej_pro_fdr <- OD$pval%in%pro_fdr$rejections_list$pvalue

        adjPval_bon <- p.adjust(OD$pval, method="BH")
        adjPval_rdw <- p.adjust(OD$pval/weight_rdw, method="BH")
        adjPval_dbn <- p.adjust(OD$pval/dbn_wgt, method="BH")
        adjPval_ihw <- adj_pvalues(ihw_fdr)

        FDR_pro <- sum(rej_pro_fdr[OD$et == 0])/max(1, sum(rej_pro_fdr))
        FDR_bh  <- sum(adjPval_bon[OD$et == 0] <= alpha)/max(1, sum(adjPval_bon <= alpha))
        FDR_rdw <- sum(adjPval_rdw[OD$et == 0] <= alpha)/max(1, sum(adjPval_rdw <= alpha))
        FDR_dbn <- sum(adjPval_dbn[OD$et == 0] <= alpha)/max(1, sum(adjPval_dbn <= alpha))
        FDR_ihw <- sum(adjPval_ihw[OD$et == 0] <= alpha)/max(1, rejections(ihw_fdr))

        FDR_POWER_pro <- sum(rej_pro_fdr[OD$et != 0])/n_alt
        FDR_POWER_bh  <- sum(adjPval_bon[OD$et != 0] <= alpha)/n_alt
        FDR_POWER_rdw <- sum(adjPval_rdw[OD$et != 0] <= alpha)/n_alt
        FDR_POWER_dbn <- sum(adjPval_dbn[OD$et != 0] <= alpha)/n_alt
        FDR_POWER_ihw <- sum(adjPval_ihw[OD$et != 0] <= alpha)/n_alt

        return(c(FWER_pro, FWER_bon, FWER_rdw, FWER_dbn, FWER_ihw,
                 POWER_pro, POWER_bon, POWER_rdw, POWER_dbn, POWER_ihw,
                 FDR_pro, FDR_bh, FDR_rdw, FDR_dbn, FDR_ihw,
                 FDR_POWER_pro, FDR_POWER_bh, FDR_POWER_rdw, FDR_POWER_dbn, FDR_POWER_ihw))
    }

    fwerPowerFdrPower_bysimu <- sapply(1:simu, fwerPowerFdrPower_simu)
    fwerPowerFdrPower <- rowMeans(fwerPowerFdrPower_bysimu, na.rm = TRUE)

    return(fwerPowerFdrPower)

}





