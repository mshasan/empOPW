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







