optimal_group_df <- function(group = 5L, pvalue, filter, h_breaks = 71L, m, m1,
           alpha = .05, mean_testEffect, effectType = c("continuous", "binary"),
           method = c("BH", "BON"))
{
    optimal_df <- function(df)
    {
        # ranks probability--------------
        ranksProb <- prob_rank_givenEffect_emp(pvalue = pvalue, filter = filter,
                        group = group, h_breaks = h_breaks, df = df,
                        effectType = effectType)

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

        return(c(df = df, n_rej = n_rejections))
    }

   message(paste("computing for group", group))
   # not necessary to use all df
   df_and_rej <- sapply(seq(2, group, 3), optimal_df)

   op_df_rej <- df_and_rej[ , which.max(df_and_rej[2,])]

   return(c(group = group, op_df_rej))
}



