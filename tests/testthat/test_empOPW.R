empOPW <- function(pvalue, filter, weight = NULL, ranksProb = NULL, mean_testEffect = NULL,
                alpha = .05, tail = 1L, delInterval = .0001, max.group = 5L, h_breaks = 71L,
                effectType = c("continuous", "binary"), method = c("BH", "BON"), ... )
{
    # formulate a data set-------------
    Data = tibble(pvalue, filter)
    data_omit_na <- Data[which(!is.na(Data$pvalue)),]
    OD <- data_omit_na[order(data_omit_na$filter, decreasing = TRUE), ]
    OD_pvalue <- OD$pvalue


    # compute the number of tests------------
    m = length(OD_pvalue)
    nullProp = qvalue(p = OD_pvalue, pi0.method = "bootstrap")$pi0
    m0 = ceiling(nullProp*m)
    m1 = m - m0

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

        # find the optimal number of groups and degrees of freedom----------
        grp_seq <- seq(5, max.group, 5)
        op_grp_df <- sapply(grp_seq, optimal_group_df, pvalue = OD$pvalue,
                        filter = OD$filter, h_breaks = h_breaks, m = m, m1 = m1,
                        alpha = alpha, mean_testEffect = mean_testEffect,
                        effectType = effectType, method = method)

        grp <- op_grp_df[1, which.max(op_grp_df[3,])]
        df <- op_grp_df[2, which.max(op_grp_df[3,])]


        message("computing ranks probabilities")
        # compute the ranks probability of the tests given the mean effect
        ranksProb <- prob_rank_givenEffect_emp(pvalue = pvalue, filter = filter,
                                    group = grp, h_breaks = h_breaks, df = df,
                                    effectType = effectType)
        message("finished computing the ranks probabilities")


        message("computing weights")
        if(effectType == "continuous"){
            wgt = weight_continuous(alpha = alpha, et = mean_testEffect,
                                    m = grp, ranksProb = ranksProb)
        } else {
            wgt = weight_binary(alpha = alpha, et = mean_testEffect, m = group,
                                m1 = m1/m*grp, ranksProb = ranksProb)
        }
        message("finished computing the weights")
    }

    grpSize <- ceiling(m/grp)
    wgt_all = rep(wgt, each = grpSize)[1:m]

    message("comparing pvalues with thresholds")
    if(method == "BH"){
        padj <- p.adjust(OD_pvalue/wgt_all, method = "BH")
        rejections_list = OD[which((padj <= alpha) == TRUE), ]
    } else {
        rejections_list = OD[which((OD_pvalue <= alpha*wgt_all/m) == TRUE), ]
    }


    # outputs--------------
    n_rejections = dim(rejections_list)[1]

    return(list(totalTests = length(pvalue), nullProp = nullProp,
                ranksProb = ranksProb, weight = wgt_all,
                rejections = n_rejections, rejections_list = rejections_list))
}





