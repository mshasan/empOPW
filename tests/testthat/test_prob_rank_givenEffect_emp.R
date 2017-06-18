prob_rank_givenEffect_emp <- function(pvalue, filter, group = 5L, h_breaks = 71L,
                                      df = 3, effectType = c("continuous", "binary"))
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
                    hist_dens <- hist(pval_perGrp, freq = FALSE,
                              breaks = seq(0, 1, length = h_breaks))$density
                    probAll = hist_dens/sum(hist_dens)
                    prob = probAll[1]
                } else {
                    prob <- 1 - propTrueNull(p = pval_perGrp, method = "lfdr")
                }

                return(prob)
            }

        probVec <- sapply(1:group, fun_prob)

        # smooting and nomalizing the ranks probability-------------
        probVec_smooth <- smooth.spline(x = 1:group, y = probVec, df = df)$y
        # if(any(probVec_smooth < 0)){
        #     neg_val <- probVec_smooth[probVec_smooth < 0]
        #     probVec_smooth <- probVec_smooth - neg_val + .000001
        # }
        probVec_smooth_norm <- probVec_smooth/sum(probVec_smooth, na.rm = TRUE)

        return(probVec_smooth_norm)
    }









