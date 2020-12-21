library(fitdistrplus)

################################################################################################
## plot a histogram and normal distribution fits of the polygenic risk scores of the individuals
##
## input: dataframe of "sample", "label", "score"
##        label ("case" or "ctrl")
################################################################################################

plot.prs.fit = function(prs, label="both") {
    ## just in case
    colnames(prs) = c("sample", "label", "score")
    ## use method `plot.fitdist`
    if (label=="both") {
        fit = fitdist(prs$score, "norm")
    } else {
        fit = fitdist(prs$score[prs$label==label], "norm")
    }
    plot(fit)
}
