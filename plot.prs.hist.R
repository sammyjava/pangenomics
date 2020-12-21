################################################################################################
## plot a histogram of the polygenic risk scores of the individuals
##
## input: dataframe of "sample", "label", "score"
################################################################################################

plot.prs.hist = function(prs, title="GWAS PRS") {

    ## just in case
    colnames(prs) = c("sample", "label", "score")
    
    col.cases = rgb(1,0,0,1/2)
    col.controls = rgb(0,0,1,1/2)
    
    nCases = nrow(prs[prs$label=="case",])
    nControls = nrow(prs[prs$label=="ctrl",])
    
    h.cases = hist(prs$score[prs$label=="case"], breaks=20, plot=FALSE)
    h.controls = hist(prs$score[prs$label=="ctrl"], breaks=20, plot=FALSE)

    ## plot density to accomodate different numbers of cases and controls
    ymax = max(h.cases$density, h.controls$density)
    plot(h.cases, freq=FALSE, col=col.cases, xlab="log10(PRS)", ylim=c(0,ymax), main=title)
    plot(h.controls, freq=FALSE, col=col.controls, add=TRUE)
    lines(c(0,0), c(0,1000), col="red")
    legend(x="topleft",
           c(paste(nCases,"cases"),paste(nControls,"controls"),paste(nCases+nControls,"both")),
           fill=c(col.cases,col.controls,rgb(.5,0,.5)), bty="n")
}
