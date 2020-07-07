##
## plot PCA loadings given a prcomp output object
##

require(stringr)
require(ggplot2)
require(RColorBrewer)

## FANCY
## theta = seq(0,2*pi,length.out = 100)
## circle = data.frame(x = cos(theta), y = sin(theta))
## p = ggplot(circle,aes(x,y)) + geom_path()
## loadings = data.frame(pca$rotation, .names=row.names(pca$rotation))
## p + geom_text(data=loadings, mapping=aes(x=PC1, y=PC2, label=.names, colour=.names)) +
##     coord_fixed(ratio=1) +
##     labs(x="PC1", y="PC2")

## NOT FANCY
num = length(colnames(pca$rotation))

## case/control labels and colors
labels = array(dim=length(rownames(pca$rotation)))
colors = array(dim=length(rownames(pca$rotation)))
cases = array(dim=length(rownames(pca$rotation)))
ctrls = array(dim=length(rownames(pca$rotation)))
for (i in 1:length(rownames(pca$rotation))) {
    casestart = str_locate_all(pattern='case', rownames(pca$rotation)[i])[[1]][1]
    ctrlstart = str_locate_all(pattern='ctrl', rownames(pca$rotation)[i])[[1]][1]
    if (!is.na(casestart)) {
        colors[i] = "darkred"
        labels[i] = substring(rownames(pca$rotation)[i], 1, casestart-2)
        cases[i] = TRUE
        ctrls[i] = FALSE
    } else if (!is.na(ctrlstart)) {
        colors[i] = "darkblue"
        labels[i] = substring(rownames(pca$rotation)[i], 1, ctrlstart-2)
        cases[i] = FALSE
        ctrls[i] = TRUE
    }
}



## these colors add to gray
ctrlcol = rgb(173,216,230, max = 255, alpha = 80, names = "lt.blue")
casecol = rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")

## the histograms
for (pc in 1:num) {
    xlabel = paste("PC",pc," ",round(summary(pca)$importance["Proportion of Variance",pc]*100,1),"% of variance", sep="")
    h.cases = hist(pca$rotation[cases,pc], breaks=33, plot=FALSE)
    h.ctrls = hist(pca$rotation[ctrls,pc], breaks=33, plot=FALSE)
    ymax = max(h.cases$density, h.ctrls$density)
    plot(h.cases, freq=FALSE, col=casecol, xlab=paste("PC",pc), ylim=c(0,ymax), main=prefix)
    plot(h.ctrls, freq=FALSE, col=ctrlcol, ylim=c(0,ymax), add=TRUE)
    ## legend(x="topleft", c(paste(nCases,"cases"),paste(nCtrls,"ctrls"),paste(nCases+nCtrls,"both")),
    ##        fill=c(col.cases,col.ctrls,rgb(.5,0,.5)), bty="n")
}
