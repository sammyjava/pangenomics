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
## loadings = data.frame(pca.pathfrs$rotation, .names=row.names(pca.pathfrs$rotation))
## p + geom_text(data=loadings, mapping=aes(x=PC1, y=PC2, label=.names, colour=.names)) +
##     coord_fixed(ratio=1) +
##     labs(x="PC1", y="PC2")

## NOT FANCY
num = length(colnames(pca.pathfrs$rotation))

## case/control labels and colors
labels = array(dim=length(rownames(pca.pathfrs$rotation)))
colors = array(dim=length(rownames(pca.pathfrs$rotation)))
cases = array(dim=length(rownames(pca.pathfrs$rotation)))
ctrls = array(dim=length(rownames(pca.pathfrs$rotation)))
for (i in 1:length(rownames(pca.pathfrs$rotation))) {
    casestart = str_locate_all(pattern='case', rownames(pca.pathfrs$rotation)[i])[[1]][1]
    ctrlstart = str_locate_all(pattern='ctrl', rownames(pca.pathfrs$rotation)[i])[[1]][1]
    if (!is.na(casestart)) {
        colors[i] = "darkred"
        labels[i] = substring(rownames(pca.pathfrs$rotation)[i], 1, casestart-2)
        cases[i] = TRUE
        ctrls[i] = FALSE
    } else if (!is.na(ctrlstart)) {
        colors[i] = "darkblue"
        labels[i] = substring(rownames(pca.pathfrs$rotation)[i], 1, ctrlstart-2)
        cases[i] = FALSE
        ctrls[i] = TRUE
    }
}



## these colors add to gray
ctrlcol = rgb(173,216,230, max = 255, alpha = 80, names = "lt.blue")
casecol = rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")

## the histograms -- first 10
for (pc in 1:10) {
    xlabel = paste("PC",pc," ",round(summary(pca.pathfrs)$importance["Proportion of Variance",pc]*100,1),"% of variance", sep="")
    h = hist(pca.pathfrs$rotation[,pc], breaks=33, plot=FALSE)
    h.cases = hist(pca.pathfrs$rotation[cases,pc], breaks=h$breaks, plot=FALSE)
    h.ctrls = hist(pca.pathfrs$rotation[ctrls,pc], breaks=h$breaks, plot=FALSE)
    ymax = max(h.cases$density, h.ctrls$density)
    plot(h.cases, freq=FALSE, col=casecol, xlab=paste("PC",pc), ylim=c(0,ymax), main=prefix)
    plot(h.ctrls, freq=FALSE, col=ctrlcol, ylim=c(0,ymax), add=TRUE)
    legend(x="topleft", c(paste(caseNum,"cases"),paste(ctrlNum,"controls"),paste(caseNum+ctrlNum,"both")), fill=c(casecol,ctrlcol,"lightgray"), bty="n")
}
