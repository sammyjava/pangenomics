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
num = length(colnames(pca.pathfrs$rotation))

## casecol = "red"
## ctrlcol = "blue"

## ctrlcol = rgb(173,216,230, max = 255, alpha = 80, names = "lt.blue")
## casecol = rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")

## ctrlcol = rgb(173,216,230, max=255)
## casecol = rgb(255,192,203, max=255)

casecol = rgb(1,0,0, maxColorValue=1, alpha = 0.15, names = "red")
ctrlcol = rgb(0,0,1, maxColorValue=1, alpha = 0.15, names = "blue")

labels = array(length(rownames(pathfrs)))
colors = array(dim=length(rownames(pca.pathfrs$rotation)))
parts = strsplit(rownames(pathfrs), ".", fixed=T)
for (i in 1:length(rownames(pathfrs))) {
    labels[i] = parts[[i]][2]
    if (labels[i]=="case") {
        colors[i] = casecol
    } else if (labels[i]=="ctrl") {
        colors[i] = ctrlcol
    } else {
        colors[i] = "black"
    }
}

## plot the first 10 PCA plots
for (i in 1:10) {
    xlabel = paste("PC",i,  " ",round(summary(pca.pathfrs)$importance["Proportion of Variance",i]*100,1),"% of variance", sep="")
    ylabel = paste("PC",i+1," ",round(summary(pca.pathfrs)$importance["Proportion of Variance",i+1]*100,1),"% of variance", sep="")
    plot(pca.pathfrs$rotation[,i], pca.pathfrs$rotation[,i+1], xlab=xlabel, ylab=ylabel, pch=19, cex=1.0, col=colors, main=prefix)
    legend(x="topleft", c(paste(caseNum,"cases"),paste(ctrlNum,"controls")), pch=19, col=c(casecol,ctrlcol), bty="n", cex=1.0)
}
