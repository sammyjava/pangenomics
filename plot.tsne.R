##
## plot tSNE results
##
plot.tsne = function(tsne, showLabels=FALSE, xmin=0, xmax=0, ymin=0, ymax=0) {
    if (xmin!=0 || xmax!=0) {
        xlim = c(xmin,xmax)
    }
    if (ymin!=0 || ymax!=0) {
        ylim = c(ymin,ymax)
    }
    if (xmin!=0 || xmax!=0 || ymin!=0 || ymax!=0) {
        plot(tsne, main=prefix, xlab="t-SNE dim 1", ylab="t-SNE dim 2", col=colors[labels], xlim=xlim, ylim=ylim, pch=1, cex=0.5);
    } else {
        plot(tsne, main=prefix, xlab="t-SNE dim 1", ylab="t-SNE dim 2", col=colors[labels], pch=1, cex=0.5);
    }
    if (showLabels) {
        text(tsne, labels=pathnames, pos=4, col=colors[labels], cex=0.8)
    }
}
