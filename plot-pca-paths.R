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
num = length(colnames(pca.paths$rotation))

## labels
##                  1         2            3           4           5
##        pheno_241.2 pheno_555 pheno_695.42 pheno_714.1 pheno_250.1
## 476296        ctrl      ctrl         case        ctrl        ctrl
## 261580        ctrl      ctrl         ctrl        unkn        ctrl

## assign colors to phenotypes
phenotypes = colnames(labels)
colors = brewer.pal(length(phenotypes), "Dark2")

## assign phenotype colors to paths; or default gray
phenotypeColor = array(dim=nrow(paths))
transparent = rgb(125, 125, 125, max=255, alpha=32)
for (i in 1:nrow(paths)) {
    phenotypeColor[i] = transparent
    for (j in 1:length(phenotypes)) {
        if (labels[i,j]=="case") {
            phenotypeColor[i] = colors[j]
        }
    }
}

## plot the first 10 PCA plots
for (i in 1:10) {
    xlabel = paste("PC",i,  " ",round(summary(pca.paths)$importance["Proportion of Variance",i]*100,1),"% of variance", sep="")
    ylabel = paste("PC",i+1," ",round(summary(pca.paths)$importance["Proportion of Variance",i+1]*100,1),"% of variance", sep="")
    plot(pca.paths$rotation[,i], pca.paths$rotation[,i+1], xlab=xlabel, ylab=ylabel, pch=19, cex=0.3, col=phenotypeColor, main=prefix)
    legend(x="topleft", phenotypes, pch=19, col=colors, bty="n", cex=1.0)
}
