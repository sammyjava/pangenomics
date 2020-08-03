##
## Calculate path-based PCA
##
library(ggfortify)
library(factoextra)

pca.paths = prcomp(t(paths), center=TRUE)

## get the results for variables (nodes)
## res.var$coord          # Coordinates
## res.var$contrib        # Contributions to the PCs
## res.var$cos2           # Quality of representation 
pca.paths.var = get_pca_var(pca.paths)

## get the results for individuals
## res.ind$coord          # Coordinates
## res.ind$contrib        # Contributions to the PCs
## res.ind$cos2           # Quality of representation
pca.paths.ind = get_pca_ind(pca.paths)
