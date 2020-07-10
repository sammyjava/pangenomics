library(factoextra)

## PCA on path FR vectors, oddly in columns so transpose pathfrs
pca.pathfrs = prcomp(t(pathfrs), center=FALSE)

## get the results for variables (nodes)
## res.var$coord          # Coordinates
## res.var$contrib        # Contributions to the PCs
## res.var$cos2           # Quality of representation 
pca.pathfrs.var = get_pca_var(pca.pathfrs)

## get the results for individuals
## res.ind$coord          # Coordinates
## res.ind$contrib        # Contributions to the PCs
## res.ind$cos2           # Quality of representation
pca.pathfrs.ind = get_pca_ind(pca.pathfrs)

