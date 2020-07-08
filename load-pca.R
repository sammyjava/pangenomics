library(factoextra)

source("load-pathfrs.R")

## PCA on path FR vectors
pca = prcomp(pathfrs, center=FALSE)
## get the results for variables (nodes)
## res.var$coord          # Coordinates
## res.var$contrib        # Contributions to the PCs
## res.var$cos2           # Quality of representation 
pca.var = get_pca_var(pca)
## get the results for individuals
## res.ind$coord          # Coordinates
## res.ind$contrib        # Contributions to the PCs
## res.ind$cos2           # Quality of representation
pca.ind = get_pca_ind(pca)

## determine which are cases and which are controls
cases = endsWith(rownames(FRs), "case")
cases.0 = endsWith(rownames(FRs), "0.case")
cases.1 = endsWith(rownames(FRs), "1.case")
controls = endsWith(rownames(FRs), "ctrl")
controls.0 = endsWith(rownames(FRs), "0.ctrl")
controls.1 = endsWith(rownames(FRs), "1.ctrl")

