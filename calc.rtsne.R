require(Rtsne)
## Rapid t-SNE on paths data
##
##                   X: matrix; Data matrix (each row is an observation, each column is a variable)
##                dims: integer; Output dimensionality (default: 2)
##        initial_dims: integer; the number of dimensions that should be retained in the initial PCA step (default: 50)
##          perplexity: numeric; Perplexity parameter (should not be bigger than 3 perplexity < nrow(X) - 1, see details for interpretation)
##               theta: numeric; Speed/accuracy trade-off (increase for less accuracy), set to 0.0 for exact TSNE (default: 0.5)
##    check_duplicates: logical; Checks whether duplicates are present. It is best to make sure there are no duplicates present and set this option to FALSE, especially for large datasets (default: TRUE)
##                 pca: logical; Whether an initial PCA step should be performed (default: TRUE)
##         partial_pca: logical; Whether truncated PCA should be used to calculate principal components (requires the irlba package). This is faster for large input matrices (default: FALSE)
##            max_iter: integer; Number of iterations (default: 1000)
##             verbose: logical; Whether progress updates should be printed (default: global "verbose" option, or FALSE if that is not set)
##         is_distance: logical; Indicate whether X is a distance matrix (experimental, default: FALSE)
##              Y_init: matrix; Initial locations of the objects. If NULL, random initialization will be used (default: NULL).
##                      Note that when using this, the initial stage with exaggerated perplexity values and a larger momentum term will be skipped.
##          pca_center: logical; Should data be centered before pca is applied? (default: TRUE)
##           pca_scale: logical; Should data be scaled before pca is applied? (default: FALSE)
##           normalize: logical; Should data be normalized internally prior to distance calculations with ‘normalize_input’? (default: TRUE)
##     stop_lying_iter: integer; Iteration after which the perplexities are no longer exaggerated (default: 250, except when Y_init is used, then 0)
##     mom_switch_iter: integer; Iteration after which the final momentum is used (default: 250, except when Y_init is used, then 0)
##            momentum: numeric; Momentum used in the first part of the optimization (default: 0.5)
##      final_momentum: numeric; Momentum used in the final part of the optimization (default: 0.8)
##                 eta: numeric; Learning rate (default: 200.0)
## exaggeration_factor: numeric; Exaggeration factor used to multiply the P matrix in the first part of the optimization (default: 12.0)
##         num_threads: integer; Number of threads to use using OpenMP, default 1. 0 corresponds to using all available cores
##               index: integer matrix; Each row contains the identity of the nearest neighbors for each observation
##            distance: numeric matrix; Each row contains the distance to the nearest neighbors in ‘index’ for each observation
calc.rtsne = function(uniquepathmatrix, prevResult=NULL) {
    ## initial values
    Y_init = NULL
    pca = FALSE
    partial_pca = TRUE
    theta = 0.5
    if (!is.null(prevResult)) {
        ## refinement values
        Y_init = prevResult$Y
        pca = FALSE
        partial_pca = FALSE
        theta = prevResult$theta/2
    }
    ## run it
    return(Rtsne(uniquepathmatrix, Y_init=Y_init,
                 num_threads=0,
                 verbose=TRUE,
                 partial_pca=partial_pca,
                 pca=pca,
                 theta=theta,
                 pca_scale=TRUE, 
                 normalize=TRUE,
                 perplexity=30,
                 momentum=0.5,
                 final_momentum=0.8,
                 eta=200.0,
                 exaggeration_factor=12.0,
                 max_iter=1000))
}
