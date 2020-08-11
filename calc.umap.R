library(umap)

## wrapper for umap
calc.umap = function(pathfrs, n_neighbors=15, n_components=2, n_epochs=200, min_dist=0.1, alpha=1.0, gamma=1.0) {

    ## umap configuration parameters
    ##            n_neighbors: 15
    ##           n_components: 2
    ##                 metric: euclidean
    ##               n_epochs: 200
    ##                  input: data
    ##                   init: spectral
    ##               min_dist: 0.1
    ##       set_op_mix_ratio: 1
    ##     local_connectivity: 1
    ##              bandwidth: 1
    ##                  alpha: 1
    ##                  gamma: 1
    ##   negative_sample_rate: 5
    ##                      a: 1.57694361269457
    ##                      b: 0.895060718151928
    ##                 spread: 1
    ##           random_state: 282099559
    ##        transform_state: NA
    ##            knn_repeats: 1
    ##                verbose: FALSE
    ##        umap_learn_args: NA
    ##                 method: naive
    ##        metric.function: [function]

    config = umap.defaults
    config$n_neighbors = n_neighbors
    config$n_components = n_components
    config$n_epochs = n_epochs
    config$min_dist = min_dist
    config$alpha = alpha
    config$gamma = gamma
    
    return(umap(pathfrs, config=config))
}
