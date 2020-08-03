##
## t-SNE on pathfrs data
##

require(stringr)
require(tsne)

colors = c("red","blue")

ecb = function(x,y){ plot(x, col=colors[labels], pch=19, cex=0.5 ) } # text(x, labels=pathnames, col=colors[labels], cex=0.75)

tsne.pathfrs = tsne(pathfrs, initial_config=NULL, k=2, min_cost=0, whiten=FALSE, initial_dims=30, perplexity=30, epoch_callback=ecb, max_iter=1000, epoch=10)
