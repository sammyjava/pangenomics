##
## t-SNE on paths data
##

require(stringr)
require(tsne)

colors = c("red","blue")
names(colors) = c("case", "ctrl")

ecb = function(x,y){ plot(x, col=colors[labels], pch=19, cex=0.5 ) } # text(x, labels=pathnames, col=colors[labels], cex=0.75)

tsne.paths = tsne(paths, initial_config=NULL, k=2, min_cost=0, whiten=FALSE, initial_dims=30, perplexity=30, epoch_callback=ecb, max_iter=1000, epoch=10)
