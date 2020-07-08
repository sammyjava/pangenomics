##
## t-SNE on FR data
##

require(stringr)
require(tsne)

colors = c("darkred","darkgreen")
names(colors) = c("case", "ctrl")

tpathfrs = t(pathfrs)

labels = array(length(rownames(tpathfrs)))
pathnames = array(length(labels))
for (i in 1:length(rownames(tpathfrs))) {
    casestart = str_locate_all(pattern='case', rownames(tpathfrs))[[i]][1]
    ctrlstart = str_locate_all(pattern='ctrl', rownames(tpathfrs))[[i]][1] 
    if (!is.na(casestart)) {
        labels[i] = "case"
        pathnames[i] = substring(rownames(tpathfrs)[i], 1, casestart-2)
    } else if (!is.na(ctrlstart)) {
        labels[i] = "ctrl"
        pathnames[i] = substring(rownames(tpathfrs)[i], 1, ctrlstart-2)
    }
}

ecb = function(x,y){ plot(x,t='n'); text(x, labels=pathnames, col=colors[labels], cex=0.75) }
tsne_frs = tsne(tpathfrs, epoch_callback=ecb, perplexity=50, initial_dims=30, max_iter=3000, epoch=100)
