source("load-libsvm.R")

col.graph = "blue"
col.fr = "darkgreen"
colors = c(col.graph, col.fr)

plot(graph.specificity, graph.sensitivity, 
     xlim=c(.40,.75), ylim=c(.40,.70),
     xlab="specificity (TNR)", ylab="sensitivity (TPR)",
     type="n"
     )

## lines first so they're under the points
for (tc in (1:19)/20) {
    x = (0:10)/10
    y = 2*tc - x
    lines(x, y, lty=2, col="gray")
    text(x, y, paste(tc*100,"%",sep=""), col="gray")
}

points(graph.specificity, graph.sensitivity, pch=19, col=col.graph)
points(fr.specificity, fr.sensitivity, pch=19, col=col.fr)

for (i in 1:length(graph)) {
    offset = (i-1)*2
    text(specificity[offset+1], sensitivity[offset+1], graph[i], pos=1, offset=1)
    arrows(specificity[offset+1], sensitivity[offset+1], specificity[offset+2], sensitivity[offset+2], length=0.1)
}


legend(x="topright", bty="n",
       c("graph path traversal", "FR path support"),
       pch=19, col=colors,
       )


## plot(totCorrect, MCC, col=colors, xlim=c(0,100), ylim=c(-1,1))
## for (i in 1:length(graph)) {
##     offset = (i-1)*3
##     text(totCorrect[offset+1], MCC[offset+1], graph[i], pos=1)
##     arrows(totCorrect[offset+1], MCC[offset+1], totCorrect[offset+2], MCC[offset+2], col="darkgreen", length=0.1)
##     arrows(totCorrect[offset+2], MCC[offset+2], totCorrect[offset+3], MCC[offset+3], col="blue", length=0.1)
## }
