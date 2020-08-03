source("load-libsvm.R")

colors = c("blue", "darkgreen")
        
plot(specificity, sensitivity, col=colors,
     xlim=c(.40,.75), ylim=c(.45,.70),
     pch=19,
     ## xlim=c(0,1), ylim=c(0,1),
     xlab="specificity (TNR)", ylab="sensitivity (TPR)"
     )

for (i in 1:length(graph)) {
    offset = (i-1)*2
    text(specificity[offset+1], sensitivity[offset+1], graph[i], pos=1)
    arrows(specificity[offset+1], sensitivity[offset+1], specificity[offset+2], sensitivity[offset+2], length=0.1)
    ## arrows(specificity[offset+2], sensitivity[offset+2], specificity[offset+3], sensitivity[offset+3], length=0.1)
}

for (tc in (1:19)/20) {
    x = (0:10)/10
    y = 2*tc - x
    lines(x, y, lty=2, col="gray")
    text(x, y, paste(tc*100,"%",sep=""), col="gray")
}

legend(x="topright", bty="n",
       c("graph", "FRs"),
       pch=19, col=colors
       )


## plot(totCorrect, MCC, col=colors, xlim=c(0,100), ylim=c(-1,1))
## for (i in 1:length(graph)) {
##     offset = (i-1)*3
##     text(totCorrect[offset+1], MCC[offset+1], graph[i], pos=1)
##     arrows(totCorrect[offset+1], MCC[offset+1], totCorrect[offset+2], MCC[offset+2], col="darkgreen", length=0.1)
##     arrows(totCorrect[offset+2], MCC[offset+2], totCorrect[offset+3], MCC[offset+3], col="blue", length=0.1)
## }
