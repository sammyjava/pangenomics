source("load-libsvm.R")

colors = c("black", "darkgreen", "blue")
        
plot(precision, recall, col=colors,
     xlim=c(.45,.70), ylim=c(.45,.70),
     xlab="precision", ylab="recall"
     )

for (i in 1:length(graph)) {
    offset = (i-1)*3
    text(precision[offset+1], recall[offset+1], graph[i], pos=1)
    arrows(precision[offset+1], recall[offset+1], precision[offset+2], recall[offset+2], col="darkgreen", length=0.1)
    arrows(precision[offset+2], recall[offset+2], precision[offset+3], recall[offset+3], col="blue", length=0.1)
}

points((1.0-.412), .502, col="black")
text((1.0-.412), .502, "HLAC", pos=1)

legend(x="topright", bty="n",
       c("Graph", "FRs", "Combined"),
       pch=1, col=colors
       )

