## ROC plot with legend
plot.ROC = function(ROC) {
    plot(ROC$FPR, ROC$TPR,
         pch=19, xlab="FPR", ylab="TPR",
         xlim=c(0.0,1.0), ylim=c(0.0,1.0),
         col=ROC$col,
         cex=0.5
         )

    ## hack: we draw arrows but with very special "arrowheads"
    arrows(ROC$FPR, ROC$TPR-ROC$TPR.err, ROC$FPR, ROC$TPR+ROC$TPR.err, length=0.05, angle=90, code=3, col=ROC$col)
    arrows(ROC$FPR-ROC$FPR.err, ROC$TPR, ROC$FPR+ROC$FPR.err, ROC$TPR, length=0.05, angle=90, code=3, col=ROC$col)

    ## dotted line for 50/50 draw
    lines(c(0,1), c(0,1), lty=2, col="gray")
}

