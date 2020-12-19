##
## plot the results of a priority scan
##

source("load-pri-scan.R")

## plot(run1.sensitivity, run1.specificity,
##      xlim=c(0,1), ylim=c(0,1),
##      xlab="Sensitivity", ylab="Specificity",
##      col="darkblue", pch=19)
## text(run1.sensitivity, run1.specificity, round(run1.fracFRs,2), pos=1)
## points(run2.sensitivity, run2.specificity,
##        xlim=c(0,1), ylim=c(0,1),
##        xlab="Sensitivity", ylab="Specificity",
##        col="darkgreen", pch=19)
## text(run2.sensitivity, run2.specificity, round(run2.fracFRs,2), pos=1)
## lines(seq(0,1,0.1)*1.0, seq(1,0,-0.1)*1.0, col="gray")
## lines(seq(0,1,0.1)*1.1, seq(1,0,-0.1)*1.1, col="gray")
## lines(seq(0,1,0.1)*1.2, seq(1,0,-0.1)*1.2, col="gray")
## legend(x="topright", col=c("darkblue","darkgreen"), pch=19, c(run1.graph,run2.graph))

xmin=1
xmax=1000
ymin=0.0
ymax=0.2
plot(run1.minPri, run1.MCC, col="darkblue", pch=19,
     log="x",
     xlim=c(xmin,xmax), ylim=c(ymin,ymax),
     xlab="Minimum Priority", ylab="MCC")
points(run2.minPri, run2.MCC, col="darkgreen", pch=19)
text(run1.minPri, run1.MCC, round(run1.fracFRs,2), pos=1)
text(run2.minPri, run2.MCC, round(run2.fracFRs,2), pos=1)
legend(x="topright", col=c("darkblue","darkgreen"), pch=19, c(run1.graph,run2.graph))
