##
## plot the results of a minSupport scan
##

source("load-support-scan.R")

##################################################################################################################
## xmin=100
## xmax=10000
## ymin=0.1
## ymax=0.2
## plot(run1.minSup, run1.MCC,
##      xlim=c(xmin,xmax), ylim=c(ymin,ymax),
##      log="x",
##      pch=19, col="darkblue",
##      xlab="Minimum Support", ylab="MCC")
## points(run2.minSup, run2.MCC,
##        pch=19, col="darkgreen")
## lines(c(xmin,xmax), c(.145,.145), col="gray")
## lines(c(xmin,xmax), c(.190,.190), col="gray")
## for (i in 0:10) {
##     lines(c(i/10,i/10), c(0,1), col="gray")
## }
## text(run1.minSup, run1.MCC, paste(round(run1.fracFRs*100,1),"%"), col="darkblue", pos=1, cex=0.5)
## text(run2.minSup, run2.MCC, paste(round(run2.fracFRs*100,1),"%"), col="darkgreen", pos=1, cex=0.5)
## legend(x="topleft", col=c("darkblue","darkgreen"), pch=19, c(paste(run1.graph,"975 FRs"),paste(run2.graph,"1173 FRs")), bty="y")
##################################################################################################################


##################################################################################################################
## plot(run1.specificity, run1.sensitivity,
##      xlim=c(.55,.75), ylim=c(.4,.6),
##      pch=19, col="darkblue",
##      xlab="Specificity", ylab="Sensitivity")
## points(run2.specificity, run2.sensitivity,
##        pch=19, col="darkgreen")
## text(run1.specificity, run1.sensitivity, paste(round(run1.fracFRs*100,1),"%"), col="darkblue", pos=1, cex=0.5)
## text(run2.specificity, run2.sensitivity, paste(round(run2.fracFRs*100,1),"%"), col="darkgreen", pos=1, cex=0.5)
## legend(x="topright", col=c("darkblue","darkgreen"), pch=19, c(paste(run1.graph,"975 FRs"),paste(run2.graph,"1173 FRs")), bty="y")
##################################################################################################################

##################################################################################################################
xmin = 0.36
xmax = 0.44
ymin = 0.50
ymax = 0.58
plot(run1.FPR, run1.TPR,
     xlim=c(xmin,xmax), ylim=c(ymin,ymax),
     pch=19, col="darkblue",
     xlab="FPR (1-specificity)", ylab="TPR (sensitivity)")
points(run2.FPR, run2.TPR,
       pch=19, col="darkgreen")
text(run1.FPR, run1.TPR, run1.minSup, col="darkblue", pos=1, cex=0.5)
text(run2.FPR, run2.TPR, run2.minSup, col="darkgreen", pos=1, cex=0.5)
legend(x="topright", col=c("darkblue","darkgreen"), pch=19, c(paste(run1.graph,"975 FRs"),paste(run2.graph,"1173 FRs")), bty="y")
for (i in 0:100) {
    tot = i/100
    lines(c(0,1), c(0,1)+2*tot-1, col="gray")
    text(xmin, round(xmin+2*tot-1,2), tot, col="gray", cex=0.5)
    text(round(ymin-2*tot+1,2), ymin, tot, col="gray", cex=0.5)
}
#################################################################################################################



