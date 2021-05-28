MMC = function(TPR,FPR) {
    // assume that we have equal cases and controls
    TNR = 1.0 - FPR
    FNR = 1.0 - TPR
    return((TPR*TNR - FPR*FNR) / sqrt((TPR+FPR)*(TPR+FNR)*(TNR+FPR)*(TNR+FNR)))
}
