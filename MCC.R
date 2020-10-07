## Mathew's Correlation Coefficient for equal cases and controls
MCC = function(TPR, FPR) {
    TNR = 1.0 - FPR
    FNR = 1.0 - TPR
    denom = (TPR+FPR)*(TPR+FNR)*(TNR+FPR)*(TNR+FNR)
    if (denom==0) {
        denom = 1
    }
    return ((TPR*TNR-FPR*FNR)/sqrt(denom))
}
