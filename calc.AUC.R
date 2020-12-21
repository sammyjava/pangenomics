#####################################################################
## calculate the AUC for a set of case/control polygenic risk scores
##
## input: dataframe of "sample", "label", "score"
#####################################################################

calc.AUC = function(prs) {

    nCase = nrow(prs[prs$label=="case",])
    nControl = nrow(prs[prs$label=="ctrl",])

    caseScore = prs$score[prs$label=="case"]
    controlScore = prs$score[prs$label=="ctrl"]

    AUC = 0
    for (i in 1:nCase) {
        for (j in 1:nControl) {
            if (caseScore[i] > controlScore[j]) AUC = AUC + 1;
        }
    }

    return(AUC/nCase/nControl)
}
