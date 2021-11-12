##
## Return a data frame containing the case/control stats per node.
## nodes are in a dataframe:  id rs contig start end genotype gf
## paths are in a dataframe:  id label path color
## paths$path is a string of the form: "[3,8,9,12,15,19,22,25,26,29,34,35,51,58,...]"
##
getNodeStats = function(nodes, paths) {
    ## cases, ctrls, caseRef, caseHet, caseHom, ctrlRef, ctrlHet, ctrlHom
    stats = data.frame(row.names=nodes$id)
    stats$cases = 0
    stats$ctrls = 0
    casenum = nrow(paths[paths$label=="case",])
    ctrlnum = nrow(paths[paths$label=="ctrl",])
    for (i in 1:nrow(paths)) {
        ## extract the path vector from the string representation
        id = paths$id[i]
        label = paths$label[i]
        pathstring = paths$path[i]
        pathstring = substr(pathstring, 2, nchar(pathstring))
        pathstring = substr(pathstring, 1, nchar(pathstring)-1)
        nodesvector = scan(text=pathstring, sep=",", quiet=TRUE)
        ## increment the node case/ctrl counts
        if (label=="case") {
            stats$cases[nodesvector] = stats$cases[nodesvector] + 1
        } else if (label=="ctrl") {
            stats$ctrls[nodesvector] = stats$ctrls[nodesvector] + 1
        }
    }
    ## prune nodes with no participation
    stats = stats[!is.null(stats$cases),]
    stats = stats[(stats$cases+stats$ctrls)>0,]
    ## spin through the nodes and compute stats
    pos = stats$cases>0 & stats$ctrls>0
    stats$OR[pos] = stats$cases[pos]/casenum / (stats$ctrls[pos]/ctrlnum)
    stats$logOR[pos] = log10(stats$OR[pos])
    ## more pruning
    stats = stats[!is.na(stats$OR),]
    ## Fisher
    for (i in 1:nrow(stats)) {
        casespos = stats$cases[i]
        ctrlspos = stats$ctrls[i]
        casesneg = casenum - casespos
        ctrlsneg = ctrlnum - ctrlspos
        stats$p[i] = fisher.test(matrix(c(casespos, casesneg, ctrlspos, ctrlsneg), nrow=2))$p.value
    }
    ## return
    return(stats)
}
