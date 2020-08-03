##
## load path FR support vectors
##

prefix = readline(prompt="FR file prefix (ex. HTT.400-0.5-1 or HTT.save): ")
pathfrs = t(read.table(file=paste(prefix,".pathfrs.txt",sep=""), stringsAsFactors=FALSE, check.names=FALSE))

## which are cases and which are controls
cases = endsWith(rownames(pathfrs), "case")
ctrls = endsWith(rownames(pathfrs), "ctrl")

caseNum = nrow(pathfrs[cases,])
ctrlNum = nrow(pathfrs[ctrls,])

labels = array(length(rownames(pathfrs)))
pathnames = array(length(labels))
parts = strsplit(rownames(pathfrs), ".", fixed=T)
for (i in 1:length(rownames(pathfrs))) {
    pathnames[i] = parts[[i]][1]
    labels[i] = parts[[i]][2]
}
