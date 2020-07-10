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
