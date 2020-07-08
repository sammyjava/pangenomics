##
## load FR path support for PCA
##

prefix = readline(prompt="FR file prefix (ex. HTT.400-0.5-1 or HTT.save): ")
pathfrs = read.table(file=paste(prefix,".pathfrs.txt",sep=""), stringsAsFactors=FALSE, check.names=FALSE)
