##
## load an FR file along with its parameters
## nodes support avgLen case ctrl size
##

## for plotting and analysis
library(dplyr)

prefix = readline(prompt="FR file prefix (ex. HTT.400-0.5-1): ")
frFilename = paste(prefix, ".frs.txt", sep="")

## read the FRs table
## nodes	size	support	case	ctrl	OR	p	pri
frs = read.table(frFilename, header=TRUE, stringsAsFactors=FALSE)
frs = frs[frs$nodes!="nodes",]
frs$size = as.numeric(frs$size)
frs$support = as.numeric(frs$support)
frs$case = as.numeric(frs$case)
frs$ctrl = as.numeric(frs$ctrl)
frs$OR = as.numeric(frs$OR)
frs$p = as.numeric(frs$p)
frs$pri = as.numeric(frs$pri)
frs = distinct(frs)

## rownames(frs) = frs$nodes
## frs$nodes = NULL

## frs$size = as.numeric(frs$size)
## frs$support = as.numeric(frs$support)
## frs$case = as.numeric(frs$case)
## frs$ctrl = as.numeric(frs$ctrl)
## frs$OR = as.numeric(frs$OR)
## frs$p = as.numeric(frs$p)
## frs$pri = as.numeric(frs$pri)

prefix.parts = strsplit(prefix, "-", fixed=TRUE)
graphPrefix = prefix.parts[[1]][1]
alpha = as.numeric(prefix.parts[[1]][2])
kappa = as.numeric(prefix.parts[[1]][3])

## load the parameters from the params.txt file
source("load-params.R")

## label counts (if exists)
labelFile = paste(graphPrefix,".labelcounts.txt",sep="")
labelsExist = file.exists(labelFile)
if (labelsExist) {
    labelCounts = read.delim(file=labelFile, header=FALSE, stringsAsFactors=FALSE, row.names=1)
    colnames(labelCounts) = c("count")
}

## nodes
## 1	rs114039523	6	29910286	29910286	T/T	0.6643768400392541
## 4	rs114039523	6	29910286	29910286	./.	8.92140244446427E-5
nodes = read.table(file=paste(graphName,"nodes","txt",sep="."), row.names=1, col.names=c("node","rs","chr","start","end","genotype","p"))

