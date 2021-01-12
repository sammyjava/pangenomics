##
## load PCA data
##

## paths data frame from pathpca file
prefix = readline(prompt="graph prefix (e.g. 3q29): ")
paths = t(read.table(paste(prefix,"pathpca.txt", sep=".")))

## remove "X" prefix we get from above
rownames(paths) = substr(rownames(paths), 2, 100)

## load path labels, but sorted differently
labs = read.table(paste(prefix,"labels.txt", sep="."), header=TRUE, row.names=1)

## sort the labels according to the paths
labels = labs[rownames(paths),]

rm(labs)

