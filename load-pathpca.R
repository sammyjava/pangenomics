##
## load PCA data and run prcomp
##

## paths data frame from pathpca file
prefix = readline(prompt="graph prefix (e.g. 3q29): ")
paths = t(read.table(paste(prefix,"pathpca.txt", sep=".")))

## which are cases and which are controls
cases = endsWith(rownames(paths), "case")
ctrls = endsWith(rownames(paths), "ctrl")

caseNum = nrow(paths[cases,])
ctrlNum = nrow(paths[ctrls,])

labels = array(length(rownames(paths)))
pathnames = array(length(labels))
parts = strsplit(rownames(paths), ".", fixed=T)
for (i in 1:length(rownames(paths))) {
    pathnames[i] = parts[[i]][1]
    labels[i] = parts[[i]][2]
}

