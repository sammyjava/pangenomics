## load the graph data, using the training paths

graph = readline(prompt="Graph (ex. HTT): ")
nodesFilename = paste(graph, ".nodes.txt", sep="")
pathsFilename = paste(graph, ".training.paths.txt", sep="")

## id  rs                  contig start    end      genotype gf
## 1   AA_A_9_30018537_FS  6      0        0        AA       0.06124219
nodes = read.table(nodesFilename, header=FALSE, stringsAsFactors=FALSE)
colnames(nodes) = c("id","rs","contig","start","end","genotype","gf")
    
## id  label  path  [color]
paths = read.table(pathsFilename, header=FALSE, stringsAsFactors=FALSE)
colnames(paths) = c("id", "label", "path")
paths$color[paths$label=="ctrl"] = "blue"
paths$color[paths$label=="case"] = "red"

## drop identifiers and get unique label+path
paths.unique = unique(paths[,2:4])

## we may have some rows with different labels but the same path, so remove those completely
paths.unique = paths.unique[!duplicated(paths.unique$path),]

## create matrix of path participation, 0/1 for not/traversing a node
pathmatrix = matrix(data=0, nrow=nrow(paths.unique), ncol=max(nodes$id))

for (i in 1:nrow(paths.unique)) {
    ## extract the path vector from the string representation
    pathstring = paths.unique$path[i]
    pathstring = substr(pathstring, 2, nchar(pathstring))
    pathstring = substr(pathstring, 1, nchar(pathstring)-1)
    nodesvector = scan(text=pathstring, sep=",", quiet=TRUE)
    for (j in nodesvector) {
        pathmatrix[i,j] = 1
    }
}
