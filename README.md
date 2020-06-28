This directory contains classes for working with pan-genomic graphs and frequented regions, based on the paper
```
Cleary, et al., "Exploring Frequented Regions in Pan-Genomic Graphs", IEEE/ACM Trans Comput Biol Bioinform. 2018 Aug 9. PMID:30106690 DOI:10.1109/TCBB.2018.2864564
``` 
## Building
The project is set up with dependencies managed with the [Gradle build tool](https://gradle.org/). To build the distribution, simply run
```
$ ./gradlew installDist
```
This will create a distribution under `build/install` that is used by the various run scripts.

### org.ncgr.pangenomics
This contains two packages with similarly-named classes:

**org.ncgr.pangenomics.allele** which contains classes for working with allele-based sequence graphs
**org.ncgr.pangenomics.genotype** which contains classes for working with genotype graphs

Basic graph-related classes, not particularly specific to frequented regions:

`PangenomicGraph` extends org.jgrapht.graph.DirectedAcyclicGraph and stores a graph with methods for reading it in from files and various output methods.
There is a `main` class for creating a graph from input data such as a GFA or VCF file.

`Node` encapsulates a node in a Graph: its ID (a long) and, for sequence graphs, its sequence.

`NodeSet` encapsulates a set of nodes in a Graph. NodeSet implements Comparable. There is a method `merge()` for merging two NodeSets.
(These are called "node clusters" in the paper above, but since I've implemented it as an extension of TreeSet, I've used "Set").

`Path` encapsulates a path through a Graph, along with its full sequence in the case of sequence graphs.

### org.ncgr.pangenomics.[allele/genotype].fr
Frequented regions-related code.

`FrequentedRegion` stores a FrequentedRegion, containing a NodeSet along with the supporting subpaths of the full set of Paths in a Graph, and lots of methods.

`FRFinder` contains a `main()` method for finding FRs based on a bunch of parameters.

`FRPair` is a utility class that contains two FRs and the result of merging them, and is used in the search loop in `FRFinder.findFRs()`.

### org.ncgr.svm
LIBSVM-based Support Vector Machine classes.

### org.ncgr.weka
Weka-based supervised classifier classes.
