package org.ncgr.pangenomics.genotype;

import java.io.File;
import java.io.IOException;

import java.util.TreeMap;
import java.util.TreeSet;

/**
 * Abstract class which defines the class variables and methods for graph importers.
 */
public abstract class Importer {
    // the Nodes we import, keyed by id
    public TreeMap<Long,Node> nodes = new TreeMap<>();

    // map each Sample to the NodeSet it traverses
    public TreeMap<Sample,NodeSet> sampleNodeSets = new TreeMap<>();
    
    // map each Node to the Set of Samples that traverse it
    public TreeMap<Node,TreeSet<Sample>> nodeSamples = new TreeMap<>();

    // load nodes from a file without restrictions
    abstract void readNodes() throws IOException;
    
    // load nodes from a file, with MAF restrictions
    abstract void readNodes(double minMAF, double maxMAF) throws IOException;

    // load paths through the desired nodes restricting to samples from a labels file
    abstract void readPaths(File labelsFile, TreeMap<Long,Node> desiredNodes) throws IOException;
}
