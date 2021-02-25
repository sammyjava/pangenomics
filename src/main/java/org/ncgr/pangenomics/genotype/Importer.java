package org.ncgr.pangenomics.genotype;

import java.io.File;
import java.io.IOException;

import java.util.TreeMap;
import java.util.TreeSet;

/**
 * Abstract class which defines the public objects and methods for graph importers.
 */
public abstract class Importer {
    // the Nodes we import, keyed by Node.id
    public TreeMap<Long,Node> nodeIdMap = new TreeMap<>();

    // the Nodes we import, keyed by Node.getKey() (optional)
    public TreeMap<String,Node> nodeKeyMap = new TreeMap<>();

    // map each Sample to the NodeSet it traverses
    public TreeMap<Sample,NodeSet> sampleNodeSets = new TreeMap<>();
    
    // map each Node to the Set of Samples that traverse it
    public TreeMap<Node,TreeSet<Sample>> nodeSamples = new TreeMap<>();

    // all Importers must enable verbosity
    public boolean verbose;

    // load nodes from a file without restrictions
    abstract void readNodes() throws IOException;
    
    // load nodes from a file, with MAF restrictions
    abstract void readNodes(double minMAF, double maxMAF) throws IOException;

    // load paths for the samples in a labels file
    abstract void readPaths(File labelsFile) throws IOException;

    // set verbosity flag
    abstract void setVerbose(boolean flag);
}
