package org.ncgr.pangenomics.genotype;

import java.util.List;
import java.util.LinkedList;
import java.util.Map;
import java.util.HashMap;
import java.util.Set;

/**
 * Abstract class which defines the class variables for graph importers.
 */
public abstract class Importer {

    // verbosity flag
    public boolean verbose = false;

    // the Nodes we import
    public List<Node> nodes = new LinkedList<>();

    // map each sample to the NodeSet it traverses
    public Map<String,NodeSet> sampleNodeSets = new HashMap<>();
    
    // map each Node to the Set of samples that traverse it
    public Map<Node,Set<String>> nodeSamples = new HashMap<>();
    
}
