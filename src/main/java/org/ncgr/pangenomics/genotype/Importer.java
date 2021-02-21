package org.ncgr.pangenomics.genotype;

import java.util.Map;
import java.util.LinkedHashMap;
import java.util.TreeSet;

/**
 * Abstract class which defines the class variables for graph importers.
 */
public abstract class Importer {

    // the Nodes we import, keyed by id
    public Map<Long,Node> nodes = new LinkedHashMap<>();

    // map each sample to the NodeSet it traverses
    public Map<String,NodeSet> sampleNodeSets = new LinkedHashMap<>();
    
    // map each Node to the Set of samples that traverse it
    public Map<Node,TreeSet<String>> nodeSamples = new LinkedHashMap<>();
}
