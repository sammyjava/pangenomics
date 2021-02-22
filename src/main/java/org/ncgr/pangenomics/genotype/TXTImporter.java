package org.ncgr.pangenomics.genotype;

import org.jgrapht.Graph;
import org.jgrapht.io.GraphImporter;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;

import java.util.List;
import java.util.LinkedList;
import java.util.Map;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.Set;
import java.util.TreeSet;

/**
 * Importer for TXT files [graph].nodes.txt and [graph].paths.txt containing Nodes and Paths.
 *
 *    abstract void readNodes(double minMAF, double maxMAF) throws IOException;
 *    abstract void readPaths(File labelsFile, TreeMap<Long,Node> desiredNodes) throws IOException;
 *
 * @author Sam Hokin
 */
public class TXTImporter extends Importer {

    boolean verbose;

    File nodesFile;
    File pathsFile;

    // construct from a nodes.txt file and paths.txt file
    public TXTImporter(File nodesFile, File pathsFile) {
	this.nodesFile = nodesFile;
	this.pathsFile = pathsFile;
    }
    
    /**
     * Import all nodes from a nodes.txt file without constraints.
     */
    @Override
    public void readNodes() throws IOException {
	readNodes(0.0, 1.0);
    }

    /**
     * Import nodes from a nodes.txt file with MAF constraints.
     */
    @Override
    public void readNodes(double minMAF, double maxMAF) throws IOException {
        String line = null;
        BufferedReader reader = new BufferedReader(new FileReader(nodesFile));
        while ((line=reader.readLine())!=null) {
	    if (line.startsWith("#")) continue;
	    Node n = new Node(line);
            nodes.put(n.id, n);
        }
        reader.close();
        if (verbose) System.err.println("Read "+nodes.size()+" nodes from "+nodesFile.getName());
    }

    /**
     * Import paths from the pathsFile and a labels file to restrict the samples.
     */
    @Override
    public void readPaths(File labelsFile, TreeMap<Long,Node> desiredNodes) throws IOException {
	if (nodes.size()==0) {
	    System.err.println("TXTImporter: readNodes() must be called before readPaths().");
	    System.exit(1);
	}
	// read in the samples
	TreeSet<Sample> desiredSamples = Sample.readSamples(labelsFile);
	// read in the paths for the desired samples
	String line = null;
        BufferedReader reader = new BufferedReader(new FileReader(pathsFile));
        while ((line=reader.readLine())!=null) {
	    Sample sample = new Sample(line);
	    if (desiredSamples.contains(sample)) {
		NodeSet nodeSet = new NodeSet(line);
		sampleNodeSets.put(sample, nodeSet);
		// add to the nodeSamples
		for (Node n : nodeSet) {
		    TreeSet<Sample> samples = nodeSamples.get(n);
		    if (samples==null) samples = new TreeSet<>(); // first sample
		    samples.add(sample);
		    nodeSamples.put(n, samples); // update
		}
	    }
	}
	// show message if not all desired samples were read
	if (sampleNodeSets.size()!=desiredSamples.size()) {
	    System.err.println("WARNING: "+labelsFile.getName()+" contains "+desiredSamples.size()+" samples, "+
			       "while only "+sampleNodeSets.size()+" were loaded from "+pathsFile.getName());
	}
	// wrap up
	if (verbose) System.err.println("Read "+sampleNodeSets.size()+" sample paths from "+labelsFile.getName()+" and "+pathsFile.getName());
    }

    /**
     * Import paths from a paths.txt file, which means populating sampleNodeSets and nodeSamples.
     * NOTE: readNodes must be called first to populate nodes map!
     * name     label   nodeset
     * 642913	case	[1,8,17,21,24,25,28,33,35,37,...]
     */
    public void readPaths(TreeMap<Long,Node> desiredNodes) throws IOException {
	if (nodes.size()==0) {
	    System.err.println("TXTImporter: readNodes() must be called before readPaths().");
	    System.exit(1);
	}
	// read in the paths
	String line = null;
        BufferedReader reader = new BufferedReader(new FileReader(pathsFile));
        while ((line=reader.readLine())!=null) {
	    Sample sample = new Sample(line);
	    NodeSet nodeSet = new NodeSet(line);
	    sampleNodeSets.put(sample, nodeSet);
	    // add to the nodeSamples
	    for (Node n : nodeSet) {
		TreeSet<Sample> samples = nodeSamples.get(n);
		if (samples==null) samples = new TreeSet<>(); // first sample
		samples.add(sample);
		nodeSamples.put(n, samples); // update
	    }
	}
	// wrap up
	if (verbose) System.err.println("Read "+sampleNodeSets.size()+" sample paths from "+pathsFile.getName());
    }
}
