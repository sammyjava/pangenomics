package org.ncgr.pangenomics.genotype;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;

import java.util.Collections;
import java.util.List;
import java.util.LinkedList;
import java.util.Arrays;
import java.util.Map;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.Set;
import java.util.TreeSet;

/**
 * Importer for PLINK --list output files.
 *
 * @author Sam Hokin
 */
public class ListImporter extends Importer {

    File listFile;

    /**
     * Construct for building nodes from a List file
     */
    public ListImporter(File listFile) {
	this.listFile = listFile;
    }

    /**
     * Construct for building paths from a List file and nodes
     */
    public ListImporter(File listFile, TreeMap<Long,Node> nodeIdMap, TreeMap<String,Node> nodeKeyMap) {
	this.listFile = listFile;
	this.nodeIdMap = nodeIdMap;
	this.nodeKeyMap = nodeKeyMap;
    }

    /**
     * Read the nodes in from a List file that have MGF within the given limits.
     * This is done by reading them all in and then pruning.
     *
     */
    @Override
    public void readNodes(double minMGF, double maxMGF) throws IOException {
	// read all nodes first
	readNodes();
	// place each locus' GFs in a set
	Map<String,TreeSet<Double>> gfMap = new HashMap<>(); // keyed by rs
	for (Node n : nodeIdMap.values()) {
	    TreeSet<Double> gfSet = gfMap.get(n.rs);
	    if (gfSet==null) gfSet = new TreeSet<>();
	    gfSet.add(n.gf);
	    gfMap.put(n.rs, gfSet);
	}
	// determine the MGF for each rs
	Map<String,Double> mgfMap = new HashMap<>();
	for (String rs : gfMap.keySet()) {
	    TreeSet<Double> gfSet = gfMap.get(rs);
	    gfSet.remove(gfSet.last());   // last was majority genotype
	    mgfMap.put(rs, gfSet.last()); // last is now largest minor genotype
	}
	// determine nodes to remove
	Set<Node> nodesToRemove = new TreeSet<>();
	for (Node n : nodeIdMap.values()) {
	    double mgf = mgfMap.get(n.rs);
	    if (mgf<minMGF || mgf>maxMGF) {
		nodesToRemove.add(n);
	    }
	}
	// remove nodes from TreeMap<Long,Node> nodes and HashMap<String,Node> nodeKeyMap
	for (Node n : nodesToRemove) {
	    if (nodeIdMap.remove(n.id)==null) {
		System.err.println("ListImporter: removal of Node "+n+" from nodes failed.");
		System.exit(1);
	    }
	    if (nodeKeyMap.remove(n.getKey())==null) {
		System.err.println("ListImporter: removal of Node "+n+" from nodeKeyMap failed.");
		System.exit(1);
	    }
	}
	if (verbose) System.err.println("Removed "+nodesToRemove.size()+" nodes with MGF<"+minMGF+" and MGF>"+maxMGF);
    }

    /**
     * Read the nodes in from a List file.
     * Each locus has four lines: HOM, HET, REF, and NC with each sample listed TWICE!!
     * contig   identifier              gtype   smp1    smp1    smp2    smp2    smp3    smp3    smp4   ...
     * 6	AA_A_9_30018537_FS	AA	174	174	509	509	1099	1099	1360   ...
     * 6	AA_A_9_30018537_FS	AP	42	42	58	58	76	76	107    ...
     * 6	AA_A_9_30018537_FS	PP	45	45	55	55	57	57	59     ...
     * 6	AA_A_9_30018537_FS	00
     */
    @Override
    public void readNodes() throws IOException {
	// maps to calculate genotype frequencies
	Map<String,Map<String,Integer>> genotypeCounts = new HashMap<>(); // counts per identifier+genotype
	Map<String,Integer> totalCounts = new HashMap<>();                // total counts per identifier
	// spin through the List records storing the loaded desiredSamples and nodes in local maps
	long nodeId = 0;
	String line = null;
	BufferedReader listReader = new BufferedReader(new FileReader(listFile));
	while ((line=listReader.readLine())!=null) {
	    String[] fields = line.split("\\t");
	    if (fields.length==3) continue; // empty no-call line
	    String contig = fields[0];
	    String rs = fields[1];
	    String genotype = fields[2];    // includes "00" if any samples
	    String nodeKey = Node.getKey(contig, 0, 0, rs, genotype);
	    // handle possible multiple lines with same node
	    Node n = nodeKeyMap.get(nodeKey);
	    if (n==null) {
		// add a new node
		nodeId++;
		// plink -list files do not supply start, end; GF will be updated below
		n = new Node(nodeId, rs, contig, 0, 0, genotype, 0.0);
		nodeIdMap.put(nodeId, n);
		nodeKeyMap.put(nodeKey, n);
	    }
	    // store this node's sample count for GF calculation
	    int sampleCount = readSamples(line).size();
	    if (!totalCounts.containsKey(rs)) {
		totalCounts.put(rs, sampleCount);
	    } else {
		int totalSampleCount = totalCounts.get(rs) + sampleCount;
		totalCounts.put(rs, totalSampleCount);
	    }
	    if (!genotypeCounts.containsKey(rs)) {
		genotypeCounts.put(rs, new HashMap<String,Integer>());
	    }
	    Map<String,Integer> counts = genotypeCounts.get(rs);
	    counts.put(genotype, sampleCount);
	}
	listReader.close();
	// update the nodes with their genotype frequencies
	for (Node n : nodeIdMap.values()) {
	    int totalCount = totalCounts.get(n.rs); // identifier is called rs from NCBI
	    int genotypeCount = genotypeCounts.get(n.rs).get(n.genotype);
	    n.gf = (double)genotypeCount / (double)totalCount;
	}
	if (verbose) {
	    System.err.println("ListImporter read "+nodeIdMap.size()+" nodes from "+listFile.getName());
	}
    }

    /**
     * Load paths through the desired nodes for samples given in a labels file
     * Each locus has four lines: HOM, HET, REF, and NC with each sample listed TWICE!!
     * contig   identifier              gtype   smp1    smp1    smp2    smp2    smp3    smp3    smp4   ...
     * 6	AA_A_9_30018537_FS	AA	174	174	509	509	1099	1099	1360   ...
     * 6	AA_A_9_30018537_FS	AP	42	42	58	58	76	76	107    ...
     * 6	AA_A_9_30018537_FS	PP	45	45	55	55	57	57	59     ...
     * 6	AA_A_9_30018537_FS	00
     */
    @Override
    public void readPaths(File labelsFile) throws IOException {
	if (nodeIdMap.size()==0) {
	    System.err.println("ERROR: ListImporter nodeIdMap is empty.");
	    System.exit(1);
	}
	if (nodeKeyMap.size()==0) {
	    System.err.println("ERROR: ListImporter nodeKeyMap is empty.");
	    System.exit(1);
	}
	// load the desired samples (with labels)
	TreeSet<Sample> desiredSamples = Sample.readSamples(labelsFile);
	// spin through the list file
	String line = null;
	BufferedReader listReader = new BufferedReader(new FileReader(listFile));
	while ((line=listReader.readLine())!=null) {
	    String[] fields = line.split("\\t");
	    if (fields.length==3) continue; // empty no-call line
	    String contig = fields[0];
	    String rs = fields[1];
	    String genotype = fields[2];    // includes "00" if any samples
	    String nodeKey = Node.getKey(contig, 0, 0, rs, genotype);
	    Node n = nodeKeyMap.get(nodeKey);
	    if (n==null) continue; // purged node
	    TreeSet<Sample> lineSamples = readSamples(line); // these do not have labels
	    for (Sample sample : desiredSamples) {
		if (!lineSamples.contains(sample)) continue; // this sample isn't in this call
		// samples per node
		TreeSet<Sample> samples = nodeSamples.get(n);
		if (samples==null) {
		    samples = new TreeSet<>();
		}
		samples.add(sample);
		nodeSamples.put(n, samples);
		// nodes per sample
		NodeSet sampleNodes = sampleNodeSets.get(sample);
		if (sampleNodes==null) {
		    sampleNodes = new NodeSet();
		}
		sampleNodes.add(n);
		sampleNodeSets.put(sample, sampleNodes);
	    }
	}
	listReader.close();
	if (verbose) System.err.println("ListImporter loaded "+sampleNodeSets.size()+" samples from "+listFile.getName());
    }

    /**
     * Utility to read the sample names from a list file line. (Sample.label is null, of course.)
     */
    TreeSet<Sample> readSamples(String line) {
	TreeSet<Sample> samples = new TreeSet<>();
	String[] fields = line.split("\\t");
	for (int i=3; i<fields.length; i++) {
	    samples.add(new Sample(fields[i], null));
	}
	return samples;
    }

    /**
     * Set verbosity
     */
    public void setVerbose(boolean flag) {
	verbose = flag;
    }
}
