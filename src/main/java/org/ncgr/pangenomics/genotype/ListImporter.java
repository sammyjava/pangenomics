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
import java.util.Set;
import java.util.TreeSet;

/**
 * Importer for PLINK --list output files.
 *
 * @author Sam Hokin
 */
public class ListImporter extends Importer {

    boolean verbose;
    
    /**
     * Read the nodes and samples in from a List file, for the samples listed in desiredSamples.
     * Each locus has four lines: HOM, HET, REF, and NC with each sample listed TWICE!
     * contig   identifier              gtype   smp1    smp1    smp2    smp2    smp3    smp3    smp4   ...
     * 6	AA_A_9_30018537_FS	AA	174	174	509	509	1099	1099	1360   ...
     * 6	AA_A_9_30018537_FS	AP	42	42	58	58	76	76	107    ...
     * 6	AA_A_9_30018537_FS	PP	45	45	55	55	57	57	59     ...
     * 6	AA_A_9_30018537_FS	00
     */
    public void read(File listFile, Set<String> desiredSamples) throws FileNotFoundException, IOException {
        if (verbose) System.err.println("Reading nodes and samples from PLINK LIST file "+listFile.getName());
	// force alpha sample ordering
	TreeSet<String> orderedSamples = new TreeSet<>(desiredSamples);
	// spin through the List records storing the loaded orderedSamples and nodes in local maps
	long nodeId = 0;
	Map<String,Node> nodesMap = new HashMap<>();
	String line = null;
	BufferedReader listReader = new BufferedReader(new FileReader(listFile));
	while ((line=listReader.readLine())!=null) {
	    String[] fields = line.split("\\t");
	    if (fields.length==3) continue; // empty no-call line
	    String contig = fields[0];
	    String id = fields[1];
	    String genotype = fields[2]; // includes "00" if any samples
	    String nodeString = id+"_"+genotype;
	    // handle multiple lines with same node
	    Node n = nodesMap.get(nodeString);
	    if (n==null) {
		// add a new node
		nodeId++;
		n = new Node(nodeId, id, contig, 0, 0, genotype, 0.0); // GF will be updated later
		nodes.put(nodeId, n);
		nodesMap.put(nodeString, n);
	    }
	    // read the samples for this line that are in orderedSamples
	    TreeSet<String> lineSamples = new TreeSet<>();
	    for (int i=3; i<fields.length; i++) {
		if (orderedSamples.contains(fields[i])) {
		    lineSamples.add(fields[i]);
		}
	    }
	    // spin through these samples, incrementing the maps
	    for (String sampleName : lineSamples) {
		// nodes per sample
		NodeSet sampleNodes = sampleNodeSets.get(sampleName);
		if (sampleNodes==null) {
		    // new sample
		    sampleNodes = new NodeSet();
		}
		sampleNodes.add(n);
		sampleNodeSets.put(sampleName, sampleNodes);
		// samples per node
		TreeSet<String> samples = nodeSamples.get(n);
		if (samples==null) samples = new TreeSet<>();
		samples.add(sampleName);
		nodeSamples.put(n, samples);
	    }
	}
	// update the nodes with their genotype (not allele) frequencies
	for (Node n : nodeSamples.keySet()) {
            Set<String> samples = nodeSamples.get(n);
            n.gf = (double)samples.size() / (double)sampleNodeSets.size();
        }
	if (verbose) {
	    System.err.println("ListImporter read "+sampleNodeSets.size()+" samples and "+nodes.size()+" nodes.");
	}
    }
}
