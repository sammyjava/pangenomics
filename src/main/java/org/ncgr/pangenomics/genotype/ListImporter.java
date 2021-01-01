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
import java.util.HashSet;

/**
 * Importer for PLINK list output files. Creates public collections that can be used to build a graph.
 *
 * @author Sam Hokin
 */
public class ListImporter {

    // verbosity flag
    public boolean verbose = false;

    // the nodes we've imported, in order
    public List<Node> nodes = new LinkedList<>();

    // map each sample to the ordered list nodes it traverses
    public Map<String,List<Node>> sampleNodesMap = new HashMap<>();
    
    // map each Node to the Set of samples that traverse it
    public Map<Node,Set<String>> nodeSamplesMap = new HashMap<>();

    /**
     * Read the nodes and paths in from a List file, for the samples listed in samples.
     * Each locus has four lines: HOM, HET, REF, and NC with each sample listed TWICE!
     * contig   identifier              gtype   smp1    smp1    smp2    smp2    smp3    smp3    smp4   ...
     * 6	AA_A_9_30018537_FS	AA	174	174	509	509	1099	1099	1360   ...
     * 6	AA_A_9_30018537_FS	AP	42	42	58	58	76	76	107    ...
     * 6	AA_A_9_30018537_FS	PP	45	45	55	55	57	57	59     ...
     * 6	AA_A_9_30018537_FS	00
     */
    public void read(File listFile, Set<String> samples) throws FileNotFoundException, IOException {
        if (verbose) System.err.println("Reading nodes from LIST file");
	// spin through the List records storing the loaded samples and nodes in local maps
	long nodeId = 0;
	Map<String,Node> nodesMap = new HashMap<>();
	BufferedReader listReader = new BufferedReader(new FileReader(listFile));
	String listLine = null;
	while ((listLine=listReader.readLine())!=null) {
	    if (verbose) System.out.print(".");
	    String[] fields = listLine.split("\\t");
	    if (fields.length==3) continue; // empty no-call line
	    String contig = fields[0];
	    String rs = fields[1];
	    String genotype = fields[2]; // includes "00" if any samples
	    String nodeString = rs+"_"+genotype;
	    Node n;
	    if (nodesMap.containsKey(nodeString)) {
		n = nodesMap.get(nodeString);
	    } else {
		// add a new node
		nodeId++;
		n = new Node(nodeId, rs, contig, 0, 0, genotype, 0.0); // GF will be updated later
		nodes.add(n);
		nodesMap.put(nodeString, n);
	    }
	    // read the samples for this node that are in the supplied samples list, using a Set for uniqueness
	    Set<String> lineSamples = new HashSet<>();
	    for (int i=3; i<fields.length; i++) {
		if (samples.contains(fields[i])) {
		    lineSamples.add(fields[i]);
		}
	    }
	    // spin through the samples, incrementing the maps
	    for (String sampleName : lineSamples) {
		// nodes per sample
		List<Node> sampleNodes;
		if (sampleNodesMap.containsKey(sampleName)) {
		    sampleNodes = sampleNodesMap.get(sampleName);
		} else {
		    sampleNodes = new LinkedList<>();
		}
		sampleNodes.add(n);
		sampleNodesMap.put(sampleName, sampleNodes);
		// samples per node
		Set<String> nodeSamples;
		if (nodeSamplesMap.containsKey(n)) {
		    nodeSamples = nodeSamplesMap.get(n);
		} else {
		    nodeSamples = new HashSet<>();
		}
		nodeSamples.add(sampleName);
		nodeSamplesMap.put(n, nodeSamples);
	    }
	}
	// update the nodes with their genotype (not allele) frequencies
	for (Node n : nodeSamplesMap.keySet()) {
            Set<String> nodeSamples = nodeSamplesMap.get(n);
            n.gf = (double)nodeSamples.size() / (double)sampleNodesMap.size();
            // DEBUG
            if (n.rs.equals("HLA_DQB1_0302")) {
                System.out.println(n.id+"\t"+n.rs+"\t"+n.genotype+"\tnodeSamples.size()="+nodeSamples.size()+"\t"+"sampleNodesMap.size()="+sampleNodesMap.size()+"\tn.gf="+n.gf);
            }
        }
	if (verbose) {
	    System.err.println("");
	    System.err.println("ListImporter read "+sampleNodesMap.size()+" samples and "+nodes.size()+" nodes.");
	}
    }
}
