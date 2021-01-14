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
import java.util.Set;
import java.util.HashSet;

/**
 * Importer for TXT files [graph].nodes.txt and [graph].paths.txt containing Nodes and Paths.
 *
 * @author Sam Hokin
 */
public class TXTImporter extends Importer {

    // map each sample name to its label
    public Map<String,String> sampleLabels = new HashMap<>();

    /**
     * Import instance values from a nodes text file and a paths text file.
     *
     * node line:
     * id     rs      contig  start   end     genotype  
     * 12764  rs12345 2       3228938 3229006 CCCACCCCTGCCCTGTCTGGGGCTGAAGTACAGTGCCACCCCTGCCCTGTCTGGGGCTGAAGGACAGTG/C
     *
     * path line:
     * name     label   nodes
     * 642913	case	[1,8,17,21,24,25,28,33,35,37,...]
     *
     * @param nodesFile the nodes File (typically [graph].nodes.txt)
     * @param pathsFile the paths File (typically [graph].paths.txt)
     */
    public void read(File nodesFile, File pathsFile) throws IOException {
        // read the nodes, storing in a map for path building
        if (verbose) System.err.println("Reading nodes from "+nodesFile.getName()+"...");
	Map<Long,Node> nodesMap = new HashMap<>();
        String line = null;
        BufferedReader nodesReader = new BufferedReader(new FileReader(nodesFile));
        while ((line=nodesReader.readLine())!=null) {
            String[] parts = line.split("\t");
            long id = Long.parseLong(parts[0]);
            String rs = parts[1];
            String contig = parts[2];
            int start = Integer.parseInt(parts[3]);
            int end = Integer.parseInt(parts[4]);
            String genotype = parts[5];
            double gf = Double.parseDouble(parts[6]);
            if (rs.equals(".")) rs = null;
            Node n = new Node(id, rs, contig, start, end, genotype, gf);
            nodesMap.put(n.id, n);
        }
        nodesReader.close();
        if (verbose) System.err.println("...read "+nodesMap.size()+" nodes.");
        // read the paths file with lines like
	// 0      1    2
	// 123ABC case [2,3,6,7,8,9,12,15,24]
        if (verbose) System.err.println("Reading paths from "+pathsFile.getName()+"...");
        BufferedReader pathsReader = new BufferedReader(new FileReader(pathsFile));
        List<String> lines = new LinkedList<String>();
        while ((line=pathsReader.readLine())!=null) {
	    String[] fields = line.split("\t");
	    String name = fields[0];
	    String label = fields[1];
	    sampleLabels.put(name, label);
	    // input NodeSet contains different Nodes from above, so grab the ones above by id identity
	    NodeSet inputNodeSet = new NodeSet(fields[2]);
	    NodeSet outputNodeSet = new NodeSet();
	    for (Node n : inputNodeSet) {
		if (nodesMap.containsKey(n.id)) {
		    outputNodeSet.add(nodesMap.get(n.id)); // id is key
		} else {
		    System.err.println("ERROR: node "+n.id+" in paths file is not present in nodes file.");
		}
	    }
	    sampleNodeSets.put(name, outputNodeSet);
	    for (Node n : outputNodeSet) {
		Set<String> samples;
		if (nodeSamples.containsKey(n)) {
		    samples = nodeSamples.get(n);
		} else {
		    samples = new HashSet<>();
		}
		samples.add(name);
		nodeSamples.put(n, samples);
	    }
	}
	// wrap up
	this.nodes = new LinkedList<Node>(nodesMap.values());
	if (verbose) System.err.println("TXTImporter read "+sampleNodeSets.size()+" sample paths and "+nodes.size()+" nodes.");
    }
}
