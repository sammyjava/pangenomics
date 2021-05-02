package org.ncgr.pangenomics.genotype;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import java.text.DecimalFormat;

import java.util.Arrays;
import java.util.Collection;
import java.util.StringJoiner;
import java.util.List;
import java.util.LinkedList;
import java.util.Set;
import java.util.HashSet;
import java.util.TreeSet;
import java.util.Map;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

/**
 * Encapsulates a set of nodes in a Graph. NodeSets are comparable based on their content.
 *
 * @author Sam Hokin
 */
public class NodeSet extends TreeSet<Node> implements Comparable {
    static DecimalFormat dec = new DecimalFormat("0.00");

    /**
     * Empty constructor.
     */
    public NodeSet() {
        super();
    }

    /**
     * Construct given a Collection of Nodes.
     */
    public NodeSet(Collection<Node> nodes) {
        this.addAll(nodes);
    }

    /**
     * Construct given a single Node.
     */
    public NodeSet(Node node) {
	this.add(node);
    }

    /**
     * Construct given a string representation but no underlying graph.
     */
    public NodeSet(String str) {
	String afterOpenBracket = str.split("\\[")[1];
	String nodesWithCommas = afterOpenBracket.split("\\]")[0];
        List<String> nodeIds = Arrays.asList(nodesWithCommas.split(","));
        for (String s : nodeIds) {
            this.add(new Node(Long.parseLong(s)));
        }
    }

    /**
     * Construct from a map of id to Nodes and a string representation, e.g. "[5,7,15,33]".
     */
    public NodeSet(Map<Long,Node> nodeMap, String str) {
        Set<String> nodeStrings = new HashSet<>(Arrays.asList(str.replace("[","").replace("]","").split(",")));
        for (String s : nodeStrings) {
            long id = Long.parseLong(s);
            if (nodeMap.containsKey(id)) this.add(nodeMap.get(id));
        }
    }

    /**
     * Return the result of merging two NodeSets.
     * NOTE: does NOT run update() on the result!
     */
    public static NodeSet merge(NodeSet ns1, NodeSet ns2) {
        NodeSet merged = new NodeSet();
        merged.addAll(ns1);
        merged.addAll(ns2);
        return merged;
    }

    /**
     * Equality if exactly the same nodes, meaning the same string.
     */
    @Override
    public boolean equals(Object o) {
	NodeSet that = (NodeSet) o;
	return this.toString().equals(that.toString());
    }

    /**
     * Must override hashCode() for Map keys.
     */
    @Override
    public int hashCode() {
	return this.toString().hashCode();
    }

    /**
     * Compare alphabetically.
     */
    public int compareTo(Object o) {
	NodeSet that = (NodeSet) o;
        return this.toString().compareTo(that.toString());
    }

    /**
     * Return a readable summary string.
     */
    @Override
    public String toString() {
        String s = "[";
        StringJoiner joiner = new StringJoiner(",");
        for (Node node : this) {
            joiner.add(String.valueOf(node.id));
        }
        s += joiner.toString();
        s += "]";
        return s;
    }

    /**
     * Return true if this NodeSet is a superset of that NodeSet.
     */
    public boolean isSupersetOf(NodeSet that) {
        return this.size()>that.size() && this.containsAll(that);
    }
    
    /**
     * Find the Levenshtein distance between this and another NodeSet.
     *
     * @param that the other Nodeset
     * @return result distance, or -1
     */
    public int distanceFrom(NodeSet that) {
        // copy the NodeSets into lists for indexed access
        List<Node> left = new LinkedList<>(this);
        List<Node> right = new LinkedList<>(that);
        int n = left.size();
        int m = right.size();
        // trivial distance
        if (n == 0) {
            return m;
        } else if (m == 0) {
            return n;
        }
        if (n>m) {
            // swap the Lists to consume less memory
            final List<Node> tmp = left;
            left = right;
            right = tmp;
            n = m;
            m = right.size();
        }
        int[] p = new int[n + 1];
        // indexes into Lists left and right
        int i; // iterates through left
        int j; // iterates through right
        int upper_left;
        int upper;
        Node rightJ; // jth Node of right
        int cost; // cost
        for (i=0; i<=n; i++) {
            p[i] = i;
        }
        for (j=1; j<=m; j++) {
            upper_left = p[0];
            rightJ = right.get(j - 1);
            p[0] = j;
            for (i=1; i<=n; i++) {
                upper = p[i];
                cost = left.get(i-1).equals(rightJ) ? 0 : 1;
                // minimum of cell to the left+1, to the top+1, diagonally left and up +cost
                p[i] = Math.min(Math.min(p[i - 1] + 1, p[i] + 1), upper_left + cost);
                upper_left = upper;
            }
        }
        return p[n];
    }

    /**
     * Return true if the nodes in this nodeset all have the same contig:start value.
     */
    public boolean haveSamePosition() {
	String contig = "";
	int start = 0;
	boolean sameStart = true;
	for (Node node : this) {
	    if (start==0) {
		start = node.start;
		contig = node.contig;
	    } else if (node.start!=start || !node.contig.equals(contig)) {
		sameStart = false;
		break;
	    }
	}
	return sameStart;
    }

    /**
     * Some handy static methods.
     */
    public static void main(String[] args) throws IOException, FileNotFoundException {
	Options options = new Options();
        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd;

	// INPUT: NodeSet in bracket notation
	Option nodeSetOption = new Option("n", "nodeset", true, "NodeSet in [34,37,38,42] notation (required)");
	nodeSetOption.setRequired(true);
	options.addOption(nodeSetOption);
	// INPUT: nodes.txt
	Option nodesFileOption = new Option("nodes", "nodesfile", true, "read graph nodes from the given nodes.txt file (required)");
	nodesFileOption.setRequired(true);
	options.addOption(nodesFileOption);
	// ACTION: spit out nodes in NodeSet
	Option printOption = new Option("p", "print", false, "print out the Nodes in the given NodeSet");
	printOption.setRequired(false);
	options.addOption(printOption);
	
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println(e.getMessage());
            formatter.printHelp("NodeSet", options);
            System.exit(1);
            return;
        }

        // spit out help and exit if nothing supplied
        if (cmd.getOptions().length==0) {
            formatter.printHelp("NodeSet", options);
            System.exit(1);
            return;
        }

	// load the Nodes
	File nodesFile = new File(cmd.getOptionValue("nodesfile"));
	TXTImporter importer = new TXTImporter(nodesFile);
        importer.readNodes();

	// read the NodeSet
	NodeSet nodeSet = new NodeSet(importer.nodeIdMap, cmd.getOptionValue("nodeset"));

	// print out the NodeSet's nodes, but using tabs as separators
	if (cmd.hasOption("print")) {
	    for (Node n : nodeSet) {
		String location = n.contig;
		if (n.start>0) location += ":"+n.start;
		if (n.end>0) location += "-"+n.end;
		String identifier = "";
		if (n.rs!=null) identifier = n.rs;
		System.out.println("["+n.id+"]\t"+location+"\t"+identifier+"\t"+n.genotype+"\t"+dec.format(n.gf));
	    }
	}
    }
}
