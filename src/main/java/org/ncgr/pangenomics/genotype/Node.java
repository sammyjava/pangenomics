package org.ncgr.pangenomics.genotype;

import java.io.Serializable;

import java.text.DecimalFormat;

/**
 * Container for a node with genotype information as well as the range on the genome.
 * An int id is used for simplicity and the calling code is expected to enforce uniqueness.
 */
public class Node implements Comparable, Serializable {
    static DecimalFormat dec = new DecimalFormat("0.00");

    public long id;
    public String contig;
    public int start;
    public int end;
    public String rs;        // NCBI rsID or other identifier
    public String genotype;
    public double gf;        // genotype (not allele) frequency
    public boolean isCalled;

    /**
     * Minimal constructor, nothing but an id.
     */
    public Node(long id) {
        this.id = id;
    }

    /**
     * Construct the full Monty.
     */
    public Node(long id, String rs, String contig, int start, int end, String genotype, double gf) {
        this.id = id;
        this.rs = rs;
	this.contig = contig;
	this.start = start;
        this.end = end;
	this.genotype = genotype;
        this.gf = gf;
	isCalled = !isNoCall();
    }

    /**
     * Construct from a nodes.txt file line.
     * 0   1           2      3        4        5        6
     * id  rs          contig start    end      genotype gf
     * 1   rs112943240 6      25726329 25726329 C/CTT    0.43
     */
    public Node(String line) {
	String[] parts = line.split("\t");
	id = Long.parseLong(parts[0]);
	rs = parts[1];
	contig = parts[2];
	start = Integer.parseInt(parts[3]);
	end = Integer.parseInt(parts[4]);
	genotype = parts[5];
	gf = Double.parseDouble(parts[6]);
	if (rs.equals(".")) rs = null;
	isCalled = !isNoCall();
    }

    /**
     * Return true if this Node is a no-call. Several no-call strings can be put here.
     */
    public boolean isNoCall() {
	if (genotype==null) {
	    System.err.println("ERROR: Node.isNoCall() called but genotype is null.");
	    System.exit(1);
	}
        return genotype.equals("./.") || genotype.equals("00");
    }

    String[] getAlleles() {
	String[] alleles = new String[0];
	if (genotype.contains("/")) {
            // G/T
            alleles = genotype.split("/"); // unphased
        } else if (genotype.contains("|")) {
            // G|T
            alleles = genotype.split("/"); // phased
        } else if (genotype.length()==2) {
            // PA
            alleles = new String[2];
            alleles[0] = String.valueOf(genotype.charAt(0));
            alleles[1] = String.valueOf(genotype.charAt(1));
        }
        return alleles;
    }

    /**
     * Return true if this Node is a homozygous call.
     */
    public boolean isHomozygous() {
        String[] alleles = getAlleles();
	if (alleles.length==2) {
	    return alleles[0].equals(alleles[1]);
	} else {
	    return false;
	}
    }

    /**
     * Return true if this Node is a heterozygous call.
     */
    public boolean isHeterozygous() {
	String[] alleles = getAlleles();
	if (alleles.length==2) {
	    return !alleles[0].equals(alleles[1]);
	} else {
	    return false;
	}
    }

    /**
     * Return a string summary of this node.
     */
    @Override
    public String toString() {
        String location = contig;
        if (start>0) location += ":"+start;
        if (end>0) location += "-"+end;
        String identifier = "";
        if (rs!=null) identifier = rs;
	return "["+id+"] "+location+" "+identifier+" "+genotype+" "+dec.format(gf);
    }

    /**
     * Return a unique key for this node based on genomic quantities.
     */
    public String getKey() {
	return getKey(contig, start, end, rs, genotype);
    }

    /**
     * Return a key for the given genomic quanties.
     */
    public static String getKey(String contig, int start, int end, String rs, String genotype) {
	return contig+":"+start+"-"+end+" "+rs+" "+genotype;
    }

    /**
     * Nodes are equal if they:
     * 1. have the same id if rs==null or genotype==null or contig==null.
     * 2. have the same id and rs and genotype if those are not null.
     * 3. have the same id, contig, start, end, and genotype if those are all present.
     */
    @Override
    public boolean equals(Object o) {
	Node that = (Node) o;
	boolean equal = (this.id==that.id);
	if (this.rs!=null && that.rs!=null && this.genotype!=null && that.genotype!=null) {
	    equal = equal && this.rs.equals(that.rs) && this.genotype.equals(that.genotype);
	}
	if (this.contig!=null && that.contig!=null && this.start!=0 && that.start!=0 && this.end!=0 && that.end!=0 && this.genotype!=null && that.genotype!=null) {
	    equal = equal && this.contig.equals(that.contig) && this.start==that.start && this.end==that.end && this.genotype.equals(that.genotype);
	}
	return equal;
    }

    /**
     * Have to override hashCode() for Map keys.
     */
    @Override
    public int hashCode() {
	return (int) this.id;
    }

    /**
     * Nodes are compared by their id.
     */
    public int compareTo(Object o) {
	Node that = (Node) o;
        return Long.compare(this.id, that.id);
    }
}
