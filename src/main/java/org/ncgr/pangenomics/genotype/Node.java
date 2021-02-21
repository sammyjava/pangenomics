package org.ncgr.pangenomics.genotype;

import java.io.Serializable;

/**
 * Container for a node with genotype information as well as the range on the genome.
 * An int id is used for simplicity and the calling code is expected to enforce uniqueness.
 */
public class Node implements Comparable, Serializable {
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
	this.contig = contig;
	this.start = start;
        this.end = end;
        this.rs = rs;
	this.genotype = genotype;
        this.gf = gf;
	isCalled = !isNoCall();
    }

    /**
     * Return true if this Node is a no-call. Several no-call strings can be put here.
     */
    public boolean isNoCall() {
	if (genotype==null) {
	    System.err.println("ERROR: Node.isNoCall() called but genotype="+genotype);
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
	return "["+id+"] "+contig+":"+start+"-"+end+" "+rs+" "+genotype+" "+gf;
    }

    /**
     * Nodes are equal if they have the same id, contig, start, end, and genotype.
     */
    @Override
    public boolean equals(Object o) {
	Node that = (Node) o;
        return this.id==that.id &&
	    this.contig.equals(that.contig) &&
	    this.start==that.start &&
	    this.end==that.end &&
	    this.genotype.equals(that.genotype);
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
        return Long.compare(this.id,that.id);
    }
}
