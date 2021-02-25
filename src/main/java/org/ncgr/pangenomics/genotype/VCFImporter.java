package org.ncgr.pangenomics.genotype;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import java.util.Collections;
import java.util.List;
import java.util.LinkedList;
import java.util.Map;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.Set;
import java.util.HashSet;
import java.util.TreeSet;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

/**
 * Importer for VCF files.
 *
 * @author Sam Hokin
 */
public class VCFImporter extends Importer {

    File vcfFile;

    /**
     * Construct for building nodes from a VCF file
     */
    public VCFImporter(File vcfFile) {
	this.vcfFile = vcfFile;
    }

    /**
     * Construct for building paths from a VCF file and nodes
     */
    public VCFImporter(File vcfFile, TreeMap<Long,Node> nodeIdMap, TreeMap<String,Node> nodeKeyMap) {
	this.vcfFile = vcfFile;
	this.nodeIdMap = nodeIdMap;
	this.nodeKeyMap = nodeKeyMap;
    }

    /**
     * Read nodes from a VCF file without any restrictions.
     */
    @Override
    public void readNodes() throws IOException {
	readNodes(0.0, 1.0);
    }

    /**
     * Read nodes from a VCF file for ALL of its samples with the given MAF restrictions.
     */
    @Override
    public void readNodes(double minMAF, double maxMAF) throws IOException {
        // create the VCF file reader
	VCFFileReader reader = new VCFFileReader(vcfFile);
	readNodes(reader, minMAF, maxMAF);
	reader.close();
       	if (verbose) System.err.println("Read "+nodeIdMap.size()+" nodes from "+vcfFile.getName());
    }

    /**
     * Read paths from a VCF file for the samples given in the labels file and the desired nodes.
     * NOTE: readNodes must be called first to populate nodeKeyMap!
     */
    @Override
    public void readPaths(File labelsFile) throws IOException {
	if (nodeKeyMap.size()==0) {
	    System.err.println("ERROR: you must call VCFImporter.readNodes() before calling VCFImporter.readPaths().");
	    System.exit(1);
	}
        // create the VCF file reader and header
	VCFFileReader reader = new VCFFileReader(vcfFile);
	VCFHeader header = reader.getFileHeader();
	// get all samples from the header
	Set<String> headerSampleNames = new TreeSet<String>(header.getGenotypeSamples());
	// get the desired samples from the labelsFile
	TreeSet<Sample> samples = Sample.readSamples(labelsFile);
	// check that the VCF contains our desired samples
	for (Sample sample : samples) {
	    if (!headerSampleNames.contains(sample.name)) {
		System.err.println("ERROR: file "+vcfFile.getName()+" does not contain sample "+sample);
		System.exit(1);
	    }
	}
	// now read in the paths
	readPaths(reader, samples);
	reader.close();
       	if (verbose) System.err.println("Read "+sampleNodeSets.size()+" sample paths from "+vcfFile.getName()+" and "+labelsFile.getName());
    }

    /**
     * Read the nodes TreeMap<Long,Node> nodes from a VCFFileReader with the given MAF restrictions.
     *
     * 1 877558 rs4372192 C T 71.55 PASS AC=1;AF=4.04e-05;AN=24736;BaseQRankSum=-1.369;CCC=24750;... GT:AD:DP:GQ:PL 0/0:7,0:7:21:0,21,281 0/0:7,0:7:21:0,21,218 ...
     *
     * NOTE: non-calls (./.) are ignored.
     */
    protected void readNodes(VCFFileReader reader, double minMAF, double maxMAF) {
	// spin through the VCF records storing the qualified nodes
	int minMAFCount = 0;
	int maxMAFCount = 0;
	long nodeId = 0;
	// load the qualified VariantContexts into a Set for parallel processing
	Set<VariantContext> vcSet = new HashSet<>();
	for (VariantContext vc : reader) {
	    // check this VariantContext's minor allele frequency, bail if outside range
	    double maf = getMAF(vc);
	    if (maf<minMAF) {
		minMAFCount++;
		continue;
	    } else if (maf>maxMAF) {
		maxMAFCount++;
		continue;
	    }
	    vcSet.add(vc);
	    // create a node for each called genotype, calculating the genotype frequency as well
	    TreeMap<String,Integer> genotypes = new TreeMap<>();
	    int totalCount = 0;
	    for (Genotype g : vc.getGenotypes()) {
		if (g.isNoCall() || g.isMixed()) continue; // only keep fully called
		String gstring = g.getGenotypeString();
		if (genotypes.containsKey(gstring)) {
		    int count = genotypes.get(gstring) + 1;
		    genotypes.put(gstring, count);
		} else {
		    genotypes.put(gstring, 1);
		}
		totalCount++;
	    }
	    for (String gstring : genotypes.keySet()) {
		String key = getNodeKey(vc, gstring);
		double gf = (double) genotypes.get(gstring) / (double) totalCount;
		Node n = new Node(++nodeId, vc.getID(), vc.getContig(), vc.getStart(), vc.getEnd(), gstring, gf);
		nodeIdMap.put(nodeId, n);
		nodeKeyMap.put(key, n);
	    }
	}
	if (verbose) {
	    System.err.println(minMAFCount+" loci were removed with MAF<"+minMAF);
	    System.err.println(maxMAFCount+" loci were removed with MAF>"+maxMAF);
	}
    }

    /**
     * Read the sample paths from a VCFFileReader for the given samples and ALL nodes.
     */
    protected void readPaths(VCFFileReader reader, TreeSet<Sample> samples) {
	for (VariantContext vc : reader) {
	    // spin through each sample to add to its NodeSet
	    for (Sample sample : samples) {
		// get the Node
		Genotype g = vc.getGenotype(sample.name);
		// bail if this is not a called genotype
		if (!g.isCalled()) continue;
		// add this Node to this sample's NodeSet if it's been loaded
		String key = getNodeKey(vc, g);
		if (nodeKeyMap.containsKey(key)) {
		    Node n = nodeKeyMap.get(key);
		    NodeSet sampleNodes = sampleNodeSets.get(sample);
		    if (sampleNodes==null) {
			sampleNodes = new NodeSet(); // first time
			sampleNodeSets.put(sample, sampleNodes);
		    }
		    sampleNodes.add(n);
		    // add to this Node's TreeSet of samples
		    TreeSet<Sample> sampleSet = nodeSamples.get(n);
		    if (sampleSet==null) {
			sampleSet = new TreeSet<>();
			nodeSamples.put(n, sampleSet);
		    }
		    sampleSet.add(sample);
		}
	    }
        }
    }
    
    /**
     * Return the minor allele frequency from a VariantContext, defined as the fraction of non-majority alleles / all alleles.
     * NOTE: this handles the common case where the REF allele is NOT the majority.
     */
    public static double getMAF(VariantContext vc) {
	int majorityCount = 0;
	Allele majorityAllele = null;
        for (Allele a : vc.getAlleles()) {
	    int count = vc.getCalledChrCount(a);
	    if (count>majorityCount) {
		majorityCount = count;
		majorityAllele = a;
	    }
        }
	int minorityCount = 0;
	for (Allele a : vc.getAlleles()) {
	    if (!a.equals(majorityAllele)) {
		minorityCount += vc.getCalledChrCount(a);
	    }
	}
	return (double)minorityCount / (double)vc.getCalledChrCount();
    }

    /**
     * Return a string uniquely representing a VariantContext and Genotype, suitable for keying a map of nodes.
     * 6:23456-23457[AA/AT]
     */
    static String getNodeKey(VariantContext vc, Genotype g) {
	return Node.getKey(vc.getContig(), vc.getStart(), vc.getEnd(), vc.getID(), g.getGenotypeString());
    }

    /**
     * Return a string uniquely representing a VariantContext and Genotype, suitable for keying a map of nodes.
     * 6:23456-23457[AA/AT]
     */
    static String getNodeKey(VariantContext vc, String gstring) {
	return Node.getKey(vc.getContig(), vc.getStart(), vc.getEnd(), vc.getID(), gstring);
    }

    /**
     * Set verbosity
     */
    public void setVerbose(boolean flag) {
	verbose = flag;
    }
}
