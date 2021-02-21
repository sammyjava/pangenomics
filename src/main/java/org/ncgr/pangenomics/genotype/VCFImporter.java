package org.ncgr.pangenomics.genotype;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import java.util.Collections;
import java.util.List;
import java.util.LinkedList;
import java.util.HashMap;
import java.util.Set;
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

    public boolean verbose = false;

    /**
     * Read the nodes and paths from a VCF file for ALL of its samples with the given MAF restrictions.
     */
    public void read(File file, double minMAF, double maxMAF) throws FileNotFoundException, IOException {
	if (verbose) System.err.println("Reading nodes and all samples from "+file.getName()+" with "+minMAF+"<=MAF<="+maxMAF);
        // create the VCF file reader
	VCFFileReader reader = new VCFFileReader(file);
	// load the VCF sample names
	VCFHeader vcfHeader = reader.getFileHeader();
	// get all samples
	Set<String> samples = new TreeSet<String>(vcfHeader.getGenotypeSamples());
	read(reader, samples, minMAF, maxMAF);
    }

    /**
     * Read the nodes and paths for the given samples in from a VCF file with the given MAF restrictions.
     *
     * 1 877558 rs4372192 C T 71.55 PASS AC=1;AF=4.04e-05;AN=24736;BaseQRankSum=-1.369;CCC=24750;... GT:AD:DP:GQ:PL 0/0:7,0:7:21:0,21,281 0/0:7,0:7:21:0,21,218 ...
     *
     * NOTE: non-calls (./.) are ignored.
     */
    public void read(File file, Set<String> samples, double minMAF, double maxMAF) throws FileNotFoundException, IOException {
        if (verbose) System.err.println("Reading nodes and "+samples.size()+" samples from "+file.getName()+" with "+minMAF+"<=MAF<="+maxMAF);
        // create the VCF file reader
	VCFFileReader reader = new VCFFileReader(file);
	read(reader, samples, minMAF, maxMAF);
    }
    
    /**
     * Read the nodes and paths in from a VCFFileReader for the desired samples and the given MAF restrictions.
     *
     * 1 877558 rs4372192 C T 71.55 PASS AC=1;AF=4.04e-05;AN=24736;BaseQRankSum=-1.369;CCC=24750;... GT:AD:DP:GQ:PL 0/0:7,0:7:21:0,21,281 0/0:7,0:7:21:0,21,218 ...
     *
     * NOTE: non-calls (./.) are ignored.
     */
    void read(VCFFileReader reader, Set<String> desiredSamples, double minMAF, double maxMAF) {
	// check that all desired samples are actually in the VCF
	VCFHeader header = reader.getFileHeader();
	List<String> genotypeSamples = header.getGenotypeSamples();
	for (String sampleName : desiredSamples) {
	    if (!genotypeSamples.contains(sampleName)) {
		System.err.println("ERROR: VCF file does not contain desired sample "+sampleName);
		System.exit(1);
	    }
	}
	// force alpha sample ordering
	TreeSet<String> orderedSamples = new TreeSet<>(desiredSamples);
	// spin through the VCF records storing the nodes and samples in maps
	int minMAFCount = 0;
	int maxMAFCount = 0;
        long nodeId = 0;
        HashMap<String,Node> nodesMap = new HashMap<>(); // key nodes by locus+genotype
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
	    // create Nodes for ALL of the distinct fully-called genotypes of this VariantContext
	    // NOTE: getGenotypes() is an ordered collection, so resulting Nodes should be numbered consistently
	    for (Genotype g : vc.getGenotypes()) {
		if (g.isCalled() && !g.isMixed()) {
		    String genotypeString = getGenotypeString(g);
		    String nodeKey = getNodeKey(vc, g);
		    if (!nodesMap.containsKey(nodeKey)) {
			nodeId++;
			Node n = new Node(nodeId, vc.getID(), vc.getContig(), vc.getStart(), vc.getEnd(), genotypeString, 0.0); // GF will be updated later
			nodes.put(nodeId, n);
			nodesMap.put(nodeKey, n);
		    }
                }
	    }
	    // add to sample paths for the desired samples and this VariantContext
	    for (String sampleName : orderedSamples) {
		// get the Node
		Genotype g = vc.getGenotype(sampleName);
		// bail if this is not a called genotype
		if (!g.isCalled()) continue;
		// add this Node to this sample's NodeSet
		Node n = nodesMap.get(getNodeKey(vc,g));
		NodeSet sampleNodes = sampleNodeSets.get(sampleName);
		if (sampleNodes==null) {
		    sampleNodes = new NodeSet(); // first time
		    sampleNodeSets.put(sampleName, sampleNodes);
		}
                sampleNodes.add(n);
		// add to this Node's TreeSet of sample names
                TreeSet<String> samples = nodeSamples.get(n);
		if (samples==null) {
		    samples = new TreeSet<>();
		    nodeSamples.put(n, samples);
		}
                samples.add(sampleName);
	    }
        }
        // update all of nodes with their genotype (not allele) frequencies WITHIN THE DESIRED SAMPLES
        for (Node n : nodeSamples.keySet()) {
            TreeSet<String> samples = nodeSamples.get(n);
            n.gf = (double)samples.size() / (double)orderedSamples.size();
        }
	if (verbose) {
	    System.err.println("VCFImporter read "+nodes.size()+" nodes and "+sampleNodeSets.size()+" samples.");
	    System.err.println(minMAFCount+" loci were removed with MAF<"+minMAF);
	    System.err.println(maxMAFCount+" loci were removed with MAF>"+maxMAF);
	}
    }

    /**
     * Return the minor allele frequency from a VariantContext, defined as the fraction of non-majority alleles / all alleles.
     * NOTE: this handles the common case where the REF allele is NOT the majority.
     */
    public static double getMAF(VariantContext vc) {
        List<Allele> alleles = vc.getAlleles();
        int highestCount = 0;
        int secondHighestCount = 0;
        for (Allele a : alleles) {
            int count = vc.getCalledChrCount(a);
            if (count>highestCount) {
                highestCount = count;
            } else if (count>secondHighestCount) {
                secondHighestCount = count;
            }
        }
	return (double)secondHighestCount / (double)vc.getCalledChrCount();
    }

    /**
     * Return a string uniquely representing a Genotype of the form AA/AT/TT
     */
    public static String getGenotypeString(Genotype g) {
	List<Allele> alleles = g.getAlleles();
	// ignore phase
	Collections.sort(alleles);
	String genotypeString = alleles.get(0).toString().replace("<","").replace(">",""); // these mess up HTML
	for (int i=1; i<alleles.size(); i++) {
	    genotypeString += "/"+alleles.get(i).toString().replace("<","").replace(">","");
	}
	return genotypeString;
    }

    /**
     * Return a string uniquely representing a VariantContext and Genotype, suitable for keying a map of nodes.
     * 6:23456-23457[AA/AT]
     */
    public static String getNodeKey(VariantContext vc, Genotype g) {
	return vc.getContig()+":"+vc.getStart()+"-"+vc.getEnd()+"["+getGenotypeString(g)+"]";
    }
}
