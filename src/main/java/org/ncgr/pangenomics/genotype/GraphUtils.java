package org.ncgr.pangenomics.genotype;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;

import java.util.LinkedList;
import java.util.Map;
import java.util.HashMap;
import java.util.Set;
import java.util.HashSet;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

/**
 * A collection of static methods to do utility stuff on a graph.
 */
public class GraphUtils {

    /**
     * The main event.
     */
    public static void main(String[] args) throws IOException, FileNotFoundException {
	Options options = new Options();
        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd;

        // REQUIRED parameters
        Option graphOption = new Option("g", "graph", true, "name of graph");
        graphOption.setRequired(true);
        options.addOption(graphOption);

        // other parameters
        Option mgfOption = new Option("mgf", "mgf", true, "minimum MGF for inclusion");
        mgfOption.setRequired(false);
        options.addOption(mgfOption);
        // samples in VCF file
        Option vcfFileOption = new Option("vcf", "vcffile", true, "load sample(s) from a VCF file");
        vcfFileOption.setRequired(false);
        options.addOption(vcfFileOption);
	// samples in LIST file
	Option listFileOption = new Option("list", "listfile", true, "load sample(s) from a PLINK list file");
	listFileOption.setRequired(false);
	options.addOption(listFileOption);

        // actions
        Option prsOption = new Option("prs", "prs", false, "compute polygenic risk scores");
        prsOption.setRequired(false);
        options.addOption(prsOption);
	//
	Option sampleFileOption = new Option("samplefile", true, "file containing samples+labels to determine paths through the given graph (-g) from a VCF (-vcf) or LIST (-list)");
	sampleFileOption.setRequired(false);
	options.addOption(sampleFileOption);

        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println(e.getMessage());
            formatter.printHelp("GraphUtils", options);
            System.exit(1);
            return;
        }
        // spit out help and exit if nothing supplied
        if (cmd.getOptions().length==0) {
            formatter.printHelp("GraphUtils", options);
            System.exit(1);
            return;
        }

        // our PangenomicGraph
        PangenomicGraph graph = new PangenomicGraph();

        // populate graph instance vars from parameters
        graph.name = cmd.getOptionValue("graph");

        // load graph from TXT files
        graph.nodesFile = new File(graph.name+".nodes.txt");
        graph.pathsFile = new File(graph.name+".paths.txt");
        graph.loadTXT();
        graph.tallyLabelCounts();
        System.err.println("Graph has "+graph.vertexSet().size()+" nodes and "+graph.paths.size()+" paths: "+graph.labelCounts.get("case")+"/"+graph.labelCounts.get("ctrl")+" cases/controls");

        // options
        boolean computePRS = cmd.hasOption("prs");
	boolean computeSamplePaths = cmd.hasOption("samplefile");

        // actions
        if (computePRS) {
            double minMGF = 1e-2;
            if (cmd.hasOption("mgf")) {
                minMGF = Double.parseDouble(cmd.getOptionValue("mgf"));
            }
            computePRS(graph, minMGF);
        } else if (computeSamplePaths) {
	    Map<String,String> samples = new HashMap<>();
	    BufferedReader reader = new BufferedReader(new FileReader(cmd.getOptionValue("samplefile")));
	    String line = null;
	    while ((line=reader.readLine())!=null) {
		if (line.startsWith("#") || line.startsWith("sample") || line.trim().length()==0) {
		    continue;
		}
		String[] fields = line.split("\t");
		String sampleName = fields[0];
		String label = fields[1];
		samples.put(sampleName, label);
	    }
	    if (cmd.hasOption("vcf")) {
		computePathFromVCF(graph, cmd.getOptionValue("vcf"), samples);
	    } else if (cmd.hasOption("list")) {
		computePathFromList(graph, cmd.getOptionValue("list"), samples);
	    } else {
		System.err.println("ERROR: you must provide a VCF file (-vcf) or LIST file (-list) that contains the given samples");
	    }			
	}
    }

    /**
     * Compute the polygenic risk scores from a graph, tossing nodes with gf<minMGF.
     */
    public static void computePRS(PangenomicGraph graph, double minMGF) {
        // need node paths
        graph.buildNodePaths();
        // calculate map of log odds ratio
        Map<Node,Double> nodeOddsRatios = new HashMap<>();
        for (Node n : graph.getNodes()) {
            // public long id;
            // public String contig;
            // public int start;
            // public int end;
            // public String rs;
            // public String genotype;
            // public double gf;
            // public boolean isCalled;
            if (!(n.gf<minMGF)) {
                nodeOddsRatios.put(n, Math.log(graph.oddsRatio(n)));
            }
        }
        // sum log odds ratio over each path, storing in a map
        Map<Path,Double> pathPRS = new HashMap<>();
        for (Path p : graph.paths) {
            double prs = 0.0;
            int num = 0;
            for (Node n : p.getNodes()) {
                if (nodeOddsRatios.containsKey(n)) {
                    num++;
                    prs += nodeOddsRatios.get(n);
                }
            }
            pathPRS.put(p, prs/num);
        }
        // output
        System.out.println("sample\tlabel\tlogPRS");
        for (Path p : pathPRS.keySet()) {
            System.out.println(p.name+"\t"+p.label+"\t"+pathPRS.get(p));
        }
    }

    /**
     * Compute the path through the graph for the given sample in the given VCF file
     */
    public static void computePathFromVCF(PangenomicGraph graph, String vcfFilename, Map<String,String> samples) {
	System.err.println("STUB: computePathFromVCF("+graph.name+","+vcfFilename+","+samples+")");
    }

    /**
     * Compute the path through the graph for the given sample in the given LIST file
     */
    public static void computePathFromList(PangenomicGraph graph, String listFilename, Map<String,String> samples) throws FileNotFoundException, IOException {
	// get the desired samples' NodeSets
	ListImporter importer = new ListImporter();
	importer.verbose = true;
	importer.read(new File(listFilename), samples.keySet());
	// output
	for (String sampleName : samples.keySet()) {
	    NodeSet nodeSet = importer.sampleNodeSets.get(sampleName);
	    Path path = new Path(graph, new LinkedList<Node>(nodeSet), sampleName, samples.get(sampleName));
	    System.out.println(path.toString());
	}
    }
}
