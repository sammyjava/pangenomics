package org.ncgr.pangenomics.genotype;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;

import java.util.List;
import java.util.LinkedList;
import java.util.Map;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.Set;
import java.util.HashSet;
import java.util.TreeSet;

import java.util.concurrent.ConcurrentSkipListSet;
import java.util.concurrent.ConcurrentSkipListMap;

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
	//
	Option minMAFOption = new Option("minmaf", "minmaf", true, "minimum MAF of loci to be included from LIST or VCF file [0.0]");
	minMAFOption.setRequired(false);
	options.addOption(minMAFOption);
	//
	Option maxMAFOption = new Option("maxmaf", "maxmaf", true, "maximum MAF of loci to be included from LIST or VCF file [1.0]");
	maxMAFOption.setRequired(false);
	options.addOption(maxMAFOption);

        // actions
        Option prsOption = new Option("prs", "prs", false, "compute polygenic risk scores");
        prsOption.setRequired(false);
        options.addOption(prsOption);
	//
	Option sampleFileOption = new Option("samplefile", true, "file containing samples+labels to determine paths through the given graph (-g) from a VCF (-vcf) or LIST (-list)");
	sampleFileOption.setRequired(false);
	options.addOption(sampleFileOption);

	// output
	Option printPcaFileOption = new Option("pca", "printpcafile", false, "print out path pca file");
	printPcaFileOption.setRequired(false);
	options.addOption(printPcaFileOption);
	//
	Option printArffFileOption = new Option("arff", "printarfffile", false, "print out an ARFF file of path node participation");
	printArffFileOption.setRequired(false);
	options.addOption(printArffFileOption);
	//
	Option printSvmFileOption = new Option("svm", "printsvmfile", false, "print out a LIBSVM-compatible file of path node participation");
	printSvmFileOption.setRequired(false);
	options.addOption(printSvmFileOption);
	//
	Option printNodePathsFileOption = new Option("nodepaths", "printnodepathsfile", false, "print out node paths file [false]");
	printNodePathsFileOption.setRequired(false);
	options.addOption(printNodePathsFileOption);


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
        PangenomicGraph graph = new PangenomicGraph(cmd.getOptionValue("graph"));
	graph.verbose = true;

        // load graph from TXT files
	graph.loadPathsFromTXT(graph.getNodesFile(), graph.getPathsFile());
        graph.tallyLabelCounts();
        System.err.println(graph.name+" has "+graph.vertexSet().size()+" nodes and "+graph.paths.size()+" paths: "+graph.labelCounts.get("case")+"/"+graph.labelCounts.get("ctrl")+" cases/controls");

        // options
        boolean computePRS = cmd.hasOption("prs");
	boolean computeSamplePaths = cmd.hasOption("samplefile");

	if (cmd.hasOption("printnodepathsfile")) {
	    String nodePathsFilename = graph.name+".nodepaths.txt";
	    // if (cmd.hasOption("pathsfile")) nodePathsFilename = cmd.getOptionValue("pathsfile")+".nodepaths.txt";
	    if (graph.verbose) System.err.println("Writing node paths file "+nodePathsFilename);
	    printNodePaths(graph, new PrintStream(nodePathsFilename));
	}

	if (cmd.hasOption("printpcafile")) {
	    String pcaFilename = graph.name+".pathpca.txt";
	    // if (cmd.hasOption("pathsfile")) pcaFilename = cmd.getOptionValue("pathsfile")+".pathpca.txt";
	    if (graph.verbose) System.err.println("Writing path PCA file "+pcaFilename);
	    printPcaData(graph, new PrintStream(pcaFilename));
	}

	if (cmd.hasOption("printarfffile")) {
	    String arffFilename = graph.name+".arff";
	    // if (cmd.hasOption("pathsfile")) arffFilename = cmd.getOptionValue("pathsfile")+".arff";
	    if (graph.verbose) System.err.println("Writing path ARFF file "+arffFilename);
	    printArffData(graph, new PrintStream(arffFilename));
	}
	
	if (cmd.hasOption("printsvmfile")) {
	    String svmFilename = graph.name+".svm.txt";
	    // if (cmd.hasOption("pathsfile")) svmFilename = cmd.getOptionValue("pathsfile")+".svm.txt";
	    if (graph.verbose) System.err.println("Writing path SVM file "+svmFilename);
	    printSvmData(graph, new PrintStream(svmFilename));
	}

	// GraphUtils???
        // build the node-paths map
        // graph.buildNodePaths();

        // actions
        if (computePRS) {
            double minMGF = 1e-2;
            if (cmd.hasOption("mgf")) {
                minMGF = Double.parseDouble(cmd.getOptionValue("mgf"));
            }
            computePRS(graph, minMGF);
        } else if (computeSamplePaths) {
	    TreeSet<Sample> samples = Sample.readSamples(new File(cmd.getOptionValue("samplefile")));
	    if (cmd.hasOption("vcf")) {
		double minMAF = 0.0;
		double maxMAF = 1.0;
		if (cmd.hasOption("minmaf")) minMAF = Double.parseDouble(cmd.getOptionValue("minmaf"));
		if (cmd.hasOption("maxmaf")) maxMAF = Double.parseDouble(cmd.getOptionValue("maxmaf"));
		computePathsFromVCF(graph, cmd.getOptionValue("vcf"), samples, minMAF, maxMAF);
	    } else if (cmd.hasOption("list")) {
		computePathsFromList(graph, cmd.getOptionValue("list"), samples);
	    } else {
		System.err.println("ERROR: you must provide a VCF file (-vcf) or LIST file (-list) that contains the given samples");
	    }			
	}
    }

    /**
     * Compute the polygenic risk scores from a graph, tossing nodes with gf<minMGF.
     */
    public static void computePRS(PangenomicGraph graph, double minMGF) {
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
            System.out.println(p.getName()+"\t"+p.getLabel()+"\t"+pathPRS.get(p));
        }
    }

    /**
     * Compute the paths through the graph for the given samples in the given VCF file
     */
    public static void computePathsFromVCF(PangenomicGraph graph, String vcfFilename, TreeSet<Sample> samples, double minMAF, double maxMAF) throws FileNotFoundException, IOException {
        VCFImporter importer = new VCFImporter(new File(vcfFilename));
	importer.verbose = true;
        importer.readNodes(minMAF, maxMAF);
	importer.readPaths(samples);
	// make sure the imported nodes match the graph's nodes for the same id (this will not be true if MAF is different)
	checkImportedNodes(graph, importer.nodes);
	// output
	outputPaths(graph, importer.sampleNodeSets);
    }

    /**
     * Compute the paths through the graph for the given samples in the given LIST file
     */
    public static void computePathsFromList(PangenomicGraph graph, String listFilename, TreeSet<Sample> samples) throws FileNotFoundException, IOException {
	// ListImporter importer = new ListImporter(new File(listFilename));
	// importer.verbose = true;
	// importer.readNodes();
	// importer.readPaths(samples);
	// // make sure the imported nodes match the graph's nodes for the same id (this will not be true if MAF is different)
	// checkImportedNodes(graph, importer.nodes);
	// // output
	// outputPaths(graph, importer.sampleNodeSets);
    }

    /**
     * Check that the imported nodes and graph nodes are identical
     */
    public static void checkImportedNodes(PangenomicGraph graph, Map<Long,Node> importedNodes) {
	for (long id : importedNodes.keySet()) {
	    Node importedNode = importedNodes.get(id);
	    Node graphNode = graph.getNode(id);
	    if (graphNode==null) {
		System.err.println("ERROR: imported node "+importedNode+" is not present in the graph.");
		System.exit(1);
	    } else if (!importedNode.equals(graphNode)) {
		System.err.println("ERROR: imported node "+importedNode+" does not equal the corresponding graph node "+graphNode);
		System.exit(1);
	    }
	}
    }

    /**
     * Output the paths
     */
    public static void outputPaths(PangenomicGraph graph, TreeMap<Sample,NodeSet> sampleNodeSets) {
    	for (Sample sample : sampleNodeSets.keySet()) {
	    NodeSet nodeSet = sampleNodeSets.get(sample);
	    Path path = new Path(graph, new LinkedList<Node>(nodeSet), sample);
	    System.out.println(path.toString());
	}
    }

    /**
     * Print ARFF of node participation by path, for Weka analysis.
     * NOTE: no-call nodes are dropped.
     *
     * @RELATION iris
     *
     * @ATTRIBUTE ID           STRING
     * @ATTRIBUTE sepallength  NUMERIC
     * @ATTRIBUTE sepalwidth   NUMERIC
     * @ATTRIBUTE petallength  NUMERIC
     * @ATTRIBUTE petalwidth   NUMERIC
     * @ATTRIBUTE class        {Iris-setosa,Iris-versicolor,Iris-virginica}
     *
     * @DATA
     * 5.1,3.5,1.4,0.2,Iris-setosa
     * 4.9,3.0,1.4,0.2,Iris-virginica
     * 4.7,3.2,1.3,0.2,Iris-versicolor
     * 4.6,3.1,1.5,0.2,Iris-setosa
     * 5.0,3.6,1.4,0.2,Iris-viginica
     */
    public static void printArffData(PangenomicGraph graph, PrintStream out) {
	out.println("@RELATION "+graph.name);
        out.println("");
	out.println("@ATTRIBUTE ID STRING"); // path ID
	// attributes: each node is labeled Nn where n is the node ID
	for (Node node : graph.getNodes()) {
	    if (!node.isNoCall()) {
		out.println("@ATTRIBUTE N"+node.id+" NUMERIC");
	    }
        }
        // add the path label class attribute
        out.println("@ATTRIBUTE class {ctrl,case}");
        out.println("");
        // data
        out.println("@DATA");
	///////////////////////////////////////////////////////////////////////////////////////
	ConcurrentSkipListSet<Path> concurrentPaths = new ConcurrentSkipListSet<>(graph.paths);
	ConcurrentSkipListSet<String> arffData = new ConcurrentSkipListSet<>();
	concurrentPaths.parallelStream().forEach(path -> {
		String arff = path.getName();
		for (Node node : graph.getNodes()) {
		    if (path.traverses(node)) {
			arff += ",1";
		    } else {
			arff += ",0";
		    }
		}
		arff += ","+path.getLabel();
		arffData.add(arff);
	    });
	///////////////////////////////////////////////////////////////////////////////////////
	for (String arff : arffData) {
	    out.println(arff);
	}
        out.close();
    }

    /**
     * Print the labeled path node participation for LIBSVM analysis.
     * NOTE: no-call nodes are dropped.
     *
     * path1 case 1:1 2:1 3:1 4:0 ... (tab separated, don't skip any indexes)
     * path2 ctrl 1:0 2:0 3:0 4:1 ...
     *
     * This is similar, but not identical to, the SVMlight format.
     */
    public static void printSvmData(PangenomicGraph graph, PrintStream out) {
	///////////////////////////////////////////////////////////////////////////////////////
	ConcurrentSkipListSet<Path> concurrentPaths = new ConcurrentSkipListSet<>(graph.paths);
	ConcurrentSkipListSet<String> svmData = new ConcurrentSkipListSet<>();
	concurrentPaths.parallelStream().forEach(path -> {
		String svm = path.getName()+"\t"+path.getLabel();
		int n = 0;
		for (Node node : graph.getNodes()) {
		    n++;
		    if (path.traverses(node)) {
			svm += "\t"+n+":+1";
		    } else {
			svm += "\t"+n+":-1";
		    }
		}
		svmData.add(svm);
	    });
	///////////////////////////////////////////////////////////////////////////////////////
	for (String svm : svmData) {
	    out.println(svm);
	}
        out.close();
    }

    /**
     * Print node participation by path, appropriate for PCA analysis.
     */
    public static void printPcaData(PangenomicGraph graph, PrintStream out) throws FileNotFoundException, IOException {
        StringBuilder headerBuilder = new StringBuilder();
        // header is paths
	List<Path> pathList = new LinkedList<>(graph.paths);
	headerBuilder.append(pathList.get(0).getName());
	for (int i=1; i<pathList.size(); i++) {
	    headerBuilder.append("\t"+pathList.get(i).getName());
	}
        out.println(headerBuilder.toString());
	////////////////////////////////////////////////////////////////////////////////////////////
	// rows are nodes and counts of path support of each node
	ConcurrentSkipListSet<Node> concurrentNodes = new ConcurrentSkipListSet<>(graph.getNodes());	
	ConcurrentSkipListMap<Node,String> pcaData = new ConcurrentSkipListMap<>();
	concurrentNodes.parallelStream().forEach(node -> {
		StringBuilder lineBuilder = new StringBuilder();
		lineBuilder.append("N"+node.id);
		// spin through every path, printing 0/1 if path doesn't/does traverse this node
		for (Path path : graph.paths) {
		    if (path.traverses(node)) {
			lineBuilder.append("\t1");
		    } else {
			lineBuilder.append("\t0");
		    }
		}
		pcaData.put(node, lineBuilder.toString());
	    });
	////////////////////////////////////////////////////////////////////////////////////////////
	for (String line : pcaData.values()) {
	    out.println(line);
	}
    }

    /**
     * Print out the node paths along with counts.
     */
    public static void printNodePaths(PangenomicGraph graph, PrintStream out) throws FileNotFoundException, IOException {
        for (Node node : graph.getNodePaths().keySet()) {
	    List<Path> pathList = graph.getNodePaths().get(node);
	    if (pathList.size()>0) {
		StringBuilder builder = new StringBuilder();
		builder.append(node.id);
		for (Path path : pathList) {
		    builder.append("\t"+path.getName());
		}
		out.println(builder.toString());
	    }
        }
    }
}
