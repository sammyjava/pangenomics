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
import java.util.concurrent.ConcurrentHashMap;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import org.jgrapht.graph.AbstractGraph;
import org.jgrapht.graph.AsSubgraph;

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
	// graph nodes in TXT file
	Option nodesFileOption = new Option("nodes", "nodesfile", true, "load graph nodes from this nodes.txt file");
	nodesFileOption.setRequired(true);
	options.addOption(nodesFileOption);
	// graph paths in TXT file
	Option pathsFileOption = new Option("paths", "pathsfile", true, "load graph paths from this paths.txt file");
	pathsFileOption.setRequired(true);
	options.addOption(pathsFileOption);
	// FILTER min MGF
	Option minMGFOption = new Option("minmgf", "minmgf", true, "minimum MGF of loci to be included from a VCF/List file [0.01]");
	minMGFOption.setRequired(false);
	options.addOption(minMGFOption);
	// FILTER max MGF
	Option maxMGFOption = new Option("maxmgf", "maxmgf", true, "maximum MGF of loci to be included from a VCF/List file [1.0]");
	maxMGFOption.setRequired(false);
	options.addOption(maxMGFOption);
        // actions
        Option prsOption = new Option("prs", "computeprs", false, "compute polygenic risk scores");
        prsOption.setRequired(false);
        options.addOption(prsOption);
	//
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
	Option printNodePathsFileOption = new Option("nodepaths", "printnodepathsfile", false, "print out node paths file");
	printNodePathsFileOption.setRequired(false);
	options.addOption(printNodePathsFileOption);
	//
	Option printSubgraphOption = new Option("sg", "subgraph", false, "print out a subgraph; requires startnode and endnode");
	printSubgraphOption.setRequired(false);
	options.addOption(printSubgraphOption);
	//
	Option startNodeOption = new Option("sn", "startnode", true, "starting node for subgraph");
	startNodeOption.setRequired(false);
	options.addOption(startNodeOption);
	//
	Option endNodeOption = new Option("en", "endnode", true, "ending node for subgraph");
	endNodeOption.setRequired(false);
	options.addOption(endNodeOption);

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

	// parameters
	double minMGF = 0.01;
	if (cmd.hasOption("minmgf")) {
	    minMGF = Double.parseDouble(cmd.getOptionValue("minmgf"));
	}
	long startNodeId = 0;
	long endNodeId = 0;
	if (cmd.hasOption("subgraph")) {
	    if (!cmd.hasOption("startnode") || !cmd.hasOption("endnode")) {
		System.err.println("ERROR: startnode and endnode are required by subgraph option.");
		System.exit(1);
	    }
	    try {
		startNodeId = Long.parseLong(cmd.getOptionValue("startnode"));
		endNodeId = Long.parseLong(cmd.getOptionValue("endnode"));
	    } catch (NumberFormatException ex) {
		System.err.println("ERROR: startnode and/or endnode are not integers.");
		System.exit(1);
	    }
	}
	
        // build our PangenomicGraph
        PangenomicGraph graph = new PangenomicGraph(cmd.getOptionValue("graph"));
	graph.verbose = true;
	graph.loadNodesFromTXT(new File(cmd.getOptionValue("nodesfile")));
	graph.loadPathsFromTXT(new File(cmd.getOptionValue("pathsfile")));
	graph.tallyLabelCounts();
	System.err.println(graph.name+" has "+graph.vertexSet().size()+" nodes and "+graph.paths.size()+" paths: "+
			   graph.labelCounts.get("case")+"/"+graph.labelCounts.get("ctrl")+" cases/controls");

        // actions
	if (cmd.hasOption("printnodepathsfile")) {
	    String nodePathsFilename = graph.name+".nodepaths.txt";
	    // if (cmd.hasOption("pathsfile")) nodePathsFilename = cmd.getOptionValue("pathsfile")+".nodepaths.txt";
	    if (graph.verbose) System.err.println("Writing node paths file "+nodePathsFilename);
	    printNodePaths(graph, new PrintStream(nodePathsFilename));
	} else if (cmd.hasOption("printpcafile")) {
	    String pcaFilename = graph.name+".pathpca.txt";
	    // if (cmd.hasOption("pathsfile")) pcaFilename = cmd.getOptionValue("pathsfile")+".pathpca.txt";
	    if (graph.verbose) System.err.println("Writing path PCA file "+pcaFilename);
	    printPcaData(graph, new PrintStream(pcaFilename));
	} else if (cmd.hasOption("printarfffile")) {
	    String arffFilename = graph.name+".arff";
	    // if (cmd.hasOption("pathsfile")) arffFilename = cmd.getOptionValue("pathsfile")+".arff";
	    if (graph.verbose) System.err.println("Writing path ARFF file "+arffFilename);
	    printArffData(graph, new PrintStream(arffFilename));
	} else if (cmd.hasOption("printsvmfile")) {
	    String svmFilename = graph.name+".svm.txt";
	    // if (cmd.hasOption("pathsfile")) svmFilename = cmd.getOptionValue("pathsfile")+".svm.txt";
	    if (graph.verbose) System.err.println("Writing path SVM file "+svmFilename);
	    printSvmData(graph, new PrintStream(svmFilename));
	} else if (cmd.hasOption("computeprs")) {
            computePRS(graph, minMGF);
	} else if (cmd.hasOption("subgraph")) {
	    System.out.println("Computing subgraph over node range "+startNodeId+"-"+endNodeId+".");
	    // form Set of included nodes
	    Set<Node> includedNodes = new TreeSet<>();
	    for (long id : graph.nodeIdMap.keySet()) {
		if (id>=startNodeId && id<=endNodeId) {
		    includedNodes.add(graph.nodeIdMap.get(id));
		}
	    }
	    // create the subgraph
	    AbstractGraph<Node,Edge> subgraph = new AsSubgraph<Node,Edge>(graph, includedNodes);
	    // output the subgraph nodes file
	    PrintStream nodeStream = new PrintStream(graph.name+".subgraph.nodes.txt");
	    for (Node n : subgraph.vertexSet()) {
		nodeStream.println(PangenomicGraph.toString(n));
	    }
	    System.out.println("Printed subgraph nodes to "+graph.name+".subgraph.nodes.txt");
	    // output the subgraph paths file
	    PrintStream pathStream = new PrintStream(graph.name+".subgraph.paths.txt");
	    for (Path p : graph.paths) {
		Path subpath = getSubpath(graph, p, startNodeId, endNodeId);
		pathStream.println(subpath.toString());
	    }
	    System.out.println("Printed subgraph paths to "+graph.name+".subgraph.paths.txt");
	}
    }

    /**
     * Compute the polygenic risk scores from a graph, tossing nodes with gf<minMGF.
     */
    public static void computePRS(PangenomicGraph graph, double minMGF) {
        // calculate map of log odds ratio
        Map<Node,Double> nodeOddsRatios = new HashMap<>();
	int numTossed = 0;
        for (Node n : graph.getNodes()) {
            if (n.gf<minMGF) {
		numTossed++;
	    } else {
		double OR = graph.oddsRatio(n);
		if (OR>Double.MIN_VALUE && OR<Double.MAX_VALUE) {
		    double logOR = Math.log(OR);
		    nodeOddsRatios.put(n, logOR);
		}
            }
        }
	System.err.println("Tossed "+numTossed+" nodes with gf<"+minMGF);
        // sum log odds ratio over each path, storing in a map
        Map<Path,Double> pathPRS = new HashMap<>();
        for (Path p : graph.paths) {
            double prs = 0.0;
            int num = 0;
            for (Node n : p.getNodes()) {
                if (nodeOddsRatios.containsKey(n)) {
                    num++;
		    double logOR = nodeOddsRatios.get(n);
                    prs += logOR;
                }
            }
	    prs = prs/num;
            pathPRS.put(p, prs);
        }
        // per-subject PRS output
        System.out.println("sample\tlabel\tlogPRS");
        for (Path p : pathPRS.keySet()) {
            System.out.println(p.getName()+"\t"+p.getLabel()+"\t"+pathPRS.get(p));
        }
	// summary stats
	int cases = 0;
	int controls = 0;
	int TP = 0;
	int TN = 0;
	int FP = 0;
	int FN = 0;
	for (Path p : pathPRS.keySet()) {
	    double prs = pathPRS.get(p);
	    if (p.getLabel().equals("case")) {
		cases++;
		if (prs>0.0) {
		    TP++;
		} else {
		    FN++;
		}
	    } else {
		controls++;
		if (prs<0.) {
		    TN++;
		} else {
		    FP++;
		}
	    }
	}
	int correct = TP + TN;
	double accuracy = (double)correct / (double)(cases+controls);
	double TPR = (double)TP / (double)cases;
	double FPR = (double)FP / (double)controls;
	double MCC = getMCC(TP, TN, FP, FN);
	System.err.println(correct+"\t"+accuracy+"\t"+TPR+"\t"+FPR+"\t"+MCC);
    }

    /**
     * Print ARFF of node participation by path, for Weka analysis.
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
	Set<String> arffData = ConcurrentHashMap.newKeySet();
	graph.paths.parallelStream().forEach(path -> {
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
     *
     * path1 case 1:1 2:1 3:1 4:0 ... (tab separated, don't skip any indexes)
     * path2 ctrl 1:0 2:0 3:0 4:1 ...
     *
     * This is similar, but not identical to, the SVMlight format.
     */
    public static void printSvmData(PangenomicGraph graph, PrintStream out) {
	///////////////////////////////////////////////////////////////////////////////////////
	Set<String> svmData = ConcurrentHashMap.newKeySet();
	graph.paths.parallelStream().forEach(path -> {
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
	Map<Node,String> pcaData = new ConcurrentHashMap<>();
	graph.getNodes().parallelStream().forEach(node -> {
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

    /**
     * Return the MCC from a given set of int counts.
     */
    public static double getMCC(int TP, int TN, int FP, int FN) {
	return (double)(TP*TN-FP*FN) / Math.sqrt((double)(TP+FP)*(double)(TP+FN)*(double)(TN+FP)*(double)(TN+FN));
    }

    /**
     * Return the MCC from a given set of rates, assuming equal cases and controls.
     */
    public static double getMCC(double TPR, double FPR) {
	double TNR = 1.0 - FPR;
	double FNR = 1.0 - TPR;
	return (TPR*TNR-FPR*FNR) / Math.sqrt((TPR+FPR)*(TPR+FNR)*(TNR+FPR)*(TNR+FNR));
    }

    /**
     * Return the subpath inclusively between the two node IDs, even if the path doesn't contain one and/or the other node.
     * @param start the starting Node id
     * @param finish the ending Node id
     * @return the subpath inclusively between start and finish
     */
    public static Path getSubpath(PangenomicGraph g, Path p, long start, long finish) {
        List<Node> subnodes = new LinkedList<>();
	for (Node node : p.getNodes()) {
	    if (node.id>=start && node.id<=finish) {
		subnodes.add(node);
	    }
        }
        return new Path(g, subnodes, p.getSample());
    }
}
