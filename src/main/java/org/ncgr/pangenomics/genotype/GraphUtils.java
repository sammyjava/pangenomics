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
	// graph nodes in TXT file
	Option nodesFileOption = new Option("nodes", "nodesfile", true, "load graph nodes from this nodes.txt file");
	nodesFileOption.setRequired(true);
	options.addOption(nodesFileOption);
	// graph paths in TXT file
	Option pathsFileOption = new Option("paths", "pathsfile", true, "load graph paths from this paths.txt file");
	pathsFileOption.setRequired(true);
	options.addOption(pathsFileOption);
	// FILTER min MAF
	Option minMAFOption = new Option("minmaf", "minmaf", true, "minimum MAF/MGF of loci to be included from a VCF/List file [0.0]");
	minMAFOption.setRequired(false);
	options.addOption(minMAFOption);
	// FILTER max MAF
	Option maxMAFOption = new Option("maxmaf", "maxmaf", true, "maximum MAF/MGF of loci to be included from a VCF/List file [1.0]");
	maxMAFOption.setRequired(false);
	options.addOption(maxMAFOption);
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
            double minMGF = 1e-2;
            if (cmd.hasOption("minmaf")) {
                minMGF = Double.parseDouble(cmd.getOptionValue("minmaf"));
            }
            computePRS(graph, minMGF);
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
