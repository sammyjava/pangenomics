package org.ncgr.pangenomics.genotype.fr;

import org.ncgr.pangenomics.genotype.*;

import org.jgrapht.GraphPath;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;

import java.text.DecimalFormat;

import java.util.List;
import java.util.LinkedList;
import java.util.Map;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.TreeMap;
import java.util.Optional;
import java.util.Properties;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentSkipListSet;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

/**
 * Static methods to read FrequentedRegions and do the other things.
 *
 * @author Sam Hokin
 */
public class FRUtils {    

    /**
     * Read FRs from the text files written by a FRFinder run, given the graph as well.
     * This constructor populates the FRs fully, with subpaths and the whole bit.
     *
     * nodes   size    support case    ctrl    OR      p       pri
     * [1341]  1       21      19      2       9.500   9.48E-3 202
     * 509678.ctrl:[18,20,21,23,24,26,27,29,30,33,34]
     * 628863.case:[18,20,21,23,24,26,27,29,30,33,34]
     * etc.
     */
    public static TreeSet<FrequentedRegion> readFrequentedRegions(String inputPrefix, PangenomicGraph graph) throws FileNotFoundException, IOException {
        // get alpha, kappa from the input prefix
        double alpha = readAlpha(inputPrefix);
        int kappa = readKappa(inputPrefix);
        // parse the priorityOption and impose defaults
        String priorityOption = readPriorityOption(inputPrefix);
        String[] priorityParts = priorityOption.split(":");
        int priorityOptionKey = Integer.parseInt(priorityParts[0]);
        String priorityOptionLabel = null;
        if (priorityParts.length>1) {
            String priorityOptionParameter = priorityParts[1];
            if (priorityParts[1].equals("case") || priorityParts[1].equals("ctrl")) {
                priorityOptionLabel = priorityParts[1];
            }
        }
        if (priorityOptionKey==1 && priorityOptionLabel==null) priorityOptionLabel = "case";
        if (priorityOptionKey==3 && priorityOptionLabel==null) priorityOptionLabel = "case";
        // get the graph from the nodes and paths files
        File nodesFile = new File(getNodesFilename(inputPrefix));
        File pathsFile = new File(getPathsFilename(inputPrefix));
        // create a node map for building subpaths
        Map<Long,Node> nodeMap = new HashMap<>();
        for (Node n : graph.getNodes()) {
            nodeMap.put(n.id, n);
        }
        // read FRs
	Map<NodeSet,FrequentedRegion> inputFRMap = new HashMap<>();
	for (FrequentedRegion fr : readFrequentedRegions(inputPrefix)) {
	    fr.graph = graph;
	    inputFRMap.put(fr.nodes, fr);
	}
	TreeSet<FrequentedRegion> sortedFRs = new TreeSet<>();
        String subpathsFilename = getFRSubpathsFilename(inputPrefix);
	File subpathsFile = new File(subpathsFilename);
	if (subpathsFile.exists()) {
	    // read the subpaths, only using subpaths matching those in the FR text file (which may have been trimmed)
	    BufferedReader reader = new BufferedReader(new FileReader(subpathsFilename));
	    String line = null;
	    while ((line=reader.readLine())!=null) {
		// 0                1       2       3           4           5     6       7
		// nodes            size    support caseSupport ctrlSupport OR    p       priority
		// [1341]           1       21      19          2           9.500 9.48E-3 202
		// 429266.case:1341]
		// 158005.case:1341]
		// ...
		String[] fields = line.split("\t");
		NodeSet nodes = graph.getNodeSet(fields[0]);
		int size = Integer.parseInt(fields[1]);
		int support = Integer.parseInt(fields[2]);
		int caseSupport = Integer.parseInt(fields[3]);
		int ctrlSupport = Integer.parseInt(fields[4]);
		double or = Double.POSITIVE_INFINITY;
		try {
		    or = Double.parseDouble(fields[5]);
		} catch (NumberFormatException e) {
		    // do nothing, it's an infinity symbol
		}
		double p = Double.parseDouble(fields[6]);
		int priority = Integer.parseInt(fields[7]);
		List<Path> subpaths = new LinkedList<>();
		for (int i=0; i<support; i++) {
		    // have to read lines to run through the file even if skipping this FR
		    line = reader.readLine();
		    if (inputFRMap.containsKey(nodes)) {
			String[] parts = line.split(":");
			String pathFull = parts[0];
			String nodeString = parts[1];
			// split out the name, label, nodes
			String[] nameParts = pathFull.split("\\.");
			String name = nameParts[0];
			String label = null;
			if (nameParts.length>1) label = nameParts[1];
			List<Node> subNodes = new LinkedList<>();
			String[] nodesAsStrings = nodeString.replace("[","").replace("]","").split(",");
			for (String nodeAsString : nodesAsStrings) {
			    long nodeId = Long.parseLong(nodeAsString);
			    subNodes.add(nodeMap.get(nodeId));
			}
			// add to the subpaths
			subpaths.add(new Path(graph, subNodes, name, label));
		    }
		}
		if (inputFRMap.containsKey(nodes)) {
		    FrequentedRegion fr = new FrequentedRegion(graph, nodes, subpaths, alpha, kappa, priorityOptionKey, priorityOptionLabel, support);
		    sortedFRs.add(fr);
		}
	    }
	} else {
	    // we don't have a subpaths file, so we have to build the subpaths for each FR, in parallel
	    System.err.println("Building FR subpaths...");
	    ConcurrentSkipListSet<FrequentedRegion> concurrentFRs = new ConcurrentSkipListSet<>();
	    concurrentFRs.addAll(inputFRMap.values());
	    ////////////////////////////////////////////////////////////////////////////////////////////
	    concurrentFRs.parallelStream().forEach(fr -> {
		    fr.updateSupport();
		});
	    ////////////////////////////////////////////////////////////////////////////////////////////
	    sortedFRs.addAll(concurrentFRs);
	}
        return sortedFRs;
    }

    /**
     * Read FRs with only parameters from the text file written by a FRFinder run.
     * This constructor populates the FRs with only their main class variables; no subpaths, for example.
     * 0       1       2       3       4       5       6       7
     * nodes   size    support case    ctrl    OR      p       pri
     * [1341]  1       21      19      2       9.500   9.48E-3 202
     */
    public static TreeSet<FrequentedRegion> readFrequentedRegions(String inputPrefix) throws FileNotFoundException, IOException {
        // get alpha, kappa from the input prefix
        double alpha = readAlpha(inputPrefix);
        int kappa = readKappa(inputPrefix);
        // read the FRs
        TreeSet<FrequentedRegion> sortedFRs = new TreeSet<>();
        String frFilename = getFRsFilename(inputPrefix);
        BufferedReader reader = new BufferedReader(new FileReader(frFilename));
        String line = null;
        while ((line=reader.readLine())!=null) {
	    FrequentedRegion fr = new FrequentedRegion(line);
            sortedFRs.add(fr);
        }
	System.err.println("FRUtils loaded "+sortedFRs.size()+" unique FRs from "+inputPrefix);
        return sortedFRs;
    }

    /**
     * Return alpha from an FRFinder run.
     */
    public static double readAlpha(String inputPrefix) throws FileNotFoundException, IOException {
        String paramsFilename = getParamsFilename(inputPrefix);
        BufferedReader reader = new BufferedReader(new FileReader(paramsFilename));
        String line = null;
        double alpha = 0;
        while ((line=reader.readLine())!=null) {
            if (line.startsWith("#alpha")) {
                String[] parts = line.split("=");
                alpha = Double.parseDouble(parts[1]);
            }
        }
        return alpha;
    }

    /**
     * Return kappa from a an FRFinder run.
     */
    public static int readKappa(String inputPrefix) throws FileNotFoundException, IOException {
        String paramsFilename = getParamsFilename(inputPrefix);
        BufferedReader reader = new BufferedReader(new FileReader(paramsFilename));
        String line = null;
        int kappa = 0;
        while ((line=reader.readLine())!=null) {
            if (line.startsWith("#kappa")) {
                String[] parts = line.split("=");
                kappa = Integer.parseInt(parts[1]);
            }
        }
        return kappa;
    }

    /**
     * Return priorityOption from an FRFinder run.
     */
    public static String readPriorityOption(String inputPrefix) throws FileNotFoundException, IOException {
        String paramsFilename = getParamsFilename(inputPrefix);
        BufferedReader reader = new BufferedReader(new FileReader(paramsFilename));
        Properties parameters = new Properties();
        parameters.load(reader);
        return parameters.getProperty("priorityOption");
    }        

    /**
     * Form the FRs output filename
     */
    public static String getFRsFilename(String prefix) {
        return prefix+".frs.txt";
    }

    /**
     * Form the FRSubpaths output filename
     */
    public static String getFRSubpathsFilename(String prefix) {
        return prefix+".subpaths.txt";
    }

    /**
     * Form the pathFRs output filename
     */
    public static String getPathFRsFilename(String prefix) {
        return prefix+".pathfrs.txt";
    }

    /**
     * Form the SVM version of the pathFRs output filename
     */
    public static String getPathFRsSVMFilename(String prefix) {
        return prefix+".svm.txt";
    }

    /**
     * Form the ARFF version of the pathFRs output filename
     */
    public static String getPathFRsARFFFilename(String prefix) {
        return prefix+".arff";
    }

    /**
     * Form the parameters output filename
     */
    public static String getParamsFilename(String prefix) {
        return prefix+".params.txt";
    }

    /**
     * Form the graph nodes filename
     * if prefix = HTT.1k-1.0-0 then filename = HTT.nodes.txt
     */
    public static String getNodesFilename(String prefix) {
        String[] parts = prefix.split("-");
        return parts[0]+".nodes.txt";
    }

    /**
     * Form the graph paths filename
     * if prefix = HTT.1k-1.0-0 then filename = HTT.paths.txt
     */
    public static String getPathsFilename(String prefix) {
        String[] parts = prefix.split("-");
        return parts[0]+".paths.txt";
    }

    /**
     * Print out the parameters.
     */
    public static void printParameters(Properties parameters, String outputPrefix, double alpha, int kappa, long clockTime) throws IOException {
        PrintStream out = new PrintStream(getParamsFilename(outputPrefix));
        String comments = "alpha="+alpha+"\n"+"kappa="+kappa+"\n"+"clocktime="+formatTime(clockTime);
        parameters.store(out, comments);
        out.close();
    }
    
    /**
     * Read the parameters from a previous run's properties file.
     */
    public static Properties readParameters(String inputPrefix) throws FileNotFoundException, IOException {
        String paramsFilename = getParamsFilename(inputPrefix);
        BufferedReader reader = new BufferedReader(new FileReader(paramsFilename));
        Properties parameters = new Properties();
        parameters.load(reader);
        parameters.setProperty("paramsFile", paramsFilename);
        parameters.setProperty("inputPrefix", inputPrefix);
        return parameters;
    }

    /**
     * Format a time duration given in milliseconds.
     */
    public static String formatTime(long millis) {
        DecimalFormat tf = new DecimalFormat("00"); // hours, minutes, seconds
        long hours = (millis / 1000) / 60 / 60;
	long minutes = (millis / 1000 / 60) % 60;
        long seconds = (millis / 1000) % 60;
	return tf.format(hours)+":"+tf.format(minutes)+":"+tf.format(seconds);
    }

    /**
     * Form an outputPrefix from inputPrefix, minSupport, and minSize.
     */
    public static String formOutputPrefix(String inputPrefix, int minSupport, int minSize) {
        return inputPrefix+"-"+minSupport+"."+minSize;
    }

    /**
     * Post-process a list of FRs read in from the inputPrefix file for given minSupport and minSize.
     * Outputs a new set of FR files with minSupport and minSize in the outputPrefix.
     */
    public static void postprocess(String inputPrefix, int minSupport, int minSize) throws FileNotFoundException, IOException {
	Properties parameters = readParameters(inputPrefix);
	parameters.setProperty("minSupport", String.valueOf(minSupport));
	parameters.setProperty("minSize", String.valueOf(minSize));
	Set<FrequentedRegion> frequentedRegions = readFrequentedRegions(inputPrefix);
	Set<FrequentedRegion> filteredFRs = new TreeSet<>();
        for (FrequentedRegion fr : frequentedRegions) {
            boolean passes = true;
            String reason = "";
            if (fr.support<minSupport) {
                passes = false;
                reason += " support";
            } else {
                reason += " SUPPORT";
            }
            if (fr.size<minSize) {
                passes = false;
                reason += " size";
            } else {
                reason += " SIZE";
            }
            if (passes) filteredFRs.add(fr);
        }
        System.out.println(filteredFRs.size()+" FRs passed minSupport="+minSupport+", minSize="+minSize);
	// output the filtered FRs
	if (filteredFRs.size()>0) {
	    printFrequentedRegions(filteredFRs, formOutputPrefix(inputPrefix, minSupport, minSize));
	}
    }

    /**
     * Print the labeled path FR support for SVM analysis. Lines are like:
     *
     * path1 case 1:1 2:1 3:1 4:0 ...
     * path2 case 1:0 2:0 3:0 4:1 ...
     * path3 ctrl 1:0 2:1 3:0 4:2 ...
     *
     * which is similar, but not identical to, the SVMlight format.
     *
     * If numCasePaths>0 or numCtrlPaths>0, randomly select numCasePaths cases and/or numCtrlPaths controls.
     * This also allows pruning on size, support, and p-value.
     */
    public static void printPathFRsSVM(String inputPrefix, int numCasePaths, int numCtrlPaths, int minSize, int minSupport, double maxPValue, int minPriority) throws IOException {
        // read the graph and FRs from files
        PangenomicGraph graph = readGraph(inputPrefix, numCasePaths, numCtrlPaths);
	TreeSet<FrequentedRegion> frequentedRegions = readFrequentedRegions(inputPrefix, graph);
	frequentedRegions = pruneFrequentedRegions(frequentedRegions, minSize, minSupport, maxPValue, minPriority);
	System.err.println("FRUtils: "+frequentedRegions.size()+" FRs to be printed to SVM file.");
	// collect the paths, cases and controls
        ConcurrentSkipListSet<Path> concurrentPaths = buildConcurrentPaths(graph, numCasePaths, numCtrlPaths);
	///////////////////////////////////////////////////////
	// build the SVM strings in parallel
	ConcurrentHashMap<String,String> pathSVM = new ConcurrentHashMap<>();
        ConcurrentSkipListSet<FrequentedRegion> concurrentFRs = new ConcurrentSkipListSet<>(frequentedRegions);
	concurrentPaths.parallelStream().forEach(path -> {
		String svm = path.label;
		int c  = 0;
		for (FrequentedRegion fr : concurrentFRs) {
		    c++;
		    svm += "\t"+c+":"+fr.countSubpathsOf(path);
		}
		pathSVM.put(path.name,svm);
	    });
	///////////////////////////////////////////////////////
	// sort by SVM string
	LinkedHashMap<String,String> sortedPathSVM = new LinkedHashMap<>();
	pathSVM.entrySet().stream().sorted(Map.Entry.comparingByValue()).forEachOrdered(x->sortedPathSVM.put(x.getKey(),x.getValue()));
	// output: no header, one path per row
        PrintStream out = new PrintStream(getPathFRsSVMFilename(inputPrefix));
        for (String name : sortedPathSVM.keySet()) {
	    String paddedName = name;
	    if (name.length()<6) paddedName = "0" + name; // special hack for 6-digit samples
            out.println(paddedName+"\t"+sortedPathSVM.get(name));
        }
        out.close();
    }
    
    /**
     * Print the (unlabeled) path FR support in ARFF format for the given number of case and control paths (0 for all).
     * The rows are sorted by case/control and then FR support vector.
     *
     * @RELATION iris
     *
     * @ATTRIBUTE ID           STRING
     * @ATTRIBUTE sepallength  NUMERIC
     * @ATTRIBUTE sepalwidth   NUMERIC
     * @ATTRIBUTE petallength  NUMERIC
     * @ATTRIBUTE petalwidth   NUMERIC
     * @ATTRIBUTE class        {Iris-setosa,Iris-versicolor,Iris-virginica}
     
     * @DATA
     * 5.1,3.5,1.4,0.2,Iris-setosa
     * 4.9,3.0,1.4,0.2,Iris-virginica
     * 4.7,3.2,1.3,0.2,Iris-versicolor
     * 4.6,3.1,1.5,0.2,Iris-setosa
     * 5.0,3.6,1.4,0.2,Iris-viginica
     *
     *
     * This also allows pruning on size, support, and p-value.
     */
    public static void printPathFRsARFF(String inputPrefix, int numCasePaths, int numCtrlPaths, int minSize, int minSupport, double maxPValue, int minPriority) throws IOException {
        // load the graph
	PangenomicGraph graph = readGraph(inputPrefix, numCasePaths, numCtrlPaths);
	// load the frequented regions and update support in parallel
	TreeSet<FrequentedRegion> frequentedRegions = readFrequentedRegions(inputPrefix, graph);
	frequentedRegions = pruneFrequentedRegions(frequentedRegions, minSize, minSupport, maxPValue, minPriority);
	System.err.println("FRUtils: "+frequentedRegions.size()+" FRs to be printed to ARFF file.");
	// collect the paths, cases and controls
        ConcurrentSkipListSet<Path> concurrentPaths = buildConcurrentPaths(graph, numCasePaths, numCtrlPaths);
	/////////////////////////////////////////////////////////
	// now load pathARFF for selected paths in parallel
	ConcurrentHashMap<String,String> pathARFF = new ConcurrentHashMap<>();
        ConcurrentSkipListSet<FrequentedRegion> concurrentFRs = new ConcurrentSkipListSet<>(frequentedRegions);
	concurrentPaths.parallelStream().forEach(path -> {
		String arff = "";
		for (FrequentedRegion fr : concurrentFRs) {
		    arff += fr.countSubpathsOf(path)+",";
		}
		arff += path.label;
		pathARFF.put(path.name, arff);
	    });
	/////////////////////////////////////////////////////////
	// sort by ARFF string
	LinkedHashMap<String,String> sortedPathARFF = new LinkedHashMap<>();
	pathARFF.entrySet().stream().sorted(Map.Entry.comparingByValue()).forEachOrdered(x->sortedPathARFF.put(x.getKey(),x.getValue()));
	// ARFF output	
        PrintStream out = new PrintStream(FRUtils.getPathFRsARFFFilename(inputPrefix));
        out.println("@RELATION "+inputPrefix);
        out.println("");
        // attributes: path ID
        out.println("@ATTRIBUTE ID STRING");
        // attributes: each FR is labeled FRn
        int c = 0;
        for (FrequentedRegion fr : frequentedRegions) {
            c++;
            String frLabel = "FR"+c;
            out.println("@ATTRIBUTE "+frLabel+" NUMERIC");
        }
        // add the class attribute
        out.println("@ATTRIBUTE class {ctrl,case}");
        out.println("");
        // data
        out.println("@DATA");
        for (String name : sortedPathARFF.keySet()) {
	    String paddedName = name;
	    if (name.length()<6) paddedName = "0" + name; // special hack for 6-digit samples
            out.println(paddedName+","+sortedPathARFF.get(name));
        }
        out.close();
    }

    /**
     * Print the (labeled) path FR support in TXT format for the given number of case and control paths (0 for all).
     * The rows are sorted by case/control and then FR support vector.
     *
     * sample1.case sample2.ctrl sample3.ctrl ... sampleN.case
     * FR1 0        3            2                1
     * ...
     *
     * This also allows pruning on size, support, and p-value.
     */
    public static void printPathFRs(String inputPrefix, int numCasePaths, int numCtrlPaths, int minSize, int minSupport, double maxPValue, int minPriority) throws IOException {
        // load the graph
	PangenomicGraph graph = readGraph(inputPrefix, numCasePaths, numCtrlPaths);
	// load the frequented regions and update support in parallel
	TreeSet<FrequentedRegion> frequentedRegions = readFrequentedRegions(inputPrefix, graph);
	frequentedRegions = pruneFrequentedRegions(frequentedRegions, minSize, minSupport, maxPValue, minPriority);
	System.err.println("FRUtils: "+frequentedRegions.size()+" FRs to be printed to pathfrs.txt file.");
	// collect the paths, cases and controls
	ConcurrentSkipListSet<Path> concurrentPaths = buildConcurrentPaths(graph, numCasePaths, numCtrlPaths);
	// number the FRs in whatever order they're in the set, into a map for parallel processing
	ConcurrentHashMap<String,FrequentedRegion> concurrentFRMap = new ConcurrentHashMap<>();
	int n = 0;
	for (FrequentedRegion fr : frequentedRegions) {
	    n++;
	    concurrentFRMap.put("FR"+n, fr);
	}
	/////////////////////////////////////////////////////////
	// now load pathFR strings for each FR in parallel
	ConcurrentSkipListSet<String> pathFRStrings = new ConcurrentSkipListSet<>();
	concurrentFRMap.entrySet().parallelStream().forEach(entry -> {
		String frName = entry.getKey();
		FrequentedRegion fr = entry.getValue();
		String pathFR = frName;
		for (Path path : concurrentPaths) {
		    pathFR += "\t"+fr.countSubpathsOf(path);
		}
		pathFRStrings.add(pathFR);
	    });
	/////////////////////////////////////////////////////////
	// output
        PrintStream out = new PrintStream(getPathFRsFilename(inputPrefix));
        boolean first = true;
        for (Path path : concurrentPaths) {
            if (first) {
                first = false;
            } else {
                out.print("\t");
            }
            out.print(path.name+"."+path.label);
        }
        out.println("");
	/////////////////////////////////////////////////////////
	// output path FR strings in parallel
	pathFRStrings.parallelStream().forEach(pathFR -> {
		out.println(pathFR);
	    });
	/////////////////////////////////////////////////////////
        out.close();
    }

    /**
     * Print out the best (last per run) FRs from a combined run file. Uses a TreeSet to make sure they're distinct.
     *
     * nodes	size	support	case	ctrl	OR	p	pri
     * [891]	1	4712	2549	2163	1.178	1.01E-14	71
     * [786,891]	2	113	10	103	0.097	8.54E-21	1012
     * [711,786,891]	3	104	2	102	0.020	3.25E-28	1707 <== print this one
     * nodes	size	support	case	ctrl	OR	p	pri
     * [892]	1	5149	2377	2772	0.858	2.50E-15	66
     * [259,892]	2	204	66	138	0.478	3.92E-7	320
     * [259,411,892]	3	101	24	77	0.312	1.02E-7	506
     * [259,411,873,892]	4	100	23	77	0.299	4.76E-8	524 <== print this one
     */
    public static void printBestFRs(String inputPrefix) throws IOException {
	TreeSet<FrequentedRegion> sortedFRs = new TreeSet<>();
        String frFilename = getFRsFilename(inputPrefix);
        BufferedReader reader = new BufferedReader(new FileReader(frFilename));
	String lastLine = null;
        String thisLine = null;
        while ((thisLine=reader.readLine())!=null) {
	    if (thisLine.startsWith("nodes")) {
		if (lastLine!=null) {
		    sortedFRs.add(new FrequentedRegion(lastLine));
		}
		lastLine = null;
	    } else {
		lastLine = thisLine;
	    }
	}
	sortedFRs.add(new FrequentedRegion(lastLine));
	// output
	System.out.println("nodes\tsize\tsupport\tcase\tctrl\tOR\tp\tpri");
	for (FrequentedRegion fr : sortedFRs) {
	    System.out.println(fr.toString());
	}
    }

     /**
      * Main class with methods for post-processing.
      */
    public static void main(String[] args) throws FileNotFoundException, IOException {
        Options options = new Options();
        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd;

        Option inputPrefixOption = new Option("i", "inputprefix", true, "prefix of input files (e.g. HLAA-0.1-Inf)");
        inputPrefixOption.setRequired(true);
        options.addOption(inputPrefixOption);
	//
	Option postprocessOption = new Option("p", "postprocess", false, "postprocess by applying new minsup, minsize to an FR set given by inputprefix");
	postprocessOption.setRequired(false);
	options.addOption(postprocessOption);
        //
        Option minSupportOption = new Option("m", "minsupport", true, "minimum number of supporting paths for an FR to be considered interesting");
        minSupportOption.setRequired(false);
        options.addOption(minSupportOption);
        //
        Option maxPValueOption = new Option("mp", "maxpvalue", true, "maximum p-value for an FR to be considered interesting");
        maxPValueOption.setRequired(false);
        options.addOption(maxPValueOption);
	//
        Option minPriorityOption = new Option("mpri", "minpriority", true, "minimum priority value for an FR to be considered interesting");
        minPriorityOption.setRequired(false);
        options.addOption(minPriorityOption);
	//
        Option minSizeOption = new Option("s", "minsize", true, "minimum number of nodes that a FR must contain to be considered interesting");
        minSizeOption.setRequired(false);
        options.addOption(minSizeOption);
	//
	Option printPathFRsSVMOption = new Option("svm", "svm", false, "print out an SVM style file from the data given by inputprefix");
	printPathFRsSVMOption.setRequired(false);
	options.addOption(printPathFRsSVMOption);
	//
	Option printPathFRsARFFOption = new Option("arff", "arff", false, "print out an ARFF style file from the data given by inputprefix");
	printPathFRsARFFOption.setRequired(false);
	options.addOption(printPathFRsARFFOption);
	//
	Option printPathFRsOption = new Option("pathfrs", "pathfrs", false, "print out a pathfrs.txt style file from the FR data given by inputprefix");
	printPathFRsOption.setRequired(false);
	options.addOption(printPathFRsOption);
        //
	Option numCasePathsOption = new Option("ncase", "numcasepaths", true, "number of case paths to include in SVM or ARFF or TXT output");
	numCasePathsOption.setRequired(false);
	options.addOption(numCasePathsOption);
	//
	Option numCtrlPathsOption = new Option("nctrl", "numcontrolpaths", true, "number of control paths to include in SVM or ARFF or TXT output");
	numCtrlPathsOption.setRequired(false);
	options.addOption(numCtrlPathsOption);
	//
	Option extractBestFRsOption = new Option("extractbest", "extractbestfrs", false, "extract the best (last) FRs from a file containing many runs, each starting with the heading line");
	extractBestFRsOption.setRequired(false);
	options.addOption(extractBestFRsOption);

        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println(e.getMessage());
            formatter.printHelp("FRUtils", options);
            System.exit(1);
            return;
        }

        if (cmd.getOptions().length==0) {
            formatter.printHelp("FRUtils", options);
            System.exit(1);
            return;
        }
	
	String inputPrefix = cmd.getOptionValue("inputprefix");

	if (cmd.hasOption("p")) {
	    int minSupport = Integer.parseInt(cmd.getOptionValue("minsupport"));
	    int minSize = Integer.parseInt(cmd.getOptionValue("minsize"));
	    postprocess(inputPrefix, minSupport, minSize);
	}

	if (cmd.hasOption("svm")) {
	    int numCasePaths = 0;
	    int numCtrlPaths = 0;
	    int minSize = 0;
	    int minSupport = 0;
	    int minPriority = 0;
	    double maxPValue = 1.0;
	    if (cmd.hasOption("ncase")) numCasePaths = Integer.parseInt(cmd.getOptionValue("ncase"));
	    if (cmd.hasOption("nctrl")) numCtrlPaths = Integer.parseInt(cmd.getOptionValue("nctrl"));
	    if (cmd.hasOption("minsize")) minSize = Integer.parseInt(cmd.getOptionValue("minsize"));
	    if (cmd.hasOption("minsupport")) minSupport = Integer.parseInt(cmd.getOptionValue("minsupport"));
	    if (cmd.hasOption("maxpvalue")) maxPValue = Double.parseDouble(cmd.getOptionValue("maxpvalue"));
	    if (cmd.hasOption("minpriority")) minPriority = Integer.parseInt(cmd.getOptionValue("minpriority"));
	    printPathFRsSVM(inputPrefix, numCasePaths, numCtrlPaths, minSize, minSupport, maxPValue, minPriority);
	}

	if (cmd.hasOption("arff")) {
	    int numCasePaths = 0;
	    int numCtrlPaths = 0;
	    int minSize = 0;
	    int minSupport = 0;
	    int minPriority = 0;
	    double maxPValue = 1.0;
	    if (cmd.hasOption("ncase")) numCasePaths = Integer.parseInt(cmd.getOptionValue("ncase"));
	    if (cmd.hasOption("nctrl")) numCtrlPaths = Integer.parseInt(cmd.getOptionValue("nctrl"));
	    if (cmd.hasOption("minsize")) minSize = Integer.parseInt(cmd.getOptionValue("minsize"));
	    if (cmd.hasOption("minsupport")) minSupport = Integer.parseInt(cmd.getOptionValue("minsupport"));
	    if (cmd.hasOption("maxpvalue")) maxPValue = Double.parseDouble(cmd.getOptionValue("maxpvalue"));
	    if (cmd.hasOption("minpriority")) minPriority = Integer.parseInt(cmd.getOptionValue("minpriority"));
	    printPathFRsARFF(inputPrefix, numCasePaths, numCtrlPaths, minSize, minSupport, maxPValue, minPriority);
	}

	if (cmd.hasOption("pathfrs")) {
	    int numCasePaths = 0;
	    int numCtrlPaths = 0;
	    int minSize = 0;
	    int minSupport = 0;
	    int minPriority = 0;
	    double maxPValue = 1.0;
	    if (cmd.hasOption("ncase")) numCasePaths = Integer.parseInt(cmd.getOptionValue("ncase"));
	    if (cmd.hasOption("nctrl")) numCtrlPaths = Integer.parseInt(cmd.getOptionValue("nctrl"));
	    if (cmd.hasOption("minsize")) minSize = Integer.parseInt(cmd.getOptionValue("minsize"));
	    if (cmd.hasOption("minsupport")) minSupport = Integer.parseInt(cmd.getOptionValue("minsupport"));
	    if (cmd.hasOption("maxpvalue")) maxPValue = Double.parseDouble(cmd.getOptionValue("maxpvalue"));
	    if (cmd.hasOption("minpriority")) minPriority = Integer.parseInt(cmd.getOptionValue("minpriority"));
	    printPathFRs(inputPrefix, numCasePaths, numCtrlPaths, minSize, minSupport, maxPValue, minPriority);
	}

	if (cmd.hasOption("extractbestfrs")) {
	    printBestFRs(inputPrefix);
	}
    }

    /**
     * Print out a set of FRs.
     */
    public static void printFrequentedRegions(Set<FrequentedRegion> frequentedRegions, String inputPrefix) throws IOException {
        if (frequentedRegions.size()==0) {
            System.err.println("NO FREQUENTED REGIONS!");
            return;
        }
        PrintStream out = new PrintStream(FRUtils.getFRsFilename(inputPrefix));
        boolean first = true;
        for (FrequentedRegion fr : frequentedRegions) {
            if (first) {
                out.println(fr.columnHeading());
                first = false;
            }
            out.println(fr.toString());
        }
        out.close();
    }

    /**
     * Load a graph from a pair of TXT files along with checks on the number of case and control paths desired.
     */
    static PangenomicGraph readGraph(String inputPrefix, int numCasePaths, int numCtrlPaths) throws FileNotFoundException, IOException {
	PangenomicGraph graph = new PangenomicGraph();
	graph.nodesFile = new File(getNodesFilename(inputPrefix));
	graph.pathsFile = new File(getPathsFilename(inputPrefix)); 
	graph.loadTXT();
	graph.tallyLabelCounts();
	// check on numCasePaths, numCtrlPaths if nonzero
	if (numCasePaths>graph.labelCounts.get("case")) {
	    System.err.println("FRUtils ERROR: numCasePaths > "+graph.labelCounts.get("case")+" cases.");
	    System.exit(1);
	}
	if (numCtrlPaths>graph.labelCounts.get("ctrl")) {
	    System.err.println("FRUtils ERROR: numCtrlPaths > "+graph.labelCounts.get("ctrl")+" controls.");
	    System.exit(1);
	}
        return graph;
    }

    /**
     * Prune a TreeSet of FrequentedRegions based on size, support, p-value and priority.
     */
    public static TreeSet<FrequentedRegion> pruneFrequentedRegions(TreeSet<FrequentedRegion> frequentedRegions, int minSize, int minSupport, double maxPValue, int minPriority) {
	frequentedRegions = removeNoCalls(frequentedRegions);
        frequentedRegions = pruneOnSize(frequentedRegions, minSize);
        frequentedRegions = pruneOnSupport(frequentedRegions, minSupport);
        frequentedRegions = pruneOnPValue(frequentedRegions, maxPValue);
        frequentedRegions = pruneOnPriority(frequentedRegions, minPriority);
	return frequentedRegions;
    }

    /**
     * Remove FRs with no-call nodes.
     */
    static TreeSet<FrequentedRegion> removeNoCalls(TreeSet<FrequentedRegion> frequentedRegions) {
        List<FrequentedRegion> frsToRemove = new LinkedList<>();
        for (FrequentedRegion fr : frequentedRegions) {
            fr.updateNodes();
            if (fr.containsNoCallNode()) frsToRemove.add(fr);
        }
        TreeSet<FrequentedRegion> prunedFRs = new TreeSet<>(frequentedRegions);
        prunedFRs.removeAll(frsToRemove);
        System.err.println("FRUtils removed "+frsToRemove.size()+" FRs with no-call nodes.");
        return prunedFRs;
    }    

    /**
     * Remove FRs with size<minSize.
     */
    static TreeSet<FrequentedRegion> pruneOnSize(TreeSet<FrequentedRegion> frequentedRegions, int minSize) {
        if (minSize<=1) return frequentedRegions;
        List<FrequentedRegion> frsToRemove = new LinkedList<>();
        for (FrequentedRegion fr : frequentedRegions) {
            if (fr.size<minSize) frsToRemove.add(fr);
        }
        TreeSet<FrequentedRegion> prunedFRs = new TreeSet<>(frequentedRegions);
        prunedFRs.removeAll(frsToRemove);
        System.err.println("FRUtils removed "+frsToRemove.size()+" FRs with size<"+minSize);
        return prunedFRs;
    }

    /**
     * Remove FRs with support<minSupport.
     */
    static TreeSet<FrequentedRegion> pruneOnSupport(TreeSet<FrequentedRegion> frequentedRegions, int minSupport) {
        if (minSupport<=1) return frequentedRegions;
        List<FrequentedRegion> frsToRemove = new LinkedList<>();
        for (FrequentedRegion fr : frequentedRegions) {
            if (fr.support<minSupport) frsToRemove.add(fr);
        }
        TreeSet<FrequentedRegion> prunedFRs = new TreeSet<>(frequentedRegions);
        prunedFRs.removeAll(frsToRemove);
        System.err.println("FRUtils removed "+frsToRemove.size()+" FRs with support<"+minSupport);
        return prunedFRs;
    }

    /**
     * Remove FRs with p>minPValue.
     */
    static TreeSet<FrequentedRegion> pruneOnPValue(TreeSet<FrequentedRegion> frequentedRegions, double maxPValue) {
    	if (maxPValue>=1.0) return frequentedRegions;
        List<FrequentedRegion> frsToRemove = new LinkedList<>();
        for (FrequentedRegion fr : frequentedRegions) {
            if (fr.pValue>maxPValue) frsToRemove.add(fr);
        }
        TreeSet<FrequentedRegion> prunedFRs = new TreeSet<>(frequentedRegions);
        prunedFRs.removeAll(frsToRemove);
        System.err.println("FRUtils removed "+frsToRemove.size()+" FRs with p>"+maxPValue);
        return prunedFRs;
    }

    /**
     * Remove FRs with priority<minPriority.
     */
    static TreeSet<FrequentedRegion> pruneOnPriority(TreeSet<FrequentedRegion> frequentedRegions, int minPriority) {
        if (minPriority<1) return frequentedRegions;
        List<FrequentedRegion> frsToRemove = new LinkedList<>();
        for (FrequentedRegion fr : frequentedRegions) {
            if (fr.priority<minPriority) frsToRemove.add(fr);
        }
        TreeSet<FrequentedRegion> prunedFRs = new TreeSet<>(frequentedRegions);
        prunedFRs.removeAll(frsToRemove);
        System.err.println("FRUtils removed "+frsToRemove.size()+" FRs with priority<"+minPriority);
        return prunedFRs;
    }

    /**
     * Build a ConcurrentSkipListSet of paths according to the given filters.
     */
    static ConcurrentSkipListSet<Path> buildConcurrentPaths(PangenomicGraph graph, int numCasePaths, int numCtrlPaths) {
	ConcurrentSkipListSet<Path> concurrentPaths = new ConcurrentSkipListSet<>();
	int nCases = 0;
	if (numCasePaths==0) {
	    // select all case paths
	    for (Path path : graph.paths) {
		if (path.isCase()) {
		    nCases++;
		    concurrentPaths.add(path);
		}
	    }
	} else {
	    // randomly select case paths
	    while (nCases<numCasePaths) {
		Optional<Path> optional = graph.paths.stream().skip((int)(graph.paths.size()*Math.random())).findFirst();
		if (optional.isPresent()) {
		    Path path = optional.get();
		    if (path.isCase() && !concurrentPaths.contains(path)) {
			nCases++;
			concurrentPaths.add(path);
		    }
		}
	    }
	}
	int nControls = 0;
	if (numCtrlPaths==0) {
	    // select all control paths
	    for (Path path : graph.paths) {
		if (path.isControl()) {
		    nControls++;
		    concurrentPaths.add(path);
		}
	    }
	} else {
	    // randomly select control paths
	    while (nControls<numCtrlPaths) {
		Optional<Path> optional = graph.paths.stream().skip((int)(graph.paths.size()*Math.random())).findFirst();
		if (optional.isPresent()) {
		    Path path = optional.get();
		    if (path.isControl() && !concurrentPaths.contains(path)) {
			nControls++;
			concurrentPaths.add(path);
		    }
		}
	    }
	}
        return concurrentPaths;
    }
}
