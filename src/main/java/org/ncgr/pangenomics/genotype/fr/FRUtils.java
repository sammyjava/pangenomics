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

import java.util.Arrays;
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
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.concurrent.ConcurrentSkipListSet;
import java.util.stream.Collectors;

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
     * This method populates the FR subpaths based on alpha, kappa, and the given priority option.
     *
     * nodes   size    support case    ctrl    OR      p       pri
     * [1341]  1       21      19      2       9.500   9.48E-3 202
     * 509678.ctrl:[18,20,21,23,24,26,27,29,30,33,34]
     * 628863.case:[18,20,21,23,24,26,27,29,30,33,34]
     * etc.
     */
    public static TreeSet<FrequentedRegion> readFrequentedRegions(PangenomicGraph graph, String frsPrefix, int priorityOptionKey, String priorityOptionLabel) throws FileNotFoundException, IOException {
        // get alpha, kappa, and graph from the input prefix
        double alpha = readAlpha(frsPrefix);
        int kappa = readKappa(frsPrefix);
        String frFilename = getFRsFilename(frsPrefix);
	// read the FRs
        TreeSet<FrequentedRegion> frequentedRegions = new TreeSet<>();
        BufferedReader reader = new BufferedReader(new FileReader(frFilename));
        String line = null;
        while ((line=reader.readLine())!=null) {
	    // this does NOT run update() and pValue and orValue must be reset to get non-rounded versions below
	    FrequentedRegion fr = new FrequentedRegion(graph, alpha, kappa, line);
	    if (fr.size>0) {
		fr.priorityOptionKey = priorityOptionKey;
		fr.priorityOptionLabel = priorityOptionLabel;
		fr.pValue = Double.NEGATIVE_INFINITY;
		fr.orValue = Double.NEGATIVE_INFINITY;
		frequentedRegions.add(fr);
	    }
	}
	System.err.println("Loaded "+frequentedRegions.size()+" FRs from "+frsPrefix);
	/////////////////////////////////////////////////////////////////////////////////
	// build the subpaths for each FR
	System.err.println("Updating FR support and priorities...");
	Set<FrequentedRegion> updatedFRs = frequentedRegions.parallelStream().map(fr -> {
		fr.update();
		return fr;
	    }).collect(Collectors.toSet());
	/////////////////////////////////////////////////////////////////////////////////
	return new TreeSet<FrequentedRegion>(updatedFRs);
    }

    /**
     * Read FRs from the text file with no filtering on size, support, p, or priority.
     */
    public static TreeSet<FrequentedRegion> readFrequentedRegions(PangenomicGraph graph, String frsPrefix) throws FileNotFoundException, IOException {
	return readFrequentedRegions(graph, frsPrefix, 1, 1, 1.0, 0);
    }

    /**
     * Read FRs from the text file written by an FR run, filtering on size, support, p-value and priority as given in the file.
     * This method populates the FR subpaths based on alpha, kappa.
     * 0       1       2       3       4       5       6       7
     * nodes   size    support case    ctrl    OR      p       pri
     * [1341]  1       21      19      2       9.500   9.48E-3 202
     * 509678.ctrl:[18,20,21,23,24,26,27,29,30,33,34]
     * 628863.case:[18,20,21,23,24,26,27,29,30,33,34]
     * etc.
     */
    public static TreeSet<FrequentedRegion> readFrequentedRegions(PangenomicGraph graph, String frsPrefix,
								  int minSize, int minSupport, double maxPValue, int minPriority) throws FileNotFoundException, IOException {
        // get alpha, kappa, and graph from the input prefix
        double alpha = readAlpha(frsPrefix);
        int kappa = readKappa(frsPrefix);
        String frFilename = getFRsFilename(frsPrefix);
	// read the FRs
        TreeSet<FrequentedRegion> frequentedRegions = new TreeSet<>();
        BufferedReader reader = new BufferedReader(new FileReader(frFilename));
        String line = null;
	int droppedSizeCount = 0;
	int droppedSupportCount = 0;
	int droppedPCount = 0;
	int droppedPriCount = 0;
        while ((line=reader.readLine())!=null) {
	    if (line.startsWith("#") || line.startsWith("nodes")) continue;
	    String[] fields = line.split("\t");
	    int size = Integer.parseInt(fields[1]);
	    int support = Integer.parseInt(fields[2]);
	    double p = Double.parseDouble(fields[6]);
	    int pri = Integer.parseInt(fields[7]);
	    if (size<minSize) {
		droppedSizeCount++;
	    } else if (support<minSupport) {
		droppedSupportCount++;
	    } else if (p>maxPValue) {
		droppedPCount++;
	    } else if (pri<minPriority) {
		droppedPriCount++;
	    } else {
		frequentedRegions.add(new FrequentedRegion(graph, alpha, kappa, line));
	    }
	}
	// DX
	System.err.println("## size<"+minSize+": "+droppedSizeCount);
	System.err.println("## support<"+minSupport+": "+droppedSupportCount);
	System.err.println("## p>"+maxPValue+": "+droppedPCount);
	System.err.println("## pri<"+minPriority+": "+droppedPriCount);
	System.err.println("## total dropped: "+(droppedSizeCount+droppedSupportCount+droppedPCount+droppedPriCount));
	// number the FRs
        int number = 0;
        for (FrequentedRegion fr : frequentedRegions) {
            fr.number = number++;
        }
	System.err.println("Loaded "+frequentedRegions.size()+" unique FRs from "+frsPrefix);
	/////////////////////////////////////////////////////////////////////////////////
	// build the subpaths for each FR
	System.err.println("Updating FR support and priorities...");
	Set<FrequentedRegion> updatedFRs = frequentedRegions.parallelStream().map(fr -> {
		fr.update();
		return fr;
	    }).collect(Collectors.toSet());
	/////////////////////////////////////////////////////////////////////////////////
	return frequentedRegions;
    }

    /**
     * Return alpha from an FRFinder run.
     */
    public static double readAlpha(String frsPrefix) throws FileNotFoundException, IOException {
        String paramsFilename = getParamsFilename(frsPrefix);
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
    public static int readKappa(String frsPrefix) throws FileNotFoundException, IOException {
        String paramsFilename = getParamsFilename(frsPrefix);
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
    public static String readPriorityOption(String frsPrefix) throws FileNotFoundException, IOException {
        String paramsFilename = getParamsFilename(frsPrefix);
        BufferedReader reader = new BufferedReader(new FileReader(paramsFilename));
        Properties parameters = new Properties();
        parameters.load(reader);
        return parameters.getProperty("priorityOption");
    }        

    /**
     * Form the FRs output filename
     */
    public static String getFRsFilename(String frsPrefix) {
        return frsPrefix+".frs.txt";
    }

    /**
     * Form the FRSubpaths output filename
     */
    public static String getFRSubpathsFilename(String frsPrefix) {
        return frsPrefix+".subpaths.txt";
    }

    /**
     * Form the pathFRs output filename
     */
    public static String getPathFRsFilename(String pathsPrefix, String frsPrefix) {
	if (pathsPrefix==null) {
	    return frsPrefix+".pathfrs.txt";
	} else {
	    return frsPrefix+"."+pathsPrefix+".pathfrs.txt";
	}
    }

    /**
     * Form the SVM version of the pathFRs output filename
     */
    public static String getPathFRsSVMFilename(String pathsPrefix, String frsPrefix) {
	
        return frsPrefix+"."+pathsPrefix+".svm.txt";
    }

    /**
     * Form the ARFF version of the pathFRs output filename
     */
    public static String getPathFRsARFFFilename(String pathsPrefix, String frsPrefix) {
        return frsPrefix+"."+pathsPrefix+".arff";
    }

    /**
     * Form the parameters output filename
     */
    public static String getParamsFilename(String prefix) {
        return prefix+".params.txt";
    }

    /**
     * Grab the graph name from an FR prefix.
     */
    public static String getGraphName(String frPrefix) {
	String[] parts = frPrefix.split("-");
	if (parts.length==1) {
	    System.err.println("ERROR: FR prefix "+frPrefix+" does not contain '-'.");
	    System.exit(1);
	}
	return parts[0];
    }

    /**
     * Form the graph nodes filename from an FR prefix.
     * if prefix = HTT.1k-1.0-0 then filename = HTT.nodes.txt
     */
    public static String getNodesFilename(String frPrefix) {
        return getGraphName(frPrefix)+".nodes.txt";
    }

    /**
     * Form the graph paths filename from an FR prefix.
     * if prefix = HTT.1k-1.0-0 then filename = HTT.paths.txt
     */
    public static String getPathsFilename(String frPrefix) {
        return getGraphName(frPrefix)+".paths.txt";
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
    public static Properties readParameters(String frsPrefix) throws FileNotFoundException, IOException {
        String paramsFilename = getParamsFilename(frsPrefix);
        BufferedReader reader = new BufferedReader(new FileReader(paramsFilename));
        Properties parameters = new Properties();
        parameters.load(reader);
        parameters.setProperty("paramsFile", paramsFilename);
        parameters.setProperty("frsPrefix", frsPrefix);
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
     * Form an outputPrefix from frsPrefix, minSupport, and minSize.
     */
    public static String formOutputPrefix(String frsPrefix, int minSupport, int minSize) {
        return frsPrefix+"-"+minSupport+"."+minSize;
    }

    /**
     * Post-process a list of FRs read in from the frsPrefix file for given minSupport and minSize.
     * Outputs a new set of FR files with minSupport and minSize in the outputPrefix.
     */
    public static void postprocess(String graphPrefix, String pathsPrefix, String frsPrefix, int minSupport, int minSize) throws FileNotFoundException, IOException {
	// read the graph
	PangenomicGraph graph = readGraph(graphPrefix, pathsPrefix);
	// parameters
	Properties parameters = readParameters(frsPrefix);
	parameters.setProperty("minSupport", String.valueOf(minSupport));
	parameters.setProperty("minSize", String.valueOf(minSize));
	Set<FrequentedRegion> frequentedRegions = readFrequentedRegions(graph, frsPrefix);
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
        System.err.println(filteredFRs.size()+" FRs passed minSupport="+minSupport+", minSize="+minSize);
	// output the filtered FRs
	if (filteredFRs.size()>0) {
	    printFrequentedRegions(filteredFRs, formOutputPrefix(frsPrefix, minSupport, minSize));
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
     */
    public static void printPathFRsSVM(String graphPrefix, String pathsPrefix, String frsPrefix,
				       int minSize, int minSupport, double maxPValue, int minPriority) throws IOException, FileNotFoundException {
	// read the graph
	final PangenomicGraph graph = readGraph(graphPrefix, pathsPrefix);
        // read the FRs from files
	final TreeSet<FrequentedRegion> frequentedRegions = readFrequentedRegions(graph, frsPrefix, minSize, minSupport, maxPValue, minPriority);
	// collect the paths, cases and controls
        final ConcurrentSkipListSet<Path> concurrentPaths = buildConcurrentPaths(graph);
	///////////////////////////////////////////////////////
	// build the SVM strings in parallel
	ConcurrentSkipListMap<String,String> pathSVM = new ConcurrentSkipListMap<>();
	concurrentPaths.parallelStream().forEach(path -> {
		String svm = path.getLabel();
		int c  = 0;
		for (FrequentedRegion fr : frequentedRegions) {
		    c++;
		    svm += "\t"+c+":"+fr.countSubpathsOf(path);
		}
		pathSVM.put(path.getName(), svm);
	    });
	///////////////////////////////////////////////////////
	// sort by SVM string
	LinkedHashMap<String,String> sortedPathSVM = new LinkedHashMap<>();
	pathSVM.entrySet().stream().sorted(Map.Entry.comparingByValue()).forEachOrdered(x->sortedPathSVM.put(x.getKey(),x.getValue()));
	// output: no header, one path per row
        PrintStream out = new PrintStream(getPathFRsSVMFilename(pathsPrefix, frsPrefix));
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
    public static void printPathFRsARFF(String graphPrefix, String pathsPrefix, String frsPrefix) throws IOException {
	// read the graph
	final PangenomicGraph graph = readGraph(graphPrefix, pathsPrefix);
    	// read the FRs from files
    	final TreeSet<FrequentedRegion> frequentedRegions = readFrequentedRegions(graph, frsPrefix);
    	// collect the paths, cases and controls
        final ConcurrentSkipListSet<Path> concurrentPaths = buildConcurrentPaths(graph);
    	/////////////////////////////////////////////////////////
    	// now load pathARFF for selected paths in parallel
    	ConcurrentSkipListMap<String,String> pathARFF = new ConcurrentSkipListMap<>();
    	concurrentPaths.parallelStream().forEach(path -> {
    		String arff = "";
    		for (FrequentedRegion fr : frequentedRegions) {
    		    arff += fr.countSubpathsOf(path)+",";
    		}
    		arff += path.getLabel();
    		pathARFF.put(path.getName(), arff);
    	    });
    	/////////////////////////////////////////////////////////
    	// sort by ARFF string
    	LinkedHashMap<String,String> sortedPathARFF = new LinkedHashMap<>();
    	pathARFF.entrySet().stream().sorted(Map.Entry.comparingByValue()).forEachOrdered(x->sortedPathARFF.put(x.getKey(),x.getValue()));
    	// ARFF output	
        PrintStream out = new PrintStream(FRUtils.getPathFRsARFFFilename(pathsPrefix, frsPrefix));
        out.println("@RELATION "+frsPrefix);
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
     */
    public static void printPathFRs(String graphPrefix, String pathsPrefix, String frsPrefix) throws IOException {
	// read the graph
	final PangenomicGraph graph = readGraph(graphPrefix, pathsPrefix);
    	// read the FRs from files
    	final TreeSet<FrequentedRegion> frequentedRegions = readFrequentedRegions(graph, frsPrefix);
    	// collect the paths, cases and controls
    	final ConcurrentSkipListSet<Path> concurrentPaths = buildConcurrentPaths(graph);
        // build a map of FR# labels to FR for row labels
    	ConcurrentSkipListMap<String,FrequentedRegion> concurrentFRMap = new ConcurrentSkipListMap<>();
    	for (FrequentedRegion fr : frequentedRegions) {
            String frLabel = "FR"+fr.number;
    	    concurrentFRMap.put(frLabel, fr);
    	}
  	/////////////////////////////////////////////////////////
    	// now load pathFR strings for each FR in parallel
    	TreeSet<String> pathFRStrings = new TreeSet<>();
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
        PrintStream out = new PrintStream(getPathFRsFilename(pathsPrefix, frsPrefix));
        boolean first = true;
        for (Path path : concurrentPaths) {
            if (first) {
                first = false;
            } else {
                out.print("\t");
            }
            out.print(path.getName()+"."+path.getLabel());
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
     * Calculate and print out polynomial risk scores per sample from FR support and odds ratio.
     */
    public static void calculatePRS(String graphPrefix, String pathsPrefix, String frsPrefix) throws IOException {
	// read in the graph
	final PangenomicGraph graph = readGraph(graphPrefix, pathsPrefix);
    	// read the FRs from files
	final ConcurrentSkipListSet<FrequentedRegion> concurrentFRs = new ConcurrentSkipListSet<>(readFrequentedRegions(graph, frsPrefix));
    	// build the paths
        final ConcurrentSkipListSet<Path> concurrentPaths = buildConcurrentPaths(graph);
	// maps
    	ConcurrentSkipListMap<Path,Integer> concurrentPathSupport = new ConcurrentSkipListMap<>(); // stores sum of support
        ConcurrentSkipListMap<Path,Double> concurrentPathPRS = new ConcurrentSkipListMap<>();      // stores PRS
  	/////////////////////////////////////////////////////////
        // initialize the result maps
        concurrentPaths.parallelStream().forEach(path -> {
                concurrentPathSupport.put(path, 0);
                concurrentPathPRS.put(path, 0.0);
            });
  	/////////////////////////////////////////////////////////
  	/////////////////////////////////////////////////////////
        // build the support and PRS values
        concurrentFRs.parallelStream().forEach(fr -> {
                double logOR = Math.log(fr.oddsRatio());
                concurrentPaths.parallelStream().forEach(path -> {
                        int prevSupport = concurrentPathSupport.get(path);
                        double prevPRS = concurrentPathPRS.get(path);
                        int support = fr.countSubpathsOf(path);
                        double deltaPRS = support*logOR;
                        concurrentPathSupport.put(path, prevSupport+support);
                        concurrentPathPRS.put(path, prevPRS+deltaPRS);
                    });
            });
    	/////////////////////////////////////////////////////////
    	/////////////////////////////////////////////////////////
        // divide the PRS sums by total support
        concurrentPaths.parallelStream().forEach(path -> {
                int totalSupport = concurrentPathSupport.get(path);
                double totalPRS = concurrentPathPRS.get(path);
                concurrentPathPRS.put(path, totalPRS/totalSupport);
            });
        /////////////////////////////////////////////////////////
    	// output
	System.out.println("sample\tlabel\tscore");
        for (Path path : concurrentPathPRS.keySet()) {
            System.out.println(path.getName()+"\t"+path.getLabel()+"\t"+concurrentPathPRS.get(path));
        }
    }

    /**
     * Print out the best (last per run) FRs from a combined run file. Uses a TreeSet to make sure they're unique.
     *
     * 0        1       2       3       4       5       6       7
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
    public static void printBestFRs(String frFilename, double alpha, int kappa) throws IOException {
	TreeSet<String> sortedLines = new TreeSet<>();
        BufferedReader reader = new BufferedReader(new FileReader(frFilename));
	String header = null;
	String lastLine = null;
        String thisLine = null;
        while ((thisLine=reader.readLine())!=null) {
	    if (lastLine==null) {
		// header
		header = thisLine;
	    } else if (thisLine.startsWith("nodes") && lastLine!=null) {
		// best FR from this round
		sortedLines.add(lastLine);
	    }
	    lastLine = thisLine;
	}
	// last best FR
	sortedLines.add(lastLine);
	// spit 'em out
	System.out.println(header);
	for (String line : sortedLines) {
	    System.out.println(line);
	}
    }

     /**
      * Main class with methods for various utility tasks
      */
    public static void main(String[] args) throws FileNotFoundException, IOException {
        Options options = new Options();
        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd;

        Option graphPrefixOption = new Option("g", "graphprefix", true, "prefix of graph file (e.g. HLAA)");
        graphPrefixOption.setRequired(false);
        options.addOption(graphPrefixOption);
        //
        Option pathsPrefixOption = new Option("p", "pathsprefix", true, "prefix of paths file (e.g. HLAA)");
        pathsPrefixOption.setRequired(false);
        options.addOption(pathsPrefixOption);
        //
        Option frsPrefixOption = new Option("f", "frsprefix", true, "prefix of FRs file (e.g. HLAA-0.0-Inf)");
        frsPrefixOption.setRequired(false);
        options.addOption(frsPrefixOption);
        //
        Option priorityOptionOption = new Option("pri", "priorityoption", true, "option for priority weighting of FRs: "+FrequentedRegion.PRIORITY_OPTIONS);
        priorityOptionOption.setRequired(false);
        options.addOption(priorityOptionOption);
	//
	Option postprocessOption = new Option("post", "postprocess", false, "postprocess by applying new minsup, minsize to an FR set given by inputprefix");
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
	Option extractBestFRsOption = new Option("extractbest", "extractbestfrs", false, "extract the best FRs from a file containing many runs, each starting with the header");
	extractBestFRsOption.setRequired(false);
	options.addOption(extractBestFRsOption);
        //
        Option prsOption = new Option("prs", "prs", false, "generate polynomial risk scores from FRFinder output");
        prsOption.setRequired(false);
        options.addOption(prsOption);

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

	// file stuff
	String graphPrefix = cmd.getOptionValue("graphprefix");
	String pathsPrefix = cmd.getOptionValue("pathsprefix");
	String frsPrefix = cmd.getOptionValue("frsprefix");

	// param defaults
	int minSize = 0;
	int minSupport = 0;
	double maxPValue = 1.0;
	int minPriority = 0;
	// supplied params
	if (cmd.hasOption("minsize")) minSize = Integer.parseInt(cmd.getOptionValue("minsize"));
	if (cmd.hasOption("minsupport")) minSupport = Integer.parseInt(cmd.getOptionValue("minsupport"));
	if (cmd.hasOption("maxpvalue")) maxPValue = Double.parseDouble(cmd.getOptionValue("maxpvalue"));
	if (cmd.hasOption("minpriority")) minPriority = Integer.parseInt(cmd.getOptionValue("minpriority"));

	// ACTIONS
	if (cmd.hasOption("post")) {
	    postprocess(graphPrefix, pathsPrefix, frsPrefix, minSupport, minSize);
	} else if (cmd.hasOption("svm")) {
	    printPathFRsSVM(graphPrefix, pathsPrefix, frsPrefix, minSize, minSupport, maxPValue, minPriority);
	} else if (cmd.hasOption("arff")) {
	    printPathFRsARFF(graphPrefix, pathsPrefix, frsPrefix);
	} else if (cmd.hasOption("pathfrs")) {
	    printPathFRs(graphPrefix, pathsPrefix, frsPrefix);
	} else if (cmd.hasOption("prs")) {
	    calculatePRS(graphPrefix, pathsPrefix, frsPrefix);
        } else if (cmd.hasOption("extractbestfrs")) {
	    double alpha = readAlpha(frsPrefix);
	    int kappa = readKappa(frsPrefix);
	    String frFilename = getFRsFilename(frsPrefix);
	    printBestFRs(frFilename, alpha, kappa);
	}
    }

    /**
     * Print out a set of FRs.
     */
    public static void printFrequentedRegions(Set<FrequentedRegion> frequentedRegions, String graphPrefix) throws IOException {
        if (frequentedRegions.size()==0) {
            System.err.println("NO FREQUENTED REGIONS!");
            return;
        }
        PrintStream out = new PrintStream(FRUtils.getFRsFilename(graphPrefix));
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
    public static PangenomicGraph readGraph(String graphPrefix, String pathsPrefix) throws IOException {
	String nodesFilename = graphPrefix+".nodes.txt";
	String pathsFilename = pathsPrefix+".paths.txt";
	System.err.println("## Reading nodes file "+nodesFilename);
	System.err.println("## Reading paths file "+pathsFilename);
	File nodesFile = new File(nodesFilename);
	File pathsFile = new File(pathsFilename);
	if (!nodesFile.isFile()) {
	    System.err.println("ERROR: file "+nodesFilename+" does not exist.");
	    System.exit(1);
	} else if (!pathsFile.isFile()) {
	    System.err.println("ERROR: file "+pathsFilename+" does not exist.");
	    System.exit(1);
	}
	PangenomicGraph graph = new PangenomicGraph(graphPrefix);
	graph.loadNodesFromTXT(nodesFile);
	graph.loadPathsFromTXT(pathsFile);
	graph.tallyLabelCounts();
	return graph;
    }
    
    /**
     * Prune a TreeSet of FrequentedRegions based on size, support, p-value and priority.
     */
    public static TreeSet<FrequentedRegion> pruneFrequentedRegions(TreeSet<FrequentedRegion> frequentedRegions,
								   int minSize, int minSupport, double maxPValue, int minPriority) {
	TreeSet<FrequentedRegion> frNoCalls = removeNoCalls(frequentedRegions);
	TreeSet<FrequentedRegion> frMinSize = pruneOnSize(frNoCalls, minSize);
        TreeSet<FrequentedRegion> frMinSupport = pruneOnSupport(frMinSize, minSupport);
        TreeSet<FrequentedRegion> frPValue = pruneOnPValue(frMinSupport, maxPValue);
        TreeSet<FrequentedRegion> frPriority = pruneOnPriority(frPValue, minPriority);
	return frPriority;
    }

    /**
     * Remove FRs with no-call nodes.
     */
    static TreeSet<FrequentedRegion> removeNoCalls(TreeSet<FrequentedRegion> frequentedRegions) {
        TreeSet<FrequentedRegion> prunedFRs = new TreeSet<>();
	int count = 0;
        for (FrequentedRegion fr : frequentedRegions) {
            fr.updateNodes();
            if (fr.containsNoCallNode()) {
		count++;
	    } else {
		prunedFRs.add(fr);
	    }
        }
        System.err.println("Removed "+count+" FRs with no-call nodes.");
        return prunedFRs;
    }    

    /**
     * Remove FRs with size<minSize.
     */
    static TreeSet<FrequentedRegion> pruneOnSize(TreeSet<FrequentedRegion> frequentedRegions, int minSize) {
        if (minSize<=1) return frequentedRegions;
	TreeSet<FrequentedRegion> prunedFRs = new TreeSet<>();
	int count = 0;
        for (FrequentedRegion fr : frequentedRegions) {
            if (fr.size<minSize) {
		count++;
	    } else {
		prunedFRs.add(fr);
	    }
        }
        System.err.println("Removed "+count+" FRs with size<"+minSize);
        return prunedFRs;
    }

    /**
     * Remove FRs with support<minSupport.
     */
    static TreeSet<FrequentedRegion> pruneOnSupport(TreeSet<FrequentedRegion> frequentedRegions, int minSupport) {
        if (minSupport<=1) return frequentedRegions;
	TreeSet<FrequentedRegion> prunedFRs = new TreeSet<>();
	int count = 0;
        for (FrequentedRegion fr : frequentedRegions) {
            if (fr.support<minSupport) {
		count++;
	    } else {
		prunedFRs.add(fr);
	    }
        }
        System.err.println("Removed "+count+" FRs with support<"+minSupport);
        return prunedFRs;
    }

    /**
     * Remove FRs with p>minPValue.
     */
    static TreeSet<FrequentedRegion> pruneOnPValue(TreeSet<FrequentedRegion> frequentedRegions, double maxPValue) {
    	if (maxPValue>=1.0) return frequentedRegions;
	TreeSet<FrequentedRegion> prunedFRs = new TreeSet<>();
	int count = 0;
        for (FrequentedRegion fr : frequentedRegions) {
            if (fr.pValue>maxPValue) {
		count++;
	    } else {
		prunedFRs.add(fr);
	    }
        }
        System.err.println("Removed "+count+" FRs with p>"+maxPValue);
        return prunedFRs;
    }

    /**
     * Remove FRs with priority<minPriority.
     */
    static TreeSet<FrequentedRegion> pruneOnPriority(TreeSet<FrequentedRegion> frequentedRegions, int minPriority) {
        if (minPriority<1) return frequentedRegions;
	TreeSet<FrequentedRegion> prunedFRs = new TreeSet<>();
	int count = 0;
        for (FrequentedRegion fr : frequentedRegions) {
            if (fr.priority<minPriority) {
		count++;
	    } else {
		prunedFRs.add(fr);
	    }
        }
        System.err.println("Removed "+count+" FRs with priority<"+minPriority);
        return prunedFRs;
    }

    /**
     * Build a TreeSet of paths from the given graph.
     */
    static ConcurrentSkipListSet<Path> buildConcurrentPaths(PangenomicGraph graph) {
	ConcurrentSkipListSet<Path> concurrentPaths = new ConcurrentSkipListSet<>();
	int nCases = 0;
	// select all case paths
	for (Path path : graph.paths) {
	    if (path.isCase()) {
		nCases++;
		concurrentPaths.add(path);
	    }
	}
	int nControls = 0;
	// select all control paths
	for (Path path : graph.paths) {
	    if (path.isControl()) {
		nControls++;
		concurrentPaths.add(path);
	    }
	}
        return concurrentPaths;
    }

    /**
     * Return true iff both FR nodesets are on same chromosome.
     */
    static boolean onSameChromosome(FrequentedRegion fr1, FrequentedRegion fr2) {
	if (!fr1.onOneChromosome() || !fr2.onOneChromosome()) {
	    return false;
	} else {
	    return (fr1.nodes.first().contig.equals(fr2.nodes.first().contig));
	}
    }
}
