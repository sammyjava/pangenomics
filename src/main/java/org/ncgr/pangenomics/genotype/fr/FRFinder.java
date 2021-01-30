package org.ncgr.pangenomics.genotype.fr;

import org.ncgr.pangenomics.genotype.*;

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
import java.util.Properties;
import java.util.Set;
import java.util.HashSet;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentSkipListSet;

import java.text.DecimalFormat;
import java.time.ZonedDateTime;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

/**
 * Finds frequented regions in a pan-genomic graph.
 *
 * See Cleary, et al., "Exploring Frequented Regions in Pan-Genomic Graphs",
 * IEEE/ACM Trans Comput Biol Bioinform. 2018 Aug 9.
 * PMID:30106690 DOI:10.1109/TCBB.2018.2864564
 *
 * @author Sam Hokin
 */
public class FRFinder {

    // save filename suffixes
    static String FREQUENTED_REGIONS_SAVE = "save.frs.txt";
    static String ALL_FREQUENTED_REGIONS_SAVE = "save.allFrequentedRegions.txt";
    static String REJECTED_NODESETS_SAVE = "save.rejectedNodeSets.txt";
    static String ACCEPTED_FRPAIRS_SAVE = "save.acceptedFRPairs.txt";

    // number formats
    static DecimalFormat percf = new DecimalFormat("0.000%");

    // a PangenomicGraph must be supplied to the constructor
    PangenomicGraph graph;

    // store paths locally since post-processing doesn't build a graph
    Set<Path> paths;

    // parameters are stored in a Properties object
    Properties parameters = new Properties();

    // the output FRs
    Map<String,FrequentedRegion> frequentedRegions;

    // priority option
    String priorityOption = "4";    // default
    int priorityOptionKey;          // 0, 1, 2, etc.
    String priorityOptionLabel;     // the current label for priority update emphasis: "case" or "ctrl"
    String priorityOptionParameter; // parameter for priority update emphasis, can be null or "alt" or "case" or "ctrl"

    // keep option
    String keepOption = "0";        // default
    int keepOptionKey;
    
    // diagnostics
    boolean verbose = false;
    boolean debug = false;
    PrintStream logOut;

    // parameters, with defaults
    boolean resume = false;
    boolean writeSaveFiles = false;
    boolean writePathFRs = false;
    boolean writeFRSubpaths = false;
    boolean requireBestNodeSet = false;
    boolean includeNoCalls = false;
    boolean requireSamePosition = false;
    boolean requireHomozygous = false;
    double minMGF = 0.01;
    double maxPVal = 1.0;
    int minSize = 1;
    int maxSize = Integer.MAX_VALUE;
    int maxRound = 0;
    int minPriority = 0;
    int minSupport = 1;
    int maxClocktime = 0;
    String graphName = null;
    String requiredNodeString = "[]";
    String excludedNodeString = "[]";
    String includedNodeString = "[]";
    String excludedPathNodeString = "[]";
    String includedPathNodeString = "[]";
    
    /**
     * Construct with a populated Graph and default parameters.
     */
    public FRFinder(PangenomicGraph graph) {
        this.graph = graph;
        this.paths = graph.paths;
        initializeParameters();
    }


    /**
     * Initialize the default parameters.
     */
    void initializeParameters() {
        parameters.setProperty("resume", String.valueOf(resume));
        parameters.setProperty("writeSaveFiles", String.valueOf(writeSaveFiles));
        parameters.setProperty("writePathFRs", String.valueOf(writePathFRs));
        parameters.setProperty("writeFRSubpaths", String.valueOf(writeFRSubpaths));
        parameters.setProperty("minSupport", String.valueOf(minSupport));
        parameters.setProperty("minSize", String.valueOf(minSize));
	parameters.setProperty("maxSize", String.valueOf(maxSize));
        parameters.setProperty("maxRound", String.valueOf(maxRound));
	parameters.setProperty("maxClocktime", String.valueOf(maxClocktime));
        parameters.setProperty("minPriority", String.valueOf(minPriority));
        parameters.setProperty("requiredNodeString", requiredNodeString);
        parameters.setProperty("excludedNodeString", excludedNodeString);
        parameters.setProperty("includedNodeString", includedNodeString);
	parameters.setProperty("excludedPathNodeString", excludedPathNodeString);
	parameters.setProperty("includedPathNodeString", includedPathNodeString);
        parameters.setProperty("priorityOption", priorityOption);
        parameters.setProperty("keepOption", keepOption);
        parameters.setProperty("keepOption", "null");
        parameters.setProperty("minMGF", String.valueOf(minMGF));
	parameters.setProperty("maxPVal", String.valueOf(maxPVal));
	parameters.setProperty("includeNoCalls", String.valueOf(includeNoCalls));
	parameters.setProperty("requireBestNodeSet", String.valueOf(requireBestNodeSet));
	parameters.setProperty("requireSamePosition", String.valueOf(requireSamePosition));
	parameters.setProperty("requireHomozygous", String.valueOf(requireHomozygous));
    }

    /**
     * Find the frequented regions in this Graph for the given alpha and kappa values.
     * double alpha = penetrance: the fraction of a supporting strain's sequence that actually supports the FR;
     *                alternatively, `1-alpha` is the fraction of inserted sequence
     * int kappa = maximum insertion: the maximum number of inserted nodes that a supporting path may have.
      */
    public void findFRs(double alpha, int kappa) throws FileNotFoundException, IOException {
	// refresh this on each run
	frequentedRegions = new HashMap<>();
	
        // optional required and excluded nodes
	NodeSet requiredNodes = graph.getNodeSet(requiredNodeString);
	NodeSet includedNodes = graph.getNodeSet(includedNodeString);
	NodeSet excludedNodes = graph.getNodeSet(excludedNodeString);
	NodeSet initialRequiredNodes = graph.getNodeSet(requiredNodeString); // this one is static for outputPrefix

        // store all scanned FRs in a Map
	ConcurrentHashMap<String,FrequentedRegion> allFrequentedRegions = new ConcurrentHashMap<>();

        // rejected NodeSets (strings), so we don't bother scanning them more than once
        ConcurrentSkipListSet<String> rejectedNodeSets = new ConcurrentSkipListSet<>();

        // accepted FRPairs so we don't merge them more than once
        ConcurrentHashMap<String,FRPair> acceptedFRPairs = new ConcurrentHashMap<>();

        // initialize log file
        logOut = new PrintStream(formOutputPrefix(alpha, kappa, initialRequiredNodes)+".log");
	printToLog("# Starting findFRs: " +
                   "alpha="+alpha+" " +
                   "kappa="+kappa+" " +
		   "priorityOption="+priorityOption+" " +
                   "minPriority="+minPriority+" " +
                   "minSupport="+minSupport+" " +
                   "minSize="+minSize+" " +
		   "maxSize="+maxSize+" " +
                   "minMGF="+minMGF+" " +
		   "maxPVal="+maxPVal+" " +
                   "maxRound="+maxRound+" " +
		   "maxClocktime="+maxClocktime+" " +
		   "keepOption="+keepOption+" " +
		   "includeNoCalls="+includeNoCalls+" " +
		   "requireBestNodeSet="+requireBestNodeSet+" " +
		   "requireSamePosition="+requireSamePosition+" " +
		   "requireHomozygous="+requireHomozygous+" " +
                   "requiredNodes="+requiredNodeString+" " +
		   "includedNodes="+includedNodeString+" " +
                   "excludedNodes="+excludedNodeString+" " +
		   "excludedPathNodes="+excludedPathNodeString+" " +
		   "includedPathNodes="+includedPathNodeString);
	

	// validate if requireSamePosition
	if (requireSamePosition && requiredNodes.size()>0 && !requiredNodes.haveSamePosition()) {
	    System.err.println("ERROR: requireSamePosition=true but required nodes do not have same position.");
	    System.exit(1);
	}

        // FR-finding round counter
        int round = 0;

        if (resume) {
            // resume from a previous run
            printToLog("# Resuming from previous run");
            String line = null;
            // allFrequentedRegions
            BufferedReader sfrReader = new BufferedReader(new FileReader(getGraphName()+"."+ALL_FREQUENTED_REGIONS_SAVE));
            while ((line=sfrReader.readLine())!=null) {
                String[] parts = line.split("\t");
                NodeSet nodes = graph.getNodeSet(parts[0]);
                allFrequentedRegions.put(nodes.toString(), new FrequentedRegion(graph, nodes, alpha, kappa, priorityOptionKey, priorityOptionLabel));
            }
	    // rejectedNodeSets
	    BufferedReader rnsReader = new BufferedReader(new FileReader(getGraphName()+"."+REJECTED_NODESETS_SAVE));
	    while ((line=rnsReader.readLine())!=null) {
		rejectedNodeSets.add(line);
	    }
            // frequentedRegions
	    // 0                                                                                                1   2         3  4       
	    // [1353,1355,1356,1357,1359,1360,1361,1363,1364,1366,1367,1368,...,1463,1464,1465,1467,1468,1469]	27  18621.00  1	 26
            BufferedReader frReader = new BufferedReader(new FileReader(getGraphName()+"."+FREQUENTED_REGIONS_SAVE));
	    line = frReader.readLine(); // header
            while ((line=frReader.readLine())!=null) {
                String[] parts = line.split("\t");
                NodeSet nodes = graph.getNodeSet(parts[0]);
                frequentedRegions.put(nodes.toString(), new FrequentedRegion(graph, nodes, alpha, kappa, priorityOptionKey, priorityOptionLabel));
                round++; // increment round so we start where we left off
            }
	    // informativationalism
            printToLog("# Loaded "+allFrequentedRegions.size()+" allFrequentedRegions.");
	    printToLog("# Loaded "+rejectedNodeSets.size()+" rejectedNodeSets.");
            printToLog("# Loaded "+frequentedRegions.size()+" frequentedRegions.");
            printToLog("# Now continuing with FR finding...");
        } else {
	    // load the single-node FRs into allFrequentedRegions
	    // - keep if gf>=minMGF
	    // - keep if p<=maxPVal
	    // - reject those with no genotype call unless includeNoCalls
	    // - keep those not in excludedNodes
	    // - exclude those with insufficient support if alpha=1.0
	    // - reject those with HET call if requireHomozygous
	    ConcurrentSkipListSet<Node> excRejects = new ConcurrentSkipListSet<>();
	    ConcurrentSkipListSet<Node> ngcRejects = new ConcurrentSkipListSet<>();
	    ConcurrentSkipListSet<Node> mgfRejects = new ConcurrentSkipListSet<>();
	    ConcurrentSkipListSet<Node> pvalRejects = new ConcurrentSkipListSet<>();
	    ConcurrentSkipListSet<Node> supportRejects = new ConcurrentSkipListSet<>();
	    ConcurrentSkipListSet<Node> homRejects = new ConcurrentSkipListSet<>();
	    ConcurrentSkipListSet<Node> nodes = new ConcurrentSkipListSet<>(graph.getNodes());
	    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	    nodes.parallelStream().forEach(node -> {
		    if (excludedNodes.contains(node)) {
			excRejects.add(node);
		    } else if (!includeNoCalls && !node.isCalled) {
			ngcRejects.add(node);
		    } else if (requireHomozygous && !node.isHomozygous()) {
			homRejects.add(node);
		    } else if (node.gf<minMGF) {
			mgfRejects.add(node);
		    } else if (maxPVal<1.0 && graph.fisherExactP(node)>maxPVal) {
			pvalRejects.add(node);
		    } else {
			FrequentedRegion fr = new FrequentedRegion(graph, new NodeSet(node), alpha, kappa, priorityOptionKey, priorityOptionLabel);
			if (alpha==1.0 && fr.support<minSupport) {
			    supportRejects.add(node);
			} else {
			    allFrequentedRegions.put(fr.nodes.toString(), fr);
			}
		    }
		});
	    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	    // DX
	    printToLog("# "+excRejects.size()+" nodes excluded because contained in excludedNodes");
	    printToLog("# "+ngcRejects.size()+" nodes excluded because not called and includeNoCalls=false");
	    printToLog("# "+mgfRejects.size()+" nodes excluded because allele frequency<"+minMGF);
	    printToLog("# "+pvalRejects.size()+" nodes excluded because p>"+maxPVal);
	    if (alpha==1.0) printToLog("# "+supportRejects.size()+" nodes excluded because FR support<"+minSupport);
	    if (requireHomozygous) printToLog("# "+homRejects.size()+" nodes excluded because not homozygous");
	    // store interesting single-node FRs in round 0, since we won't hit them in the loop
	    for (FrequentedRegion fr : allFrequentedRegions.values()) {
		if (isInteresting(fr)) {
		    if (requiredNodes.size()==0 && includedNodes.size()==0) {
			frequentedRegions.put(fr.nodes.toString(), fr);
		    } else {
			for (Node n : requiredNodes) {
			    if (fr.nodes.contains(n)) frequentedRegions.put(fr.nodes.toString(), fr);
			}
			for (Node n : includedNodes) {
			    if (fr.nodes.contains(n)) frequentedRegions.put(fr.nodes.toString(), fr);
			}
		    }
		}
	    }
        }

        // add the full required nodes FR to allFrequentedRegions, and frequentedRegions if interesting
	if (requiredNodes.size()>0) {
	    FrequentedRegion requiredFR = new FrequentedRegion(graph, requiredNodes, alpha, kappa, priorityOptionKey, priorityOptionLabel);
	    allFrequentedRegions.put(requiredFR.nodes.toString(), requiredFR);
	    if (isInteresting(requiredFR)) {
		frequentedRegions.put(requiredFR.nodes.toString(), requiredFR);
	    }
	}

	// add the full included nodes FR to allFrequentedRegions, and frequentedRegions if interesting
	if (includedNodes.size()>0) {
	    FrequentedRegion includedFR = new FrequentedRegion(graph, includedNodes, alpha, kappa, priorityOptionKey, priorityOptionLabel);
	    allFrequentedRegions.put(includedFR.nodes.toString(), includedFR);
	    if (isInteresting(includedFR)) {
		frequentedRegions.put(includedFR.nodes.toString(), includedFR);
	    }
	}

        // DX: dump out the interesting single-node FRs sorted by priority
        printToLog("# "+allFrequentedRegions.size()+" single-node FRs will be used to initiate search.");
        TreeSet<FrequentedRegion> sortedFRs = new TreeSet<>(frequentedRegions.values());
        for (FrequentedRegion fr : sortedFRs) {
            printToLog("0:"+fr.toString());
        }

	// used if requireBestNodeSet
	FrequentedRegion bestFR = null;

        // build the FRs round by round
	long startTime = System.currentTimeMillis();
	boolean clocktimeExceeded = false;
        boolean added = true;
        while (added && (maxRound==0 || round<maxRound) && !clocktimeExceeded) {
	    round++;
            added = false;
	    long roundStartTime = System.currentTimeMillis();
            // store FRPairs in a map keyed by merged nodes in THIS round for parallel operation and sorting
            ConcurrentSkipListSet<FRPair> interestingFRPairs = new ConcurrentSkipListSet<>();
	    // requiredNodes and bestFR need to be final for the parallel stream
	    final NodeSet finalRequiredNodes = requiredNodes;
	    final FrequentedRegion finalBestFR = bestFR;
	    // NOTE: the fr1 loop need not be parallel since the fr2 loop will consume all the CPUs anyway;
	    // by running the fr1 loop in series we ensure that the fr2 loop cycles use the same fr1 loop data.
	    for (FrequentedRegion fr1 : allFrequentedRegions.values()) {
		if (finalBestFR==null) {
		    if (finalRequiredNodes.size()>0) {
			boolean contains = true;
			for (Node n : finalRequiredNodes) {
			    if (!fr1.nodes.contains(n)) {
				contains = false;
				break;
			    }
			}
			if (!contains) continue;
		    }
		} else if (!fr1.equals(finalBestFR)) {
		    continue;
		}
		// abort if max clock time exceeded
		long clocktime = System.currentTimeMillis() - startTime;
		if (maxClocktime!=0 && clocktime>maxClocktime*60000) {
		    clocktimeExceeded = true;
		    System.err.println("Maximum clock time of "+maxClocktime+" minutes exceeded.");
		    break;
		}
		////////////////////////////////////////////////////////////////////////////////////////////////
		// start fr2 parallel stream
		allFrequentedRegions.entrySet().parallelStream().forEach(entry2 -> {
			FrequentedRegion fr2 = entry2.getValue();
			FRPair frpair = new FRPair(fr1, fr2, priorityOptionKey, priorityOptionLabel);
			String nodesKey = frpair.nodes.toString();
			boolean rejected = false;
			if (finalBestFR!=null && fr2.size>1) {
			    // we just append one node at a time in this mode
			    rejected = true;
			} else if (frequentedRegions.containsKey(nodesKey)) {
			    // if already found, bail
			    rejected = true;
			} else if (acceptedFRPairs.containsKey(nodesKey)) {
			    // get stored FRPair since already accepted and merged in a previous round
			    frpair = acceptedFRPairs.get(nodesKey);
			} else if (rejectedNodeSets.contains(nodesKey)) {
			    // already rejected, bail
			    rejected = true;
			} else if (requireSamePosition && !frpair.nodes.haveSamePosition()) {
			    // can't have an FR with nodes at different positions
			    rejected = true;
			}
			// reject if not all the required nodes are present
			if (!rejected && finalRequiredNodes.size()>0 && finalBestFR==null) {
			    for (Node n : finalRequiredNodes) {
				if (!frpair.nodes.contains(n)) {
				    rejected = true; // lacks a required node
				    break;
				}
			    }
			}
			// reject if one of the excluded nodes is present
			if (!rejected && excludedNodes.size()>0) {
			    for (Node n : excludedNodes) {
				if (frpair.nodes.contains(n)) {
				    rejected = true; // contains an excluded node
				    break;
				}
			    }
			}
			// reject if not any of the included nodes is present
			if (!rejected && includedNodes.size()>0) {
			    rejected = true;
			    for (Node n : includedNodes) {
				if (frpair.nodes.contains(n)) {
				    rejected = false; // contains an included node
				    break;
				}
			    }
			}
			// NOW WE MERGE
			if (rejected) {
			    // add this rejected NodeSet to the rejected list
			    rejectedNodeSets.add(nodesKey);
			} else {
			    // merge this FRPair if not already merged
			    if (frpair.merged==null) {
				frpair.merge();
			    }
			    // should we keep this merged FR according to keepOption?
			    boolean keep = true; // default if keepOption not set
			    if (keepOption.startsWith("subset") && frpair.merged.size>=keepOptionKey) {
				// keep FRs that are subsets of others or have higher priority
				for (FrequentedRegion frOld : frequentedRegions.values()) {
				    if (frpair.merged.nodes.isSupersetOf(frOld.nodes) && frpair.merged.priority<=frOld.priority) {
					keep = false;
					break;
				    }
				}
			    } else if (keepOption.startsWith("distance")) {
				// keep FRs that are at least a distance of keepOptionKey away from others or have higher priority
				for (FrequentedRegion frOld : frequentedRegions.values()) {
				    if (frpair.merged.nodes.distanceFrom(frOld.nodes)<keepOptionKey && frpair.merged.priority<=frOld.priority) {
					keep = false;
					break;
				    }
				}
			    }
			    if (keep) {
				// add this keeper pair to acceptedFRPairs
				acceptedFRPairs.put(nodesKey, frpair);
				// add this pair to the current interestingFRPairs if interesting
				if (isInteresting(frpair.merged)) {
				    interestingFRPairs.add(frpair);
				    if (debug) System.err.println("int:"+frpair.merged.toString()+"|"+(System.currentTimeMillis()-roundStartTime)/1000);
				}
			    } else {
				// add this to the rejected list
				rejectedNodeSets.add(nodesKey);
			    }
			}
		    });
		// end fr2 parallelStream
		////////////////////////////////////////////////////////////////////////////////////////////////
	    }
	    // end of fr1 loop
	    // add our new best merged FR from this round and store the interesting merged FRs
            if (interestingFRPairs.size()>0) {
                added = true;
		// add all of this round's interesting FRs to allFrequentedRegions OUTSIDE of the fr2 loop
		for (FRPair pair : interestingFRPairs) {
		    allFrequentedRegions.put(pair.merged.nodes.toString(), pair.merged);
		}
		// add the best FR to the output FRs map
                FrequentedRegion fr = interestingFRPairs.last().merged;
                frequentedRegions.put(fr.nodes.toString(), fr);
                // toggle priorityOptionLabel for next round if so desired
                if (priorityOptionParameter!=null && priorityOptionParameter.equals("alt")) togglePriorityOptionLabel();
                // output the best FR from this round
                printToLog(round+":"+fr.toString());
		// update requiredNodes if requireBestNodeSet
		if (requireBestNodeSet) {
		    bestFR = fr;
		    requiredNodes = fr.nodes;
		}
            } else {
                // show the top remaining FR that wasn't added
                TreeSet<FRPair> remainingFRPairs = new TreeSet<>(acceptedFRPairs.values());
		for (FRPair pair : remainingFRPairs.descendingSet()) {
		    if (!frequentedRegions.containsKey(pair.merged.nodes.toString())) {
			printToLog("-------------------------------------------------------------------------------------------------------");
			printToLog("TR:"+pair.merged.toString());
			break;
		    }
		}
            }
            // output current state for continuation if aborted
            if (frequentedRegions.size()>0 && writeSaveFiles()) {
                // params with current clock time
                FRUtils.printParameters(parameters, getGraphName()+".save", alpha, kappa, System.currentTimeMillis()-startTime);
                // allFrequentedRegions
                PrintStream sfrOut = new PrintStream(getGraphName()+"."+ALL_FREQUENTED_REGIONS_SAVE);
                for (FrequentedRegion fr : allFrequentedRegions.values()) {
                    sfrOut.println(fr.toString());
                }
                sfrOut.close();
                // acceptedFRPairs
                PrintStream afrpOut = new PrintStream(getGraphName()+"."+ACCEPTED_FRPAIRS_SAVE);
                for (FRPair frpair : acceptedFRPairs.values()) {
                    afrpOut.println(frpair.toString());
                }
                // rejectedNodeSets
                PrintStream rnsOut = new PrintStream(getGraphName()+"."+REJECTED_NODESETS_SAVE);
                for (String nodesKey : rejectedNodeSets) {
                    rnsOut.println(nodesKey);
                }
                rnsOut.close();
                // frequentedRegions
                PrintStream frOut = new PrintStream(getGraphName()+"."+FREQUENTED_REGIONS_SAVE);
                boolean first = true;
                for (FrequentedRegion fr : frequentedRegions.values()) {
                    if (first) {
                        frOut.println(fr.columnHeading()); // header
                        first = false;
                    }
                    frOut.println(fr.toString());
                }
                frOut.close();
            }
        }

        // timing
	long elapsedTime = System.currentTimeMillis() - startTime;
        printToLog("Found "+frequentedRegions.size()+" FRs.");
	printToLog("Clock time: "+FRUtils.formatTime(elapsedTime));
        
	// final output
	if (frequentedRegions.size()>0) {
            FRUtils.printParameters(parameters, formOutputPrefix(alpha, kappa, initialRequiredNodes), alpha, kappa, elapsedTime);
            printFrequentedRegions(formOutputPrefix(alpha, kappa, initialRequiredNodes));
            if (writePathFRs) {
		System.err.println("Writing path FRs file...");
		printPathFRs(formOutputPrefix(alpha, kappa, initialRequiredNodes));
	    }
            if (writeFRSubpaths) {
		System.err.println("Writing FR subpaths file...");
		printFRSubpaths(formOutputPrefix(alpha, kappa, initialRequiredNodes));
	    }
	}
    }

    // parameters file value getters
    public boolean writeSaveFiles() {
        return Boolean.parseBoolean(parameters.getProperty("writeSaveFiles"));
    }
    public boolean writePathFRs() {
        return Boolean.parseBoolean(parameters.getProperty("writePathFRs"));
    }
    public boolean writeFRSubpaths() {
        return Boolean.parseBoolean(parameters.getProperty("writeFRSubpaths"));
    }
    public boolean resume() {
        return Boolean.parseBoolean(parameters.getProperty("resume"));
    }
    public int getMinSupport() {
        return Integer.parseInt(parameters.getProperty("minSupport"));
    }
    public int getMinSize() {
        return Integer.parseInt(parameters.getProperty("minSize"));
    }
    public int getMaxSize() {
	return Integer.parseInt(parameters.getProperty("maxSize"));
    }
    public String getPriorityOption() {
        return parameters.getProperty("priorityOption");
    }
    public String getGraphName() {
        return parameters.getProperty("graphName");
    }
    public int getMaxRound() {
        return Integer.parseInt(parameters.getProperty("maxRound"));
    }
    public int getMaxClocktime() {
	return Integer.parseInt(parameters.getProperty("maxClocktime"));
    }
    public int getMinPriority() {
        return Integer.parseInt(parameters.getProperty("minPriority"));
    }
    public String getRequiredNodes() {
        return parameters.getProperty("requiredNodeString");
    }
    public String getIncludedNodes() {
	return parameters.getProperty("includedNodeString");
    }
    public String getExcludedNodes() {
        return parameters.getProperty("excludedNodeString");
    }
    public String getExcludedPathNodes() {
	return parameters.getProperty("excludedPathNodeString");
    }
    public String getIncludedPathNodes() {
	return parameters.getProperty("includedPathNodeString");
    }
    public String getKeepOption() {
        return parameters.getProperty("keepOption");
    }
    public double getMinMGF() {
        return Double.parseDouble(parameters.getProperty("minMGF"));
    }
    public double getMaxPVal() {
	return Double.parseDouble(parameters.getProperty("maxPVal"));
    }
    public boolean getRequireBestNodeSet() {
	return Boolean.parseBoolean(parameters.getProperty("requireBestNodeSet"));
    }
    public boolean getIncludeNoCalls() {
	return Boolean.parseBoolean(parameters.getProperty("includeNoCalls"));
    }
    public boolean getRequireSamePosition() {
	return Boolean.parseBoolean(parameters.getProperty("requireSamePosition"));
    }
    
    // parameter setters - set instance vars as well as value in parameters
    public void setPriorityOption(String priorityOption) {
	this.priorityOption = priorityOption;
        parsePriorityOption(priorityOption);
        parameters.setProperty("priorityOption", priorityOption);
    }
    public void setResume() {
	this.resume = true;
        parameters.setProperty("resume", "true");
    }
    public void setWriteSaveFiles() {
	this.writeSaveFiles = true;
        parameters.setProperty("writeSaveFiles", "true");
    }
    public void setWritePathFRs() {
	this.writePathFRs = true;
        parameters.setProperty("writePathFRs", "true");
    }
    public void setWriteFRSubpaths() {
	this.writeFRSubpaths = true;
        parameters.setProperty("writeFRSubpaths", "true");
    }
    public void setMinSupport(int minSupport) {
	this.minSupport = minSupport;
        parameters.setProperty("minSupport", String.valueOf(minSupport));
    }
    public void setMinSize(int minSize) {
	this.minSize = minSize;
        parameters.setProperty("minSize", String.valueOf(minSize));
    }
    public void setMaxSize(int maxSize) {
	this.maxSize = maxSize;
	parameters.setProperty("maxSize", String.valueOf(maxSize));
    }
    public void setMaxRound(int maxRound) {
	this.maxRound = maxRound;
        parameters.setProperty("maxRound", String.valueOf(maxRound));
    }
    public void setMaxClocktime(int maxClocktime) {
	this.maxClocktime = maxClocktime;
	parameters.setProperty("maxClocktime", String.valueOf(maxClocktime));
    }
    public void setMinPriority(int minPriority) {
	this.minPriority = minPriority;
        parameters.setProperty("minPriority", String.valueOf(minPriority));
    }
    public void setRequiredNodes(String requiredNodeString) {
	this.requiredNodeString = requiredNodeString;
        parameters.setProperty("requiredNodeString", requiredNodeString);
    }
    public void setIncludedNodes(String includedNodeString) {
	this.includedNodeString = includedNodeString;
	parameters.setProperty("includedNodeString", includedNodeString);
    }
    public void setExcludedNodes(String excludedNodeString) {
	this.excludedNodeString = excludedNodeString;
        parameters.setProperty("excludedNodeString", excludedNodeString);
    }
    public void setExcludedPathNodes(String excludedPathNodeString) {
	this.excludedPathNodeString = excludedPathNodeString;
	parameters.setProperty("excludedPathNodeString", excludedPathNodeString);
    }
    public void setIncludedPathNodes(String includedPathNodeString) {
	this.includedPathNodeString = includedPathNodeString;
	parameters.setProperty("includedPathNodeString", includedPathNodeString);
    }
    public void setKeepOption(String keepOption) {
	this.keepOption = keepOption;
        if (keepOption.startsWith("subset") || keepOption.startsWith("distance")) {
            String[] parts = keepOption.split(":");
            if (parts.length==2) {
                keepOptionKey = Integer.parseInt(parts[1]);
            } else if (keepOption.equals("subset")) {
                keepOptionKey = 3; // default
            } else if (keepOption.equals("distance")) {
                keepOptionKey = 2; // default
            }
	    parameters.setProperty("keepOption", keepOption);
        } else {
            System.err.println("ERROR: allowed keepoption values are: subset[:N]|distance[:N]");
            System.exit(1);
        }
    }
    public void setGraphName(String graphName) {
	this.graphName = graphName;
        parameters.setProperty("graphName", graphName);
    }
    public void setMinMGF(double minMGF) {
	this.minMGF = minMGF;
        parameters.setProperty("minMGF", String.valueOf(minMGF));
    }
    public void setMaxPVal(double maxPVal) {
	this.maxPVal = maxPVal;
	parameters.setProperty("maxPVal", String.valueOf(maxPVal));
    }
    public void setRequireBestNodeSet() {
	this.requireBestNodeSet = true;
	parameters.setProperty("requireBestNodeSet", "true");
    }
    public void setIncludeNoCalls() {
	this.includeNoCalls = true;
	parameters.setProperty("includeNoCalls", "true");
    }
    public void setRequireSamePosition() {
	this.requireSamePosition = true;
	parameters.setProperty("requireSamePosition", "true");
    }
    public void setRequireHomozygous() {
	this.requireHomozygous = true;
	parameters.setProperty("requireHomozygous", "true");
    }
    
    /**
     * Print to both the console and the log file
     */
    void printToLog(String text) {
        System.out.println(text);
        logOut.println(text);
    }

    /**
     * Print a delineating heading, for general use.
     */
    static void printHeading(String heading) {
        for (int i=0; i<heading.length(); i++) System.out.print("="); System.out.println("");
        System.out.println(heading);
        for (int i=0; i<heading.length(); i++) System.out.print("="); System.out.println("");
    }

    /**
     * Print the path names and the count of subpaths for each FR.
     * This can be used as input to a classification routine.
     */
    void printPathFRs(String outputPrefix) throws IOException {
        PrintStream out = new PrintStream(FRUtils.getPathFRsFilename(outputPrefix));
        // columns are paths
        boolean first = true;
        for (Path path : paths) {
            if (first) {
                first = false;
            } else {
                out.print("\t");
            }
            out.print(path.name+"."+path.label);
        }
        out.println("");
        // rows are FRs
        int c = 1;
        for (FrequentedRegion fr : frequentedRegions.values()) {
            out.print("FR"+(c++));
            for (Path path : paths) {
                out.print("\t"+fr.countSubpathsOf(path));
            }
            out.println("");
        }
        out.close();
    }

    /**
     * Print out the FRs, in order of priority.
     */
    void printFrequentedRegions(String outputPrefix) throws IOException {
        if (frequentedRegions.size()==0) {
            System.err.println("NO FREQUENTED REGIONS!");
            return;
        }
        PrintStream out = new PrintStream(FRUtils.getFRsFilename(outputPrefix));
        boolean first = true;
        TreeSet<FrequentedRegion> sortedFRs = new TreeSet<>(frequentedRegions.values());
        for (FrequentedRegion fr : sortedFRs) {
            if (first) {
                out.println(fr.columnHeading());
                first = false;
            }
            out.println(fr.toString());
        }
        out.close();
    }
    
    /**
     * Print out the FRs along with their subpaths, strictly to an output file.
     */
    void printFRSubpaths(String outputPrefix) throws IOException {
        if (frequentedRegions.size()==0) {
            System.err.println("NO FREQUENTED REGIONS!");
            return;
        }
        PrintStream out = new PrintStream(FRUtils.getFRSubpathsFilename(outputPrefix));
        TreeSet<FrequentedRegion> sortedFRs = new TreeSet<>(frequentedRegions.values());
        for (FrequentedRegion fr : sortedFRs) {
            out.println(fr.toString());
            out.print(fr.subpathsString()); // contains \n at end of every line
        }
        out.close();
    }

    /**
     * Print a crude histogram of FR node sizes.
     */
    public void printFRHistogram() {
        Map<Integer,Integer> countMap = new TreeMap<>();
        for (FrequentedRegion fr : frequentedRegions.values()) {
            if (countMap.containsKey(fr.size)) {
                int count = countMap.get(fr.size);
                count++;
                countMap.put(fr.size, count);
            } else {
                countMap.put(fr.size, 1);
            }
        }
        for (int num : countMap.keySet()) {
            System.out.println("FR node size (#):"+num+" ("+countMap.get(num)+")");
        }
    }

    /**
     * Form an outputPrefix with given alpha and kappa and single required node, if required.
     */
    String formOutputPrefix(double alpha, int kappa, NodeSet requiredNodes) {
        DecimalFormat alphaFormat = new DecimalFormat("0.0");
	String prefix = "";
        if (kappa==Integer.MAX_VALUE) {
            prefix = getGraphName()+"-"+alphaFormat.format(alpha)+"-Inf";
        } else {
            prefix = getGraphName()+"-"+alphaFormat.format(alpha)+"-"+kappa;
        }
	if (requiredNodes.size()==1) {
	    prefix += "-"+requiredNodes.first();
	}
	return prefix;
    }

    /**
     * Return true if the given FR is "interesting". (Note that the alpha, kappa requirements are enforced by fr.support>=minSupport.)
     * This also uses priorityOptionLabel for the O.R- and p-based priorities.
     */
    boolean isInteresting(FrequentedRegion fr) {
        boolean interesting = fr.support>=minSupport;
	interesting = interesting && fr.size>=minSize;
	interesting = interesting && fr.size<=maxSize;
	if (minPriority!=0) interesting = interesting && fr.priority>=minPriority;
        // if (priorityOptionKey==3 || priorityOptionKey==4) {
        //     if (priorityOptionLabel!=null && priorityOptionLabel.equals("case")) {
        //         interesting = interesting && fr.oddsRatio()>1.0;
        //     } else if (priorityOptionLabel!=null && priorityOptionLabel.equals("ctrl")) {
        //         interesting = interesting && fr.oddsRatio()<1.0;
        //     }
        // }
        return interesting;   
    }

    /**
     * Toggle the current priority option label between "case" and "ctrl". (Used with "alt" priority option parameter.)
     */
    void togglePriorityOptionLabel() {
        if (priorityOptionLabel==null) {
            priorityOptionLabel = "case";
        } else if (priorityOptionLabel.equals("ctrl")) {
            priorityOptionLabel = "case";
        } else {
            priorityOptionLabel = "ctrl";
        }
    }

    /**
     * Parse out the priority option key and parameter like "case", "ctrl", "alt", plus set the label if needed.
     */
    void parsePriorityOption(String priorityOption) {
        String[] parts = priorityOption.split(":");
        priorityOptionKey = Integer.parseInt(parts[0]);
        if (parts.length>1) {
            priorityOptionParameter = parts[1];
            if (priorityOptionParameter.equals("case") || priorityOptionParameter.equals("ctrl")) {
                priorityOptionLabel = priorityOptionParameter;
            }
        } else {
            priorityOptionParameter = null;
            priorityOptionLabel = null;
        }
        // impose defaults
        if (priorityOptionKey==1 && priorityOptionLabel==null) priorityOptionLabel = "case";
    }

    /**
     * Command-line utility
     */
    public static void main(String[] args) throws FileNotFoundException, IOException {
        Options options = new Options();
        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd;

        Option alphaOption = new Option("a", "alpha", true, "value of alpha (required)");
        alphaOption.setRequired(true);
        options.addOption(alphaOption);
        //
        Option kappaOption = new Option("k", "kappa", true, "value of kappa; -1=infinity (required)");
        kappaOption.setRequired(true);
        options.addOption(kappaOption);
        //
        Option graphOption = new Option("graph", "graph", true, "graph name");
        graphOption.setRequired(true);
        options.addOption(graphOption);
        //
        Option txtOption = new Option("txt", "txt", false, "load from [graph].nodes.txt and [graph].paths.txt");
        txtOption.setRequired(false);
        options.addOption(txtOption);
        //
        Option minSupportOption = new Option("m", "minsupport", true, "minimum number of supporting paths for an FR to be considered interesting [1]");
        minSupportOption.setRequired(false);
        options.addOption(minSupportOption);
        //
        Option minSizeOption = new Option("s", "minsize", true, "minimum number of nodes that a FR must contain to be considered interesting [1]");
        minSizeOption.setRequired(false);
        options.addOption(minSizeOption);
	//
        Option maxSizeOption = new Option("maxs", "maxsize", true, "maximum number of nodes that a FR can contain to be considered interesting [unlimited]");
        maxSizeOption.setRequired(false);
        options.addOption(maxSizeOption);
	//
        Option labelsOption = new Option("p", "pathlabels", true, "tab-delimited file with pathname<tab>label");
        labelsOption.setRequired(false);
        options.addOption(labelsOption);
        //
        Option verboseOption = new Option("v", "verbose", false, "verbose output [false]");
        verboseOption.setRequired(false);
        options.addOption(verboseOption);
        //
        Option debugOption = new Option("do", "debug", false, "debug output [false]");
        debugOption.setRequired(false);
        options.addOption(debugOption);
        //
        Option priorityOptionOption = new Option("pri", "priorityoption", true, "option for priority weighting of FRs: "+FrequentedRegion.PRIORITY_OPTIONS);
        priorityOptionOption.setRequired(true);
        options.addOption(priorityOptionOption);
        //
        Option resumeOption = new Option("r", "resume", false, "resume from a previous run [false]");
        resumeOption.setRequired(false);
        options.addOption(resumeOption);
        //
        Option maxRoundOption = new Option("mr", "maxround", true, "maximum FR-finding round to run [0=unlimited]");
        maxRoundOption.setRequired(false);
        options.addOption(maxRoundOption);
        //
        Option minPriorityOption = new Option("mp", "minpriority", true, "minimum priority for an FR to be considered interesting [0=disabled]");
        minPriorityOption.setRequired(false);
        options.addOption(minPriorityOption);
        //
        Option writeSaveFilesOption = new Option("wsf", "writesavefiles", false, "write save files after each FR is found [false]");
        writeSaveFilesOption.setRequired(false);
        options.addOption(writeSaveFilesOption);
        //
        Option writePathFRsOption = new Option("wpfr", "writepathfrs", false, "write the path FR file at end [false]");
        writePathFRsOption.setRequired(false);
        options.addOption(writePathFRsOption);
        //
        Option writeFRSubpathsOption = new Option("wfrs", "writefrsubpaths", false, "write the FR subpaths file at end [false]");
        writeFRSubpathsOption.setRequired(false);
        options.addOption(writeFRSubpathsOption);
        //
        Option requiredNodesOption = new Option("rn", "requirednodes", true, "require that interesting FRs contain the given nodes []");
        requiredNodesOption.setRequired(false);
        options.addOption(requiredNodesOption);
        //
	Option requiredNodeScanOption = new Option("rns", "requirednodescan", true, "starting and ending node for a single required node scan, e.g. 1-12345");
	requiredNodeScanOption.setRequired(false);
	options.addOption(requiredNodeScanOption);
	//
        Option includedNodesOption = new Option("in", "includednodes", true, "require that interesting FRs contain at least one of the given nodes []");
        includedNodesOption.setRequired(false);
        options.addOption(includedNodesOption);
        //
        Option excludedNodesOption = new Option("en", "excludednodes", true, "require that interesting FRs NOT contain the given nodes []");
        excludedNodesOption.setRequired(false);
        options.addOption(excludedNodesOption);
        //
        Option excludedPathNodesOption = new Option("ep", "excludedpathnodes", true, "exclude paths that include any of the given nodes []");
        excludedPathNodesOption.setRequired(false);
        options.addOption(excludedPathNodesOption);
        //
        Option includedPathNodesOption = new Option("ip", "includedpathnodes", true, "include only paths that include at least one of the given nodes []");
        includedPathNodesOption.setRequired(false);
        options.addOption(includedPathNodesOption);
        //
        Option keepOptionOption = new Option("keep", "keepoption", true, "option for keeping FRs in finder loop: subset[:N]|distance[:N] [keep all]");
        keepOptionOption.setRequired(false);
        options.addOption(keepOptionOption);
        //
        Option minMGFOption = new Option("minmgf", "minmgf", true, "minimum MGF for nodes included in search [0.01]");
        minMGFOption.setRequired(false);
        options.addOption(minMGFOption);
	//
        Option maxPValOption = new Option("maxp", "maxpval", true, "maximum p-value for nodes included in search [1.0]");
        maxPValOption.setRequired(false);
        options.addOption(maxPValOption);
	//
	Option requireBestNodeSetOption = new Option("rbns", "requirebestnodeset", false, "require the best NodeSet from the previous round in the next round [false]");
	requireBestNodeSetOption.setRequired(false);
	options.addOption(requireBestNodeSetOption);
	//
	Option includeNoCallsOption = new Option("inc", "includenocalls", false, "include nodes which do not have genotype calls (./.) [false]");
	includeNoCallsOption.setRequired(false);
	options.addOption(includeNoCallsOption);
	//
	Option requireSamePositionOption = new Option("rsp", "requiresameposition", false, "require that interesting FRs consist of nodes at same genomic position [false]");
	requireSamePositionOption.setRequired(false);
	options.addOption(requireSamePositionOption);
	//
	Option maxClocktimeOption = new Option("maxct", "maxclocktime", true, "limit the computation to the given clock time in minutes [0=unlimited]");
	maxClocktimeOption.setRequired(false);
	options.addOption(maxClocktimeOption);
	//
	Option requireHomozygousOption = new Option("rh", "requirehomozygous", false, "require that FR nodes be homozygous genotypes");
	requireHomozygousOption.setRequired(false);
	options.addOption(requireHomozygousOption);
	
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println(e.getMessage());
            formatter.printHelp("FRFinder", options);
            System.exit(1);
            return;
        }

        if (cmd.getOptions().length==0) {
            formatter.printHelp("FRFinder", options);
            System.exit(1);
            return;
        }

        // path labels file
        File labelsFile = null;
        if (cmd.hasOption("pathlabels")) labelsFile = new File(cmd.getOptionValue("pathlabels"));

        // required run parameters
        double alpha = Double.parseDouble(cmd.getOptionValue("alpha"));
        int kappa = Integer.parseInt(cmd.getOptionValue("kappa"));
	if (kappa==-1) kappa = Integer.MAX_VALUE; // effectively infinity

	// check syntax of requirednodescan if present
	int requiredNodeStart = 0;
	int requiredNodeEnd = 0;
	if (cmd.hasOption("requirednodescan")) {
	    if (cmd.getOptionValue("requirednodescan").contains("-")) {
		String[] parts = cmd.getOptionValue("requirednodescan").split("-");
		requiredNodeStart = Integer.parseInt(parts[0]);
		requiredNodeEnd = Integer.parseInt(parts[1]);
	    } else {
		System.err.println("ERROR: --requirednodescan requires two numbers separated by dash, like 23-45.");
		System.exit(1);
	    }
	}
	// load graph from a pair of TXT files
	String graphName = cmd.getOptionValue("graph");
	PangenomicGraph pg = new PangenomicGraph();
	pg.name = graphName;
	pg.nodesFile = new File(graphName+".nodes.txt");
	pg.pathsFile = new File(graphName+".paths.txt");
	if (cmd.hasOption("verbose")) pg.verbose = true;
	pg.loadTXT();
	// remove paths that contain an excluded path node, if there are any
	String excludedPathNodeString = "[]";
	if (cmd.hasOption("excludedpathnodes")) {
	    excludedPathNodeString = cmd.getOptionValue("excludedpathnodes");
	}
	NodeSet excludedPathNodes = pg.getNodeSet(excludedPathNodeString);
	if (excludedPathNodes.size()>0) {
	    List<Path> pathsToRemove = new LinkedList<>();
	    for (Path path : pg.paths) {
		for (Node node : excludedPathNodes) {
		    if (path.getNodes().contains(node)) {
			pathsToRemove.add(path);
			break;
		    }
		}
	    }
	    pg.paths.removeAll(pathsToRemove);
	    System.out.println("# Graph has had "+pathsToRemove.size()+" paths removed which contained excluded nodes.");
	}
	// limit to paths that contain an included node, if given
	String includedPathNodeString = "[]";
	if (cmd.hasOption("includedpathnodes")) {
	    includedPathNodeString = cmd.getOptionValue("includedpathnodes");
	}
	NodeSet includedPathNodes = pg.getNodeSet(includedPathNodeString);
	int formerPathCount = pg.paths.size();
	if (includedPathNodes.size()>0) {
	    Set<Path> pathsToKeep = new HashSet<>();
	    for (Path path : pg.paths) {
		for (Node node : includedPathNodes) {
		    if (path.contains(node)) {
			pathsToKeep.add(path);
			break;
		    }
		}
	    }
	    pg.paths = pathsToKeep;
	    int removedCount = formerPathCount - pg.paths.size();
	    System.out.println("# Graph has had "+removedCount+" paths removed which did not contain one of the included nodes.");
	}
	// other stuff
	pg.buildNodePaths();
	pg.tallyLabelCounts();
	System.out.println("# Graph has "+pg.vertexSet().size()+" nodes and "+pg.edgeSet().size()+" edges with "+pg.paths.size()+" paths.");
	System.out.println("# Graph has "+pg.labelCounts.get("case")+" case paths and "+pg.labelCounts.get("ctrl")+" ctrl paths.");
	// instantiate the FRFinder with this PangenomicGraph
	FRFinder frf = new FRFinder(pg);
	frf.setGraphName(graphName);
	// set optional FRFinder parameters
	if (cmd.hasOption("priorityoption")) frf.setPriorityOption(cmd.getOptionValue("priorityoption"));
	if (cmd.hasOption("minsupport")) frf.setMinSupport(Integer.parseInt(cmd.getOptionValue("minsupport")));
	if (cmd.hasOption("minsize")) frf.setMinSize(Integer.parseInt(cmd.getOptionValue("minsize")));
	if (cmd.hasOption("maxsize")) frf.setMaxSize(Integer.parseInt(cmd.getOptionValue("maxsize")));
	if (cmd.hasOption("maxround")) frf.setMaxRound(Integer.parseInt(cmd.getOptionValue("maxround")));
	if (cmd.hasOption("maxclocktime")) frf.setMaxClocktime(Integer.parseInt(cmd.getOptionValue("maxclocktime")));
	if (cmd.hasOption("minpriority")) frf.setMinPriority(Integer.parseInt(cmd.getOptionValue("minpriority")));
	if (cmd.hasOption("requirednodes")) frf.setRequiredNodes(cmd.getOptionValue("requirednodes"));
	if (cmd.hasOption("includednodes")) frf.setIncludedNodes(cmd.getOptionValue("includednodes"));
	if (cmd.hasOption("excludednodes")) frf.setExcludedNodes(cmd.getOptionValue("excludednodes"));
	if (cmd.hasOption("excludedpathnodes")) frf.setExcludedPathNodes(cmd.getOptionValue("excludedpathnodes"));
	if (cmd.hasOption("includedpathnodes")) frf.setIncludedPathNodes(cmd.getOptionValue("includedpathnodes"));
	if (cmd.hasOption("keepoption")) frf.setKeepOption(cmd.getOptionValue("keepoption"));
	if (cmd.hasOption("resume")) frf.setResume();
	if (cmd.hasOption("writesavefiles")) frf.setWriteSaveFiles();
	if (cmd.hasOption("writepathfrs")) frf.setWritePathFRs();
	if (cmd.hasOption("writefrsubpaths")) frf.setWriteFRSubpaths();
	if (cmd.hasOption("minmgf")) frf.setMinMGF(Double.parseDouble(cmd.getOptionValue("minmgf")));
	if (cmd.hasOption("maxpval")) frf.setMaxPVal(Double.parseDouble(cmd.getOptionValue("maxpval")));
	if (cmd.hasOption("requirebestnodeset")) frf.setRequireBestNodeSet();
	if (cmd.hasOption("includenocalls")) frf.setIncludeNoCalls();
	if (cmd.hasOption("requiresameposition")) frf.setRequireSamePosition();
	if (cmd.hasOption("requirehomozygous")) frf.setRequireHomozygous();
	// these are not stored in parameters
	if (cmd.hasOption("verbose")) frf.verbose = true;
	if (cmd.hasOption("debug")) frf.debug = true;
	// run the requested job(s)
	if (requiredNodeStart>0 && requiredNodeEnd>0) {
	    for (long id=requiredNodeStart; id<=requiredNodeEnd; id++) {
		if (pg.nodeIdMap.containsKey(id)) {
		    Node n = pg.getNode(id);
		    if (!n.isNoCall()) {
			frf.setRequiredNodes("["+id+"]");
			frf.findFRs(alpha, kappa);
		    }
		}
	    }
	} else {
	    frf.findFRs(alpha, kappa);
	}
    }
}
