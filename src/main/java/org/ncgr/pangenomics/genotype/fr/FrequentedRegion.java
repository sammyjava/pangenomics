package org.ncgr.pangenomics.genotype.fr;

import org.ncgr.pangenomics.genotype.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import java.text.DecimalFormat;

import java.util.List;
import java.util.LinkedList;
import java.util.Map;
import java.util.TreeMap;
import java.util.Set;
import java.util.HashSet;
import java.util.TreeSet;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

/**
 * Represents a cluster of nodes along with the supporting subpaths of the full set of strain/subject/subspecies paths.
 *
 * @author Sam Hokin
 */
class FrequentedRegion implements Comparable {

    static String PRIORITY_OPTIONS =
        "0:label=total support or label support [null,case,ctrl], " +
        "1:label=(label support-other support) [case,ctrl,alt], " +
        "2=|case support-control support|, " +
        "3:label=O.R. in label's favor [null,case,ctrl], " +
        "4=Fisher's exact test two-tailed p value";

    // static utility stuff
    static DecimalFormat df = new DecimalFormat("0.00");
    static DecimalFormat pf = new DecimalFormat("0.00E0");
    static DecimalFormat orf = new DecimalFormat("0.000");
    
    // the PangenomicGraph that this FrequentedRegion belongs to
    PangenomicGraph graph;

    // the set of Nodes that encompass this FR and the size
    NodeSet nodes;
    int size;
    
    // the subpaths, identified by their originating path name and label, that start and end on this FR's nodes
    List<Path> subpaths;
    
    // the subpath support of this FR
    int support = 0;

    // the case and control subpath support of this FR
    int caseSubpathSupport = 0;
    int ctrlSubpathSupport = 0;

    // a subpath must satisfy the requirement that it traverses at least a fraction alpha of this.nodes
    double alpha;

    // a subpath must satisfy the requirement that its contiguous nodes that do NOT belong in this.nodes have number no larger than kappa
    int kappa;

    // the priority and priority option for comparison
    int priority;                // the priority metric
    int priorityOptionKey;       // 0, 1, 2, etc.
    String priorityOptionLabel;  // label for priority update emphasis, can be null

    // the calculated p value and odds ratio (saved when methods called first time)
    double pValue = Double.NEGATIVE_INFINITY;
    double orValue = Double.NEGATIVE_INFINITY;

    // a number used to label FRs for easy usage, assigned when needed
    int number = 0;

    /**
     * Construct given a PangenomicGraph, NodeSet and alpha and kappa parameters.
     */
    FrequentedRegion(PangenomicGraph graph, NodeSet nodes, double alpha, int kappa, int priorityOptionKey, String priorityOptionLabel) {
        this.graph = graph;
        this.nodes = nodes;
        this.alpha = alpha;
        this.kappa = kappa;
        this.priorityOptionKey = priorityOptionKey;
        this.priorityOptionLabel = priorityOptionLabel;
	size = nodes.size();
        update();
    }

    /**
     * Construct given a PangenomicGraph, string representation of nodes, and alpha and kappa and priorityOption parameters.
     * 0                               1       2               3       4       5       6       7
     * [7,9,14,19,103,132,174] 7       3030    21105.00        1582    1448    1554    1.373   2.89E-16
     */
    FrequentedRegion(PangenomicGraph graph, String frString, double alpha, int kappa, int priorityOptionKey, String priorityOptionLabel) {
        String[] parts = frString.split("\t");
        String nodeString = parts[0];
        support = Integer.parseInt(parts[1]);
        if (parts.length>3) {
            caseSubpathSupport = Integer.parseInt(parts[3]);
            ctrlSubpathSupport = Integer.parseInt(parts[4]);
        }
        this.graph = graph;
        this.nodes = new NodeSet(nodeString);
        this.alpha = alpha;
        this.kappa = kappa;
        this.priorityOptionKey = priorityOptionKey;
        this.priorityOptionLabel = priorityOptionLabel;
	size = nodes.size();
        update();
    }

    /**
     * Construct given a PangenomicGraph, NodeSet and Subpaths
     */
    FrequentedRegion(PangenomicGraph graph, NodeSet nodes, List<Path> subpaths, double alpha, int kappa, int priorityOptionKey, String priorityOptionLabel) {
        this.graph = graph;
        this.nodes = nodes;
        this.subpaths = subpaths;
        this.alpha = alpha;
        this.kappa = kappa;
        this.priorityOptionKey = priorityOptionKey;
        this.priorityOptionLabel = priorityOptionLabel;
	size = nodes.size();
        update();
    }

    /**
     * Construct given a PangenomicGraph, NodeSet and Subpaths and already known support 
     */
    FrequentedRegion(PangenomicGraph graph, NodeSet nodes, List<Path> subpaths, double alpha, int kappa, int priorityOptionKey, String priorityOptionLabel, int support) {
        this.graph = graph;
        this.nodes = nodes;
        this.subpaths = subpaths;
        this.alpha = alpha;
        this.kappa = kappa;
        this.support = support;
        this.priorityOptionKey = priorityOptionKey;
        this.priorityOptionLabel = priorityOptionLabel;
	size = nodes.size();
        update();
    }

    /**
     * Construct given only basic information, used for post-processing. NO GRAPH.
     */
    FrequentedRegion(NodeSet nodes, List<Path> subpaths, double alpha, int kappa, int priorityOptionKey, String priorityOptionLabel, int support) {
        this.nodes = nodes;
        this.subpaths = subpaths;
        this.alpha = alpha;
        this.kappa = kappa;
        this.support = support;
        this.priorityOptionKey = priorityOptionKey;
        this.priorityOptionLabel = priorityOptionLabel;
	size = nodes.size();
        update();
    }

    /**
     * Construct given only the pieces in the FR.toString() output (and alpha, kappa).
     * NOTE: does not update.
     */
    FrequentedRegion(NodeSet nodes, double alpha, int kappa, int support, int caseSubpathSupport, int ctrlSubpathSupport, double orValue, double pValue, int priority) {
        this.nodes = nodes;
        this.alpha = alpha;
        this.kappa = kappa;
        this.support = support;
        this.caseSubpathSupport = caseSubpathSupport;
        this.ctrlSubpathSupport = ctrlSubpathSupport;
        this.orValue = orValue;
        this.pValue = pValue;
        this.priority = priority;
	size = nodes.size();
    }

    /**
     * Construct given a graph, alpha and kappa, and a line from an frs.txt file.
     * Note: does NOT update(); uses p-value, OR value and priority from file input. (These are rounded from true values.)
     *
     * 0nodes         1size 2support 3case 4ctrl 5OR    6p        7pri
     * [711,786,891]  3     104	     2     102   0.020  3.25E-28  1707
     */
    FrequentedRegion(PangenomicGraph graph, double alpha, int kappa, String line) {
	if (line.startsWith("[")) {
	    this.graph = graph;
	    this.alpha = alpha;
	    this.kappa = kappa;
	    // the line values
	    String[] fields = line.split("\t");
	    this.nodes = new NodeSet(fields[0]);
	    this.size = Integer.parseInt(fields[1]);
	    this.support = Integer.parseInt(fields[2]);
	    this.caseSubpathSupport = Integer.parseInt(fields[3]);
	    this.ctrlSubpathSupport = Integer.parseInt(fields[4]);
	    try {
		this.orValue = Double.parseDouble(fields[5]);
	    } catch (NumberFormatException e) {
		// do nothing, it's an infinity symbol
	    }
	    this.pValue = Double.parseDouble(fields[6]);
	    this.priority = Integer.parseInt(fields[7]);
	}
    }

    /**
     * Construct given a line in an frs.txt output file. Runs update() to get true double-valued p and OR.
     *
     * 0nodes         1size 2support 3case 4ctrl 5OR    6p        7pri
     * [711,786,891]  3     104	     2     102   0.020  3.25E-28  1707
     */
    FrequentedRegion(PangenomicGraph graph, double alpha, int kappa, int priorityOptionKey, String priorityOptionLabel, String line) {
	if (line.startsWith("[")) {
	    this.graph = graph;
	    this.alpha = alpha;
	    this.kappa = kappa;
	    // the line values
	    String[] fields = line.split("\t");
	    this.nodes = new NodeSet(fields[0]);
	    this.size = Integer.parseInt(fields[1]);
	    this.support = Integer.parseInt(fields[2]);
	    this.caseSubpathSupport = Integer.parseInt(fields[3]);
	    this.ctrlSubpathSupport = Integer.parseInt(fields[4]);
	    // set the priority option
	    this.priorityOptionKey = priorityOptionKey;
	    this.priorityOptionLabel = priorityOptionLabel;
	    // run update() to get exact p-value and OR values without roundoff from reading in.
	    update();
	}
    }

    /**
     * Construct given an output line from toString(), without setting the other values. Useful for sorting.
     * 0      1     2        3     4     5   6  7
     * nodes  size  support  case  ctrl  OR  p  pri
     */
    FrequentedRegion(String line) {
	if (line.startsWith("[")) {
	    String[] fields = line.split("\t");
	    this.nodes = new NodeSet(fields[0]);
	    this.size = Integer.parseInt(fields[1]);
	    this.support = Integer.parseInt(fields[2]);
	    this.caseSubpathSupport = Integer.parseInt(fields[3]);
	    this.ctrlSubpathSupport = Integer.parseInt(fields[4]);
	    this.orValue = Double.parseDouble(fields[5]);
	    this.pValue = Double.parseDouble(fields[6]);
	    this.priority = Integer.parseInt(fields[7]);
	}
    }

    /**
     * Update the subpaths, priority, etc.
     */
    void update() {
        updateSupport();
        updatePriority();
    }

    /**
     * Equality is based on nodes.
     */
    @Override
    public boolean equals(Object o) {
	FrequentedRegion that = (FrequentedRegion) o;
        return this.nodes.equals(that.nodes);
    }

    /**
     * Must override hashCode() for Map keys.
     */
    @Override
    public int hashCode() {
	return this.nodes.toString().hashCode();
    }

    /**
     * Comparison is based on higher priority, then higher support, then smaller size, then the nodes themselves.
     */
    public int compareTo(Object o) {
	FrequentedRegion that = (FrequentedRegion) o;
        if (this.priority!=that.priority) {
            return this.priority - that.priority;
        } else if (this.support!=that.support) {
            return this.support - that.support;
        } else if (this.size!=that.size) {
            return -(this.size - that.size);
        } else {
            return this.nodes.compareTo(that.nodes);
        }
    }

    /**
     * Update the nodes in this FR so they contain node attributes from the graph rather than just an id.
     */
    public void updateNodes() {
	for (Node n : nodes) {
	    n = graph.getNode(n.id);
	}
    }

    /**
     * Update the subpaths and support from the graph paths for the current alpha and kappa values.
     */
    void updateSupport() {
        subpaths = new LinkedList<>();
        for (Path p : graph.paths) {
            List<Path> supportPaths = computeSupport(p);
            subpaths.addAll(supportPaths);
        }
        support = subpaths.size();
        caseSubpathSupport = getSubpathSupport("case");
        ctrlSubpathSupport = getSubpathSupport("ctrl");
    }

    /**
     * Set a new alpha value.
     */
    public void setAlpha(double alpha) {
        this.alpha = alpha;
    }

    /**
     * Set a new kappa value.
     */
    public void setKappa(int kappa) {
        this.kappa = kappa;
    }

    /**
     * Return the column heading for the toString() fields
     * nodes size support case ctrl OR p pri
     */
    public static String columnHeading() {
	return "nodes\tsize\tsupport\tcase\tctrl\tOR\tp\tpri";
    }

    /**
     * Return the support associated with paths with the given label (not subpaths).
     */
    public int getPathSupport(String label) {
        int count = 0;
        List<String> countedPaths = new LinkedList<>();
        for (Path subpath : subpaths) {
            if (!countedPaths.contains(subpath.getName()) && subpath.getLabel()!=null && subpath.getLabel().equals(label)) {
                countedPaths.add(subpath.getName());
                count++;
            }
        }
        return count;
    }

    /**
     * Return the count of subpaths labeled with the given label. Here multiple subpaths of a path count individually.
     */
    public int getSubpathSupport(String label) {
        int count = 0;
        for (Path subpath : subpaths) {
            if (subpath.getLabel()!=null && subpath.getLabel().equals(label)) {
                count++;
            }
        }
        return count;
    }

    /**
     * Update the integer priority metric used for ordering FRs.
     */
    void updatePriority() {
        priority = 0;
	pValue = fisherExactP();
	orValue = oddsRatio();
        if (priorityOptionKey==0) {
	    // support
	    if (priorityOptionLabel==null) {
		priority = support;
	    } else if (priorityOptionLabel.equals("case")) {
		priority = caseSubpathSupport;
	    } else if (priorityOptionLabel.equals("ctrl")) {
		priority = ctrlSubpathSupport;
	    } else {
		throwPriorityOptionError();
	    }
        } else if (priorityOptionKey==1) {
	    // support difference
            if (priorityOptionLabel.equals("case")) {
                priority = caseSubpathSupport - ctrlSubpathSupport;
            } else if (priorityOptionLabel.equals("ctrl")) {
                priority = ctrlSubpathSupport - caseSubpathSupport;
            } else {
		throwPriorityOptionError();
            }
        } else if (priorityOptionKey==2) {
	    // abs support difference
            priority = Math.abs(caseSubpathSupport - ctrlSubpathSupport);
        } else if (priorityOptionKey==3) {
	    // odds ratio
	    int mult = 1000; // multiplier to convert to integer
            if (caseSubpathSupport>0 && ctrlSubpathSupport==0) {
		priority = 2*mult;                     // set OR=100
            } else if (caseSubpathSupport==0 && ctrlSubpathSupport>0) {
		priority = -2*mult;                    // set OR=1/100
            } else {
                double log10OR = Math.log10(orValue);
		priority = (int)(log10OR*mult);
	    }
	    if (priorityOptionLabel==null) {
		priority = Math.abs(priority);       // positive definite
	    } else if (priorityOptionLabel.equals("ctrl")) {
		priority = -priority;                // flip to favor ctrl
	    }
	} else if (priorityOptionKey==4) {
	    // p-value
	    double mlog10p = -Math.log10(pValue);
            priority = (int)(mlog10p*100);
        } else {
	    throwPriorityOptionError();
        }
    }

    /**
     * Exit with an error message if an inappropriate priority option label has been given.
     */
    void throwPriorityOptionError() {
	System.err.println("ERROR: priority "+priorityOptionKey+":"+priorityOptionLabel+" is not supported by FrequentedRegion.updatePriority().");
	System.exit(1);
    }

    /**
     * Return the Fisher's exact test p value for case graph paths vs control graph paths.
     * NOTE1: this does NOT depend on subpaths support!
     * NOTE2: uses graph.fisherExact, which is computed once, since sum(cells)=sum(paths).
     *      | support         | non-support        |
     *      |--------------------------------------|
     * case | casePathSupport | casePathNonsupport |
     * ctrl | ctrlPathSupport | ctrlPathNonsupport |
     * NOTE3: this is calculated only ONCE and stored in pValue.
     */
    public double fisherExactP() {
        if (pValue==Double.NEGATIVE_INFINITY) {
            int casePaths = graph.getPathCount("case");
            int ctrlPaths = graph.getPathCount("ctrl");
            int casePathSupport = getPathSupport("case");
            int ctrlPathSupport = getPathSupport("ctrl");
            int casePathNonsupport = casePaths - casePathSupport;
            int ctrlPathNonsupport = ctrlPaths - ctrlPathSupport;
            pValue = graph.fisherExact.getTwoTailedP(casePathSupport, casePathNonsupport, ctrlPathSupport, ctrlPathNonsupport);
        }
        return pValue;
    }

    /**
     * Return the odds ratio for cases vs controls in terms of supporting subpaths vs. graph paths.
     * 0 = zero case subpath support, POSITIVE_INFINITY = zero control subpath support.
     * NOTE: this is calculated ONCE and stored in orValue.
     */
    public double oddsRatio() {
        if (orValue==Double.NEGATIVE_INFINITY) {
            int casePaths = graph.getPathCount("case");
            int ctrlPaths = graph.getPathCount("ctrl");
            if (ctrlSubpathSupport>0) {
                orValue = (double)caseSubpathSupport * (double)ctrlPaths / ( (double)ctrlSubpathSupport * (double)casePaths );
            } else {
                orValue = Double.POSITIVE_INFINITY;
            }
        }
        return orValue;
    }

    /**
     * Return a string summary of this frequented region.
     * 0      1     2        3     4     5   6  7
     * nodes  size  support  case  ctrl  OR  p  pri
     */
    @Override
    public String toString() {
        return nodes.toString()+"\t"+size+"\t"+support+"\t"+caseSubpathSupport+"\t"+ctrlSubpathSupport+
	    "\t"+orf.format(orValue)+"\t"+pf.format(pValue)+"\t"+priority;
    }

    /**
     * Return a string with the subpaths.
     */
    public String subpathsString() {
        StringBuilder sb = new StringBuilder();
        for (Path sp : subpaths) {
            sb.append(sp.toString());
	    sb.append("\n");
        }
	return sb.toString();
    }

    /**
     * Return true if this FR contains a subpath which belongs to the given Path.
     */
    public boolean containsSubpathOf(Path path) {
        if (subpaths!=null) {
            for (Path sp : subpaths) {
                if (sp.equals(path)) return true;
            }
        }
        return false;
    }

    /**
     * Return a count of subpaths of FR that belong to the given Path.
     */
    public int countSubpathsOf(Path path) {
	if (subpaths==null || subpaths.size()==0) {
	    return 0;
	} else {
	    int count = 0;
            for (Path sp : subpaths) {
                if (sp.getName().equals(path.getName())) count++;
            }
	    return count;
        }
    }

    /**
     * Return true if the nodes in this FR are a subset of the nodes in the given FR (but they are not equal!).
     */
    public boolean isSubsetOf(FrequentedRegion fr) {
        if (this.equals(fr)) {
            return false;
        } else {
            return this.nodes.equals(fr.nodes);
        }
    }
    
    /**
     * Return the count of subpaths that have the given label.
     */
    public int labelCount(String label) {
        int count = 0;
        for (Path sp : subpaths) {
            if (sp.getLabel().equals(label)) count++;
        }
        return count;
    }

    /**
     * Return the NodeSet associated with this FR.
     */
    public NodeSet getNodes() {
        return nodes;
    }

    /**
     * Return true if this FR contains the given node.
     */
    public boolean containsNode(Node n) {
        return nodes.contains(n);
    }

    /**
     * Return true if one of the nodes in this FR is a no call.
     */
    public boolean containsNoCallNode() {
	for (Node n : nodes) {
	    if (n.genotype==null) n = graph.getNode(n.id);
	    if (n.isNoCall()) {
		return true;
	    }
	}
	return false;
    }

    /**
     * Return the count of case subpaths traversing the given node.
     */
    public int getCaseCount(Node n) {
        int count = 0;
        for (Path p : subpaths) {
            if (p.isCase()) {
                for (Node node : p.getNodes()) {
                    if (node.equals(n)) {
                        count++;
                        break;
                    }
                }
            }
        }
        return count;
    }

    /**
     * Return the count of control subpaths traversing the given node.
     */
    public int getControlCount(Node n) {
        int count = 0;
        for (Path p : subpaths) {
            if (p.isControl()) {
                for (Node node : p.getNodes()) {
                    if (node.equals(n)) {
                        count++;
                        break;
                    }
                }
            }
        }
        return count;
    }

    /**
     * Command-line utility gives results for an input cluster of nodes and alpha, kappa and graph.
     */
    public static void main(String[] args) throws FileNotFoundException, IOException {
        Options options = new Options();
 
        Option alphaOption = new Option("a", "alpha", true, "starting value of alpha for a scan (can equal alphaend)");
        alphaOption.setRequired(true);
        options.addOption(alphaOption);
        //
        Option kappaOption = new Option("k", "kappa", true, "starting value of kappa for a scan (can equal kappaend)");
        kappaOption.setRequired(true);
        options.addOption(kappaOption);
        //
        Option priorityOptionOption = new Option("pri", "priorityoption", true, "option for priority weighting of FRs: "+FrequentedRegion.PRIORITY_OPTIONS);
        priorityOptionOption.setRequired(true);
        options.addOption(priorityOptionOption);
        //
        Option graphOption = new Option("graph", "graph", true, "graph name");
        graphOption.setRequired(true);
        options.addOption(graphOption);
        //
        Option gfaOption = new Option("gfa", "gfa", false, "load from [graph].paths.gfa");
        gfaOption.setRequired(false);
        options.addOption(gfaOption);
        //
        Option txtOption = new Option("txt", "txt", false, "load from [graph].nodes.txt and [graph].paths.txt");
        txtOption.setRequired(false);
        options.addOption(txtOption);
        //
        Option nodesOption = new Option("n", "nodes", true, "set of nodes to calculate FR e.g. [1,2,3,4,5]");
        nodesOption.setRequired(true);
        options.addOption(nodesOption);
	//
	Option pathsOption = new Option("p", "pathsfile", true, "paths.txt file");
	pathsOption.setRequired(true);
	options.addOption(pathsOption);
        //
        Option excludedPathNodesOption = new Option("ep", "excludedpathnodes", true, "exclude paths that include any of the given nodes []");
        excludedPathNodesOption.setRequired(false);
        options.addOption(excludedPathNodesOption);
        //
        Option includedPathNodesOption = new Option("ip", "includedpathnodes", true, "include only paths that include at least one of the given nodes []");
        includedPathNodesOption.setRequired(false);
        options.addOption(includedPathNodesOption);
        
        CommandLine cmd;
        HelpFormatter formatter = new HelpFormatter();
        try {
            CommandLineParser parser = new DefaultParser();
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println(e.getMessage());
            formatter.printHelp("FrequentedRegion", options);
            System.exit(1);
            return;
        }

        if (cmd.getOptions().length==0) {
            formatter.printHelp("FrequentedRegion", options);
            System.exit(1);
            return;
        }
        
        // alpha, kappa
        double alpha = Double.parseDouble(cmd.getOptionValue("alpha"));
        int kappa = Integer.parseInt(cmd.getOptionValue("kappa"));
        if (kappa==-1) kappa = Integer.MAX_VALUE; // effectively infinity

        // parse the priorityOption and impose defaults
        String[] parts = cmd.getOptionValue("priorityoption").split(":");
        int priorityOptionKey = Integer.parseInt(parts[0]);
        String priorityOptionLabel = null;
        if (parts.length>1) {
            String priorityOptionParameter = parts[1];
            if (priorityOptionParameter.equals("case") || priorityOptionParameter.equals("ctrl")) {
                priorityOptionLabel = priorityOptionParameter;
            }
        }
        // impose defaults
        if (priorityOptionKey==1 && priorityOptionLabel==null) priorityOptionLabel = "case";

        // import the PangenomicGraph from a pair of TXT files
        PangenomicGraph pg = new PangenomicGraph(cmd.getOptionValue("graph"));
	pg.loadNodesFromTXT(pg.getNodesFile());
	pg.loadPathsFromTXT(new File(cmd.getOptionValue("pathsfile")));
	// remove paths that contain an excluded node, if given
	String excludedPathNodeString = "[]";
	if (cmd.hasOption("excludedpathnodes")) {
	    excludedPathNodeString = cmd.getOptionValue("excludedpathnodes");
	}
	NodeSet excludedPathNodes = pg.getNodeSet(excludedPathNodeString);
	if (excludedPathNodes.size()>0) {
	    Set<Path> pathsToRemove = new HashSet<>();
	    for (Path path : pg.paths) {
		for (Node node : excludedPathNodes) {
		    if (path.traverses(node)) {
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
	    TreeSet<Path> pathsToKeep = new TreeSet<>();
	    for (Path path : pg.paths) {
		for (Node node : includedPathNodes) {
		    if (path.traverses(node)) {
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
        pg.tallyLabelCounts();
	System.out.println("# Graph has "+pg.vertexSet().size()+" nodes and "+pg.edgeSet().size()+" edges with "+pg.paths.size()+" paths.");
        System.out.println("# Graph has "+pg.labelCounts.get("case")+" case paths and "+pg.labelCounts.get("ctrl")+" ctrl paths.");

        // create the FrequentedRegion with this PangenomicGraph
        NodeSet nodes = pg.getNodeSet(cmd.getOptionValue("nodes"));
        FrequentedRegion fr = new FrequentedRegion(pg, nodes, alpha, kappa, priorityOptionKey, priorityOptionLabel);

        // print it out
        System.out.println(fr.columnHeading());
        System.out.println(fr.toString());
    }

    /**
     * Algorithm 1 from Cleary, et al. generates the supporting path segments of this path for the given NodeSet and alpha and kappa parameters.
     *
     * @param nodes the NodeSet, or cluster C as it's called in Algorithm 1
     * @param alpha the penetrance parameter = minimum fraction of nodes in C that are in subpath
     * @param kappa the insertion parameter = maximum inserted number of nodes
     * @return the set of supporting path segments
     */
    public List<Path> computeSupport(Path p) {
        // s = the supporting subpaths
        List<Path> s = new LinkedList<>();
        // m = the list of the path's nodes that are in C=nodes
        List<Node> m = new LinkedList<>();
        for (Node n : p.getNodes()) {
            if (nodes.contains(n)) m.add(n);
        }
        // find maximal subpaths
        for (int i=0; i<m.size(); i++) {
            Node nl = m.get(i);
            Node nr = null;
            int num = 0;
            for (int j=i; j<m.size(); j++) {
		if (kappa<Integer.MAX_VALUE) {
		    // kappa test
		    Path subpath = p.subpath(nl, m.get(j));
		    int maxInsertion = 0; // max insertion
		    int insertion = 0; // continguous insertion
		    for (Node n : subpath.getNodes()) {
			if (nodes.contains(n)) {
			    // reset and save previous insertion if large
			    if (insertion>maxInsertion) maxInsertion = insertion;
			    insertion = 0;
			} else {
			    insertion += 1;
			}
		    }
		    if (maxInsertion>kappa) break;
		}
                // we're good, set nr from this cycle
                nr = m.get(j);
                num = j - i + 1; // number of this path's nodes in nodes collection
            }
            // is this a subpath of an already counted subpath? (maximality test)
            Path subpath = p.subpath(nl,nr);
            boolean ignore = false;
            for (Path checkpath : s) {
                if (checkpath.contains(subpath)) {
                    ignore = true;
                    break;
                }
            }
            // sanity check
            if (subpath.getNodes().size()==0) {
                System.err.println("ERROR: subpath.getNodes().size()=0; path="+this.toString()+" nl="+nl+" nr="+nr);
                ignore = true;
            }
            // alpha test on maximal subpath; use num>0 to allow alpha=0
            if (!ignore && num>0 && num>=alpha*size) s.add(subpath);
        }
        return s;
    }

    /**
     * Return true iff all nodes are on one chromosome.
     */
    public boolean onOneChromosome() {
	Node n1 = nodes.first();
	for (Node n : nodes) {
	    if (!n.contig.equals(n1.contig)) {
		return false;
	    }
	}
	return true;
    }
}
