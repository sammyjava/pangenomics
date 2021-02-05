package org.ncgr.pangenomics.genotype;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;

import java.util.Arrays;
import java.util.List;
import java.util.LinkedList;
import java.util.Map;
import java.util.HashMap;
import java.util.Optional;
import java.util.TreeMap;
import java.util.Set;
import java.util.HashSet;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentSkipListSet;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import org.jgrapht.graph.DirectedAcyclicGraph;

import org.mskcc.cbio.portal.stats.FisherExact;

public class PangenomicGraph extends DirectedAcyclicGraph<Node,Edge> {

    // output verbosity
    public boolean verbose = false;

    // a name for this graph
    public String name;

    // the TXT files holding the nodes and paths
    public File nodesFile;
    public File pathsFile;
    
    // Set of Paths that traverse this graph
    public Set<Path> paths = new HashSet<>();

    // Map of Path names to their Path (convenience)
    public Map<String,Path> pathNameMap = new HashMap<>();

    // Map of Node ids to their Nodes (convenience)
    public Map<Long,Node> nodeIdMap = new HashMap<>();

    // maps a Node to the Paths that traverse it
    public Map<Node,List<Path>> nodePaths = new HashMap<>();

    // associate samples with labels like "case" and "ctrl"
    public Map<String,String> sampleLabels = new HashMap<>();
    
    // maps a path label to a count of paths that have that label
    public Map<String,Integer> labelCounts; // keyed by label, construct in method

    // computed once to save time
    public FisherExact fisherExact;

    // drop no-call nodes and paths traversing them
    public boolean dropNoCallPaths = false;

    // make balanced graph, equal cases and controls
    public boolean equalizeCasesControls = false;

    // set the maximum number of cases
    public int maxCases;

    /**
     * Basic constructor.
     */
    public PangenomicGraph() {
        super(Edge.class);
    }

    /**
     * Build this graph from the provided Nodes and and maps. sampleLabels must already be populated.
     */
    public void buildGraph(List<Node> nodes, Map<String,NodeSet> sampleNodeSets, Map<Node,Set<String>> nodeSamples) {
	// validation
	if (sampleLabels.size()==0) {
	    System.err.println("ERROR in PangenomicGraph.buildGraph: sampleLabels has not been populated.");
	    System.exit(1);
	}
	if (nodes.size()==0) {
	    System.err.println("ERROR in PangenomicGraph.buildGraph: nodes list is empty.");
	    System.exit(1);
	}
	if (sampleNodeSets.size()==0) {
	    System.err.println("ERROR in PangenomicGraph.buildGraph: sampleNodeSets is empty.");
	    System.exit(1);
	}
	if (nodeSamples.size()==0) {
	    System.err.println("ERROR in PangenomicGraph.buildGraph: nodeSamples is empty.");
	    System.exit(1);
	}
        // remove sampleLabels entries that aren't in sampleNodeSets
	int removeCount = 0;
        for (String sampleName : sampleNodeSets.keySet()) {
            if (!sampleLabels.containsKey(sampleName)) {
		sampleLabels.remove(sampleName);
		removeCount++;
	    }
        }
	if (verbose && removeCount>0) {
	    System.err.println("Removed "+removeCount+" samples that were not in sampleNodeSets.");
	    System.err.println(sampleLabels.size()+" samples remain.");
	}
        // optionally drop no-call nodes and all paths that traverse them
        if (dropNoCallPaths) {
            Set<Node> nodesToDrop = new HashSet<>();
            Set<String> samplesToDrop = new HashSet<>();
            for (Node n : nodes) {
                if (n.isNoCall()) {
                    nodesToDrop.add(n);
		    Set<String> sampleNames = nodeSamples.get(n);
                    samplesToDrop.addAll(sampleNames);
                }
            }
            for (Node n : nodesToDrop) {
                nodes.remove(n);
                nodeSamples.remove(n);
            }
	    for (String sampleName : samplesToDrop) {
		sampleNodeSets.remove(sampleName);
		sampleLabels.remove(sampleName);
	    }
            System.err.println("Removed "+nodesToDrop.size()+" no-call nodes and "+samplesToDrop.size()+" samples containing no-call nodes.");
        }
	// optionally limit maximum cases
	if (maxCases>0) {
	    int nCases = 0;
	    int nControls = 0;
            for (String label : sampleLabels.values()) {
                if (label.equals("case")) {
                    nCases++;
                } else if (label.equals("ctrl")) {
                    nControls++;
                }
            }
            Set<String> samplesToDrop = new HashSet<>();
            int numberToRemove = nCases - maxCases;
	    int count = 0;
	    // randomly select samples to drop
	    while (count<numberToRemove) {
		Optional<String> optional = sampleLabels.keySet().stream().skip((int)(sampleLabels.size()*Math.random())).findFirst();
		if (optional.isPresent()) {
		    String sampleName = optional.get();
		    if (count<numberToRemove && sampleLabels.get(sampleName).equals("case") && !samplesToDrop.contains(sampleName)) {
			samplesToDrop.add(sampleName);
			count++;
		    }
                }
	    }
	    // remove dropped samples
	    for (String sampleName : samplesToDrop) {
		sampleLabels.remove(sampleName);
		sampleNodeSets.remove(sampleName);
		for (Node n : nodeSamples.keySet()) {
		    Set<String> sampleNames = nodeSamples.get(n);
		    sampleNames.remove(sampleName);
		}
            }
            // drop nodes that no longer have any samples
            Set<Node> nodesToDrop = new HashSet<>();
            for (Node n : nodeSamples.keySet()) {
                Set<String> sampleNames = nodeSamples.get(n);
                if (sampleNames.size()==0) nodesToDrop.add(n);
            }
	    for (Node n : nodesToDrop) {
		nodes.remove(n);
		nodeSamples.remove(n);
	    }
            System.err.println("Removed "+samplesToDrop.size()+" samples and "+nodesToDrop.size()+" orphaned nodes to set number of cases = "+maxCases+".");
	}
	// optionally enforce cases=controls
        if (equalizeCasesControls) {
            int nCases = 0;
            int nControls = 0;
            for (String label : sampleLabels.values()) {
                if (label.equals("case")) {
                    nCases++;
                } else if (label.equals("ctrl")) {
                    nControls++;
                }
            }
            Set<String> samplesToDrop = new HashSet<>();
            int numberToRemove = Math.abs(nCases - nControls);
            int count = 0;
            // randomly select samples to drop
            while (count<numberToRemove) {
                Optional<String> optional = sampleLabels.keySet().stream().skip((int)(sampleLabels.size()*Math.random())).findFirst();
                if (optional.isPresent()) {
                    String sampleName = optional.get();
                    if (nCases>nControls && sampleLabels.get(sampleName).equals("case") && !samplesToDrop.contains(sampleName)) {
                        samplesToDrop.add(sampleName);
                        count++;
                    } else if (nControls>nCases && sampleLabels.get(sampleName).equals("ctrl") && !samplesToDrop.contains(sampleName)) {
                        samplesToDrop.add(sampleName);
                        count++;
                    }
                }
	    }
	    // remove dropped samples
	    for (String sampleName : samplesToDrop) {
		sampleLabels.remove(sampleName);
		sampleNodeSets.remove(sampleName);
		for (Node n : nodeSamples.keySet()) {
		    Set<String> sampleNames = nodeSamples.get(n);
		    sampleNames.remove(sampleName);
		}
            }
            // drop nodes that no longer have any samples
            Set<Node> nodesToDrop = new HashSet<>();
            for (Node n : nodeSamples.keySet()) {
                Set<String> sampleNames = nodeSamples.get(n);
                if (sampleNames.size()==0) nodesToDrop.add(n);
            }
	    for (Node n : nodesToDrop) {
		nodes.remove(n);
		nodeSamples.remove(n);
	    }
            System.err.println("Removed "+samplesToDrop.size()+" samples and "+nodesToDrop.size()+" orphaned nodes to equalize cases and controls.");
        }
        // add the nodes as graph vertices
        if (verbose) System.err.println("Adding nodes to graph vertices...");
        for (Node n : nodes) {
            addVertex(n);
	    nodeIdMap.put(n.id, n);
        }
        // build the paths and path-labeled graph edges from the sampleLabels and sampleNodeSets map
	// NOTE: sampleLabels may contain some samples not in sampleNodeSets and vice-versa
        if (verbose) System.err.println("Creating paths and adding edges to graph...");
        for (String sampleName : sampleLabels.keySet()) {
	    if (sampleNodeSets.containsKey(sampleName)) {
		NodeSet sampleNodes = sampleNodeSets.get(sampleName);
		Path path = new Path(this, new LinkedList<Node>(sampleNodes), sampleName, sampleLabels.get(sampleName));
		paths.add(path);
		pathNameMap.put(sampleName, path);
		// add edges
		Node lastNode = null;
		for (Node node : sampleNodes) {
		    if (lastNode!=null) {
			if (!containsEdge(lastNode, node)) {
			    try {
				addEdge(lastNode, node);
			    } catch (Exception e) {
				System.err.println("ERROR adding edge from "+lastNode+" to "+node);
				System.err.println(e);
				System.exit(1);
			    }
			}
		    }
		    lastNode = node;
		}
	    }
        }
	// initialize FisherExact for later use
        fisherExact = new FisherExact(paths.size());
    }

    /**
     * Build this graph from the provided Nodes and Paths.
     */
    public void buildGraph(List<Node> nodes, Set<Path> paths) {
        this.paths = paths;
        pathNameMap = new HashMap<>();
        // add the nodes as graph vertices
        if (verbose) System.err.println("Adding nodes to graph vertices...");
        for (Node n : nodes) {
            addVertex(n);
        }
        // build the path-labeled graph edges
        if (verbose) System.err.println("Adding edges to graph from paths...");
        for (Path path : paths) {
            pathNameMap.put(path.name, path);
            // add edges
            Node lastNode = null;
            for (Object n : path.getNodes()) {
                Node node = (Node) n;
                if (lastNode!=null) {
                    if (!containsEdge(lastNode, node)) {
                        try {
                            addEdge(lastNode, node);
                        } catch (Exception e) {
                            System.err.println("ERROR adding edge from "+lastNode+" to "+node);
                            System.err.println(e);
                            System.exit(1);
                        }
                    }
                }
                lastNode = node;
            }
        }
        fisherExact = new FisherExact(paths.size());
    }

    /**
     * Read sample labels from the instance tab-delimited file. Comment lines start with #. Header starts with "sample".
     * Disregard samples not labeled "case" or "ctrl".
     *
     * sample	pheno_241.2 pheno_555 pheno_695.42 pheno_714.1 pheno_250.1
     * 476296	case        ctrl      ctrl         ctrl        ctrl
     * 261580	ctrl        ctrl      ctrl         ctrl        ctrl
     * 273564	ctrl        ctrl      ctrl         ctrl        unkn
     */
    public void readSampleLabels(File labelsFile) throws FileNotFoundException, IOException {
	if (verbose) System.err.println("Reading samples and labels from "+labelsFile.getName());
        String line = null;
        BufferedReader reader = new BufferedReader(new FileReader(labelsFile));
        while ((line=reader.readLine())!=null) {
            if (line.startsWith("#")) continue;
            String[] fields = line.split("\t");
	    if (fields[0].equals("sample")) continue;
            if (fields.length>=2) {
		// only keep case or control samples
		if (fields[1].equals("case") || fields[1].equals("ctrl")) {
		    sampleLabels.put(fields[0], fields[1]); // name, label
		}
            } else {
		System.err.println("Found a line in "+labelsFile+" with less than two fields:");
		System.err.println(line);
		System.exit(1);
	    }
        }
        reader.close();
    }

    /**
     * Tally the label counts for the paths in this graph (can be a subset of all samples/paths).
     */
    public void tallyLabelCounts() {
	labelCounts = new HashMap<>();
        for (Path path : paths) {
            if (labelCounts.containsKey(path.label)) {
                int count = labelCounts.get(path.label);
                labelCounts.put(path.label, count+1);
            } else {
                labelCounts.put(path.label, 1);
            }
        }
    }

    /**
     * Build the node paths: the set of paths that run through each node.
     */
    public void buildNodePaths() {
	if (verbose) System.err.println("Building node paths...");
        // initialize empty paths for each node
        for (Node n : getNodes()) {
            nodePaths.put(n, new LinkedList<Path>());
        }
        // now load the paths
	for (Path path : paths) {
	    List<Node> nodes = path.getNodes();
	    for (Node n : nodes) {
		nodePaths.get(n).add(path);
	    }
	}
	// validation
	for (Node n : nodePaths.keySet()) {
	    List<Path> pathList = nodePaths.get(n);
	    for (Path p : pathList) {
		if (p==null) {
		    System.err.println("ERROR: p==null for node "+n.toString());
		    System.exit(1);
		}
	    }
	}
    }

    /**
     * Return the node with the given id, else null.
     */
    public Node getNode(long id) {
        return nodeIdMap.get(id);
    }

    /**
     * Return true if this graph contains the path with the given name.
     */
    public boolean hasPath(String name) {
        return pathNameMap.containsKey(name);
    }
    
    /**
     * Return the path with the given name.
     */
    public Path getPath(String name) {
        return pathNameMap.get(name);
    }

    /**
     * Return the number of this graph's paths.
     */
    public int getPathCount() {
        return paths.size();
    }

    /**
     * Get the count of paths with the given label.
     */
    public int getPathCount(String label) {
        if (labelCounts.containsKey(label)) {
            return labelCounts.get(label);
        } else {
            return 0;
        }
    }

    /**
     * Return the number of paths that traverse the given node.
     */
    public int getPathCount(Node n) {
        return getPaths(n).size();
    }

    /**
     * Get the total count of paths that follow the given Edge.
     */
    public int getPathCount(Edge e) {
        int count = 0;
        for (Path p : paths) {
            List<Edge> edges = p.getEdges();
            if (edges.contains(e)) {
                count++;
            }
        }
        return count;
    }

    /**
     * Get the count of paths with the given label that follow the given Edge.
     */
    public int getPathCount(Edge e, String label) {
        int count = 0;
        for (Path p : paths) {
            if (p.label.equals(label)) {
                List<Edge> edges = p.getEdges();
                if (edges.contains(e)) {
                    count++;
                }
            }
        }
        return count;
    }

    /**
     * Return the paths that traverse the given node.
     */
    public List<Path> getPaths(Node n) {
        return nodePaths.get(n);
    }

    /**
     * Just a synonym for vertexSet(), but puts the Nodes in a List.
     */
    public List<Node> getNodes() {
        return new LinkedList<Node>(vertexSet());
    }

    /**
     * Construct a NodeSet from a string representation, e.g. [1350,1352,1353,1465,1467,1468,1469].
     */
    public NodeSet getNodeSet(String str) {
        NodeSet nodes = new NodeSet();
        List<Node> allNodes = getNodes();
        Map<Long,Node> allNodeMap = new HashMap<>();
        for (Node n : allNodes) allNodeMap.put(n.id, n);
        List<String> nodeStrings = Arrays.asList(str.replace("[","").replace("]","").split(","));
        for (String s : nodeStrings) {
            if (s.length()>0) {
                long id = Long.parseLong(s);
                if (allNodeMap.containsKey(id)) {
                    nodes.add(allNodeMap.get(id));
                } else {
                    // bail, we're asked for a node that is not in the graph
		    System.err.println("ERROR: Graph does not contain node "+id);
		    System.exit(1);
		}
            }
        }
        return nodes;
    }

    /**
     * Get the label counts map for paths that traverse the given node.
     */
    public Map<String,Integer> getLabelCounts(Node n) {
        Map<String,Integer> map = new HashMap<>();
        for (Path p : nodePaths.get(n)) {
            if (map.containsKey(p.label)) {
                // increment count
                int count = map.get(p.label) + 1;
                map.put(p.label, count);
            } else {
                // initialize count
                map.put(p.label, 1);
            }
        }
        return map;
    }

    /**
     * Get the label counts for paths that follow the given Edge.
     * TODO: restrict path loop to those that are on Nodes connected by Edge. May not be worth it.
     */
    public Map<String,Integer> getLabelCounts(Edge e) {
        Map<String,Integer> map = new HashMap<>();
        for (Path p : paths) {
            List<Edge> edges = p.getEdges();
            if (edges.contains(e)) {
                if (map.containsKey(p.label)) {
                    // increment count
                    int count = map.get(p.label) + 1;
                    map.put(p.label, count);
                } else {
                    // initialize count
                    map.put(p.label, 1);
                }
            }
        }
        return map;
    }

    /**
     * Return the path names as a string array.
     */
    public String[] getPathNames() {
        String[] names = new String[paths.size()];
	int i = 0;
	for (Path p : paths) {
	    names[i++] = p.name;
	}
        return names;
    }

    /**
     * Return the odds ratio for case/control paths that traverse the given node vs. case/control paths that don't.
     * 0                 = all control paths traverse node
     * POSITIVE_INFINITY = all case paths traverse node
     */
    public double oddsRatio(Node n) {
        Map<String,Integer> countsMap = getLabelCounts(n);
        int nodeCasePaths = 0; if (countsMap.containsKey("case")) nodeCasePaths = countsMap.get("case");
        int nodeCtrlPaths = 0; if (countsMap.containsKey("ctrl")) nodeCtrlPaths = countsMap.get("ctrl");
        int otherCasePaths = labelCounts.get("case") - nodeCasePaths;
        int otherCtrlPaths = labelCounts.get("ctrl") - nodeCtrlPaths;
        return (double)(nodeCasePaths*otherCtrlPaths) / (double)(nodeCtrlPaths*otherCasePaths);
    }

    /**
     * Return the Fisher's exact test two-tailed p value for case/control paths that traverse the given node
     * vs. case/control paths that don't.
     *
     *      | traverse node | don't trv node |
     *      |--------------------------------|
     * case | nodeCasePaths | otherCasePaths |
     * ctrl | nodeCtrlPaths | otherCtrlPaths |
     */
    public double fisherExactP(Node n) {
        Map<String,Integer> countsMap = getLabelCounts(n);
        int nodeCasePaths = 0; if (countsMap.containsKey("case")) nodeCasePaths = countsMap.get("case");
        int nodeCtrlPaths = 0; if (countsMap.containsKey("ctrl")) nodeCtrlPaths = countsMap.get("ctrl");
        int otherCasePaths = labelCounts.get("case") - nodeCasePaths;
        int otherCtrlPaths = labelCounts.get("ctrl") - nodeCtrlPaths;
        return fisherExact.getTwoTailedP(nodeCasePaths, otherCasePaths, nodeCtrlPaths, otherCtrlPaths);
    }

    /**
     * Print out the nodes
     * public long id;
     * public String rs;
     * public String contig;
     * public int start;
     * public int end;
     * public String genotype;
     */
    public void printNodes(PrintStream out) {
        if (out==System.out) printHeading("NODES");
        for (Node n : getNodes()) {
            out.println(n.id+"\t"+n.rs+"\t"+n.contig+"\t"+n.start+"\t"+n.end+"\t"+n.genotype+"\t"+n.gf);
        }
    }

    /**
     * Print the paths, labeled by Path.name.
     */
    public void printPaths(PrintStream out) {
        if (out==System.out) printHeading("PATHS");
        for (Path path : paths) {
	    out.println(path.toString());
        }
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
     * Print out the node paths along with counts.
     */
    public void printNodePaths(PrintStream out) {
        if (out==System.out) printHeading("NODE PATHS");
        for (Node node : nodePaths.keySet()) {
            List<Path> pathList = nodePaths.get(node);
            StringBuilder builder = new StringBuilder();
            builder.append(node.id);
            for (Path path : pathList) {
                builder.append("\t"+path.name);
            }
            out.println(builder.toString());
        }
    }

    /**
     * Print node participation by path, appropriate for PCA analysis.
     * NOTE: no-call nodes are dropped.
     */
    public void printPcaData(PrintStream out) throws FileNotFoundException, IOException {
        StringBuilder headerBuilder = new StringBuilder();
        // header is paths, which are alpha sorted for mergability
	TreeSet<Path> sortedPaths = new TreeSet<>(paths);
        boolean first = true;
        for (Path path : sortedPaths) {
            if (first) {
                headerBuilder.append(path.name);
                first = false;
            } else {
                headerBuilder.append("\t"+path.name);
            }
        }
        out.println(headerBuilder.toString());
	/////////////////////////////////////////////////////////////////////////////
        // rows are nodes and counts of path support of each node
	ConcurrentSkipListSet<Node> concurrentNodes = new ConcurrentSkipListSet<>(getNodes());
	ConcurrentSkipListSet<String> lines = new ConcurrentSkipListSet<>();
	concurrentNodes.parallelStream().forEach(node -> {
		if (!node.isNoCall()) {
		    StringBuilder lineBuilder = new StringBuilder();
		    lineBuilder.append("N"+node.id);
		    List<Path> nPaths = nodePaths.get(node.id);
		    // spin through every path, printing 0/1 if path doesn't/does traverse this node
		    for (Path path : sortedPaths) {
			List<Node> nodeList = path.getNodes();
			if (nodeList.contains(node)) {
			    lineBuilder.append("\t1");
			} else {
			    lineBuilder.append("\t0");
			}		    
		    }
		    lines.add(lineBuilder.toString());
		}
	    });
	/////////////////////////////////////////////////////////////////////////////
	for (String line : lines) {
	    out.println(line);
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
    public void printArffData(PrintStream out) {
	out.println("@RELATION "+name);      // graph name
        out.println("");
	out.println("@ATTRIBUTE ID STRING"); // path ID
	// attributes: each node is labeled Nn where n is the node ID
	for (Node node : getNodes()) {
	    if (!node.isNoCall()) {
		out.println("@ATTRIBUTE N"+node.id+" NUMERIC");
	    }
        }
        // add the path label class attribute
        out.println("@ATTRIBUTE class {ctrl,case}");
        out.println("");
        // data
        out.println("@DATA");
	/////////////////////////////////////////////////////////////////////////////
	ConcurrentSkipListSet<Path> concurrentPaths = new ConcurrentSkipListSet<>();
	ConcurrentSkipListSet<String> arffData = new ConcurrentSkipListSet<>();
	concurrentPaths.addAll(paths);
	concurrentPaths.parallelStream().forEach(path -> {
		String arff = path.name;
		for (Node node : getNodes()) {
		    if (!node.isNoCall()) {
			if (path.contains(node)) {
			    arff += ",1";
			} else {
			    arff += ",0";
			}
		    }
		}
		arff += ","+path.label;
		arffData.add(arff);
	    });
	/////////////////////////////////////////////////////////////////////////////
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
    public void printSvmData(PrintStream out) {
	/////////////////////////////////////////////////////////////////////////////
	ConcurrentSkipListSet<Path> concurrentPaths = new ConcurrentSkipListSet<>();
	ConcurrentSkipListSet<String> svmData = new ConcurrentSkipListSet<>();
	concurrentPaths.addAll(paths);
	concurrentPaths.parallelStream().forEach(path -> {
		String svm = path.name+"\t"+path.label;
		int n = 0;
		for (Node node : getNodes()) {
		    if (!node.isNoCall()) {
			n++;
			if (path.contains(node)) {
			    svm += "\t"+n+":+1";
			} else {
			    svm += "\t"+n+":-1";
			}
		    }
		}
		svmData.add(svm);
	    });
	/////////////////////////////////////////////////////////////////////////////
	for (String svm : svmData) {
	    out.println(svm);
	}
        out.close();
    }

    /**
     * Load the graph from a VCF file.
     */
    public void loadVCF(File vcfFile) throws IOException {
        VCFImporter vcfImporter = new VCFImporter();
	vcfImporter.verbose = true;
        vcfImporter.read(vcfFile, true); // ignore phasing
        buildGraph(vcfImporter.nodes, vcfImporter.sampleNodeSets, vcfImporter.nodeSamples);
    }

    /**
     * Load the graph from a PLINK list output file for the samples listed in sampleLabels.
     */
    public void loadList(File listFile) throws IOException {
        ListImporter listImporter = new ListImporter();
	listImporter.verbose = true;
	if (verbose) System.err.println("Loading LIST data for "+sampleLabels.size()+" samples.");
        listImporter.read(listFile, sampleLabels.keySet());
        buildGraph(listImporter.nodes, listImporter.sampleNodeSets, listImporter.nodeSamples);
    }

    /**
     * Load the graph from a pair of TXT files.
     * nodesFile and pathsFile must already be set.
     */
    public void loadTXT() throws FileNotFoundException, IOException {
        if (nodesFile==null || pathsFile==null) {
            System.err.println("ERROR: graph.nodesFile and/or graph.pathsFile has not been set!");
            System.exit(1);
        }
        TXTImporter txtImporter = new TXTImporter();
	txtImporter.verbose = true;
        txtImporter.read(nodesFile, pathsFile);
	this.sampleLabels = txtImporter.sampleLabels;
        buildGraph(txtImporter.nodes, txtImporter.sampleNodeSets, txtImporter.nodeSamples);
    }

    /**
     * The main event.
     */
    public static void main(String[] args) throws IOException, FileNotFoundException {
	Options options = new Options();
        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd;

        Option verboseOption = new Option("v", "verbose", false, "verbose output (false)");
        verboseOption.setRequired(false);
        options.addOption(verboseOption);
        // REQUIRED parameters
        Option graphOption = new Option("g", "graph", true, "name of graph");
        graphOption.setRequired(true);
        options.addOption(graphOption);
        // if vcf then labelfile is required
        Option vcfFileOption = new Option("vcf", "vcffile", true, "load graph from <graph>.vcf.gz");
        vcfFileOption.setRequired(false);
        options.addOption(vcfFileOption);
	// if list then labelfile is required
	Option listFileOption = new Option("list", "listfile", true, "load graph from a PLINK list file");
	listFileOption.setRequired(false);
	options.addOption(listFileOption);
        // txt does not require labelfile
        Option txtFileOption = new Option("txt", "txtfile", true, "root name of a pair of TXT files containing nodes and paths (e.g. FOO means FOO.nodes.txt and FOO.paths.txt)");
        txtFileOption.setRequired(false);
        options.addOption(txtFileOption);
	// 
        Option labelFileOption = new Option("l", "labelfile", true, "tab-delimited file containing one sample<tab>label per line");
        labelFileOption.setRequired(false);
        options.addOption(labelFileOption);
        // drop no-call nodes and paths
        Option dropNoCallPathsOption = new Option("dncp", "dropnocallpaths", false, "drop no-call nodes and paths that traverse them [false]");
        dropNoCallPathsOption.setRequired(false);
        options.addOption(dropNoCallPathsOption);
        // reduce cases or controls so they are equal in number
        Option equalizeCasesControlsOption = new Option("ecc", "equalizecasescontrols", false, "reduce cases or controls to make equal in number [false]");
        equalizeCasesControlsOption.setRequired(false);
        options.addOption(equalizeCasesControlsOption);
	// set max cases
	Option maxCasesOption = new Option("maxc", "maxcases", true, "maxiumum number of cases in the graph");
	maxCasesOption.setRequired(false);
	options.addOption(maxCasesOption);
	// optional output
	Option printPcaFileOption = new Option("pca", "printpcafile", false, "print out path pca file [false]");
	printPcaFileOption.setRequired(false);
	options.addOption(printPcaFileOption);
	//
	Option printArffFileOption = new Option("arff", "printarfffile", false, "print out an ARFF file of path node participation [false]");
	printArffFileOption.setRequired(false);
	options.addOption(printArffFileOption);
	//
	Option printSvmFileOption = new Option("svm", "printsvmfile", false, "print out a LIBSVM-compatible file of path node participation [false]");
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
            formatter.printHelp("PangenomicGraph", options);
            System.exit(1);
            return;
        }
        // spit out help and exit if nothing supplied
        if (cmd.getOptions().length==0) {
            formatter.printHelp("PangenomicGraph", options);
            System.exit(1);
            return;
        }
        // validation
        boolean loadVCF = cmd.hasOption("vcffile");
	boolean loadList = cmd.hasOption("listfile");
        boolean loadTXT = cmd.hasOption("txtfile");
        if (!loadVCF && !loadList && !loadTXT) {
            System.err.println("ERROR: You must specify loading from a VCF with --vcffile [filename], a LIST file with --listfile, or a pair of TXT files with --txtfile.");
            System.exit(1);
        }
	if ((loadVCF || loadList) && !cmd.hasOption("labelfile")) {
	    System.err.println("ERROR: You must specify a label file with --labelfile if you are loading a graph from a VCF or LIST file.");
	    System.exit(1);
	}

        // our PangenomicGraph
        PangenomicGraph graph = new PangenomicGraph();
        if (cmd.hasOption("verbose")) graph.verbose = true;
        if (cmd.hasOption("dropnocallpaths")) graph.dropNoCallPaths = true;
        if (cmd.hasOption("equalizecasescontrols")) graph.equalizeCasesControls = true;
	if (cmd.hasOption("maxcases")) graph.maxCases = Integer.parseInt(cmd.getOptionValue("maxcases"));

        // load the graph from a TXT, VCF, or LIST file
        graph.name = cmd.getOptionValue("graph");
        if (loadTXT) {
            // TXT load pulls sample labels from paths file
	    // NOTE: this overrides VCF or LIST loading
            graph.nodesFile = new File(graph.name+".nodes.txt");
            graph.pathsFile = new File(graph.name+".paths.txt");
            graph.loadTXT();
            graph.tallyLabelCounts();
	    System.err.println("Graph has "+graph.vertexSet().size()+" nodes and "+graph.paths.size()+" paths: "+graph.labelCounts.get("case")+"/"+graph.labelCounts.get("ctrl")+" cases/controls");
        } else if (loadVCF) {
            // VCF load needs separate sample labels load
	    // NOTE: this overrides LIST loading
	    graph.readSampleLabels(new File(cmd.getOptionValue("labelfile")));
            graph.loadVCF(new File(cmd.getOptionValue("vcffile")));
	    graph.tallyLabelCounts();
	    System.err.println("Graph has "+graph.vertexSet().size()+" nodes, "+graph.edgeSet().size()+" edges, and "+graph.paths.size()+
                               " paths: "+graph.labelCounts.get("case")+"/"+graph.labelCounts.get("ctrl")+" cases/controls");
        } else if (loadList) {
            // LIST load needs separate sample labels load
	    graph.readSampleLabels(new File(cmd.getOptionValue("labelfile")));
	    graph.loadList(new File(cmd.getOptionValue("listfile")));
	    graph.tallyLabelCounts();
	    System.err.println("Graph has "+graph.vertexSet().size()+" nodes, "+graph.edgeSet().size()+" edges, and "+graph.paths.size()+
                               " paths: "+graph.labelCounts.get("case")+"/"+graph.labelCounts.get("ctrl")+" cases/controls");
	}

        // build the node-paths map
        graph.buildNodePaths();

        // output
	if (!cmd.hasOption("txtfile")) {
	    if (graph.verbose) System.err.println("Writing nodes file...");
	    PrintStream nodesOut = new PrintStream(graph.name+".nodes.txt");
	    graph.printNodes(nodesOut);
	    if (graph.verbose) System.err.println("Writing paths file...");
	    PrintStream pathsOut = new PrintStream(graph.name+".paths.txt");
	    graph.printPaths(pathsOut);
	}

	if (cmd.hasOption("printnodepathsfile")) {
	    if (graph.verbose) System.err.println("Writing node paths file...");
	    PrintStream out = new PrintStream(graph.name+".nodepaths.txt");
	    graph.printNodePaths(out);
	}

	if (cmd.hasOption("printpcafile")) {
	    if (graph.verbose) System.err.println("Writing path PCA file...");
	    PrintStream out = new PrintStream(graph.name+".pathpca.txt");
	    graph.printPcaData(out);
	}

	if (cmd.hasOption("printarfffile")) {
	    if (graph.verbose) System.err.println("Writing path ARFF file...");
	    PrintStream out = new PrintStream(graph.name+".arff");
	    graph.printArffData(out);
	}
	
	if (cmd.hasOption("printsvmfile")) {
	    if (graph.verbose) System.err.println("Writing path SVM file...");
	    PrintStream out = new PrintStream(graph.name+".svm.txt");
	    graph.printSvmData(out);
	}
    }
}
