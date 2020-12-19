package org.ncgr.pangenomics.genotype;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;

import java.util.Map;
import java.util.HashMap;

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
        Option mafOption = new Option("maf", "maf", true, "minimum MAF for inclusion");
        mafOption.setRequired(false);
        options.addOption(mafOption);

        // actions
        Option prsOption = new Option("prs", "prs", false, "compute polygenic risk scores");
        prsOption.setRequired(false);
        options.addOption(prsOption);

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

        // TXT load pulls sample labels from paths file
        graph.nodesFile = new File(graph.name+".nodes.txt");
        graph.pathsFile = new File(graph.name+".paths.txt");
        graph.loadTXT();
        graph.tallyLabelCounts();
        System.err.println("Graph has "+graph.vertexSet().size()+" nodes and "+graph.paths.size()+" paths: "+graph.labelCounts.get("case")+"/"+graph.labelCounts.get("ctrl")+" cases/controls");

        // options
        boolean prs = cmd.hasOption("prs");

        // actions
        if (prs) {
            double minMAF = 1e-2;
            if (cmd.hasOption("maf")) {
                minMAF = Double.parseDouble(cmd.getOptionValue("maf"));
            }
            computePRS(graph, minMAF);
        }
    }

    /**
     * Compute the polygenic risk scores from a graph, tossing nodes with af<minMAF.
     */
    public static void computePRS(PangenomicGraph graph, double minMAF) {
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
            // public double af;
            // public boolean isCalled;
            if (!(n.af<minMAF)) {
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

}
