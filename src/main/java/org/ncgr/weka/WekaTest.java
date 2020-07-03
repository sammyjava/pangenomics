package org.ncgr.weka;

import java.util.Random;
import java.util.TreeMap;
import java.util.concurrent.ConcurrentHashMap;

import java.text.DecimalFormat;

import weka.classifiers.Classifier;
import weka.classifiers.Evaluation;
import weka.classifiers.functions.*;
import weka.classifiers.rules.*;
import weka.classifiers.trees.*;
import weka.core.Instances;
import weka.core.converters.ConverterUtils.DataSource;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
 
/**
 * Run a selection of classifiers on an ARFF file.
 */
public class WekaTest {
    public static int KFOLD = 10;

    public static void main(String[] args) throws Exception {
        Options options = new Options();
        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd = null;
	//
	Option inputFileOption = new Option("i", "inputfile", true, "input file containing feature vectors in ARFF format");
	inputFileOption.setRequired(true);
	options.addOption(inputFileOption);
	//
	Option nRunsOption = new Option("nruns", true, "number of parallel runs [1]");
	nRunsOption.setRequired(false);
	options.addOption(nRunsOption);
 
        if (args.length==0) {
            formatter.printHelp("WekaTest [options]", options);
            System.exit(1);
        }
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println(e.getMessage());
            formatter.printHelp("WekaTest [options]", options);
            System.exit(1);
        }

        String arffFile = cmd.getOptionValue("inputfile");

	int nRuns = 1;
	if (cmd.hasOption("nruns")) nRuns = Integer.parseInt(cmd.getOptionValue("nruns"));

	// Read all the instances in the file (ARFF, CSV, XRFF, ...)
	DataSource source = new DataSource(arffFile);
 	Instances data = source.getDataSet();
        // remove the ID attribute
        data.deleteAttributeAt(0);
        // set the class attribute index
        data.setClassIndex(data.numAttributes() - 1);

	// load Evaluations into Map for parallel ops
	ConcurrentHashMap<String,Evaluation> evaluations = new ConcurrentHashMap<>();
	// run each model nRuns times
	for (int i=0; i<nRuns; i++) {
	    // trees
	    evaluations.put("ds_"+i, new Evaluation(data));
	    evaluations.put("rf_"+i, new Evaluation(data));
	    evaluations.put("rep_"+i, new Evaluation(data));
	    // evaluations.put("ht_"+i, new Evaluation(data));
	    // evaluations.put("j48_"+i, new Evaluation(data));
	    // evaluations.put("lmt_"+i, new Evaluation(data));
	    // evaluations.put("rt_"+i, new Evaluation(data));
	    // rules
	    // evaluations.put("oneR_"+i, new Evaluation(data));
	    // evaluations.put("jRip_"+i, new Evaluation(data));
	    // evaluations.put("part_"+i, new Evaluation(data));
	}
	System.err.println("Initialized "+evaluations.size()+" evaluations.");
    
	///////////////////////////////////////////////////////////////////////////////////////
	// cross-validate each model nRun times in a parallel streams
	evaluations.entrySet().parallelStream().forEach(entry -> {
		String key = entry.getKey();
		String[] parts = key.split("_");
		int seed = Integer.parseInt(parts[1]);
		Evaluation evaluation = entry.getValue();
		Classifier model = null;
		if (key.startsWith("ds")) {
		    model = new DecisionStump();
		} else if (key.startsWith("oneR_")) {
		    model = new OneR();
		} else if (key.startsWith("ht_")) {
		    model = new HoeffdingTree();
		} else if (key.startsWith("lmt_")) {
		    model = new LMT();
		} else if (key.startsWith("part_")) {
		    model = new PART();
		} else if (key.startsWith("j48_")) {
		    model = new J48();
		} else if (key.startsWith("rf_")) {
		    RandomForest rf = new RandomForest();
		    rf.setSeed(seed);
		    model = rf;
		} else if (key.startsWith("rt_")) {
		    RandomTree rt = new RandomTree();
		    rt.setSeed(seed);
		    model = rt;
		} else if (key.startsWith("rep_")) {
		    REPTree rep = new REPTree();
		    rep.setSeed(seed);
		    model = rep;
		} else if (key.startsWith("jRip_")) {
		    JRip jRip = new JRip();
		    jRip.setSeed((long)seed);
		    model = jRip;
		} else {
		    System.err.println("ERROR: cannot match model to key "+key);
		    System.exit(1);
		}
		try {
		    evaluation.crossValidateModel(model, data, KFOLD, new Random(seed));
		} catch (Exception e) {
		    System.err.println(e);
		    System.exit(1);
		}
	    });
	////////////////////////////////////////////////////////////////////////////////////////
	// output
	for (String key : evaluations.keySet()) {
	    Evaluation evaluation = evaluations.get(key);
	    System.out.println("===== ["+key+"] ===== ");
	    System.out.println(evaluation.toSummaryString());
	    System.out.println(evaluation.toClassDetailsString());
	}
    }
}
