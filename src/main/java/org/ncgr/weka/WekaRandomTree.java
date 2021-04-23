package org.ncgr.weka;

import java.text.DecimalFormat;
import java.util.List;
import java.util.LinkedList;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.ConcurrentHashMap;

import weka.classifiers.Evaluation;
import weka.classifiers.trees.RandomTree;
import weka.core.Instances;
import weka.core.converters.ConverterUtils.DataSource;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import org.apache.commons.math3.stat.StatUtils;
 
/**
 * Run Weka RandomTree cross-validation on an ARFF file.
 */
public class WekaRandomTree {
    static DecimalFormat decf = new DecimalFormat("0.000");
    static DecimalFormat intf = new DecimalFormat("0");
    static DecimalFormat perf = new DecimalFormat("0.0");
    static String pm = "\u00B1";
    
    public static void main(String[] args) throws Exception {
	Options options = new Options();
        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd;

        Option arffFileOption = new Option("a", "arfffile", true, "ARFF input file (required)");
        arffFileOption.setRequired(true);
        options.addOption(arffFileOption);
	//
	Option kfoldOption = new Option("k", "kfold", true, "cross-validation k-fold (required)");
	kfoldOption.setRequired(true);
	options.addOption(kfoldOption);
	
	try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println(e.getMessage());
            formatter.printHelp("WekaRandomTree", options);
            System.exit(1);
            return;
        }

        String arffFile = cmd.getOptionValue("arfffile");

	// cross-validation options
	final int kfold = Integer.parseInt(cmd.getOptionValue("kfold"));

	// Read all the instances in the file (ARFF, CSV, XRFF, ...)
	Instances data = Util.rearrange(new DataSource(arffFile).getDataSet());
 
	// populate a map with 10 Evaluations for parallel runs
	Map<Integer,Evaluation> evaluations = new ConcurrentHashMap<>();
	for (int i=0; i<10; i++) {
	    evaluations.put(i, new Evaluation(data));
	}
	
	/////////////////////////////////////////////////////////////////////////////////
	// run the cross-validations
	evaluations.entrySet().parallelStream().forEach(entry -> {
		int key = entry.getKey();
		Evaluation evaluation = entry.getValue();
		try {
		    RandomTree model = new RandomTree();
		    evaluation.crossValidateModel(model, data, kfold, new Random(key));
		} catch (Exception e) {
		    System.err.println(e);
		    System.exit(1);
		}
	    });
	////////////////////////////////////////////////////////////////////////////////////////

	// arrays for summary stats
	double[] corrects = new double[evaluations.size()];
	double[] pctCorrects = new double[evaluations.size()];
	double[] TPRs = new double[evaluations.size()];
	double[] FPRs = new double[evaluations.size()];
	double[] MCCs = new double[evaluations.size()];
	double[] precisions = new double[evaluations.size()];
	double[] recalls = new double[evaluations.size()];
	// output
	for (int key : evaluations.keySet()) {
	    Evaluation evaluation = evaluations.get(key);
	    System.out.println("===== ["+key+"] ===== ");
	    System.out.println(evaluation.toSummaryString());
	    System.out.println(evaluation.toClassDetailsString());
	    // assume case=1
	    corrects[key] = evaluation.correct();
	    pctCorrects[key] = evaluation.pctCorrect();
	    TPRs[key] = evaluation.truePositiveRate(1);
	    FPRs[key] = evaluation.falsePositiveRate(1);
	    MCCs[key] = evaluation.matthewsCorrelationCoefficient(1);
	    precisions[key] = evaluation.precision(1);
	    recalls[key] = evaluation.recall(1);
	}
	// summary stats
	System.out.println("-------------------------------------------------------------------------------------------------------------");
	System.out.println("Correct\t\tperCorrect\tTruePosRate\tFalsePosRate\tMathCorrC\tPrecision\tRecall");
	System.out.println(intf.format(StatUtils.mean(corrects))+pm+intf.format(Math.sqrt(StatUtils.variance(corrects)))+"\t\t"+
			   perf.format(StatUtils.mean(pctCorrects))+pm+perf.format(Math.sqrt(StatUtils.variance(pctCorrects)))+"\t"+
			   decf.format(StatUtils.mean(TPRs))+pm+decf.format(Math.sqrt(StatUtils.variance(TPRs)))+"\t"+
			   decf.format(StatUtils.mean(FPRs))+pm+decf.format(Math.sqrt(StatUtils.variance(FPRs)))+"\t"+
			   decf.format(StatUtils.mean(MCCs))+pm+decf.format(Math.sqrt(StatUtils.variance(MCCs)))+"\t"+
			   decf.format(StatUtils.mean(precisions))+pm+decf.format(Math.sqrt(StatUtils.variance(precisions)))+"\t"+
			   decf.format(StatUtils.mean(recalls))+pm+decf.format(Math.sqrt(StatUtils.variance(recalls))));
	System.out.println("-------------------------------------------------------------------------------------------------------------");
    }
}
