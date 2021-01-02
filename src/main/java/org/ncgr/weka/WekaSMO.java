package org.ncgr.weka;

import java.text.DecimalFormat;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.ConcurrentHashMap;

import weka.classifiers.Evaluation;
import weka.classifiers.functions.SMO;
import weka.core.Instances;
import weka.core.SelectedTag;
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
 * Run Weka SMO cross-validation on an ARFF file.
 */
public class WekaSMO {
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
	//
	Option nrunsOption = new Option("nruns", "nruns", true, "number of parallel runs [10]");
	nrunsOption.setRequired(false);
	options.addOption(nrunsOption);
	//
	Option normalizationOption = new Option("N", "normalization", true, "0=normalize, 1=standardize, 2=neither [0]");
	normalizationOption.setRequired(false);
	options.addOption(normalizationOption);
	//
	Option complexityOption = new Option("C", "complexity", true, "complexity constant C [1.0]");
	complexityOption.setRequired(false);
	options.addOption(complexityOption);
	//
	Option toleranceOption = new Option("L", "tolerance", true, "tolerance parameter [1.0e-3]");
	toleranceOption.setRequired(false);
	options.addOption(toleranceOption);
	//
	Option epsilonOption = new Option("P", "epsilon", true, "epsilon for round-off error [1.0e-12]");
	epsilonOption.setRequired(false);
	options.addOption(epsilonOption);
	//
	Option fitCalibrationModelsOption = new Option("M", "fitcalibrationmodels", false, "fit calibration models to SVM outputs");
	fitCalibrationModelsOption.setRequired(false);
	options.addOption(fitCalibrationModelsOption);
	
	// Option kernelOption = new Option("K", "kernel", true, "the kernel to use [weka.classifiers.functions.supportVector.PolyKernel]");
	// kernelOption.setRequired(false);
	// options.addOption(kernelOption);
	// Options specific to kernel weka.classifiers.functions.supportVector.PolyKernel:
	// -E <num>
	//  The Exponent to use.
	//  (default: 1.0)
	// -L
	//  Use lower-order terms.
	//  (default: no)
	// -C <num>
	//  The size of the cache (a prime number), 0 for full cache and 
	//  -1 to turn it off.
	//  (default: 250007)
	//
	// -calibrator <scheme specification>
	//  Full name of calibration model, followed by options.
	//  (default: "weka.classifiers.functions.Logistic")
	// Options specific to calibrator weka.classifiers.functions.Logistic:
	// -C
	//  Use conjugate gradient descent rather than BFGS updates.
	// -R <ridge>
	//  Set the ridge in the log-likelihood.
	// -M <number>
	//  Set the maximum number of iterations (default -1, until convergence).
	
	try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println(e.getMessage());
            formatter.printHelp("WekaSMO", options);
            System.exit(1);
            return;
        }

	// our ARFF file
        final String arffFile = cmd.getOptionValue("arfffile");

	// cross-validation kfold
	final int kfold = Integer.parseInt(cmd.getOptionValue("kfold"));

	// number of runs
	int nruns = 10; 
	if (cmd.hasOption("nruns")) nruns = Integer.parseInt(cmd.getOptionValue("nruns"));
	
	// Read all the instances in the file (ARFF, CSV, XRFF, ...)
	DataSource source = new DataSource(arffFile);
 	Instances data = source.getDataSet();
        // remove the ID attribute
        data.deleteAttributeAt(0);
        // set the class attribute index
        data.setClassIndex(data.numAttributes() - 1);
 
	// populate a map with nrun Evaluations for parallel runs
	// use key to seed random number generator
	ConcurrentHashMap<Integer,Evaluation> evaluations = new ConcurrentHashMap<>();
	ConcurrentHashMap<Integer,SMO> models = new ConcurrentHashMap<>();
	for (int i=0; i<nruns; i++) {
	    evaluations.put(i, new Evaluation(data));
	    SMO model = new SMO();
	    models.put(i, model);
	    model.turnChecksOn();
	    // set model options
	    if (cmd.hasOption("N")) {
		switch (Integer.parseInt(cmd.getOptionValue("N"))) {
		case 1: model.setFilterType(new SelectedTag(SMO.FILTER_STANDARDIZE, SMO.TAGS_FILTER));
		    break;
		case 2: model.setFilterType(new SelectedTag(SMO.FILTER_NONE, SMO.TAGS_FILTER));
		    break;
		default: model.setFilterType(new SelectedTag(SMO.FILTER_NORMALIZE, SMO.TAGS_FILTER));
		    break;
		}
	    }
	    if (cmd.hasOption("C")) {
		model.setC(Double.parseDouble(cmd.getOptionValue("C")));
	    }
	    if (cmd.hasOption("L")) {
		model.setToleranceParameter(Double.parseDouble(cmd.getOptionValue("L")));
	    }
	    if (cmd.hasOption("P")) {
		model.setEpsilon(Double.parseDouble(cmd.getOptionValue("P")));
	    }
	    if (cmd.hasOption("M")) {
		model.setBuildCalibrationModels(true);
	    }
	}
	
	/////////////////////////////////////////////////////////////////////////////////
	// run the cross-validations
	evaluations.entrySet().parallelStream().forEach(entry -> {
		int key = entry.getKey();
		Evaluation evaluation = entry.getValue();
		SMO model = models.get(key);
		try {
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
	    System.out.println("===== ["+(key+1)+"] ===== ");
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
	System.out.println("Correct\t% correct\tTPR\t\tFPR\t\tMCC\t\tPrecision\tRecall");
	System.out.println(intf.format(StatUtils.mean(corrects))+pm+intf.format(Math.sqrt(StatUtils.variance(corrects)))+"\t"+
			   perf.format(StatUtils.mean(pctCorrects))+pm+perf.format(Math.sqrt(StatUtils.variance(pctCorrects)))+"\t"+
			   decf.format(StatUtils.mean(TPRs))+pm+decf.format(Math.sqrt(StatUtils.variance(TPRs)))+"\t"+
			   decf.format(StatUtils.mean(FPRs))+pm+decf.format(Math.sqrt(StatUtils.variance(FPRs)))+"\t"+
			   decf.format(StatUtils.mean(MCCs))+pm+decf.format(Math.sqrt(StatUtils.variance(MCCs)))+"\t"+
			   decf.format(StatUtils.mean(precisions))+pm+decf.format(Math.sqrt(StatUtils.variance(precisions)))+"\t"+
			   decf.format(StatUtils.mean(recalls))+pm+decf.format(Math.sqrt(StatUtils.variance(recalls))));
	System.out.println("-------------------------------------------------------------------------------------------------------------");
    }
}
