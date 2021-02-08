 package org.ncgr.weka;

import java.text.DecimalFormat;

import java.util.List;
import java.util.LinkedList;
import java.util.Map;
import java.util.LinkedHashMap;
import java.util.Map.Entry;
import java.util.Random;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.concurrent.ConcurrentSkipListSet;
import java.util.concurrent.ThreadLocalRandom;

import weka.classifiers.Evaluation;
import weka.classifiers.functions.SMO;
import weka.classifiers.functions.supportVector.PolyKernel;
import weka.classifiers.functions.supportVector.RBFKernel;
import weka.core.Instance;
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
 * Run Weka SMO operations on an ARFF file.
 */
public class WekaSMO {
    // grid search ranges and steps
    static final int C_POWER_BEGIN = -5;
    static final int C_POWER_END = 15;
    static final int C_POWER_STEP = 1;
    static final int G_POWER_BEGIN = -25;
    static final int G_POWER_END = 0;
    static final int G_POWER_STEP = 1;

    // number formats
    static final DecimalFormat sci = new DecimalFormat("0.00E00");
    static final DecimalFormat perc = new DecimalFormat("00.0%");
    static final DecimalFormat rate = new DecimalFormat("0.000");
    static final DecimalFormat round = new DecimalFormat("0");
    
    public static void main(String[] args) throws Exception {
	Options options = new Options();
        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd;

	// REQUIRED
        Option arffFileOption = new Option("a", "arfffile", true, "input data file in ARFF format (required)");
        arffFileOption.setRequired(true);
        options.addOption(arffFileOption);
	// OPTIONAL
	Option verboseOption = new Option("v", "verbose", false, "verbosity toggle (false)");
	verboseOption.setRequired(false);
	options.addOption(verboseOption);
	//
	Option testFileOption = new Option("tf", "testfile", true, "test data file in ARFF format");
	testFileOption.setRequired(false);
	options.addOption(testFileOption);
	// ACTIONS
	Option gridSearchOption = new Option("gs", "gridsearch", false, "run a grid search over parameters");
	gridSearchOption.setRequired(false);
	options.addOption(gridSearchOption);
	//
	Option crossValidateOption = new Option("cv", "crossvalidate", false, "run cross-validation");
	crossValidateOption.setRequired(false);
	options.addOption(crossValidateOption);
	//
	Option testOption = new Option("test", "test", false, "run training/testing on datasets with given SMO and kernel parameters");
	testOption.setRequired(false);
	options.addOption(testOption);
	// PARAMETERS
	Option kfoldOption = new Option("k", "kfold", true, "cross-validation k-fold (required)");
	kfoldOption.setRequired(false);
	options.addOption(kfoldOption);
	//
	Option complexityOption = new Option("C", "complexity", true, "SMO complexity constant C");
	complexityOption.setRequired(false);
	options.addOption(complexityOption);
	//
	Option gammaOption = new Option("g", "gamma", true, "RBF kernel gamma parameter");
	gammaOption.setRequired(false);
	options.addOption(gammaOption);
	//
	Option normalizationOption = new Option("N", "normalization", true, "0=normalize, 1=standardize, 2=neither [0]");
	normalizationOption.setRequired(false);
	options.addOption(normalizationOption);
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
	//
	Option nGridSearchOption = new Option("ng", "ngridsearch", true, "set number of cases+controls (equal) to use in grid search [0=all]");
	nGridSearchOption.setRequired(false);
	options.addOption(nGridSearchOption);

	

	try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println(e.getMessage());
            formatter.printHelp("WekaSMO", options);
            System.exit(1);
            return;
        }

	boolean verbose = cmd.hasOption("verbose");

	// actions
	final boolean crossValidate = cmd.hasOption("crossvalidate");
	final boolean gridSearch = cmd.hasOption("gridsearch");
	final boolean test = cmd.hasOption("test");

	// cross-validation and grid search kfold
	int kfold = 0;
	if (crossValidate || gridSearch) {
	    kfold = Integer.parseInt(cmd.getOptionValue("kfold"));
	}

	// Read instances from the input ARFF file
 	Instances data = new DataSource(cmd.getOptionValue("arfffile")).getDataSet();
        // remove the ID attribute
        data.deleteAttributeAt(0);
        // set the class attribute index at the end (0/1)
        data.setClassIndex(data.numAttributes() - 1);

	// testing, cross-validation parameters
	double C = 0.0;
	double gamma = 0.0;
	if (test || crossValidate) {
	    C = Double.parseDouble(cmd.getOptionValue("C"));
	    gamma = Double.parseDouble(cmd.getOptionValue("gamma"));
	}

	// do work
	if (crossValidate) {
	    double correct = crossValidate(data, kfold, C, gamma);
	} else if (gridSearch) {
	    int nCases = 0;
	    int nControls = 0;
	    if (cmd.hasOption("ngridsearch")) {
		nCases = Integer.parseInt(cmd.getOptionValue("ngridsearch"))/2;
		nControls = nCases;
	    }
	    Instances searchData = reduceData(data, nCases, nControls);
	    searchGrid(searchData, kfold);
	} else if (test) {
	    Instances testData = new DataSource(cmd.getOptionValue("testfile")).getDataSet();
	    // remove the ID attribute
	    testData.deleteAttributeAt(0);
	    // set the class attribute index at the end (0/1)
	    testData.setClassIndex(data.numAttributes() - 1);
	    // test
	    test(data, testData, C, gamma);
	}
    }

    /**
     * Test a model on the given data.
     */
    public static void test(Instances trainingData, Instances testingData, double C, double gamma) throws Exception {
	// SMO model
	SMO model = new SMO();
	model.setC(C);
	RBFKernel kernel = new RBFKernel();
	kernel.setGamma(gamma);
	model.setKernel(kernel);
	// training
	model.buildClassifier(trainingData);
	Evaluation trainingEvaluation = new Evaluation(trainingData);
	System.err.println("Training SMO classifier on "+trainingData.size()+" instances with C="+sci.format(C)+" and gamma="+sci.format(gamma));
	trainingEvaluation.evaluateModel(model, trainingData);
	printResults(trainingEvaluation);
	// test
	System.err.println("Evaluating SMO classifier on "+testingData.size()+" instances.");
	Evaluation testingEvaluation = new Evaluation(testingData);
	testingEvaluation.evaluateModel(model, testingData);
	// output
	printResults(testingEvaluation);
    }

    /**
     * Print evaluation results.
     */
    public static void printResults(Evaluation evaluation) {
	double correct = evaluation.correct();
	double pctCorrect = evaluation.pctCorrect()/100.0;
	double TPR = evaluation.truePositiveRate(1);
	double FPR = evaluation.falsePositiveRate(1);
	System.out.println("correct\tpct\tTPR\tFPR");
	System.out.println(round.format(correct)+"\t"+pctCorrect+"\t"+rate.format(TPR)+"\t"+rate.format(FPR));
    }

    /**
     * Reduce the number of cases and controls from those in the given data.
     */
    public static Instances reduceData(Instances data, int nCases, int nControls) {
	if (nCases==0 || nControls==0) return data;
	System.err.println("Reducing data to "+nCases+" cases and "+nControls+" controls...");
	Instances reducedData = new Instances(data, (nCases+nControls));
	Random r = new Random(ThreadLocalRandom.current().nextInt());
	List<Integer> usedIndexes = new LinkedList<>();
	double caseValue = 1.0;
	double ctrlValue = 0.0;
	// load cases
	int caseCount = 0;
	while (caseCount<nCases) {
	    int i = r.nextInt(data.size());
	    if (usedIndexes.contains(i)) continue;
	    Instance instance = data.get(i);
	    if (instance.classValue()==caseValue) {
		caseCount++;
		usedIndexes.add(i);
		reducedData.add(instance);
	    }
	}
	// load controls
	int ctrlCount = 0;
	while (ctrlCount<nControls) {
	    int i = r.nextInt(data.size());
	    if (usedIndexes.contains(i)) continue;
	    Instance instance = data.get(i);
	    if (instance.classValue()==ctrlValue) {
		ctrlCount++;
		usedIndexes.add(i);
		reducedData.add(instance);
	    }
	}
	return reducedData;
    }

    /**
     * Run a grid search over C and gamma.
     */
    public static void searchGrid(Instances data, int kfold) {
	// parameter search grid
	ConcurrentSkipListSet<Double> CSet = new ConcurrentSkipListSet<>();
	for (int n=C_POWER_BEGIN; n<=C_POWER_END; n+=C_POWER_STEP) {
	    CSet.add(Math.pow(2.0, n));
	}
	ConcurrentSkipListSet<Double> gammaSet = new ConcurrentSkipListSet<>();
	for (int n=G_POWER_BEGIN; n<=G_POWER_END; n+=G_POWER_STEP) {
	    gammaSet.add(Math.pow(2.0, n));
	}
	// these maps are keyed by C:gamma
	ConcurrentSkipListMap<String,Double> totalCorrectMap = new ConcurrentSkipListMap<>();
	////////////////////////////////////////////////////////////////////////////////////////////////
	// C loop
	CSet.parallelStream().forEach(C -> {
		////////////////////////////////////////////////////////////////////////////////////////////////
		// gamma loop
		gammaSet.parallelStream().forEach(gamma -> {
			String key = C+":"+gamma;
			try {
			    double correct = crossValidate(data, kfold, C, gamma);
			    totalCorrectMap.put(key, correct);
			} catch (Exception ex) {
			    System.err.println(ex);
			    System.exit(1);
			}
		    });
		////////////////////////////////////////////////////////////////////////////////////////////////
	    });            
	////////////////////////////////////////////////////////////////////////////////////////////////
	// sort the totalCorrect map
	List<Entry<String,Double>> totalCorrectList = new LinkedList<>(totalCorrectMap.entrySet());
	totalCorrectList.sort(Entry.comparingByValue());
	Map<String,Double> sortedCorrectMap = new LinkedHashMap<>();
	for (Entry<String,Double> entry : totalCorrectList) {
	    String key = entry.getKey();
	    double totalCorrect = entry.getValue();
	    sortedCorrectMap.put(key, totalCorrect);
	}
	// output
	System.out.println("C\tgamma\tcorrect\tperc");
	for (String key : sortedCorrectMap.keySet()) {
	    String[] fields = key.split(":");
	    double C = Double.parseDouble(fields[0]);
	    double gamma = Double.parseDouble(fields[1]);
	    int totalCorrect = sortedCorrectMap.get(key).intValue();
	    System.out.println(sci.format(C)+"\t"+sci.format(gamma)+"\t"+totalCorrect+"\t"+perc.format((double)totalCorrect/(double)data.size()));
	}
    }
    
    /**
     * Run a cross-validation on the data for a given kfold, C and gamma.
     */
    public static double crossValidate(Instances data, int kfold, double C, double gamma) throws Exception {
	// // set model options
	// if (cmd.hasOption("N")) {
	//     switch (Integer.parseInt(cmd.getOptionValue("N"))) {
	//     case 1: model.setFilterType(new SelectedTag(SMO.FILTER_STANDARDIZE, SMO.TAGS_FILTER));
	// 	break;
	//     case 2: model.setFilterType(new SelectedTag(SMO.FILTER_NONE, SMO.TAGS_FILTER));
	// 	break;
	//     default: model.setFilterType(new SelectedTag(SMO.FILTER_NORMALIZE, SMO.TAGS_FILTER));
	// 	break;
	//     }
	// }
	// if (cmd.hasOption("L")) {
	//     model.setToleranceParameter(Double.parseDouble(cmd.getOptionValue("L")));
	// }
	// if (cmd.hasOption("P")) {
	//     model.setEpsilon(Double.parseDouble(cmd.getOptionValue("P")));
	// }
	// if (cmd.hasOption("M")) {
	//     model.setBuildCalibrationModels(true);
	// }
	// run the cross-validation
	Evaluation evaluation = new Evaluation(data);
	SMO model = new SMO();
	model.turnChecksOn();
	// complexity
	model.setC(C);
	// kernel
	RBFKernel kernel = new RBFKernel();
	kernel.setGamma(gamma);
	model.setKernel(kernel);
	evaluation.crossValidateModel(model, data, kfold, new Random(ThreadLocalRandom.current().nextInt(0, 100)));
	// assume case=1
	double correct = evaluation.correct();
	double pctCorrect = evaluation.pctCorrect();
	double TPR = evaluation.truePositiveRate(1);
	double FPR = evaluation.falsePositiveRate(1);
	double MCC = evaluation.matthewsCorrelationCoefficient(1);
	double precision = evaluation.precision(1);
	double recall = evaluation.recall(1);
	// DEBUG
	System.err.println(evaluation.toSummaryString());
	//
	return correct;
    }
}

