package org.ncgr.weka;

import java.text.DecimalFormat;

import java.util.List;
import java.util.LinkedList;
import java.util.Map;
import java.util.LinkedHashMap;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.Random;
import java.util.Set;
import java.util.HashSet;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ThreadLocalRandom;

import weka.classifiers.Classifier;
import weka.classifiers.Evaluation;
import weka.classifiers.functions.SMO;
import weka.classifiers.functions.supportVector.Kernel;
import weka.classifiers.functions.supportVector.RBFKernel;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.SelectedTag;
import weka.core.SerializationHelper;
import weka.core.converters.ConverterUtils.DataSource;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
 
/**
 * Run Weka SMO operations on an ARFF file. RBF kernel is hardcoded.
 */
public class WekaSMO {
    // number formats
    static final DecimalFormat sci = new DecimalFormat("0.00E00");
    static final DecimalFormat perc = new DecimalFormat("0.0%");
    static final DecimalFormat rate = new DecimalFormat("0.000");
    
    /**
     * Main class to run static methods.
     */
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
	Option testFileOption = new Option("tf", "testfile", true, "test data file in ARFF format");
	testFileOption.setRequired(false);
	options.addOption(testFileOption);
	// ACTIONS
	Option gridSearchOption = new Option("gs", "gridsearch", false, "run a grid search over (C,gamma); test data will be predicted for each pair if testfile is provided");
	gridSearchOption.setRequired(false);
	options.addOption(gridSearchOption);
	Option crossValidateOption = new Option("cv", "crossvalidate", false, "run cross-validation");
	crossValidateOption.setRequired(false);
	options.addOption(crossValidateOption);
	Option testOption = new Option("test", "test", false, "run training/testing on datasets with given SMO and RBF kernel parameters");
	testOption.setRequired(false);
	options.addOption(testOption);
	Option optimizeOption = new Option("op", "optimize", false, "optimize model by running noptimize training/testing runs and choosing highest MCC");
	optimizeOption.setRequired(false);
	options.addOption(optimizeOption);
	// PARAMETERS
	Option kfoldOption = new Option("k", "kfold", true, "cross-validation k-fold (required)");
	kfoldOption.setRequired(false);
	options.addOption(kfoldOption);
	Option complexityOption = new Option("C", "complexity", true, "SMO complexity constant C");
	complexityOption.setRequired(false);
	options.addOption(complexityOption);
	Option gammaOption = new Option("gamma", "gamma", true, "RBF kernel parameter gamma");
	gammaOption.setRequired(false);
	options.addOption(gammaOption);
	Option normalizationOption = new Option("N", "normalization", true, "0=normalize, 1=standardize, 2=neither [0]");
	normalizationOption.setRequired(false);
	options.addOption(normalizationOption);
	Option toleranceOption = new Option("L", "tolerance", true, "tolerance parameter [1.0e-3]");
	toleranceOption.setRequired(false);
	options.addOption(toleranceOption);
	Option epsilonOption = new Option("P", "epsilon", true, "epsilon for round-off error [1.0e-12]");
	epsilonOption.setRequired(false);
	options.addOption(epsilonOption);
	Option fitCalibrationModelsOption = new Option("M", "fitcalibrationmodels", false, "fit calibration models to SVM outputs");
	fitCalibrationModelsOption.setRequired(false);
	options.addOption(fitCalibrationModelsOption);
	Option nGridSearchOption = new Option("ng", "ngridsearch", true, "set number of cases+controls (equal) to use in grid search [0=all]");
	nGridSearchOption.setRequired(false);
	options.addOption(nGridSearchOption);
	Option numOptimizeOption = new Option("no", "numoptimize", true, "the number of training/testing runs to optimize the model [10]");
	numOptimizeOption.setRequired(false);
	options.addOption(numOptimizeOption);
	Option debugOption = new Option("D", "debug", false, "toggle on debug mode [off]");
	debugOption.setRequired(false);
	options.addOption(debugOption);
        Option CPowerBeginOption = new Option("cb", "cpowerbegin", true, "begining 2^n exponent in C scan [5]");
        CPowerBeginOption.setRequired(false);
        options.addOption(CPowerBeginOption);
        Option CPowerEndOption = new Option("ce", "cpowerend", true, "ending 2^n exponent in C scan [11]");
        CPowerEndOption.setRequired(false);
        options.addOption(CPowerEndOption);
        Option CPowerStepOption = new Option("cs", "cpowerstep", true, "step for 2^n exponent in C scan [1]");
        CPowerStepOption.setRequired(false);
        options.addOption(CPowerStepOption);
        Option gammaPowerBeginOption = new Option("gb", "gammapowerbegin", true, "begining 2^n exponent in RBF gamma scan [-30]");
        gammaPowerBeginOption.setRequired(false);
        options.addOption(gammaPowerBeginOption);
        Option gammaPowerEndOption = new Option("ge", "gammapowerend", true, "ending 2^n exponent in RBF gamma scan [-10]");
        gammaPowerEndOption.setRequired(false);
        options.addOption(gammaPowerEndOption);
        Option gammaPowerStepOption = new Option("gs", "gammapowerstep", true, "step for 2^n exponent in RBF gamma scan [1]");
        gammaPowerStepOption.setRequired(false);
        options.addOption(gammaPowerStepOption);

	try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println(e.getMessage());
            formatter.printHelp("WekaSMO", options);
            System.exit(1);
            return;
        }

	// actions
	final boolean crossValidate = cmd.hasOption("crossvalidate");
	final boolean gridSearch = cmd.hasOption("gridsearch");
	final boolean test = cmd.hasOption("test");
	final boolean optimize = cmd.hasOption("optimize");

	// parameters
	int kfold = 10;
	double C = 0.0;
	double gamma = 0.0;
	int numOptimize = 10;
	if (cmd.hasOption("complexity")) C = Double.parseDouble(cmd.getOptionValue("C"));
	if (cmd.hasOption("gamma")) gamma = Double.parseDouble(cmd.getOptionValue("gamma"));
	if (cmd.hasOption("kfold")) kfold = Integer.parseInt(cmd.getOptionValue("kfold"));
	if (cmd.hasOption("numoptimize")) numOptimize = Integer.parseInt(cmd.getOptionValue("numoptimize"));
	boolean debug = cmd.hasOption("debug");

	int CPowerBegin = 5;
	int CPowerEnd = 11;
	int CPowerStep = 1;
	int gammaPowerBegin = -30;
	int gammaPowerEnd = -10;
	int gammaPowerStep = 1;
        if (cmd.hasOption("cpowerbegin")) CPowerBegin = Integer.parseInt(cmd.getOptionValue("cpowerbegin"));
        if (cmd.hasOption("cpowerend")) CPowerEnd = Integer.parseInt(cmd.getOptionValue("cpowerend"));
        if (cmd.hasOption("cpowerstep")) CPowerStep = Integer.parseInt(cmd.getOptionValue("cpowerstep"));
        if (cmd.hasOption("gammapowerbegin")) gammaPowerBegin = Integer.parseInt(cmd.getOptionValue("gammapowerbegin"));
        if (cmd.hasOption("gammapowerend")) gammaPowerEnd = Integer.parseInt(cmd.getOptionValue("gammapowerend"));
        if (cmd.hasOption("gammapowerstep")) gammaPowerStep = Integer.parseInt(cmd.getOptionValue("gammapowerstep"));
	
	// Read training instances from the input ARFF files
 	Instances trainingData = Util.rearrange(new DataSource(cmd.getOptionValue("arfffile")).getDataSet());

	// do work
	if (crossValidate) {
	    double mcc = crossValidate(trainingData, kfold, C, gamma, debug);
	} else if (gridSearch) {
	    int nCases = 0;
	    int nControls = 0;
	    if (cmd.hasOption("ngridsearch")) {
		nCases = Integer.parseInt(cmd.getOptionValue("ngridsearch"))/2;
		nControls = nCases;
	    }
	    Instances searchData = Util.reduceData(trainingData, nCases, nControls);
	    if (cmd.hasOption("testfile")) {
		Instances testData = Util.rearrange(new DataSource(cmd.getOptionValue("testfile")).getDataSet());
		searchGrid(searchData, testData, kfold, CPowerBegin, CPowerEnd, CPowerStep, gammaPowerBegin, gammaPowerEnd, gammaPowerStep);
	    } else {
		searchGrid(searchData, kfold, CPowerBegin, CPowerEnd, CPowerStep, gammaPowerBegin, gammaPowerEnd, gammaPowerStep);
	    }
	} else if (test) {
	    Instances testData = Util.rearrange(new DataSource(cmd.getOptionValue("testfile")).getDataSet());
	    trainTest(trainingData, testData, C, gamma, true);
	} else if (optimize) {
	    Instances testData = Util.rearrange(new DataSource(cmd.getOptionValue("testfile")).getDataSet());
	    optimize(trainingData, testData, C, gamma, numOptimize);
	}
    }

    /**
     * Run a cross-validation on the data for a given kfold, C and kernel parameter, returning the MCC.
     */
    public static double crossValidate(Instances trainingData, int kfold, double C, double gamma, boolean printOutput) throws Exception {
	SMO smo = getSMO(C, gamma);
	Evaluation evaluation = new Evaluation(trainingData);
	evaluation.crossValidateModel(smo, trainingData, kfold, new Random(ThreadLocalRandom.current().nextInt()));
	if (printOutput) {
	    // output
	    System.err.println("# training result on "+trainingData.size()+" instances:");
	    Util.printResults(evaluation, C, gamma);
	}
	return evaluation.matthewsCorrelationCoefficient(1);
    }

    /**
     * Test a model on the given training and testing data, returning the MCC.
     */
    public static double trainTest(Instances trainingData, Instances testingData, double C, double gamma, boolean printOutput) throws Exception {
	SMO model = getSMO(C, gamma);
	Evaluation trainingEvaluation = new Evaluation(trainingData);
	Evaluation testingEvaluation = new Evaluation(testingData);
	if (printOutput) {
	    System.err.println("# testing/training SMO classifier with C="+sci.format(C)+" RBF gamma="+sci.format(gamma));
	}
	model.buildClassifier(trainingData);	                // train
	trainingEvaluation.evaluateModel(model, trainingData);  // validate
	testingEvaluation.evaluateModel(model, testingData);	// test
	if (printOutput) {
	    // output
	    System.err.println("# training result on "+trainingData.size()+" instances:");
	    Util.printResults(trainingEvaluation, C, gamma);
	    System.err.println("# testing result on "+testingData.size()+" instances:"); 
	    Util.printResults(testingEvaluation, C, gamma);
	}
	return testingEvaluation.matthewsCorrelationCoefficient(1);
    }

    /**
     * Run a k-fold cross-validation grid search over C and kernel param. Defines "best" as largest training MCC.
     */
    public static void searchGrid(Instances trainingData, int kfold,
				  int CPowerBegin, int CPowerEnd, int CPowerStep, int gammaPowerBegin, int gammaPowerEnd, int gammaPowerStep) {
	System.err.println("# searching (C,param) grid with "+kfold+"-fold cross-validation on "+trainingData.size()+" training instances");
	// parameter search grid - use final sets just because
	Set<Double> CTempSet = new HashSet<>();
	for (int n=CPowerBegin; n<=CPowerEnd; n+=CPowerStep) {
	    CTempSet.add(Math.pow(2.0, n));
	}
	final Set<Double> CSet = new HashSet<>(CTempSet);
	Set<Double> gammaTempSet = new HashSet<>();
	for (int n=gammaPowerBegin; n<=gammaPowerEnd; n+=gammaPowerStep) {
	    gammaTempSet.add(Math.pow(2.0, n));
	}
	final Set<Double> gammaSet = new HashSet<>(gammaTempSet);
	// this map is keyed by C:gamma
	Map<String,Double> mccMap = new ConcurrentHashMap<>();
	// C loop in series
	for (double C : CSet) {
	    ////////////////////////////////////////////////////////////////////////////////////////////////
	    // kernel parameter loop in parallel
	    gammaSet.parallelStream().forEach(gamma -> {
		    final String key = C+":"+gamma;
		    try {
			double mcc = crossValidate(trainingData, kfold, C, gamma, false);
			mccMap.put(key, mcc);
		    } catch (Exception ex) {
			System.err.println(ex);
			System.exit(1);
		    }
		});
	    ////////////////////////////////////////////////////////////////////////////////////////////////
	}
	// sort the mcc map
	List<Entry<String,Double>> mccList = new LinkedList<>(mccMap.entrySet());
	mccList.sort(Entry.comparingByValue());
	// output
	System.out.println("#-------------------------------------------------------------------------");
	for (Entry<String,Double> entry : mccList) {
	    String key = entry.getKey();
	    double mcc = entry.getValue();
	    String[] fields = key.split(":");
	    double C = Double.parseDouble(fields[0]);
	    double param = Double.parseDouble(fields[1]);
	    System.out.println(sci.format(C)+"\t"+sci.format(param)+"\t"+rate.format(mcc));
	}
	System.out.println("#-------------------------------------------------------------------------");
    }

    /**
     * Run a k-fold cross-validation grid search over C and kernel param, running prediction on a testing set for each. Defines "best" as largest test MCC.
     * NOTE: this is NOT a legit training/testing run since you must optimize training and live with the testing result.
     */
    public static void searchGrid(Instances trainingData, Instances testingData, int kfold,
				  int CPowerBegin, int CPowerEnd, int CPowerStep, int gammaPowerBegin, int gammaPowerEnd, int gammaPowerStep) {
	System.err.println("# searching (C,param) grid with "+kfold+"-fold cross-validation on "+trainingData.size()+" and "+testingData.size()+" testing instances.");
	// parameter search grid - use final sets just because
	Set<Double> CTempSet = new HashSet<>();
	for (int n=CPowerBegin; n<=CPowerEnd; n+=CPowerStep) {
	    CTempSet.add(Math.pow(2.0, n));
	}
	final Set<Double> CSet = new HashSet<>(CTempSet);
	Set<Double> gammaTempSet = new HashSet<>();
	for (int n=gammaPowerBegin; n<=gammaPowerEnd; n+=gammaPowerStep) {
	    gammaTempSet.add(Math.pow(2.0, n));
	}
	final Set<Double> gammaSet = new HashSet<>(gammaTempSet);
	// this map is keyed by C:gamma
	Map<String,Double> trainingMCCMap = new ConcurrentHashMap<>();
	Map<String,Double> testingMCCMap = new ConcurrentHashMap<>();
	// C loop in series
	for (double C : CSet) {
	    ////////////////////////////////////////////////////////////////////////////////////////////////
	    // kernel parameter loop in parallel
	    gammaSet.parallelStream().forEach(gamma -> {
		    final String key = C+":"+gamma;
		    try {
			double trainingMCC = crossValidate(trainingData, kfold, C, gamma, false);
			double testingMCC = trainTest(trainingData, testingData, C, gamma, false);
			trainingMCCMap.put(key, trainingMCC);
			testingMCCMap.put(key, testingMCC);
		    } catch (Exception ex) {
			System.err.println(ex);
			System.exit(1);
		    }
		});
	    ////////////////////////////////////////////////////////////////////////////////////////////////
	}
	// sort the mcc map
	List<Entry<String,Double>> mccList = new LinkedList<>(testingMCCMap.entrySet());
	mccList.sort(Entry.comparingByValue());
	// output
	System.out.println("#-------------------------------------------------------------------------");
	for (Entry<String,Double> entry : mccList) {
	    String key = entry.getKey();
	    double testingMCC = entry.getValue();
	    String[] fields = key.split(":");
	    double C = Double.parseDouble(fields[0]);
	    double param = Double.parseDouble(fields[1]);
	    double trainingMCC = trainingMCCMap.get(key);
	    System.out.println(sci.format(C)+"\t"+sci.format(param)+"\t"+rate.format(trainingMCC)+"\t"+rate.format(testingMCC));
	}
	System.out.println("#-------------------------------------------------------------------------");
    }

    /**
     * Optimize the model by training/testing numOptimize times and saving the best one (highest mcc).
     *
     * serialize model
     * weka.core.SerializationHelper.write("/some/where/j48.model", cls);
     *
     * deserialize model
     * Classifier cls = (Classifier) weka.core.SerializationHelper.read("/some/where/j48.model");
     */
    public static void optimize(Instances trainingData, Instances testingData, double C, double gamma, int numOptimize) throws Exception {
	TreeMap<Double,Classifier> models = new TreeMap<>();
	TreeMap<Double,Evaluation> evaluations = new TreeMap<>();
	System.err.println("# optimizing SMO with C="+sci.format(C)+" RBF gamma="+sci.format(gamma));
    	for (int round=1; round<=numOptimize; round++) {
	    SMO model = getSMO(C, gamma);
	    Evaluation testingEvaluation = new Evaluation(testingData);
	    model.buildClassifier(trainingData);                  // train
	    testingEvaluation.evaluateModel(model, testingData);  // test
	    Util.printResults(testingEvaluation);                 // output
	    // store
	    evaluations.put(testingEvaluation.matthewsCorrelationCoefficient(1), testingEvaluation);
	    models.put(testingEvaluation.matthewsCorrelationCoefficient(1), model);
	}
	// output
	Util.printStats(evaluations);
	// save
	System.out.println("# saving best model to file wekasmo.model");
	SerializationHelper.write("wekasmo.model", models.lastEntry().getValue());
    }

    /**
     * Instantiate an SMO classifier with the given RBF parameters.
     */
    public static SMO getSMO(double C, double gamma) {
	SMO model = new SMO();
	int modelseed = ThreadLocalRandom.current().nextInt();
	model.setRandomSeed(modelseed);
	model.setNumDecimalPlaces(3);
	model.setC(C);
	RBFKernel kernel = new RBFKernel();
	kernel.setGamma(gamma);
	model.setKernel(kernel);
	return model;
    }
}
