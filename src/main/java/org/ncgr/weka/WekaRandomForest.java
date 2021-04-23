package org.ncgr.weka;

import java.util.TreeMap;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;

import weka.classifiers.Classifier;
import weka.classifiers.Evaluation;
import weka.classifiers.trees.RandomForest;
import weka.classifiers.meta.AdaBoostM1;
import weka.classifiers.meta.LogitBoost;
import weka.core.Instance;
import weka.core.Instances;
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
 * Run Weka random forest operations on an ARFF file.
 *
 * setBagSizePercent(int newBagSizePercent)                           Sets the size of each bag, as a percentage of the training set size.
 * setBatchSize(String size)                                          Set the preferred batch size for batch prediction.
 * setBreakTiesRandomly(boolean newBreakTiesRandomly)                 Set whether to break ties randomly.
 * setComputeAttributeImportance(boolean computeAttributeImportance)  Set whether to compute and output attribute importance scores.
 * setDebug(boolean debug)                                            Set debugging mode.
 * setMaxDepth(int value)                                             Set the maximum depth of the tree, 0 for unlimited.
 * setNumDecimalPlaces(int num)                                       Set the number of decimal places.
 * setNumFeatures(int newNumFeatures)                                 Set the number of features to use in random selection.
 * setNumIterations(int numIterations)                                Sets the number of bagging iterations
 * setSeed(int s)                                                     Sets the seed for the random number generator.
 */
public class WekaRandomForest {
    /**
     * Main class for running all static methods.
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
	Option boostOption = new Option("bo", "boostoption", true, "boost method: Ada or Logit (none)");
	boostOption.setRequired(false);
	options.addOption(boostOption);
	// ACTIONS
	Option crossValidateOption = new Option("cv", "crossvalidate", false, "run k-fold cross-validation");
	crossValidateOption.setRequired(false);
	options.addOption(crossValidateOption);
	Option testOption = new Option("test", "test", false, "run training/testing on arfffile/testfile");
	testOption.setRequired(false);
	options.addOption(testOption);
	Option optimizeOption = new Option("op", "optimize", false, "optimize model by running noptimize training/testing runs and choosing best");
	optimizeOption.setRequired(false);
	options.addOption(optimizeOption);
	// PARAMETERS
	Option kfoldOption = new Option("k", "kfold", true, "cross-validation k-fold [10]");
	kfoldOption.setRequired(false);
	options.addOption(kfoldOption);
	Option bagSizePercentOption = new Option("bsp", "bagsizepercent", true, "the size of each bag, as a percentage of the training set size [100]");
	bagSizePercentOption.setRequired(false);
	options.addOption(bagSizePercentOption);
	Option breakTiesRandomlyOption = new Option("btr", "breaktiesrandomly", false, "break ties randomly");
	breakTiesRandomlyOption.setRequired(false);
	options.addOption(breakTiesRandomlyOption);
	Option debugOption = new Option("D", "debug", false, "turn on debugging mode");
	debugOption.setRequired(false);
	options.addOption(debugOption);
	Option maxDepthOption = new Option("md", "maxdepth", true, "maximum depth of the tree [0=unlimited]");
	maxDepthOption.setRequired(false);
	options.addOption(maxDepthOption);
	Option numFeaturesOption = new Option("nf", "numfeatures", true, "number of features to use in random selection [0=automatic]");
	numFeaturesOption.setRequired(false);
	options.addOption(numFeaturesOption);
	Option numIterationsOption = new Option("ni", "numiterations", true, "the number of bagging iterations [100]");
	numIterationsOption.setRequired(false);
	options.addOption(numIterationsOption);
	Option numOptimizeOption = new Option("no", "numoptimize", true, "the number of training/testing runs to optimize the model [10]");
	numOptimizeOption.setRequired(false);
	options.addOption(numOptimizeOption);

	try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println(e.getMessage());
            formatter.printHelp("WekaRandomForest", options);
            System.exit(1);
            return;
        }

	// actions
	final boolean crossValidate = cmd.hasOption("crossvalidate");
	final boolean test = cmd.hasOption("test");
	final boolean optimize = cmd.hasOption("optimize");

	// boosting
	String boostMethod = ""; // default = no boosting
	if (cmd.hasOption("boostoption")) boostMethod = cmd.getOptionValue("boostoption");

	// parameter defaults
	int kfold = 10;
	double bagSizePercent = 100.0;
	int maxDepth = 0;
	int numFeatures = 0;
	int numIterations = 100;
	int numOptimize = 10;
	// parameters
	if (cmd.hasOption("kfold")) kfold = Integer.parseInt(cmd.getOptionValue("kfold"));
	if (cmd.hasOption("bagsizepercent")) bagSizePercent = Double.parseDouble(cmd.getOptionValue("bagsizepercent"));
	if (cmd.hasOption("maxdepth")) maxDepth = Integer.parseInt(cmd.getOptionValue("maxdepth"));
	if (cmd.hasOption("numfeatures")) numFeatures = Integer.parseInt(cmd.getOptionValue("numfeatures"));
	if (cmd.hasOption("numiterations")) numIterations = Integer.parseInt(cmd.getOptionValue("numiterations"));
	if (cmd.hasOption("numoptimize")) numOptimize = Integer.parseInt(cmd.getOptionValue("numoptimize"));
	boolean breakTiesRandomly = cmd.hasOption("breaktiesrandomly");
	boolean debug = cmd.hasOption("debug");

	// Read instances from the input ARFF file
 	Instances data = Util.rearrange(new DataSource(cmd.getOptionValue("arfffile")).getDataSet());

	// do work
	if (crossValidate) {
	    // no boosting with cross-validation
	    double correct = crossValidate(data, kfold, numIterations, numFeatures, debug);
	} else if (test) {
	    Instances testData = Util.rearrange(new DataSource(cmd.getOptionValue("testfile")).getDataSet());
	    test(data, testData, boostMethod, numIterations, numFeatures, debug);
	} else if (optimize) {
	    Instances testData = Util.rearrange(new DataSource(cmd.getOptionValue("testfile")).getDataSet());
	    optimize(data, testData, boostMethod, numIterations, numFeatures, numOptimize, debug);
	}
    }

    /**
     * Run a cross-validation on the data for a given kfold, returning MCC.
     */
    public static double crossValidate(Instances data, int kfold, int numIterations, int numFeatures, boolean debug) throws Exception {
	Evaluation evaluation = new Evaluation(data);
	RandomForest rf = getRandomForest(numIterations, numFeatures);
	rf.setDebug(debug);
	System.out.println("# "+kfold+"-fold cross-validating (unboosted) RandomForest");
	Random random = new Random(ThreadLocalRandom.current().nextInt());
	evaluation.crossValidateModel(rf, data, kfold, random);
	// output
	Util.printResults(evaluation);
	return evaluation.matthewsCorrelationCoefficient(1);
    }

    /**
     * Train and test a model on the given data.
     */
    public static void test(Instances trainingData, Instances testingData, String boostMethod, int numIterations, int numFeatures, boolean debug) throws Exception {
	RandomForest model = getRandomForest(numIterations, numFeatures);
	Evaluation trainingEvaluation = new Evaluation(trainingData);
	Evaluation testingEvaluation = new Evaluation(testingData);
	int boostingRounds = 0;
	if (boostMethod.equals("Ada")) {
	    System.out.println("# training/testing AdaBoosted RandomForest with numIterations="+numIterations+" numFeatures="+numFeatures);
	    AdaBoostM1 boost = new AdaBoostM1();
	    boost.setDebug(debug);
	    boost.setClassifier(model);
	    boost.setUseResampling(false);
	    boost.setSeed(ThreadLocalRandom.current().nextInt());
	    boost.initializeClassifier(trainingData);
	    boostingRounds = Util.runBoosting(boost);               // boost
	    trainingEvaluation.evaluateModel(boost, trainingData);  // validate
	    testingEvaluation.evaluateModel(boost, testingData);    // test
	} else if (boostMethod.equals("Logit")) {
	    System.out.println("# training/testing LogitBoosted RandomForest with numIterations="+numIterations+" numFeatures="+numFeatures);
	    LogitBoost boost = new LogitBoost();
	    boost.setDebug(debug);
	    boost.setClassifier(model);
	    boost.setUseResampling(false);
	    boost.setSeed(ThreadLocalRandom.current().nextInt());
	    boost.initializeClassifier(trainingData);
	    boostingRounds = Util.runBoosting(boost);               // boost
	    trainingEvaluation.evaluateModel(boost, trainingData);  // validate
	    testingEvaluation.evaluateModel(boost, testingData);    // test
	} else {
	    System.out.println("# training/testing (unboosted) RandomForest with numIterations="+numIterations+" numFeatures="+numFeatures);
	    model.buildClassifier(trainingData);	                    // train
	    trainingEvaluation.evaluateModel(model, trainingData);	    // validate
	    testingEvaluation.evaluateModel(model, testingData);	    // test
	}
	// output
	System.out.println("# training result on "+trainingData.size()+" instances with "+boostingRounds+" boosting rounds");
	Util.printResults(trainingEvaluation);
	System.out.println("# testing result on "+testingData.size()+" instances with "+boostingRounds+" boosting rounds"); 
	Util.printResults(testingEvaluation);
    }

    /**
     * Optimize the model by training/testing numOptimize times and saving the best one (highest MCC).
     *
     * serialize model
     * weka.core.SerializationHelper.write("/some/where/j48.model", cls);
     *
     * deserialize model
     * Classifier cls = (Classifier) weka.core.SerializationHelper.read("/some/where/j48.model");
     */
    public static void optimize(Instances trainingData, Instances testingData, String boostMethod, int numIterations, int numFeatures, int numOptimize, boolean debug) throws Exception {
	TreeMap<Double,Classifier> models = new TreeMap<>();
	TreeMap<Double,Evaluation> evaluations = new TreeMap<>();
	if (boostMethod.equals("Ada")) {
	    System.out.println("# optimizing AdaBoosted RandomForest with numIterations="+numIterations+" numFeatures="+numFeatures+" numOptimize="+numOptimize);
	} else if (boostMethod.equals("Logit")) {
	    System.out.println("# optimizing LogitBoosted RandomForest with numIterations="+numIterations+" numFeatures="+numFeatures+" numOptimize="+numOptimize);
	} else {
	    System.out.println("# optimizing (unboosted) RandomForest with numIterations="+numIterations+" numFeatures="+numFeatures+" numOptimize="+numOptimize);
	}
	for (int round=1; round<=numOptimize; round++) {
	    RandomForest model = getRandomForest(numIterations, numFeatures);
	    Evaluation testingEvaluation = new Evaluation(testingData);
	    if (boostMethod.equals("Ada")) {
		AdaBoostM1 boost = new AdaBoostM1();
		boost.setDebug(debug);
		boost.setClassifier(model);
		boost.setUseResampling(false);
		boost.setSeed(ThreadLocalRandom.current().nextInt());
		boost.initializeClassifier(trainingData);
		int boostingRounds = Util.runBoosting(boost);             // boost
		testingEvaluation.evaluateModel(boost, testingData);  // test
		System.err.println("# "+boostingRounds+" boosting rounds");
		Util.printResults(testingEvaluation);                 // output
		// store
		evaluations.put(testingEvaluation.matthewsCorrelationCoefficient(1), testingEvaluation);
		models.put(testingEvaluation.matthewsCorrelationCoefficient(1), boost);
	    } else if (boostMethod.equals("Logit")) {
		LogitBoost boost = new LogitBoost();
		boost.setDebug(debug);
		boost.setClassifier(model);
		boost.setUseResampling(false);
		boost.setSeed(ThreadLocalRandom.current().nextInt());
		boost.initializeClassifier(trainingData);
		int boostingRounds = Util.runBoosting(boost);              // boost
		testingEvaluation.evaluateModel(boost, testingData);   // test
		System.err.println("# "+boostingRounds+" boosting rounds");
		Util.printResults(testingEvaluation);                  // output
		// store
		evaluations.put(testingEvaluation.matthewsCorrelationCoefficient(1), testingEvaluation);
		models.put(testingEvaluation.matthewsCorrelationCoefficient(1), boost);
	    } else {
		model.buildClassifier(trainingData);                     // train
		testingEvaluation.evaluateModel(model, testingData);     // test
		Util.printResults(testingEvaluation);                 // output
		// store
		evaluations.put(testingEvaluation.matthewsCorrelationCoefficient(1), testingEvaluation);
		models.put(testingEvaluation.matthewsCorrelationCoefficient(1), model);
	    }
	}
	// output
	Util.printStats(evaluations);
	// save
	System.out.println("# saving best model to file wekarf.model");
	SerializationHelper.write("wekarf.model", models.lastEntry().getValue());
    }

    /**
     * Create a RandomForest classifier with the given parameters.
     */
    public static RandomForest getRandomForest(int numIterations, int numFeatures) {
	RandomForest model = new RandomForest();
	model.setSeed(ThreadLocalRandom.current().nextInt());
	model.setNumDecimalPlaces(3);
	model.setNumExecutionSlots(0);
	model.setNumIterations(numIterations);
	model.setNumFeatures(numFeatures);
	return model;
    }
}
