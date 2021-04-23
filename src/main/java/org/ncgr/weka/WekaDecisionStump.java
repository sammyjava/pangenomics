package org.ncgr.weka;

import java.util.Random;
import java.util.TreeMap;
import java.util.concurrent.ThreadLocalRandom;

import weka.classifiers.Classifier;
import weka.classifiers.Evaluation;
import weka.classifiers.meta.AdaBoostM1;
import weka.classifiers.meta.LogitBoost;
import weka.classifiers.trees.DecisionStump;
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
 * Run Weka DecisionStump on ARFF files.
 */
public class WekaDecisionStump {
    /**
     * Main class runs static methods.
     */
    public static void main(String[] args) throws Exception {
	Options options = new Options();
        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd;

	// REQUIRED
        Option arffFileOption = new Option("a", "arfffile", true, "ARFF input file (required)");
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
	Option crossValidateOption = new Option("cv", "crossvalidate", false, "run cross-validation");
	crossValidateOption.setRequired(false);
	options.addOption(crossValidateOption);
	Option testOption = new Option("test", "test", false, "run training/testing on datasets");
	testOption.setRequired(false);
	options.addOption(testOption);
	Option optimizeOption = new Option("op", "optimize", false, "optimize model by running noptimize training/testing runs and choosing best");
	optimizeOption.setRequired(false);
	options.addOption(optimizeOption);
	// PARAMETERS
	Option kfoldOption = new Option("k", "kfold", true, "cross-validation k-fold [10]");
	kfoldOption.setRequired(false);
	options.addOption(kfoldOption);
	Option numOptimizeOption = new Option("no", "numoptimize", true, "the number of training/testing runs to optimize the model [10]");
	numOptimizeOption.setRequired(false);
	options.addOption(numOptimizeOption);
	Option debugOption = new Option("D", "debug", false, "toggle on debug mode [false]");
	debugOption.setRequired(false);
	options.addOption(debugOption);
	
	try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println(e.getMessage());
            formatter.printHelp("WekaDecisionStump", options);
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
	
	// parameters
	int kfold = 10;
	int numOptimize = 10;
	if (cmd.hasOption("kfold")) kfold = Integer.parseInt(cmd.getOptionValue("kfold"));
	if (cmd.hasOption("numoptimize")) numOptimize = Integer.parseInt(cmd.getOptionValue("numoptimize"));
	boolean debug = cmd.hasOption("debug");

	// Read instances from the input ARFF files
 	Instances data = Util.rearrange(new DataSource(cmd.getOptionValue("arfffile")).getDataSet());

	// do work
	if (crossValidate) {
	    // no boosting with cross-validation
	    crossValidate(data, kfold, debug);
	} else if (test) {
	    Instances testData = Util.rearrange(new DataSource(cmd.getOptionValue("testfile")).getDataSet());
	    test(data, testData, boostMethod, debug);
	} else if (optimize) {
	    Instances testData = Util.rearrange(new DataSource(cmd.getOptionValue("testfile")).getDataSet());
	    optimize(data, testData, boostMethod, numOptimize, debug);
	}
    }

    /**
     * Run a cross-validation on the data.
     */
    public static void crossValidate(Instances data, int kfold, boolean debug) throws Exception {
	// run the cross-validation
	DecisionStump model = getDecisionStump();
	model.setDebug(debug);
	Evaluation evaluation = new Evaluation(data);
	System.err.println("# cross-validating DecisionStump with kfold="+kfold);
	Random random = new Random(ThreadLocalRandom.current().nextInt());
	evaluation.crossValidateModel(model, data, kfold, random);
	// output
	Util.printResults(evaluation);
    }

    /**
     * Train/Test a model on the given data.
     */
    public static void test(Instances trainingData, Instances testingData, String boostMethod, boolean debug) throws Exception {
	DecisionStump model = getDecisionStump();
	model.setDebug(debug);
	Evaluation trainingEvaluation = new Evaluation(trainingData);
	Evaluation testingEvaluation = new Evaluation(testingData);
	int boostingRounds = 0;
	if (boostMethod.equals("Ada")) {
	    System.err.println("# testing/training AdaBoosted DecisionStump");
	    AdaBoostM1 boost = new AdaBoostM1();
	    boost.setDebug(debug);
	    boost.setClassifier(model);
	    boost.setUseResampling(false);
	    boost.setSeed(ThreadLocalRandom.current().nextInt());
	    boost.initializeClassifier(trainingData);
	    boostingRounds = Util.runBoosting(boost);              // boost
	    trainingEvaluation.evaluateModel(boost, trainingData); // validate
	    testingEvaluation.evaluateModel(boost, testingData);   // test
	} else if (boostMethod.equals("Logit")) {
	    System.err.println("# testing/training LogitBoosted DecisionStump");
	    LogitBoost boost = new LogitBoost();
	    boost.setDebug(debug);
	    boost.setClassifier(model);
	    boost.setUseResampling(false);
	    boost.setSeed(ThreadLocalRandom.current().nextInt());
	    boost.initializeClassifier(trainingData);
	    boostingRounds = Util.runBoosting(boost);              // boost
	    trainingEvaluation.evaluateModel(boost, trainingData); // validate
	    testingEvaluation.evaluateModel(boost, testingData);   // test
	} else {
	    System.err.println("# testing/training (unboosted) DecisionStump");
	    model.buildClassifier(trainingData);	                   // train
	    trainingEvaluation.evaluateModel(model, trainingData);	   // validate
	    testingEvaluation.evaluateModel(model, testingData);	   // test
	}
	// output
	System.err.println("# training result on "+trainingData.size()+" instances after "+boostingRounds+" boosting rounds");
	Util.printResults(trainingEvaluation);
	System.err.println("# testing result on "+testingData.size()+" instances after "+boostingRounds+" boosting rounds");
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
    public static void optimize(Instances trainingData, Instances testingData, String boostMethod, int numOptimize, boolean debug) throws Exception {
	TreeMap<Double,Classifier> models = new TreeMap<>();
	TreeMap<Double,Evaluation> evaluations = new TreeMap<>();
	if (boostMethod.equals("Ada")) {
	    System.out.println("# optimizing AdaBoosted DecisionStump with numOptimize="+numOptimize);
	} else if (boostMethod.equals("Logit")) {
	    System.out.println("# optimizing LogitBoosted DecisionStump with numOptimize="+numOptimize);
	} else {
	    System.out.println("# optimizing (unboosted) DecisionStump with numOptimize="+numOptimize);
	}
	for (int round=1; round<=numOptimize; round++) {
	    DecisionStump model = getDecisionStump();
	    model.setDebug(debug);
	    Evaluation testingEvaluation = new Evaluation(testingData);
	    if (boostMethod.equals("Ada")) {
		AdaBoostM1 boost = new AdaBoostM1();
		boost.setDebug(debug);
		boost.setClassifier(model);
		boost.setUseResampling(false);
		boost.setSeed(ThreadLocalRandom.current().nextInt());
		boost.initializeClassifier(trainingData);
		int boostingRounds = Util.runBoosting(boost);         // boost
		testingEvaluation.evaluateModel(boost, testingData);  // test
		System.err.println("# "+boostingRounds+" AdaBoosting rounds");
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
		int boostingRounds = Util.runBoosting(boost);         // boost
		testingEvaluation.evaluateModel(boost, testingData);  // test
		System.err.println("# "+boostingRounds+" LogitBoosting rounds");
		Util.printResults(testingEvaluation);                 // output
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
	System.out.println("# saving best model to file wekads.model");
	SerializationHelper.write("wekads.model", models.lastEntry().getValue());
    }

    /**
     * Create a DecisionStump classifier with standard parameters.
     */
    public static DecisionStump getDecisionStump() {
	DecisionStump model = new DecisionStump();
	model.setNumDecimalPlaces(3);
	return model;
    }
}
