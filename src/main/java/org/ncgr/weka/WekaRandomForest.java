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
import weka.classifiers.trees.RandomForest;
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
 * Run Weka random forest operations on an ARFF file.
 */
public class WekaRandomForest {
    // grid search complexity range and step
    static final int C_POWER_BEGIN = -5;
    static final int C_POWER_END = 15;
    static final int C_POWER_STEP = 1;
    
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
	Option testOption = new Option("test", "test", false, "run training/testing on datasets");
	testOption.setRequired(false);
	options.addOption(testOption);
	// PARAMETERS
	Option kfoldOption = new Option("k", "kfold", true, "cross-validation k-fold (required)");
	kfoldOption.setRequired(false);
	options.addOption(kfoldOption);
	//
	Option iterationsOption = new Option("I", "iterations", true, "some parameter to be decided");
	iterationsOption.setRequired(false);
	options.addOption(iterationsOption);
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
            formatter.printHelp("WekaRandomForest", options);
            System.exit(1);
            return;
        }

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

	// do work
	int iterations = Integer.parseInt(cmd.getOptionValue("iterations"));
	if (crossValidate) {
	    double correct = crossValidate(data, kfold, iterations);
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
	    test(data, testData, iterations);
	}
    }

    /**
     * Test a model on the given data.
     */
    public static void test(Instances trainingData, Instances testingData, int iterations) throws Exception {
	// RandomForest classifier
	RandomForest model = new RandomForest();
	model.setNumDecimalPlaces(3);
	model.setNumExecutionSlots(0);
	model.setNumIterations(iterations);
	// training
	model.buildClassifier(trainingData);
	Evaluation trainingEvaluation = new Evaluation(trainingData);
	System.err.println("Training RandomForest classifier on "+trainingData.size()+" instances with "+iterations+" iterations");
	trainingEvaluation.evaluateModel(model, trainingData);
	printResults(trainingEvaluation);
	// test
	System.err.println("Evaluating RandomForest classifier on "+testingData.size()+" instances.");
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
     * Run a grid search over C
     */
    public static void searchGrid(Instances data, int kfold) {
	// // parameter search grid
	// ConcurrentSkipListSet<Double> CSet = new ConcurrentSkipListSet<>();
	// for (int n=C_POWER_BEGIN; n<=C_POWER_END; n+=C_POWER_STEP) {
	//     CSet.add(Math.pow(2.0, n));
	// }
	// // these maps are keyed by C
	// ConcurrentSkipListMap<String,Double> totalCorrectMap = new ConcurrentSkipListMap<>();
	// ////////////////////////////////////////////////////////////////////////////////////////////////
	// // C loop
	// CSet.parallelStream().forEach(C -> {
	// 	final String key = String.valueOf(C);
	// 	try {
	// 	    final double correct = crossValidate(data, kfold, C);
	// 	    totalCorrectMap.put(key, correct);
	// 	} catch (Exception ex) {
	// 	    System.err.println(ex);
	// 	    System.exit(1);
	// 	}
	//     });
	// ////////////////////////////////////////////////////////////////////////////////////////////////
	// // sort the totalCorrect map
	// List<Entry<String,Double>> totalCorrectList = new LinkedList<>(totalCorrectMap.entrySet());
	// totalCorrectList.sort(Entry.comparingByValue());
	// Map<String,Double> sortedCorrectMap = new LinkedHashMap<>();
	// for (Entry<String,Double> entry : totalCorrectList) {
	//     String key = entry.getKey();
	//     double totalCorrect = entry.getValue();
	//     sortedCorrectMap.put(key, totalCorrect);
	// }
	// // output
	// System.out.println("C\tparam\tcorrect\tperc");
	// for (String key : sortedCorrectMap.keySet()) {
	//     String[] fields = key.split(":");
	//     double C = Double.parseDouble(fields[0]);
	//     double param = Double.parseDouble(fields[1]);
	//     int totalCorrect = sortedCorrectMap.get(key).intValue();
	//     System.out.println(sci.format(C)+"\t"+sci.format(param)+"\t"+totalCorrect+"\t"+perc.format((double)totalCorrect/(double)data.size()));
	// }
    }
    
    /**
     * Run a cross-validation on the data for a given kfold, and C.
     */
    public static double crossValidate(Instances data, int kfold, int iterations) throws Exception {
	// run the cross-validation
	Evaluation evaluation = new Evaluation(data);
	RandomForest model = new RandomForest();
	model.setNumDecimalPlaces(3);
	model.setNumExecutionSlots(0);
	// iterations
	model.setNumIterations(iterations);
	// cross-validate
	evaluation.crossValidateModel(model, data, kfold, new Random(ThreadLocalRandom.current().nextInt(0, 100)));
	// DEBUG
	System.err.println("kfold="+kfold+" iterations="+iterations);
	System.err.println(evaluation.toSummaryString());
	//
	printResults(evaluation);
	return  evaluation.correct();
    }
}
	
	
