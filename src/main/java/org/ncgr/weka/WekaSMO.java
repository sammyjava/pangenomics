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
import weka.classifiers.functions.supportVector.Kernel;
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
    // grid search complexity range and step
    static final int C_POWER_BEGIN = -5;
    static final int C_POWER_END = 15;
    static final int C_POWER_STEP = 1;

    // RBF kernel gamma range and step
    static final int G_POWER_BEGIN = -25;
    static final int G_POWER_END = 0;
    static final int G_POWER_STEP = 1;

    // Poly kernel exponent range and step
    static final double E_BEGIN = -2.0;
    static final double E_END = +2.0;
    static final double E_STEP = 0.5;
    
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
	Option kernelOption = new Option("kernel", "kernel", true, "kernel function: RBF or Poly (required)");
	kernelOption.setRequired(true);
	options.addOption(kernelOption);
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
	Option kernelParamOption = new Option("kp", "kernelparam", true, "SMO kernel parameter (RBF:gamma, Poly:exponent)");
	kernelParamOption.setRequired(false);
	options.addOption(kernelParamOption);
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

	// actions
	final boolean crossValidate = cmd.hasOption("crossvalidate");
	final boolean gridSearch = cmd.hasOption("gridsearch");
	final boolean test = cmd.hasOption("test");

	// cross-validation and grid search kfold
	int kfold = 0;
	if (crossValidate || gridSearch) {
	    kfold = Integer.parseInt(cmd.getOptionValue("kfold"));
	}

	// the kernel and parameters
	double C = 0.0;
	double kernelParam = 0.0;
	String kernelName = cmd.getOptionValue("kernel");
	if (kernelName.equals("RBF")) {
	    if (test || crossValidate) {
		C = Double.parseDouble(cmd.getOptionValue("C"));
		kernelParam = Double.parseDouble(cmd.getOptionValue("kernelparam"));
	    }
	} else if (kernelName.equals("Poly")) {
	    if (test || crossValidate) {
		C = Double.parseDouble(cmd.getOptionValue("C"));
		kernelParam = Double.parseDouble(cmd.getOptionValue("kernelparam"));
	    }
	} else {
	    System.err.println("ERROR: you must specify either RBF or Poly kernel");
	    System.exit(1);
	}

	// Read instances from the input ARFF file
 	Instances data = new DataSource(cmd.getOptionValue("arfffile")).getDataSet();
        // remove the ID attribute
        data.deleteAttributeAt(0);
        // set the class attribute index at the end (0/1)
        data.setClassIndex(data.numAttributes() - 1);

	// do work
	if (crossValidate) {
	    double correct = crossValidate(data, kernelName, kfold, C, kernelParam);
	} else if (gridSearch) {
	    int nCases = 0;
	    int nControls = 0;
	    if (cmd.hasOption("ngridsearch")) {
		nCases = Integer.parseInt(cmd.getOptionValue("ngridsearch"))/2;
		nControls = nCases;
	    }
	    Instances searchData = reduceData(data, nCases, nControls);
	    searchGrid(searchData, kernelName, kfold);
	} else if (test) {
	    Instances testData = new DataSource(cmd.getOptionValue("testfile")).getDataSet();
	    // remove the ID attribute
	    testData.deleteAttributeAt(0);
	    // set the class attribute index at the end (0/1)
	    testData.setClassIndex(data.numAttributes() - 1);
	    // test
	    test(data, testData, kernelName, C, kernelParam);
	}
    }

    /**
     * Test a model on the given data.
     */
    public static void test(Instances trainingData, Instances testingData, String kernelName, double C, double kernelParam) throws Exception {
	// SMO model
	SMO model = new SMO();
	model.setNumDecimalPlaces(3);
	model.setC(C);
	if (kernelName.equals("RBF")) {
	    model.setKernel(getRBFKernel(kernelParam));
	} else if (kernelName.equals("Poly")) {
	    model.setKernel(getPolyKernel(kernelParam));
	} else {
	    throw new Exception("ERROR: only RBF and Poly kernels are supported at this time.");
	}
	// training
	model.buildClassifier(trainingData);
	Evaluation trainingEvaluation = new Evaluation(trainingData);
	System.err.println("Training SMO classifier on "+trainingData.size()+" instances with C="+sci.format(C)+" and kernelParam="+sci.format(kernelParam));
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
     * Run a grid search over C and kernel param.
     */
    public static void searchGrid(Instances data, String kernelName, int kfold) {
	// parameter search grid
	ConcurrentSkipListSet<Double> CSet = new ConcurrentSkipListSet<>();
	for (int n=C_POWER_BEGIN; n<=C_POWER_END; n+=C_POWER_STEP) {
	    CSet.add(Math.pow(2.0, n));
	}
	ConcurrentSkipListSet<Double> kernelParamSet = new ConcurrentSkipListSet<>();
	if (kernelName.equals("RBF")) {
	    for (int n=G_POWER_BEGIN; n<=G_POWER_END; n+=G_POWER_STEP) {
		kernelParamSet.add(Math.pow(2.0, n));
	    }
	} else  if (kernelName.equals("Poly")) {
	    for (double e=E_BEGIN; e<=E_END; e+=E_STEP) {
		kernelParamSet.add(e);
	    }
	}
	// these maps are keyed by C:kernelParam
	ConcurrentSkipListMap<String,Double> totalCorrectMap = new ConcurrentSkipListMap<>();
	////////////////////////////////////////////////////////////////////////////////////////////////
	// C loop
	CSet.parallelStream().forEach(C -> {
		////////////////////////////////////////////////////////////////////////////////////////////////
		// kernel parameter loop
		kernelParamSet.parallelStream().forEach(kernelParam -> {
			final String key = C+":"+kernelParam;
			try {
			    final double correct = crossValidate(data, kernelName, kfold, C, kernelParam);
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
	System.out.println("C\tparam\tcorrect\tperc");
	for (String key : sortedCorrectMap.keySet()) {
	    String[] fields = key.split(":");
	    double C = Double.parseDouble(fields[0]);
	    double param = Double.parseDouble(fields[1]);
	    int totalCorrect = sortedCorrectMap.get(key).intValue();
	    System.out.println(sci.format(C)+"\t"+sci.format(param)+"\t"+totalCorrect+"\t"+perc.format((double)totalCorrect/(double)data.size()));
	}
    }
    
    /**
     * Run a cross-validation on the data for a given kfold, C and kernel parameter.
     */
    public static double crossValidate(Instances data, String kernelName, int kfold, double C, double kernelParam) throws Exception {
	// run the cross-validation
	Evaluation evaluation = new Evaluation(data);
	SMO model = new SMO();
	model.setNumDecimalPlaces(3);
	model.turnChecksOn();
	// complexity
	model.setC(C);
	// kernel
	if (kernelName.equals("RBF")) {
	    // kernelParam = gamma
	    model.setKernel(getRBFKernel(kernelParam));
	} else if (kernelName.equals("Poly")) {
	    // kernelParam = exponent
	    model.setKernel(getPolyKernel(kernelParam));
	} else {
	    throw new Exception("ERROR: only RBF and Poly kernels are supported at this time.");
	}
	// cross-validate
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
	System.err.println("kfold="+kfold+" C="+C+" kernelParam="+kernelParam);
	System.err.println(evaluation.toSummaryString());
	//
	return correct;
    }

    /**
     * Instantiate an RBF kernel with the given gamma value.
     */
    public static RBFKernel getRBFKernel(double gamma) {
	RBFKernel kernel = new RBFKernel();
	kernel.setGamma(gamma);
	return kernel;
    }

    /**
     * Instantiate a Poly kernel with the given exponent.
     */
    public static PolyKernel getPolyKernel(double exponent) {
	PolyKernel kernel = new PolyKernel();
	kernel.setExponent(exponent);
	return kernel;
    }
}
	
	
