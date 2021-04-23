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
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ThreadLocalRandom;

import weka.classifiers.Classifier;
import weka.classifiers.Evaluation;
import weka.classifiers.functions.SMO;
import weka.classifiers.functions.supportVector.Kernel;
import weka.classifiers.functions.supportVector.PolyKernel;
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
 * Run Weka SMO operations on an ARFF file.
 */
public class WekaSMO {
    // grid search complexity range and step
    static final int C_POWER_BEGIN = 5;
    static final int C_POWER_END = 11;
    static final int C_POWER_STEP = 1;

    // RBF kernel gamma range and step
    static final int G_POWER_BEGIN = -30;
    static final int G_POWER_END = -10;
    static final int G_POWER_STEP = 1;

    // Poly kernel exponent range and step
    static final double E_BEGIN = -2.0;
    static final double E_END = +2.0;
    static final double E_STEP = 0.5;
    
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
	Option crossValidateOption = new Option("cv", "crossvalidate", false, "run cross-validation");
	crossValidateOption.setRequired(false);
	options.addOption(crossValidateOption);
	Option testOption = new Option("test", "test", false, "run training/testing on datasets with given SMO and kernel parameters");
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
	Option kernelParamOption = new Option("kp", "kernelparam", true, "SMO kernel parameter (RBF:gamma, Poly:exponent)");
	kernelParamOption.setRequired(false);
	options.addOption(kernelParamOption);
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
	double kernelParam = 0.0;
	int numOptimize = 10;
	String kernelName = cmd.getOptionValue("kernel");
	if (cmd.hasOption("complexity")) C = Double.parseDouble(cmd.getOptionValue("C"));
	if (cmd.hasOption("kernelparam")) kernelParam = Double.parseDouble(cmd.getOptionValue("kernelparam"));
	if (cmd.hasOption("kfold")) kfold = Integer.parseInt(cmd.getOptionValue("kfold"));
	if (cmd.hasOption("numoptimize")) numOptimize = Integer.parseInt(cmd.getOptionValue("numoptimize"));
	boolean debug = cmd.hasOption("debug");

	// Read instances from the input ARFF files
 	Instances data = Util.rearrange(new DataSource(cmd.getOptionValue("arfffile")).getDataSet());

	// do work
	if (crossValidate) {
	    double mcc = crossValidate(data, kernelName, kfold, C, kernelParam, debug);
	} else if (gridSearch) {
	    int nCases = 0;
	    int nControls = 0;
	    if (cmd.hasOption("ngridsearch")) {
		nCases = Integer.parseInt(cmd.getOptionValue("ngridsearch"))/2;
		nControls = nCases;
	    }
	    Instances searchData = Util.reduceData(data, nCases, nControls);
	    searchGrid(searchData, kernelName, kfold);
	} else if (test) {
	    Instances testData = Util.rearrange(new DataSource(cmd.getOptionValue("testfile")).getDataSet());
	    test(data, testData, kernelName, C, kernelParam);
	} else if (optimize) {
	    Instances testData = Util.rearrange(new DataSource(cmd.getOptionValue("testfile")).getDataSet());
	    optimize(data, testData, kernelName, C, kernelParam, numOptimize);
	}
    }

    /**
     * Run a cross-validation on the data for a given kfold, C and kernel parameter, returning the MCC.
     */
    public static double crossValidate(Instances data, String kernelName, int kfold, double C, double kernelParam, boolean debug) throws Exception {
	SMO smo = getSMO(kernelName, C, kernelParam);
	Evaluation evaluation = new Evaluation(data);
	if (debug) System.err.println("# cross-validating SMO with kfold="+kfold+" C="+sci.format(C)+" exp/gamma="+sci.format(kernelParam));
	evaluation.crossValidateModel(smo, data, kfold, new Random(ThreadLocalRandom.current().nextInt()));
	// output
	Util.printResults(evaluation, C, kernelParam);
	return evaluation.matthewsCorrelationCoefficient(1);
    }

    /**
     * Test a model on the given data.
     */
    public static void test(Instances trainingData, Instances testingData, String kernelName, double C, double kernelParam) throws Exception {
	SMO model = getSMO(kernelName, C, kernelParam);
	Evaluation trainingEvaluation = new Evaluation(trainingData);
	Evaluation testingEvaluation = new Evaluation(testingData);
	System.err.println("# testing/training SMO classifier with C="+sci.format(C)+" kernelParam="+sci.format(kernelParam));
	model.buildClassifier(trainingData);	                  // train
	trainingEvaluation.evaluateModel(model, trainingData);  // validate
	testingEvaluation.evaluateModel(model, testingData);	  // test
	// output
	System.err.println("# training result on "+trainingData.size()+" instances:");
	Util.printResults(trainingEvaluation, C, kernelParam);
	System.err.println("# testing result on "+testingData.size()+" instances:"); 
	Util.printResults(testingEvaluation, C, kernelParam);
    }

    /**
     * Run a k-fold cross-validation grid search over C and kernel param. Defines "best" as largest MCC.
     */
    public static void searchGrid(Instances data, String kernelName, int kfold) {
	System.err.println("# searching (C,param) grid with "+kfold+"-fold cross-validation on "+data.size()+" instances");
	// parameter search grid
	Set<Double> CSet = ConcurrentHashMap.newKeySet();
	for (int n=C_POWER_BEGIN; n<=C_POWER_END; n+=C_POWER_STEP) {
	    CSet.add(Math.pow(2.0, n));
	}
	Set<Double> kernelParamSet = ConcurrentHashMap.newKeySet();
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
	Map<String,Double> mccMap = new ConcurrentHashMap<>();
	////////////////////////////////////////////////////////////////////////////////////////////////
	// C loop
	CSet.parallelStream().forEach(C -> {
		////////////////////////////////////////////////////////////////////////////////////////////////
		// kernel parameter loop
		kernelParamSet.parallelStream().forEach(kernelParam -> {
			final String key = C+":"+kernelParam;
			try {
			    double mcc = crossValidate(data, kernelName, kfold, C, kernelParam, false);
			    mccMap.put(key, mcc);
			} catch (Exception ex) {
			    System.err.println(ex);
			    System.exit(1);
			}
		    });
		////////////////////////////////////////////////////////////////////////////////////////////////
	    });            
	////////////////////////////////////////////////////////////////////////////////////////////////
	// sort the mcc map
	List<Entry<String,Double>> mccList = new LinkedList<>(mccMap.entrySet());
	mccList.sort(Entry.comparingByValue());
	Map<String,Double> sortedMCCMap = new LinkedHashMap<>();
	for (Entry<String,Double> entry : mccList) {
	    String key = entry.getKey();
	    double mcc = entry.getValue();
	    sortedMCCMap.put(key, mcc);
	}
	// output
	System.out.println("C\tparam\t\tMCC");
	for (String key : sortedMCCMap.keySet()) {
	    String[] fields = key.split(":");
	    double C = Double.parseDouble(fields[0]);
	    double param = Double.parseDouble(fields[1]);
	    double mcc = sortedMCCMap.get(key);
	    System.out.println(sci.format(C)+"\t"+sci.format(param)+"\t"+rate.format(mcc));
	}
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
    public static void optimize(Instances trainingData, Instances testingData, String kernelName, double C, double kernelParam, int numOptimize) throws Exception {
	TreeMap<Double,Classifier> models = new TreeMap<>();
	TreeMap<Double,Evaluation> evaluations = new TreeMap<>();
	System.err.println("# optimizing SMO with "+kernelName+" kernel C="+sci.format(C)+" exp/gamma="+sci.format(kernelParam));
    	for (int round=1; round<=numOptimize; round++) {
	    SMO model = getSMO(kernelName, C, kernelParam);
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
     * Instantiate an SMO classifier with the given parameters.
     */
    public static SMO getSMO(String kernelName, double C, double kernelParam) {
	SMO model = new SMO();
	int modelseed = ThreadLocalRandom.current().nextInt();
	model.setRandomSeed(modelseed);
	model.setNumDecimalPlaces(3);
	model.setC(C);
	if (kernelName.equals("RBF")) {
	    RBFKernel kernel = new RBFKernel();
	    kernel.setGamma(kernelParam);
	    model.setKernel(kernel);
	} else if (kernelName.equals("Poly")) {
	    PolyKernel kernel = new PolyKernel();
	    kernel.setExponent(kernelParam);
	    model.setKernel(kernel);
	} else {
	    System.err.println("ERROR: only RBF and Poly kernels are supported at this time.");
	    System.exit(1);
	}
	return model;
    }
}
