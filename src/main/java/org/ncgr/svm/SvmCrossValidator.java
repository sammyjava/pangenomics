package org.ncgr.svm;

import java.io.FileNotFoundException;
import java.io.IOException;

import java.text.DecimalFormat;

import java.util.List;
import java.util.ArrayList;
import java.util.Vector;
import java.util.concurrent.ConcurrentHashMap;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import libsvm.svm;
import libsvm.svm_node;
import libsvm.svm_parameter;
import libsvm.svm_problem;

/**
 * Use svm to run a stratified k-fold cross-validation of data.
 */
public class SvmCrossValidator {
    static DecimalFormat pf = new DecimalFormat("0.0%");
    static DecimalFormat df = new DecimalFormat("0.000");
    
    svm_parameter param;
    svm_problem prob;
    int nrFold = SvmUtil.NRFOLD;

    // regression results
    public double totalError = 0.0;
    public double meanSquaredError = 0.0;
    public double squaredCorrCoeff = 0.0;

    // classification results
    public int plusFails = 0;
    public int minusFails = 0;
    public int plusTotal = 0;
    public int minusTotal = 0;
    public double percentCorrect = 0.0;
    public double percentPlusCorrect = 0.0;
    public double percentMinusCorrect = 0.0;
    public int totalSamples = 0;
    public int totalCorrect = 0;
    public double accuracy = 0.0;
    public int maxIndex = 0;
    public Vector<Integer> errorIndices;

    // the input filename
    public String inputFilename;
    
    // the input samples
    public List<Sample> samples;
    // for alignment with libsvm vectors
    public Vector<String> sampleNames;
    
    /**
     * Construct with default svm_parameter object, default-fold cross-validation and the given input file name.
     */
    public SvmCrossValidator(String inputFilename, int nCases, int nControls) throws FileNotFoundException, IOException {
	this.inputFilename = inputFilename;
        param = SvmUtil.getDefaultParam();
	if (nCases==0 && nControls==0) {
	    samples = SvmUtil.readSamples(inputFilename);
	} else {
	    samples = SvmUtil.reduceSamples(SvmUtil.readSamples(inputFilename), nCases, nControls);
	}
        createProblem();
    }

    /**
     * Construct given a populated svm_parameter object, n-fold number and an input file name.
     */
    public SvmCrossValidator(svm_parameter param, int nrFold, String inputFilename, int nCases, int nControls) throws IOException {
        this.param = param;
        this.nrFold = nrFold;
	this.inputFilename = inputFilename;
        // validate nrFold
        if (nrFold<2) {
            System.err.println("Error: n-fold cross validation requires n>=2.");
            System.exit(1);
        }
	if (nCases==0 && nControls==0) {
	    samples = SvmUtil.readSamples(inputFilename);
	} else {
	    samples = SvmUtil.reduceSamples(SvmUtil.readSamples(inputFilename), nCases, nControls);
	}
        createProblem();
    }

    /**
     * Construct given a populated svm_parameter object, n-fold number and populated labels and svm_nodes.
     */
    public SvmCrossValidator(svm_parameter param, int nrFold, Vector<Double> vy, Vector<svm_node[]> vx) {
        this.param = param;
        this.nrFold = nrFold;
        // validate nrFold
        if (nrFold<2) {
            System.err.println("Error: n-fold cross validation requires n>=2.");
            System.exit(1);
        }
        createProblem(vy, vx);
    }

    /**
     * Construct given a populated svm_problem, nr-fold number and populated svm_parameter.
     */
    public SvmCrossValidator(svm_parameter param, int nrFold, svm_problem prob) {
        this.param = param;
        this.nrFold = nrFold;
        this.prob = prob;
        // validate nrFold
        if (nrFold<2) {
            System.err.println("Error: n-fold cross validation requires n>=2.");
            System.exit(1);
        }
    }

    /**
     * Perform the cross validation.
     */
    public void run() {
        double[] target = new double[prob.l];
        double sumv = 0.0, sumy = 0.0, sumvv = 0.0, sumyy = 0.0, sumvy = 0.0;
        errorIndices = new Vector<>();
        // run it
	try {
	    svm.svm_cross_validation(prob, param, nrFold, target);
	} catch (Exception e) {
	    System.err.println(e.getMessage());
	    System.exit(1);
	}
        if (param.svm_type == svm_parameter.EPSILON_SVR || param.svm_type == svm_parameter.NU_SVR) {
            // regression results
            for (int i=0; i<prob.l; i++) {
                double y = prob.y[i];
                double v = target[i];
                totalError += (v-y)*(v-y);
                sumv += v;
                sumy += y;
                sumvv += v*v;
                sumyy += y*y;
                sumvy += v*y;
            }
            meanSquaredError = totalError/prob.l;
            squaredCorrCoeff = ((prob.l*sumvy-sumv*sumy)*(prob.l*sumvy-sumv*sumy))/((prob.l*sumvv-sumv*sumv)*(prob.l*sumyy-sumy*sumy));
        } else {
            // classification results
            for (int i=0; i<prob.l; i++) {
                if (target[i] == prob.y[i]) {
                    ++totalCorrect;
                } else {
                    errorIndices.add(i);
                }
            }
            accuracy = (double)totalCorrect/(double)prob.l;
            totalSamples = prob.l;
        }
	// calculate summary stats
	plusFails = 0;
	minusFails = 0;
	for (int i : errorIndices) {
	    if (prob.y[i]==+1) {
		plusFails++;
	    } else if (prob.y[i]==-1) {
		minusFails++;
	    }
	}
	plusTotal = 0;
	minusTotal = 0;
	for (double y : prob.y) {
	    if (y==+1) plusTotal++;
	    if (y==-1) minusTotal++;
	}
	percentCorrect = (double)(plusTotal+minusTotal-plusFails-minusFails)/(double)(plusTotal+minusTotal);
	percentPlusCorrect = (double)(plusTotal-plusFails)/(double)plusTotal;
	percentMinusCorrect = (double)(minusTotal-minusFails)/(double)minusTotal;
    }

    /**
     * Print out the summary of results to System.err.
     */
    public void printSummary() {
	if (param.svm_type==svm_parameter.EPSILON_SVR || param.svm_type== svm_parameter.NU_SVR) {
	    System.err.println("Cross Validation Mean squared error = "+meanSquaredError);
	    System.err.println("Cross Validation Squared correlation coefficient = "+squaredCorrCoeff);
	    return;
	}
	// summary line
	int total = plusTotal + minusTotal;
	int TP = plusTotal - plusFails;
	int TN = minusTotal - minusFails;
	int FP = minusFails;
	int FN = plusFails;
	int totCorrect = TP + TN;
	double fracCorrect = (double)totCorrect / (double)total;
	double fracPlusCorrect = (double)TP / (double)plusTotal;
	double fracMinusCorrect = (double)TN / (double)minusTotal;
	double TPR = (double)TP / (double)plusTotal;
	double FPR = (double)FP / (double)minusTotal;
	double MCC = (double)(TP*TN-FP*FN) / Math.sqrt((double)(TP+FP)*(double)(TP+FN)*(double)(TN+FP)*(double)(TN+FN));
	System.err.println("Total Corr.\tCase Corr.\tControl Corr.\tTotal\tCase\tControl\tTPR\tFPR\tMCC");
	System.err.println(totCorrect+"/"+total+"\t"+
			   TP+"/"+plusTotal+"\t"+
			   TN+"/"+minusTotal+"\t"+
			   pf.format(fracCorrect)+"\t"+
			   pf.format(fracPlusCorrect)+"\t"+
			   pf.format(fracMinusCorrect)+"\t"+
			   df.format(TPR)+"\t"+
			   df.format(FPR)+"\t"+
			   df.format(MCC));
	return;
    }

    /**
     * Print output to System.out for further processing.
     */
    public void printOutput() {
	// output for further processing
	System.out.println(inputFilename+"\t"+maxIndex+"\t"+param.C+"\t"+param.gamma+"\t"+nrFold+"\t"+plusTotal+"\t"+minusTotal+"\t"+plusFails+"\t"+minusFails);
	// predictions by sample
	for (int i=0; i<prob.l; i++) {
	    String label = "";
	    System.out.print(samples.get(i).name+"\t"+samples.get(i).label);
	    if (errorIndices.contains(i)) {
		System.out.println("\tfalse");
	    } else {
		System.out.println("\ttrue");
	    }
	}
    }

    /**
     * Build label vector vy and nodes vector vx from a list of samples and then create the svm_problem.
     */
    public void createProblem() {
        maxIndex = 0;
        sampleNames = new Vector<>();
        Vector<Double> vy = new Vector<>();
        Vector<svm_node[]> vx = new Vector<>();
        for (Sample sample : samples) {
            sampleNames.addElement(sample.name);
	    // cases are "plus", controls are "minus"
            double dlabel = 0;
	    if (SvmUtil.isCase(sample)) {
                dlabel = 1.0;
            } else if (SvmUtil.isControl(sample)) {
                dlabel = -1.0;
            } else {
		System.err.println("ERROR: sample "+sample+" is neither case nor control. Exiting.");
		System.exit(1);
	    }
            vy.addElement(dlabel);
            svm_node[] x = new svm_node[sample.values.size()];
            int j = 0;
            for (int index : sample.values.keySet()) {
                double value = sample.values.get(index);
                x[j] = new svm_node();
                x[j].index = index;
                x[j].value = value;
                maxIndex = Math.max(maxIndex, x[j].index);
                j++;
            }
            vx.addElement(x);
        }
        createProblem(vy, vx);
    }

    /**
     * Create and populate the svm_problem from the label vector vy and nodes vector vx.
     * Sets instance vars prob and param.
     */
    void createProblem(Vector<Double> vy, Vector<svm_node[]> vx) {
        // find the maximum index value
        int max_index = 0;
        for (svm_node[] x : vx) {
            for (svm_node node : x) {
                max_index = Math.max(max_index, node.index);
            }
        }
        // create the problem round town
        prob = new svm_problem();
        prob.l = vy.size();
        prob.x = new svm_node[prob.l][];
        for (int i=0;i<prob.l;i++) {
            prob.x[i] = vx.elementAt(i);
        }
        prob.y = new double[prob.l];
        for (int i=0;i<prob.l;i++) {
            prob.y[i] = vy.elementAt(i);
        }
        // set param.gamma = 1/N if zero
        if (param.gamma==0 && max_index>0) param.gamma = 1.0/max_index;
        // validation
        if (param.kernel_type == svm_parameter.PRECOMPUTED) {
            for (int i=0;i<prob.l;i++) {
                if (prob.x[i][0].index != 0) {
                    System.err.println("Wrong kernel matrix: first column must be 0:sample_serial_number");
                    System.exit(1);
                }
                if ((int)prob.x[i][0].value <= 0 || (int)prob.x[i][0].value > max_index) {
                    System.err.println("Wrong input format: sample_serial_number out of range");
                    System.exit(1);
                }
            }
        }
        String errorMsg = svm.svm_check_parameter(prob,param);
        if (errorMsg!=null) {
            System.err.println("ERROR: "+errorMsg);
            System.exit(1);
        }
    }

    /**
     * Command line version.
     */
    public static void main(String[] args) throws IOException {
        Options options = new Options();
        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd = null;
	//
	Option inputFileOption = new Option("i", "inputfile", true, "input file containing feature vectors in SVM format");
	inputFileOption.setRequired(true);
	options.addOption(inputFileOption);
        //
        Option svmTypeOption = new Option("s", "svm-type", true, "type of SVM: 0=C-SVC 1=nu-SVC 2=one-class SVM 3=epsilon-SVR 4=nu-SVR [0]");
        svmTypeOption.setRequired(false);
        options.addOption(svmTypeOption);
        //
        Option kernelTypeOption = new Option("t", "kernel-type", true, "type of kernel function: 0=linear 1=polynomial 2=radial basis function 3=sigmoid 4=precomputed kernel [2]");
        kernelTypeOption.setRequired(false);
        options.addOption(kernelTypeOption);
        //
        Option kernelDegreeOption = new Option("d", "kernel-degree", true, "degree parameter in kernel function [3]");
        kernelDegreeOption.setRequired(false);
        options.addOption(kernelDegreeOption);
        //
        Option kernelGammaOption = new Option("gamma", "kernel-gamma", true, "gamma parameter in kernel function [1/#features]");
        kernelGammaOption.setRequired(false);
        options.addOption(kernelGammaOption);
        //
        Option kernelCoef0Option = new Option("coef0", "kernel-coef0", true, "coef0 parameter in kernel function [0]");
        kernelCoef0Option.setRequired(false);
        options.addOption(kernelCoef0Option);
        //
        Option costOption = new Option("C", "cost", true, "cost parameter C in C-SVC, epsilon-SVR and nu-SVR [1]");
        costOption.setRequired(false);
        options.addOption(costOption);
        //
        Option nuOption = new Option("nu", "nu", true, "nu parameter of nu-SVC, one-class SVM, and nu-SVR [0.5]");
        nuOption.setRequired(false);
        options.addOption(nuOption);
        //
        Option epsilonLossOption = new Option("ep", "epsilon-loss", true, "epsilon value in loss function of epsilon-SVR [0.1]");
        epsilonLossOption.setRequired(false);
        options.addOption(epsilonLossOption);
        //
        Option cacheSizeOption = new Option("m", "cachesize", true, "set cache memory size in MB [100]");
        cacheSizeOption.setRequired(false);
        options.addOption(cacheSizeOption);
        //
        Option epsilonOption = new Option("e", "epsilon", true, "set tolerance of termination criterion [0.001]");
        epsilonOption.setRequired(false);
        options.addOption(epsilonOption);
        //
        Option shrinkingOption = new Option("h", "shrinking", true, "toggle whether to use the shrinking heuristics (0/1) [1]");
        shrinkingOption.setRequired(false);
        options.addOption(shrinkingOption);
        //
        Option probabilityEstimatesOption = new Option("prob", "probability-estimates", true, "toggle whether to train a SVC or SVR model for probability estimates (0/1) [0]");
        probabilityEstimatesOption.setRequired(false);
        options.addOption(probabilityEstimatesOption);
        //
        Option weightOption = new Option("w", "weight", true, "multiply parameter C of class i by weight, for C-SVC [1]");
        weightOption.setRequired(false);
        options.addOption(weightOption);
        //
        Option nrFoldOption = new Option("k", "nrfold", true, "k value for k-fold cross-validation ["+SvmUtil.NRFOLD+"]");
        nrFoldOption.setRequired(false);
        options.addOption(nrFoldOption);
        //
        Option verboseOption = new Option("v", "verbose", false, "verbose SVMLIB output");
        verboseOption.setRequired(false);
        options.addOption(verboseOption);
	//
	Option nCasesOption = new Option("ncases", true, "set number of cases to use in search [0=all]");
	nCasesOption.setRequired(false);
	options.addOption(nCasesOption);
	//
	Option nControlsOption = new Option("ncontrols", true, "set number of controls to use in search [0=all]");
	nControlsOption.setRequired(false);
	options.addOption(nControlsOption);
	//
	Option nRunsOption = new Option("nruns", true, "number of parallel runs [1]");
	nRunsOption.setRequired(false);
	options.addOption(nRunsOption);

        if (args.length==0) {
            formatter.printHelp("SvmCrossValidator [options]", options);
            System.exit(1);
        }
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println(e.getMessage());
            formatter.printHelp("SvmCrossValidator [options]", options);
            System.exit(1);
        }

        // start with default param values
        svm_parameter param = SvmUtil.getDefaultParam();

        // update svm_parameter values from options
        if (cmd.hasOption("svm-type")) {
            param.svm_type = Integer.parseInt(cmd.getOptionValue("svm-type"));
        }
        if (cmd.hasOption("kernel-type")) {
            param.kernel_type = Integer.parseInt(cmd.getOptionValue("kernel-type"));
        }
        if (cmd.hasOption("kernel-degree")) {
            param.degree = Integer.parseInt(cmd.getOptionValue("kernel-degree"));
        }
        if (cmd.hasOption("kernel-gamma")) {
            param.gamma = Double.parseDouble(cmd.getOptionValue("kernel-gamma"));
        }
        if (cmd.hasOption("nu")) {
            param.nu = Double.parseDouble(cmd.getOptionValue("nu"));
        }
        if (cmd.hasOption("cachesize")) {
            param.cache_size = Double.parseDouble(cmd.getOptionValue("cachesize"));
        }
        if (cmd.hasOption("cost")) {
            param.C = Double.parseDouble(cmd.getOptionValue("cost"));
        }
        if (cmd.hasOption("epsilon")) {
            param.eps = Double.parseDouble(cmd.getOptionValue("epsilon"));
        }
        if (cmd.hasOption("epsilon-loss")) {
            param.p = Double.parseDouble(cmd.getOptionValue("epsilon-loss"));
        }
        if (cmd.hasOption("shrinking")) {
            param.shrinking = Integer.parseInt(cmd.getOptionValue("shrinking"));
        }
        if (cmd.hasOption("probability-estimates")) {
            param.probability = Integer.parseInt(cmd.getOptionValue("probability-estimates"));
        }
	
        int nrFold = SvmUtil.NRFOLD;
        if (cmd.hasOption("nrfold")) {
            nrFold = Integer.parseInt(cmd.getOptionValue("nrfold"));
        }

	String inputFilename = cmd.getOptionValue("i");

	int nCases = 0;
	int nControls = 0;
	if (cmd.hasOption("ncases")) nCases = Integer.parseInt(cmd.getOptionValue("ncases"));
	if (cmd.hasOption("ncontrols")) nControls = Integer.parseInt(cmd.getOptionValue("ncontrols"));

	int nRuns = 1;
	if (cmd.hasOption("nruns")) nRuns = Integer.parseInt(cmd.getOptionValue("nruns"));

        // this is weird, setting a static function in svm
        if (cmd.hasOption("v")) {
            SvmUtil.setVerbose();
        } else {
            SvmUtil.setQuiet();
        }

        // instantiate and run nRuns times in a parallelStream
	ConcurrentHashMap<Integer,SvmCrossValidator> scvMap = new ConcurrentHashMap<>();
	for (int i=0; i<nRuns; i++) {
	    scvMap.put(i, new SvmCrossValidator(param, nrFold, inputFilename, nCases, nControls));
	}
	////////////////////////////////////////////////////////////////////////////////////////
	// parallel stream
	scvMap.entrySet().parallelStream().forEach(entry -> {
		int i = entry.getKey();
		SvmCrossValidator svm = entry.getValue();
		svm.run();
	    });
	////////////////////////////////////////////////////////////////////////////////////////
	// output
	for (int i : scvMap.keySet()) {
	    SvmCrossValidator svm = scvMap.get(i);
	    System.err.println("===== ["+i+"] =====");
	    svm.printSummary();
	    svm.printOutput();
	}
    }
}
