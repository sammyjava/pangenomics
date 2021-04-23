package org.ncgr.libsvm;

import java.io.FileNotFoundException;
import java.io.IOException;

import java.text.DecimalFormat;

import java.util.List;
import java.util.LinkedList;
import java.util.Map;
import java.util.Map.Entry;
import java.util.LinkedHashMap;
import java.util.Set;
import java.util.Vector;
import java.util.concurrent.ConcurrentHashMap;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import libsvm.svm_node;
import libsvm.svm_parameter;

/**
 * Java port of the libsvm grid.py utility for searching a (C,gamma) grid to find optimal values.
 *
 * @author Sam Hokin
 */
public class GridSearcher {

    static DecimalFormat sci = new DecimalFormat("0.00E0");
    static DecimalFormat perc = new DecimalFormat("00.0%");
    static DecimalFormat rate = new DecimalFormat("0.000");

    boolean verbose = false;

    // SVM kernel default
    static final String KERNEL_STRING = "RBF";
    static final int KERNEL_TYPE = svm_parameter.RBF;

    // cross-validation default
    static final int NRFOLD = 10;

    // grid searching defaults
    static final int C_POWER_BEGIN = -5;
    static final int C_POWER_END = 15;
    static final int C_POWER_STEP = 1;
    static final int GAMMA_POWER_BEGIN = -25;
    static final int GAMMA_POWER_END = 0;
    static final int GAMMA_POWER_STEP = 1;

    // kernel params
    String kernelString = KERNEL_STRING;
    int kernel_type = KERNEL_TYPE;

    // grid search params
    int CPowerBegin = C_POWER_BEGIN;
    int CPowerEnd = C_POWER_END;
    int CPowerStep = C_POWER_STEP;
    int gammaPowerBegin = GAMMA_POWER_BEGIN;
    int gammaPowerEnd = GAMMA_POWER_END;
    int gammaPowerStep = GAMMA_POWER_STEP;

    // cross-validation fold
    int nrFold = NRFOLD;

    // quantities that only depend on samples
    Vector<Double> vy = new Vector<>();
    Vector<svm_node[]> vx = new Vector<>();

    // input data
    String inputFilename;
    public List<Sample> samples;
    public Vector<String> sampleNames;

    /**
     * Load the samples from a data file and set default parameters.
     *
     * @param inputFilename the data file in SVM format
     * @param nCases the number of randomly-chosen cases to use in the search (0=all)
     * @param nControls the number of randomly-chosen controls to use in the search (0=all)
     */
    public GridSearcher(String inputFilename, int nCases, int nControls) throws FileNotFoundException, IOException {
	// load all the samples
	this.inputFilename = inputFilename;
	if (nCases==0 && nControls==0) {
	    samples = Util.readSamples(inputFilename);
	} else {
	    samples = Util.reduceSamples(Util.readSamples(inputFilename), nCases, nControls);
	}
	// build the SVMLIB objects
	sampleNames = new Vector<>();
        int maxIndex = 0;
        for (Sample sample : samples) {
            sampleNames.addElement(sample.name);
            double dlabel = 0;
            if (Util.isCase(sample)) {
                dlabel = 1.0;
            } else if (Util.isControl(sample)) {
                dlabel = -1.0;
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
        // bail if we've got nothing
        if (maxIndex==0) {
	    System.err.println("GridSearcher ERROR: maxIndex==0, exiting.");
	    System.exit(1);
	}
	// bizarre setting of static function
	Util.setQuiet();
    }

    /**
     * Run the scan over (C,gamma), then sort on MCC.
     */
    public void run() {
	if (verbose) {
	    System.err.println(inputFilename);
	    System.err.println(samples.size()+" samples");
            System.err.println(nrFold+"-fold cross-validation");
	    System.err.println(kernelString+" kernel");
	    System.err.println("C=2^n from "+CPowerBegin+"("+Math.pow(2.0,CPowerBegin)+") to "+CPowerEnd+"("+Math.pow(2.0,CPowerEnd)+") with n-step="+CPowerStep);
	    if (kernel_type==svm_parameter.RBF) {
		System.err.println("gamma=2^n from "+gammaPowerBegin+"("+Math.pow(2.0,gammaPowerBegin)+") to "+gammaPowerEnd+"("+Math.pow(2.0,gammaPowerEnd)+") with n-step="+gammaPowerStep);
	    }
	}
	Set<Double> CSet = ConcurrentHashMap.newKeySet();
        for (int n=CPowerBegin; n<=CPowerEnd; n+=CPowerStep) {
            CSet.add(Math.pow(2.0, n));
        }
	Set<Double> gammaSet = ConcurrentHashMap.newKeySet();
	for (int n=gammaPowerBegin; n<=gammaPowerEnd; n+=gammaPowerStep) {
            gammaSet.add(Math.pow(2.0, n));
	}
	// these maps are keyed by the svm_parameter
	Map<svm_parameter,Double> accuracyMap = new ConcurrentHashMap<>();
	Map<svm_parameter,Double> tprMap = new ConcurrentHashMap<>();
	Map<svm_parameter,Double> fprMap = new ConcurrentHashMap<>();
	Map<svm_parameter,Double> mccMap = new ConcurrentHashMap<>();
        ////////////////////////////////////////////////////////////////////////////////////////////////
        // C loop
        CSet.parallelStream().forEach(C -> {
		if (kernel_type==svm_parameter.RBF) {
		    ////////////////////////////////////////////////////////////////////////////////////////////////
		    // gamma loop
		    gammaSet.parallelStream().forEach(gamma -> {
			    svm_parameter param = Util.getDefaultParam();
			    param.svm_type = svm_parameter.C_SVC;   // C_SVC, NU_SVC, ONE_CLASS, EPSILON_SVR, NU_SVR
			    param.C = C;
			    param.kernel_type = kernel_type;        // LINEAR, POLY, RBF, SIGMOID, PRECOMPUTED
			    param.gamma = gamma;
			    CrossValidator svc = new CrossValidator(param, nrFold, vy, vx);
			    svc.samples = samples;
			    svc.inputFilename = inputFilename;
			    svc.run();
			    if (verbose) {
				System.err.println(kernelString+": C="+sci.format(param.C)+" gamma="+sci.format(param.gamma)+
						   " accuracy="+perc.format(svc.accuracy)+
						   " TPR="+rate.format(svc.TPR)+
						   " FPR="+rate.format(svc.FPR)+
						   " MCC="+rate.format(svc.MCC));
			    }
			    // store the results
			    accuracyMap.put(param, svc.accuracy);
			    tprMap.put(param, svc.TPR);
			    fprMap.put(param, svc.FPR);
			    mccMap.put(param, svc.MCC);
			});
		    ////////////////////////////////////////////////////////////////////////////////////////////////
		} else {
		    svm_parameter param = Util.getDefaultParam();
		    param.svm_type = svm_parameter.C_SVC;   // C_SVC, NU_SVC, ONE_CLASS, EPSILON_SVR, NU_SVR
		    param.C = C;
		    param.kernel_type = kernel_type;        // LINEAR, POLY, RBF, SIGMOID, PRECOMPUTED
		    CrossValidator svc = new CrossValidator(param, nrFold, vy, vx);
		    svc.samples = samples;
		    svc.inputFilename = inputFilename;
		    svc.run();
		    if (verbose) {
			System.err.println(kernelString+": C="+sci.format(param.C)+
					   " accuracy="+perc.format(svc.accuracy)+
					   " TPR="+rate.format(svc.TPR)+
					   " FPR="+rate.format(svc.FPR)+
					   " MCC="+rate.format(svc.MCC));
		    }
		    // store the results
		    accuracyMap.put(param, svc.accuracy);
		    tprMap.put(param, svc.TPR);
		    fprMap.put(param, svc.FPR);
		    mccMap.put(param, svc.MCC);
		}
            });            
        ////////////////////////////////////////////////////////////////////////////////////////////////
	// sort the results by MCC
        List<Entry<svm_parameter,Double>> mccList = new LinkedList<>(mccMap.entrySet());
        mccList.sort(Entry.comparingByValue());
	Map<svm_parameter,Double> sortedAccuracyMap = new LinkedHashMap<>();
	Map<svm_parameter,Double> sortedTPRMap = new LinkedHashMap<>();
	Map<svm_parameter,Double> sortedFPRMap = new LinkedHashMap<>();
	Map<svm_parameter,Double> sortedMCCMap = new LinkedHashMap<>();
	for (Entry<svm_parameter,Double> entry : mccList) {
	    svm_parameter param = entry.getKey();
	    sortedMCCMap.put(param, entry.getValue());
	    sortedAccuracyMap.put(param, accuracyMap.get(param));
	    sortedTPRMap.put(param, tprMap.get(param));
	    sortedFPRMap.put(param, fprMap.get(param));
	}
	// output
	if (kernel_type==svm_parameter.RBF) {
	    System.out.println("C\tgamma\taccur\tTPR\tFPR\tMCC");
	} else {
	    System.out.println("C\taccur\tTPR\tFPR\tMCC");
	}
	for (svm_parameter param : sortedMCCMap.keySet()) {
	    if (kernel_type==svm_parameter.RBF) {
		System.out.println(sci.format(param.C)+"\t"+sci.format(param.gamma)+
				   "\t"+perc.format(sortedAccuracyMap.get(param))+
				   "\t"+perc.format(sortedTPRMap.get(param))+
				   "\t"+perc.format(sortedFPRMap.get(param))+
				   "\t"+rate.format(sortedMCCMap.get(param)));
	    } else {
		System.out.println(sci.format(param.C)+
				   "\t"+perc.format(sortedAccuracyMap.get(param))+
				   "\t"+perc.format(sortedTPRMap.get(param))+
				   "\t"+perc.format(sortedFPRMap.get(param))+
				   "\t"+rate.format(sortedMCCMap.get(param)));
	    }
	}
    }

    /**
     * Set verbose flag to true.
     */
    public void setVerbose() {
        verbose = true;
    }

    /**
     * Set the kernel based on the string representation; exit(1) if incorrect string.
     */
    public void setKernel(String kernelString) {
	if (kernelString.equals("LINEAR")) {
	    kernel_type = svm_parameter.LINEAR;
	} else if (kernelString.equals("POLY")) {
	    kernel_type = svm_parameter.POLY;
	} else if (kernelString.equals("RBF")) {
	    kernel_type = svm_parameter.RBF;
	} else if (kernelString.equals("SIGMOID")) {
	    kernel_type = svm_parameter.SIGMOID;
	} else if (kernelString.equals("PRECOMPUTED")) {
	    kernel_type = svm_parameter.PRECOMPUTED;
	} else {
	    System.err.println("ERROR: the accepted values for kernel are: LINEAR, POLY, RBF, SIGMOID, or PRECOMPUTED");
	    System.exit(1);
	}
	this.kernelString = kernelString;
    }

    /**
     * Command-line operation.
     */
    public static void main(String[] args) throws IOException {
        Options options = new Options();
        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd = null;

	Option dataFileOption = new Option("d", "datafile", true, "input data file in SVM format [required]");
	dataFileOption.setRequired(true);
	options.addOption(dataFileOption);
	//
	Option kernelOption = new Option("k", "kernel", true, "choose the SVM kernel: LINEAR, POLY, RBF, SIGMOID, PRECOMPUTED ["+KERNEL_STRING+"]");
	kernelOption.setRequired(false);
	options.addOption(kernelOption);
	//
        Option CPowerBeginOption = new Option("cb", "cpowerbegin", true, "begining 2^n exponent in C scan ["+C_POWER_BEGIN+"]");
        CPowerBeginOption.setRequired(false);
        options.addOption(CPowerBeginOption);
        //
        Option CPowerEndOption = new Option("ce", "cpowerend", true, "ending 2^n exponent in C scan ["+C_POWER_END+"]");
        CPowerEndOption.setRequired(false);
        options.addOption(CPowerEndOption);
        //
        Option CPowerStepOption = new Option("cs", "cpowerstep", true, "step for 2^n exponent in C scan ["+C_POWER_STEP+"]");
        CPowerStepOption.setRequired(false);
        options.addOption(CPowerStepOption);
	//
        Option gammaPowerBeginOption = new Option("gb", "gammapowerbegin", true, "begining 2^n exponent in RBF gamma scan ["+GAMMA_POWER_BEGIN+"]");
        gammaPowerBeginOption.setRequired(false);
        options.addOption(gammaPowerBeginOption);
        //
        Option gammaPowerEndOption = new Option("ge", "gammapowerend", true, "ending 2^n exponent in RBF gamma scan ["+GAMMA_POWER_END+"]");
        gammaPowerEndOption.setRequired(false);
        options.addOption(gammaPowerEndOption);
        //
        Option gammaPowerStepOption = new Option("gs", "gammapowerstep", true, "step for 2^n exponent in RBF gamma scan ["+GAMMA_POWER_STEP+"]");
        gammaPowerStepOption.setRequired(false);
        options.addOption(gammaPowerStepOption);
	//
        Option nFoldOption = new Option("nr", "nrfold", true, "nr-fold for cross validation ["+NRFOLD+"]");
        nFoldOption.setRequired(false);
        options.addOption(nFoldOption);
	//
        Option vOption = new Option("v", "verbose", false, "toggle verbose output");
        vOption.setRequired(false);
        options.addOption(vOption);
	//
	Option nCasesOption = new Option("ncases", "ncases", true, "set number of cases to use in search [0=all]");
	nCasesOption.setRequired(false);
	options.addOption(nCasesOption);
	//
	Option nControlsOption = new Option("ncontrols", "ncontrols", true, "set number of controls to use in search [0=all]");
	nControlsOption.setRequired(false);
	options.addOption(nControlsOption);

        if (args.length<2) {
            formatter.printHelp("GridSearcher [options]", options);
            System.exit(1);
        }
	
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println(e.getMessage());
            formatter.printHelp("GridSearcher", options);
            System.exit(1);
        }

        String datafile = cmd.getOptionValue("datafile");

	// number of cases/controls to run
	int nCases = 0;
	int nControls = 0;
	if (cmd.hasOption("ncases")) nCases = Integer.parseInt(cmd.getOptionValue("ncases"));
	if (cmd.hasOption("ncontrols")) nControls = Integer.parseInt(cmd.getOptionValue("ncontrols"));

        // initialize
	GridSearcher gs = new GridSearcher(datafile, nCases, nControls);
        if (cmd.hasOption("verbose")) {
            gs.setVerbose();
        }

	// kernel type
	if (cmd.hasOption("kernel")) gs.setKernel(cmd.getOptionValue("kernel"));

	// cross-validation k-fold
	if (cmd.hasOption("nrfold")) gs.nrFold = Integer.parseInt(cmd.getOptionValue("nrfold"));
	
	// optionally set grid search parameters
        if (cmd.hasOption("cpowerbegin")) gs.CPowerBegin = Integer.parseInt(cmd.getOptionValue("cpowerbegin"));
        if (cmd.hasOption("cpowerend")) gs.CPowerEnd = Integer.parseInt(cmd.getOptionValue("cpowerend"));
        if (cmd.hasOption("cpowerstep")) gs.CPowerStep = Integer.parseInt(cmd.getOptionValue("cpowerstep"));
        if (cmd.hasOption("gammapowerbegin")) gs.gammaPowerBegin = Integer.parseInt(cmd.getOptionValue("gammapowerbegin"));
        if (cmd.hasOption("gammapowerend")) gs.gammaPowerEnd = Integer.parseInt(cmd.getOptionValue("gammapowerend"));
        if (cmd.hasOption("gammapowerstep")) gs.gammaPowerStep = Integer.parseInt(cmd.getOptionValue("gammapowerstep"));

        // run the search
        gs.run();
    }
}
