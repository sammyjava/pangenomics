package org.ncgr.svm;

import java.io.FileNotFoundException;
import java.io.IOException;

import java.text.DecimalFormat;

import java.util.List;
import java.util.LinkedList;
import java.util.Map;
import java.util.HashMap;
import java.util.Vector;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.concurrent.ConcurrentSkipListSet;

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

    static DecimalFormat df = new DecimalFormat("0.00E0");
    static DecimalFormat pf = new DecimalFormat("00.0%");

    boolean verbose = false;

    // grid defaults
    static final int C_POWER_BEGIN = -10;
    static final int C_POWER_END = 10;
    static final int C_POWER_STEP = 1;
    static final int G_POWER_BEGIN = -20;
    static final int G_POWER_END = 0;
    static final int G_POWER_STEP = 1;

    // grid sparsity
    int CPowerStep = C_POWER_STEP;
    int GPowerStep = G_POWER_STEP;

    // SVM kernel
    int kernel_type = svm_parameter.RBF;


    // quantities that only depend on samples
    Vector<Double> vy = new Vector<>();
    Vector<svm_node[]> vx = new Vector<>();

    // cross-validation
    int nrFold = SvmUtil.NRFOLD;

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
	    samples = SvmUtil.readSamples(inputFilename);
	} else {
	    samples = SvmUtil.reduceSamples(SvmUtil.readSamples(inputFilename), nCases, nControls);
	}
	// build the SVMLIB objects
	sampleNames = new Vector<>();
        int maxIndex = 0;
        for (Sample sample : samples) {
            sampleNames.addElement(sample.name);
            double dlabel = 0;
            if (SvmUtil.isCase(sample)) {
                dlabel = 1.0;
            } else if (SvmUtil.isControl(sample)) {
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
	SvmUtil.setQuiet();
    }

    /**
     * Run the search for each C and gamma value
     */
    public void run() {
	if (verbose) {
	    System.err.println("GridSearcher: "+inputFilename);
	    System.err.println("GridSearcher: "+samples.size()+" samples");
            System.err.println("GridSearcher: "+nrFold+"-fold cross-validation");
	    System.err.println("GridSearcher: C=2^n from n="+C_POWER_BEGIN+" to "+C_POWER_END+" in steps of "+CPowerStep);
	    System.err.println("GridSearcher: gamma=2^n from n="+G_POWER_BEGIN+" to "+G_POWER_END+" in steps of "+GPowerStep);
	}
	ConcurrentSkipListSet<Double> CSet = new ConcurrentSkipListSet<>();
        for (int n=C_POWER_BEGIN; n<=C_POWER_END; n+=CPowerStep) {
            CSet.add(Math.pow(2.0, n));
        }
	ConcurrentSkipListSet<Double> gammaSet = new ConcurrentSkipListSet<>();
	for (int n=G_POWER_BEGIN; n<=G_POWER_END; n+=GPowerStep) {
            gammaSet.add(Math.pow(2.0, n));
	}
	// these maps are keyed by a string representing C and gamma
	ConcurrentSkipListMap<String,svm_parameter> paramMap = new ConcurrentSkipListMap<>();
	ConcurrentSkipListMap<String,Integer> totalCorrectMap = new ConcurrentSkipListMap<>();
	ConcurrentSkipListMap<String,Double> accuracyMap = new ConcurrentSkipListMap<>();
        ////////////////////////////////////////////////////////////////////////////////////////////////
        // C loop
        CSet.parallelStream().forEach(C -> {
                ////////////////////////////////////////////////////////////////////////////////////////////////
                // gamma loop
                gammaSet.parallelStream().forEach(gamma -> {
                        svm_parameter param = SvmUtil.getDefaultParam();
                        param.C = C;
                        param.gamma = gamma;
			param.kernel_type = kernel_type;        // LINEAR, POLY, RBF, SIGMOID, PRECOMPUTED
			param.svm_type = svm_parameter.C_SVC;   // C_SVC, NU_SVC, ONE_CLASS, EPSILON_SVR, NU_SVR
                        String key = param.C+":"+param.gamma;
                        SvmCrossValidator svc = new SvmCrossValidator(param, nrFold, vy, vx);
                        svc.samples = samples;
                        svc.inputFilename = inputFilename;
                        svc.run();
                        if (verbose) {
                            System.err.println("GridSearcher: C="+df.format(param.C)+" gamma="+df.format(param.gamma)+" totalCorrect="+svc.totalCorrect+" accuracy="+pf.format(svc.accuracy));
                        }
                        // store the results
                        paramMap.put(key, param);
                        totalCorrectMap.put(key, svc.totalCorrect);
                        accuracyMap.put(key, svc.accuracy);
                    });
                ////////////////////////////////////////////////////////////////////////////////////////////////
            });              
        ////////////////////////////////////////////////////////////////////////////////////////////////
	// scan the results for best results
        svm_parameter bestParam = null;
	int bestTotalCorrect = 0;
	double bestAccuracy = 0.0;
	for (String key : paramMap.keySet()) {
	    svm_parameter param = paramMap.get(key);
	    int totalCorrect = totalCorrectMap.get(key);
	    double accuracy = accuracyMap.get(key);
	    if (totalCorrect>bestTotalCorrect) {
		bestParam = param;
		bestTotalCorrect = totalCorrect;
		bestAccuracy = accuracy;
	    }
	}
        // now get all parameters that match the best results
        Map<svm_parameter,Integer> bestCorrectMap = new HashMap<>();
        Map<svm_parameter,Double> bestAccuracyMap = new HashMap<>();
        for (String key : paramMap.keySet()) {
	    svm_parameter param = paramMap.get(key);
	    int totalCorrect = totalCorrectMap.get(key);
	    double accuracy = accuracyMap.get(key);
            if (totalCorrect==bestTotalCorrect) {
                bestCorrectMap.put(param, totalCorrect);
                bestAccuracyMap.put(param, accuracy);
            }
        }
	// output
        System.out.println("C\tgamma\tcorrect\taccuracy");
        for (svm_parameter param : bestCorrectMap.keySet()) {
            int totalCorrect = bestCorrectMap.get(param);
            double accuracy = bestAccuracyMap.get(param);
            System.out.println(df.format(param.C)+"\t"+df.format(param.gamma)+"\t"+totalCorrect+"\t"+pf.format(accuracy));
        }
    }

    /**
     * Set verbose flag to true.
     */
    public void setVerbose() {
        verbose = true;
    }

    /**
     * Command-line operation.
     */
    public static void main(String[] args) throws IOException {
        Options options = new Options();
        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd = null;

	Option dataFileOption = new Option("datafile", true, "input data file in SVM format [required]");
	dataFileOption.setRequired(true);
	options.addOption(dataFileOption);
        //
        Option CStepOption = new Option("Cstep", true, "step for n in C=2^n [1]");
        CStepOption.setRequired(false);
        options.addOption(CStepOption);
        //
        Option GStepOption = new Option("gammastep", true, "step for n in gamma=2^n [1]");
        GStepOption.setRequired(false);
        options.addOption(GStepOption);
	//
        Option nFoldOption = new Option("k", true, "k-fold for cross validation ["+SvmUtil.NRFOLD+"]");
        nFoldOption.setRequired(false);
        options.addOption(nFoldOption);
	//
        Option vOption = new Option("v", false, "toggle verbose output");
        vOption.setRequired(false);
        options.addOption(vOption);
	//
	Option nCasesOption = new Option("ncases", true, "set number of cases to use in search [0=all]");
	nCasesOption.setRequired(false);
	options.addOption(nCasesOption);
	//
	Option nControlsOption = new Option("ncontrols", true, "set number of controls to use in search [0=all]");
	nControlsOption.setRequired(false);
	options.addOption(nControlsOption);
	//
	Option kernelOption = new Option("kernel", true, "choose the SVM kernel: LINEAR, POLY, RBF, SIGMOID, PRECOMPUTED [RBF]");
	kernelOption.setRequired(false);
	options.addOption(kernelOption);

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

	int nCases = 0;
	int nControls = 0;
	if (cmd.hasOption("ncases")) nCases = Integer.parseInt(cmd.getOptionValue("ncases"));
	if (cmd.hasOption("ncontrols")) nControls = Integer.parseInt(cmd.getOptionValue("ncontrols"));


        // initialize
	GridSearcher gs = new GridSearcher(datafile, nCases, nControls);
        if (cmd.hasOption("v")) {
            gs.setVerbose();
        }
        if (cmd.hasOption("Cstep")) {
            gs.CPowerStep = Integer.parseInt(cmd.getOptionValue("Cstep"));
        }
        if (cmd.hasOption("gammastep")) {
            gs.GPowerStep = Integer.parseInt(cmd.getOptionValue("gammastep"));
        }
        if (cmd.hasOption("k")) {
            gs.nrFold = Integer.parseInt(cmd.getOptionValue("k"));
        }
	if (cmd.hasOption("kernel")) {
	    String kernelString = cmd.getOptionValue("kernel");
	    if (kernelString.equals("LINEAR")) {
		gs.kernel_type = svm_parameter.LINEAR;
	    } else if (kernelString.equals("POLY")) {
		gs.kernel_type = svm_parameter.POLY;
	    } else if (kernelString.equals("RBF")) {
		gs.kernel_type = svm_parameter.RBF;
	    } else if (kernelString.equals("SIGMOID")) {
		gs.kernel_type = svm_parameter.SIGMOID;
	    } else if (kernelString.equals("PRECOMPUTED")) {
		gs.kernel_type = svm_parameter.PRECOMPUTED;
	    } else {
		System.err.println("ERROR: the accepted values for kernel are: LINEAR, POLY, RBF, SIGMOID, or PRECOMPUTED");
		System.exit(1);
	    }
	}

        // run the search
        gs.run();
    }
}
