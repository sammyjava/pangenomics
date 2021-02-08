package org.ncgr.svm;

import java.io.FileNotFoundException;
import java.io.IOException;

import java.text.DecimalFormat;

import java.util.List;
import java.util.LinkedList;
import java.util.Map;
import java.util.Map.Entry;
import java.util.LinkedHashMap;
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
public class SvmGridSearcher {

    static DecimalFormat df = new DecimalFormat("0.00E0");
    static DecimalFormat pf = new DecimalFormat("00.0%");

    boolean verbose = false;

    // grid defaults
    static final int C_POWER_BEGIN = -5;
    static final int C_POWER_END = 15;
    static final int C_POWER_STEP = 1;
    static final int G_POWER_BEGIN = -25;
    static final int G_POWER_END = 0;
    static final int G_POWER_STEP = 1;

    // grid sparsity
    int CPowerStep = C_POWER_STEP;
    int GPowerStep = G_POWER_STEP;

    // SVM kernel default
    String kernelString = "RBF";
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
    public SvmGridSearcher(String inputFilename, int nCases, int nControls) throws FileNotFoundException, IOException {
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
	    System.err.println("SvmGridSearcher ERROR: maxIndex==0, exiting.");
	    System.exit(1);
	}
	// bizarre setting of static function
	SvmUtil.setQuiet();
    }

    /**
     * Run the search for each C (and gamma)
     */
    public void run() {
	if (verbose) {
	    System.err.println("SvmGridSearcher: "+inputFilename);
	    System.err.println("SvmGridSearcher: "+samples.size()+" samples");
            System.err.println("SvmGridSearcher: "+nrFold+"-fold cross-validation");
	    System.err.println("SvmGridSearcher: "+kernelString+" kernel");
	    System.err.println("SvmGridSearcher: C=2^n from "+Math.pow(2.0,C_POWER_BEGIN)+" to "+Math.pow(2.0,C_POWER_END)+" with n-step="+CPowerStep);
	    if (kernel_type==svm_parameter.RBF) {
		System.err.println("SvmGridSearcher: gamma=2^n from "+Math.pow(2.0,G_POWER_BEGIN)+" to "+Math.pow(2.0,G_POWER_END)+" with n-step="+GPowerStep);
	    }
	}
	ConcurrentSkipListSet<Double> CSet = new ConcurrentSkipListSet<>();
        for (int n=C_POWER_BEGIN; n<=C_POWER_END; n+=CPowerStep) {
            CSet.add(Math.pow(2.0, n));
        }
	ConcurrentSkipListSet<Double> gammaSet = new ConcurrentSkipListSet<>();
	for (int n=G_POWER_BEGIN; n<=G_POWER_END; n+=GPowerStep) {
            gammaSet.add(Math.pow(2.0, n));
	}
	// these maps are keyed by the svm_parameter
	ConcurrentSkipListMap<svm_parameter,Integer> totalCorrectMap = new ConcurrentSkipListMap<>();
	ConcurrentSkipListMap<svm_parameter,Double> accuracyMap = new ConcurrentSkipListMap<>();
        ////////////////////////////////////////////////////////////////////////////////////////////////
        // C loop
        CSet.parallelStream().forEach(C -> {
		if (kernel_type==svm_parameter.RBF) {
		    ////////////////////////////////////////////////////////////////////////////////////////////////
		    // gamma loop
		    gammaSet.parallelStream().forEach(gamma -> {
			    svm_parameter param = SvmUtil.getDefaultParam();
			    param.svm_type = svm_parameter.C_SVC;   // C_SVC, NU_SVC, ONE_CLASS, EPSILON_SVR, NU_SVR
			    param.C = C;
			    param.kernel_type = kernel_type;        // LINEAR, POLY, RBF, SIGMOID, PRECOMPUTED
			    param.gamma = gamma;
			    SvmCrossValidator svc = new SvmCrossValidator(param, nrFold, vy, vx);
			    svc.samples = samples;
			    svc.inputFilename = inputFilename;
			    svc.run();
			    if (verbose) {
				System.err.println("SvmGridSearcher: C="+df.format(param.C)+" gamma="+df.format(param.gamma)+" totalCorrect="+svc.totalCorrect+" accuracy="+pf.format(svc.accuracy));
			    }
			    // store the results
			    totalCorrectMap.put(param, svc.totalCorrect);
			    accuracyMap.put(param, svc.accuracy);
			});
		    ////////////////////////////////////////////////////////////////////////////////////////////////
		} else {
		    svm_parameter param = SvmUtil.getDefaultParam();
		    param.svm_type = svm_parameter.C_SVC;   // C_SVC, NU_SVC, ONE_CLASS, EPSILON_SVR, NU_SVR
		    param.C = C;
		    param.kernel_type = kernel_type;        // LINEAR, POLY, RBF, SIGMOID, PRECOMPUTED
		    SvmCrossValidator svc = new SvmCrossValidator(param, nrFold, vy, vx);
		    svc.samples = samples;
		    svc.inputFilename = inputFilename;
		    svc.run();
		    if (verbose) {
			System.err.println("SvmGridSearcher: C="+df.format(param.C)+" totalCorrect="+svc.totalCorrect+" accuracy="+pf.format(svc.accuracy));
		    }
		    // store the results
		    totalCorrectMap.put(param, svc.totalCorrect);
		    accuracyMap.put(param, svc.accuracy);
		}
            });            
        ////////////////////////////////////////////////////////////////////////////////////////////////
	// sort the results by total correct
        List<Entry<svm_parameter,Integer>> totalCorrectList = new LinkedList<>(totalCorrectMap.entrySet());
        totalCorrectList.sort(Entry.comparingByValue());
        Map<svm_parameter,Integer> sortedCorrectMap = new LinkedHashMap<>();
	Map<svm_parameter,Double> sortedAccuracyMap = new LinkedHashMap<>();
	for (Entry<svm_parameter,Integer> entry : totalCorrectList) {
	    svm_parameter param = entry.getKey();
	    int totalCorrect = entry.getValue();
	    double accuracy = accuracyMap.get(param);
            sortedCorrectMap.put(param, totalCorrect);
	    sortedAccuracyMap.put(param, accuracy);
	}
	// output
	if (kernel_type==svm_parameter.RBF) {
	    System.out.println("C\tgamma\tcorrect\taccuracy");
	    System.err.println("C\tgamma\tcorrect\taccuracy");
	} else {
	    System.out.println("C\tcorrect\taccuracy");
	    System.err.println("C\tcorrect\taccuracy");
	}
	for (svm_parameter param : sortedCorrectMap.keySet()) {
	    int totalCorrect = sortedCorrectMap.get(param);
	    double accuracy = sortedAccuracyMap.get(param);
	    if (kernel_type==svm_parameter.RBF) {
		System.out.println(df.format(param.C)+"\t"+df.format(param.gamma)+"\t"+totalCorrect+"\t"+pf.format(accuracy));
		System.err.println(df.format(param.C)+"\t"+df.format(param.gamma)+"\t"+totalCorrect+"\t"+pf.format(accuracy));
	    } else {
		System.out.println(df.format(param.C)+"\t"+totalCorrect+"\t"+pf.format(accuracy));
		System.err.println(df.format(param.C)+"\t"+totalCorrect+"\t"+pf.format(accuracy));
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
            formatter.printHelp("SvmGridSearcher [options]", options);
            System.exit(1);
        }
	
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println(e.getMessage());
            formatter.printHelp("SvmGridSearcher", options);
            System.exit(1);
        }

        String datafile = cmd.getOptionValue("datafile");

	int nCases = 0;
	int nControls = 0;
	if (cmd.hasOption("ncases")) nCases = Integer.parseInt(cmd.getOptionValue("ncases"));
	if (cmd.hasOption("ncontrols")) nControls = Integer.parseInt(cmd.getOptionValue("ncontrols"));


        // initialize
	SvmGridSearcher gs = new SvmGridSearcher(datafile, nCases, nControls);
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
	    gs.setKernel(cmd.getOptionValue("kernel"));
	}

        // run the search
        gs.run();
    }
}
