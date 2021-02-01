package org.ncgr.svm;

import java.io.BufferedReader;
import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;

import java.util.StringTokenizer;
import java.util.Vector;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import libsvm.svm;
import libsvm.svm_model;
import libsvm.svm_node;
import libsvm.svm_parameter;
import libsvm.svm_problem;

public class SvmTrainer {
    
    svm_parameter param;
    svm_problem prob;
    svm_model model;

    /**
     * Instantiate with given param values.
     */
    public SvmTrainer(svm_parameter param) {
        this.param = param;
    }
    
    /**
     * Train the SVM.
     */
    void run() throws IOException {
        model = svm.svm_train(prob, param);
    }

    /**
     * Save the model to a file.
     */
    void save(String filename) throws IOException {
        svm.svm_save_model(filename, model);
    }

    /**
     * Read a problem in from a file
     * sample   label   feature feature feature feature ...
     * 92550	case	1:-1.0	2:1.0	3:-1.0	4:-1.0  ...
     */
    void readProblem(String filename) throws IOException {
        BufferedReader fp = new BufferedReader(new FileReader(filename));
        Vector<Double> vy = new Vector<Double>();
        Vector<svm_node[]> vx = new Vector<svm_node[]>();
        int max_index = 0;
        // 92550	case	1:-1.0	2:1.0	3:-1.0	4:-1.0
	String line = null;
        while ((line=fp.readLine())!=null) {
            StringTokenizer st = new StringTokenizer(line," \t\n\r\f:");
            String sample = st.nextToken();
            String label = st.nextToken();
            if (label.equals("ctrl")) {
                vy.addElement(-1.0);
            } else if (label.equals("case")) {
                vy.addElement(+1.0);
            } else {
                System.err.println("ERROR: label values in "+filename+" must be 'case' or 'ctrl'.");
                System.exit(1);
            }
            int m = st.countTokens()/2 - 1;
            svm_node[] x = new svm_node[m];
            for (int j=0;j<m;j++) {
                x[j] = new svm_node();
                x[j].index = atoi(st.nextToken());
                x[j].value = atof(st.nextToken());
            }
            if (m>0) max_index = Math.max(max_index, x[m-1].index);
            vx.addElement(x);
        }

        prob = new svm_problem();
        prob.l = vy.size();
        prob.x = new svm_node[prob.l][];
        for (int i=0;i<prob.l;i++)
            prob.x[i] = vx.elementAt(i);
        prob.y = new double[prob.l];
        for (int i=0;i<prob.l;i++)
            prob.y[i] = vy.elementAt(i);

        if (param.gamma == 0 && max_index > 0)
            param.gamma = 1.0/max_index;

        if (param.kernel_type == svm_parameter.PRECOMPUTED) {
            for (int i=0;i<prob.l;i++) {
                if (prob.x[i][0].index != 0) {
                    System.err.print("Wrong kernel matrix: first column must be 0:sample_serial_number\n");
                    System.exit(1);
                }
                if ((int)prob.x[i][0].value <= 0 || (int)prob.x[i][0].value > max_index) {
                    System.err.print("Wrong input format: sample_serial_number out of range\n");
                    System.exit(1);
                }
            }
        }

        fp.close();

        // check the parameters, bail if a problem
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

	Option dataFileOption = new Option("datafile", true, "input data file in SVM format [required]");
	dataFileOption.setRequired(true);
	options.addOption(dataFileOption);
	//
	Option modelFileOption = new Option("modelfile", true, "output model file [required]");
        modelFileOption.setRequired(true);
        options.addOption(modelFileOption);
        //
        Option verboseOption = new Option("v", "verbose", false, "toggle verbose output");
        verboseOption.setRequired(false);
        options.addOption(verboseOption);
        //
	Option kernelOption = new Option("t", "kernel-type", true, "SVM kernel function: LINEAR, POLY, RBF, SIGMOID [RBF]");
	kernelOption.setRequired(false);
	options.addOption(kernelOption);
        //
        Option svmTypeOption = new Option("s", "svm-type", true, "type of SVM: C-SVC, nu-SVC, one-class-SVM, epsilon-SVR, nu-SVR [C-SVC]");
        svmTypeOption.setRequired(false);
        options.addOption(svmTypeOption);
        // 
        Option kernelDegreeOption = new Option("d", "kernel-degree", true, "set degree in kernel function [3]");
        kernelDegreeOption.setRequired(false);
        options.addOption(kernelDegreeOption);
        //
        Option kernelGammaOption = new Option("g", "kernel-gamma", true, "set gamma in kernel function [1/#features]");
        kernelGammaOption.setRequired(false);
        options.addOption(kernelGammaOption);
        //
        Option kernelCoef0Option = new Option("r", "kernel-coef0", true, "set coef0 in kernel function [0]");
        kernelCoef0Option.setRequired(false);
        options.addOption(kernelCoef0Option);
        //
        Option costOption = new Option("c", "cost", true, "set the cost parameter C in C-SVC, epsilon-SVR and nu-SVR [1]");
        costOption.setRequired(false);
        options.addOption(costOption);
        //
        Option nuOption = new Option("n", "nu", true, "set the parameter nu of nu-SVC, one-class SVM, and nu-SVR [0.5]");
        nuOption.setRequired(false);
        options.addOption(nuOption);
        //
        Option epsilonLossOption = new Option("p", "epsilon-loss", true, "set the epsilon value in loss function of epsilon-SVR [0.1]");
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
        Option shrinkingOption = new Option("h", "shrinking", true, "0/1 toggle whether to use the shrinking heuristics [1]");
        shrinkingOption.setRequired(false);
        options.addOption(shrinkingOption);
        //
        Option probabilityEstimatesOption = new Option("b", "probability-estimates", true, "0/1 toggle whether to train a SVC or SVR model for probability estimates [0]");
        probabilityEstimatesOption.setRequired(false);
        options.addOption(probabilityEstimatesOption);
        //
        Option weightOption = new Option("w", "weight", true, "set the parameter C of class i to weight*C, for C-SVC [1]");
        weightOption.setRequired(false);
        options.addOption(weightOption);

        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println(e.getMessage());
            formatter.printHelp("SvmTrainer", options);
            System.exit(1);
        }

        // required options
        String datafilename = cmd.getOptionValue("datafile");
        String modelfilename = cmd.getOptionValue("modelfile");

        // toggles
        boolean verbose = cmd.hasOption("verbose");
	
        // svm_type
        // C_SVC = 0;
        // NU_SVC = 1;
        // ONE_CLASS = 2;
        // EPSILON_SVR = 3;
        // NU_SVR = 4;
        int svmTypeCode = 0;
        if (cmd.hasOption("svm-type")) {
            String svmType = cmd.getOptionValue("svm-type");
            if (svmType.equals("C-SVC")) {
                svmTypeCode = svm_parameter.C_SVC;
            } else if (svmType.equals("nu-SVC")) {
                svmTypeCode = svm_parameter.NU_SVC;
            } else if (svmType.equals("one-class-SVM")) {
                svmTypeCode = svm_parameter.ONE_CLASS;
            } else if (svmType.equals("epsilon-SVR")) {
                svmTypeCode = svm_parameter.EPSILON_SVR;
            } else if (svmType.equals("nu-SVR")) {
                svmTypeCode = svm_parameter.NU_SVR;
            } else {
                System.err.println("ERROR: incorrect svm-type supplied: "+svmType);
                System.exit(1);
            }
        }

	// kernel_type
        // LINEAR = 0;
        // POLY = 1;
        // RBF = 2;
        // SIGMOID = 3;
        int kernelTypeCode = 2;
        if (cmd.hasOption("kernel-type")) {
            String kernelType = cmd.getOptionValue("kernel-type");
            if (kernelType.equals("LINEAR")) {
                kernelTypeCode = svm_parameter.LINEAR;
            } else if (kernelType.equals("POLYNOMIAL")) {
                kernelTypeCode = svm_parameter.POLY;
            } else if (kernelType.equals("RBF")) {
                kernelTypeCode = svm_parameter.RBF;
            } else if (kernelType.equals("SIGMOID")) {
                kernelTypeCode = svm_parameter.SIGMOID;
            } else {
                System.err.println("ERROR: incorrect kernel-type supplied: "+kernelType);
                System.exit(1);
            }
        }

        // -d degree : set degree in kernel function (default 3)
        // -g gamma : set gamma in kernel function (default 1/num_features)
        // -r coef0 : set coef0 in kernel function (default 0)
        // -c cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)
        // -n nu : set the parameter nu of nu-SVC, one-class SVM, and nu-SVR (default 0.5)
        // -p epsilon : set the epsilon in loss function of epsilon-SVR (default 0.1)
        // -m cachesize : set cache memory size in MB (default 100)
        // -e epsilon : set tolerance of termination criterion (default 0.001)
        // -h shrinking: whether to use the shrinking heuristics, 0 or 1 (default 1)
        // -b probability_estimates: whether to train a SVC or SVR model for probability estimates, 0 or 1 (default 0)
        // -wi weight: set the parameter C of class i to weight*C, for C-SVC (default 1)
        // The k in the -g option means the number of attributes in the input data.

        // default SVM parameter values
        svm_parameter param = new svm_parameter();
        param.degree = 3;
        param.gamma = 0;	// = 1/num_features
        param.coef0 = 0;
        param.nu = 0.5;
        param.cache_size = 100;
        param.C = 1;
        param.eps = 1e-3;
        param.p = 0.1;
        param.shrinking = 1;
        param.probability = 0;
        param.nr_weight = 0;
        param.weight_label = new int[0];
        param.weight = new double[0];

        // input values
        param.svm_type = svmTypeCode;
        param.kernel_type = kernelTypeCode;
        if (cmd.hasOption("d")) {
            param.degree = Integer.parseInt(cmd.getOptionValue("d"));
        }
        if (cmd.hasOption("g")) {
            param.gamma = Double.parseDouble(cmd.getOptionValue("g"));
        }
        if (cmd.hasOption("n")) {
            param.nu = Double.parseDouble(cmd.getOptionValue("n"));
        }
        if (cmd.hasOption("m")) {
            param.cache_size = Double.parseDouble(cmd.getOptionValue("m"));
        }
        if (cmd.hasOption("c")) {
            param.C = Double.parseDouble(cmd.getOptionValue("c"));
        }
        if (cmd.hasOption("e")) {
            param.eps = Double.parseDouble(cmd.getOptionValue("e"));
        }
        if (cmd.hasOption("p")) {
            param.p = Double.parseDouble(cmd.getOptionValue("p"));
        }
        if (cmd.hasOption("h")) {
            param.shrinking = Integer.parseInt(cmd.getOptionValue("h"));
        }
        if (cmd.hasOption("b")) {
            param.probability = Integer.parseInt(cmd.getOptionValue("b"));
        }

        // this is weird, setting a static function in svm
        if (verbose) {
            SvmUtil.setVerbose();
        } else {
            SvmUtil.setQuiet();
        }

        // instantiate with this param object
        SvmTrainer st = new SvmTrainer(param);

        // load the problem
        st.readProblem(datafilename);

        // run it
        st.run();

        // save the model file
        st.save(modelfilename);
    }

    /**
     * Python-to-Java noise.
     */
    static double atof(String s) {
        double d = Double.valueOf(s).doubleValue();
        if (Double.isNaN(d) || Double.isInfinite(d)) {
            System.err.print("NaN or Infinity in input\n");
            System.exit(1);
        }
        return(d);
    }

    /**
     * Python-to-Java noise.
     */
    static int atoi(String s) {
        return Integer.parseInt(s);
    }
}
