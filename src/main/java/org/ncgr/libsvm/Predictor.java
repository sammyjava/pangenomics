package org.ncgr.libsvm;

import java.io.BufferedReader;
import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FileNotFoundException;
import java.io.IOException;

import java.util.Formatter;
import java.util.StringTokenizer;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import libsvm.svm;
import libsvm.svm_node;
import libsvm.svm_model;
import libsvm.svm_parameter;
import libsvm.svm_print_interface;

public class Predictor {

    int svm_type;
    int nr_class;
    boolean predictProbability;
    
    // files
    BufferedReader input;
    svm_model model;
    DataOutputStream output;
    
    // classification output
    public int total;
    public int correct;
    public double accuracy;

    // regression output
    public double error;
    public double meanSquareError;
    public double squaredCorrCoeff;

    // file names
    String inputFilename;
    String modelFilename;
    String outputFilename;

    // the chosen output
    static svm_print_interface svm_print_string;

    /**
     * Construct given file names and the predictProbability flag
     */
    public Predictor(String inputFilename, String modelFilename, String outputFilename, boolean predictProbability) throws FileNotFoundException, IOException {
        this.inputFilename = inputFilename;
        this.modelFilename = modelFilename;
        this.outputFilename = outputFilename;
        
        this.input = new BufferedReader(new FileReader(inputFilename));
        this.model = svm.svm_load_model(modelFilename);
        this.output = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFilename)));
        this.predictProbability = predictProbability;
        
        // validation
        if (model == null) {
            System.err.println("Can't open model file:"+modelFilename);
            System.exit(1);
        }
        if (predictProbability) {
            if (svm.svm_check_probability_model(model)==0) {
                System.err.println("Model does not support probability estimates");
                System.exit(1);
            }
        } else if (svm.svm_check_probability_model(model)!=0) {
            System.out.println("Model supports probability estimates, but disabled in prediction.");
        }
    }

    /**
     * Run the prediction on the input file.
     */
    void run() throws IOException {
        double sumv = 0, sumy = 0, sumvv = 0, sumyy = 0, sumvy = 0;

        svm_type = svm.svm_get_svm_type(model);
        nr_class = svm.svm_get_nr_class(model);

        double[] prob_estimates = null;
        
        if (predictProbability) {
            if (svm_type == svm_parameter.EPSILON_SVR || svm_type == svm_parameter.NU_SVR) {
                System.out.println("Prob. model for test data: target value = predicted value + z,\nz: Laplace distribution e^(-|z|/sigma)/(2sigma),sigma="+svm.svm_get_svr_probability(model));
            } else {
                int[] labels=new int[nr_class];
                svm.svm_get_labels(model,labels);
                prob_estimates = new double[nr_class];
                output.writeBytes("labels");
                for (int j=0;j<nr_class;j++) {
                    output.writeBytes(" "+labels[j]);
                }
                output.writeBytes("\n");
            }
        }
        
        String line;
        while ((line=input.readLine())!=null) {
            StringTokenizer st = new StringTokenizer(line," \t\n\r\f:");
            String sample = st.nextToken();
            String label = st.nextToken();
            double target = 0;
            if (label.equals("ctrl")) {
                target = -1;
            } else if (label.equals("case")) {
                target = +1;
            } else {
                System.err.println("ERROR: label values in "+inputFilename+" must be 'case' or 'ctrl'.");
                System.exit(1);
            }
            int m = st.countTokens()/2 - 1;
            svm_node[] x = new svm_node[m];
            for (int j=0; j<m; j++) {
                x[j] = new svm_node();
                x[j].index = atoi(st.nextToken());
                x[j].value = atof(st.nextToken());
            }
            
            double v;
            if (predictProbability && (svm_type==svm_parameter.C_SVC || svm_type==svm_parameter.NU_SVC)) {
                v = svm.svm_predict_probability(model, x, prob_estimates);
                output.writeBytes(v+" ");
                for (int j=0;j<nr_class;j++) {
                    output.writeBytes(prob_estimates[j]+"\t");
                }
                output.writeBytes("\n");
            } else {
                v = svm.svm_predict(model,x);
                output.writeBytes((int)v+"\t");
            }

            if (v == target) ++correct;
            error += (v-target)*(v-target);
            sumv += v;
            sumy += target;
            sumvv += v*v;
            sumyy += target*target;
            sumvy += v*target;
            ++total;
        }
        input.close();
        output.close();
        
        if (svm_type == svm_parameter.EPSILON_SVR || svm_type == svm_parameter.NU_SVR) {
            meanSquareError = error/total;
            squaredCorrCoeff = ((total*sumvy-sumv*sumy)*(total*sumvy-sumv*sumy))/((total*sumvv-sumv*sumv)*(total*sumyy-sumy*sumy));
        } else {
            accuracy = (double)correct/(double)total;
        }
    }

    public static void main(String[] args) throws IOException {
        Options options = new Options();
        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd;
        
	Option dataFileOption = new Option("datafile", true, "input data file in SVM format [required]");
	dataFileOption.setRequired(true);
	options.addOption(dataFileOption);
	//
	Option modelFileOption = new Option("modelfile", true, "output model file [required]");
        modelFileOption.setRequired(true);
        options.addOption(modelFileOption);
        //
        Option outputFileOption = new Option("outputfile", true, "output prediction file [required]");
        outputFileOption.setRequired(true);
        options.addOption(outputFileOption);
        //
        Option predictProbabilityOption = new Option("b", "predict-probability", false, "toggle prediction of probability estimates [false]");
        predictProbabilityOption.setRequired(false);
        options.addOption(predictProbabilityOption);
        // 
        Option verboseOption = new Option("v", "verbose", false, "verbose output");
        verboseOption.setRequired(false);
        options.addOption(verboseOption);

        if (args.length==0) {
            formatter.printHelp("Predictor [options]", options);
            System.exit(1);
        }
        
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.out.println(e.getMessage());
            formatter.printHelp("Predictor [options]", options);
            System.exit(1);
            return;
        }

        // set values based on options
        boolean predictProbability = cmd.hasOption("b");
        
        // odd way to toggle printing with a static var
        if (cmd.hasOption("v")) {
            Util.setVerbose();
        } else {
            Util.setVerbose();
        }

        // file names
        String datafilename = cmd.getOptionValue("datafile");
        String modelfilename = cmd.getOptionValue("modelfile");
        String outputfilename = cmd.getOptionValue("outputfile");

        // create the predictor
        Predictor sp = new Predictor(datafilename, modelfilename, outputfilename, predictProbability);
        
        // run it
        sp.run();

        // output
        if (sp.svm_type==svm_parameter.EPSILON_SVR || sp.svm_type==svm_parameter.NU_SVR) {
            System.out.println("Mean squared error = "+sp.meanSquareError+" (regression)");
            System.out.println("Squared correlation coefficient = "+sp.squaredCorrCoeff+" (regression)");
        } else {
            System.out.println("Accuracy = "+sp.accuracy*100+"% ("+sp.correct+"/"+sp.total+") (classification)");
        }
    }

    /**
     * Silly Python to Java function.
     */
    static double atof(String s) {
        return Double.valueOf(s).doubleValue();
    }
    
    /**
     * Silly Python to Java function.
     */
    static int atoi(String s) {
        return Integer.parseInt(s);
    }

}
