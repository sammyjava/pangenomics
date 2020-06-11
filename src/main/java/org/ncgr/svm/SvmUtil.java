package org.ncgr.svm;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileNotFoundException;
import java.io.IOException;

import java.util.List;;
import java.util.LinkedList;
import java.util.Optional;

import libsvm.svm;
import libsvm.svm_parameter;
import libsvm.svm_print_interface;

/**
 * Some utility static methods.
 */
public class SvmUtil {

    // cross-validation k-fold default value
    public static int NRFOLD = 10;

    // null output
    static svm_print_interface svm_print_null = new svm_print_interface() { public void print(String s) {} };

    // stdout
    static svm_print_interface svm_print_stdout = new svm_print_interface() { public void print(String s) { System.err.print(s); } };
    
    public static svm_parameter getDefaultParam() {
        svm_parameter param = new svm_parameter();
        param.svm_type = svm_parameter.C_SVC;
        param.kernel_type = svm_parameter.RBF;
        param.degree = 3;
        param.gamma = 0;	// set to 1/num_features in readProblem
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
        return param;
    }

    public static void setQuiet() {
        svm.svm_set_print_string_function(svm_print_null);
    }
    public static void setVerbose() {
        svm.svm_set_print_string_function(svm_print_stdout);
    }

    /**
     * Load a list of samples from a data file.
     */
    public static List<Sample> readSamples(String dataFilename) throws FileNotFoundException, IOException {
        List<Sample> samples = new LinkedList<>();
        BufferedReader reader = new BufferedReader(new FileReader(dataFilename));
        String line = null;
        while ((line=reader.readLine())!=null) {
            Sample sample = new Sample(line);
            samples.add(sample);
        }
        reader.close();
        return samples;
    }

    /**
     * Reduce the number of samples to nCases and nControls
     */
    public static List<Sample> reduceSamples(List<Sample> samples, int nCases, int nControls) {
	List<Sample> chosenSamples = new LinkedList<>();
	// randomly select case and control paths
	if (nCases>0) {
	    int n = 0;
	    while (n<nCases) {
		Optional<Sample> optional = samples.stream().skip((int)(samples.size()*Math.random())).findFirst();
		if (optional.isPresent()) {
		    Sample sample = optional.get();
		    if (isCase(sample) && !chosenSamples.contains(sample)) {
			n++;
			chosenSamples.add(sample);
		    }
		}
	    }
	}
	if (nControls>0) {
	    int n = 0;
	    while (n<nControls) {
		Optional<Sample> optional = samples.stream().skip((int)(samples.size()*Math.random())).findFirst();
		if (optional.isPresent()) {
		    Sample sample = optional.get();
		    if (isControl(sample) && !chosenSamples.contains(sample)) {
			n++;
			chosenSamples.add(sample);
		    }
		}
	    }
	}
	return chosenSamples;
    }
    
    /**
     * Return true if Sample is labeled a case
     */
    public static boolean isCase(Sample sample) {
	return sample.label.equals("case") || sample.label.equals("1") || sample.label.equals("+1");
    }

    /**
     * Return true if Sample is labeled a control
     */
    public static boolean isControl(Sample sample) {
	return sample.label.equals("ctrl") || sample.label.equals("-1");
    }
}
