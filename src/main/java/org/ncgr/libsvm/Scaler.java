package org.ncgr.libsvm;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;

import java.text.DecimalFormat;

import java.util.Formatter;
import java.util.StringTokenizer;
import java.util.List;
import java.util.ArrayList;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

/**
 * Utility to scale data to a given interval like [-1,1].
 */
public class Scaler {
    
    // defaults
    double lower = -1.0;
    double upper = 1.0;
    double yLower;
    double yUpper;
    boolean yScaling = false;

    double[] featureMax;
    double[] featureMin;
    double yMax = +1; // case=+1
    double yMin = -1; // ctrl=-1
    int max_index = 0;

    // our samples
    List<Sample> samples;

    // optional restore file
    String restoreFilename;

    /**
     * Construct by reading samples in from an SVM-format file.
     */
    public Scaler(String svmFilename) throws FileNotFoundException, IOException {
        this.samples = Util.readSamples(svmFilename);
    }

    /**
     * Construct given a List of samples.
     */
    public Scaler(List<Sample> samples) {
        this.samples = samples;
    }

    /**
     * Pass 1: find out max index of attributes from the samples or the restore file.
     * Assumption: min index of attributes is 1.
     */
    void findMaxIndex() throws FileNotFoundException, IOException {
        max_index = 0;
	if (restoreFilename==null) {
	    // samples data
	    for (Sample sample : samples) {
		for (int i : sample.values.keySet()) {
		    max_index = Math.max(max_index, i);
		}
	    }
	} else {
	    // x
	    // -1.000000000000000 1.000000000000000
	    // 0 -1.000000000000000 1.000000000000000
	    int idx, c;
	    BufferedReader restoreReader = new BufferedReader(new FileReader(restoreFilename));
	    // skip some lines
	    if ((c=restoreReader.read())=='y') {
		restoreReader.readLine();
		restoreReader.readLine();
		restoreReader.readLine();
	    }
	    restoreReader.readLine();
	    restoreReader.readLine();
	    // read the range lines
	    String line = null;
	    while ((line=restoreReader.readLine())!=null) {
		StringTokenizer st = new StringTokenizer(line);
		idx = Integer.parseInt(st.nextToken());
 		max_index = Math.max(max_index, idx);
 	    }
	    restoreReader.close();
	}
    }
    
    /**
     * Pass 2: find out min/max value of each feature.
     */
    void findMinMaxValues() throws FileNotFoundException, IOException {
	// get featureMin/featureMax from samples
        featureMax = new double[max_index+1];
        featureMin = new double[max_index+1];
        for (int i=0; i<max_index; i++) {
            featureMax[i] = -Double.MAX_VALUE;
            featureMin[i] =  Double.MAX_VALUE;
        }
        for (Sample sample : samples) {
            for (int index : sample.values.keySet()) {
                double value = sample.values.get(index);
                int i = index - 1;
                featureMax[i] = Math.max(featureMax[i], value);
                featureMin[i] = Math.min(featureMin[i], value);
            }
        }
	// Pass 2.5: get featureMin/featureMax from restore file.
	if (restoreFilename!=null) {
	    int idx, c;
	    double fmin, fmax;
	    BufferedReader restoreReader = new BufferedReader(new FileReader(restoreFilename));
	    c = restoreReader.read();
	    if (c=='y') {
		// NOTHING HAPPENS HERE!
		System.err.println("'y' found, so nothing happens!");
		System.exit(1);
		// restoreReader.readLine();
		// StringTokenizer st = new StringTokenizer(restoreReader.readLine());
		// double y_lower = Double.parseDouble(st.nextToken());
		// double y_upper = Double.parseDouble(st.nextToken());
		// st = new StringTokenizer(restoreReader.readLine());
		// double y_min = Double.parseDouble(st.nextToken());
		// double y_max = Double.parseDouble(st.nextToken());
		// boolean y_scaling = true;
	    } else if (c=='x') {
		restoreReader.readLine();		// pass the '\n' after 'x'
		StringTokenizer st = new StringTokenizer(restoreReader.readLine());
		lower = Double.parseDouble(st.nextToken());
		upper = Double.parseDouble(st.nextToken());
		String line = null;
		while ((line=restoreReader.readLine())!=null) {
		    StringTokenizer st2 = new StringTokenizer(line);
		    idx = Integer.parseInt(st2.nextToken());
		    fmin = Double.parseDouble(st2.nextToken());
		    fmax = Double.parseDouble(st2.nextToken());
		    if (idx<=max_index) {
			featureMin[idx] = fmin;
			featureMax[idx] = fmax;
		    }
		}
	    }
	    restoreReader.close();
	}
    }

    /**
     * Pass 3: scale to an output file.
     * TODO: just scale to new list of samples; move output to another method.
     */
    void scaleToOutput(PrintStream out) throws IOException {
        for (Sample sample : samples) {
            out.print(sample.name+"\t"+sample.label);
            int next_index = 1;
            for (int index : sample.values.keySet()) {
                for (int i=next_index; i<index; i++) output(out, i, 0); // ???
                output(out, index, sample.values.get(index));
                next_index = index + 1;
            }
            for (int i=next_index; i<=max_index; i++) output(out, i, 0); // ???
            out.println("");
        }
    }
    void output(PrintStream out, int index, double value) {
        int i = index - 1;
        // skip single-valued attribute
        if (featureMax[i]==featureMin[i]) return;
        if (value==featureMin[i]) {
            value = lower;
        } else if (value==featureMax[i]) {
            value = upper;
        } else {
            value = lower + (upper-lower) *  (value-featureMin[i])/(featureMax[i]-featureMin[i]);
        }
        if (value!=0) {
            out.print("\t"+index+":"+value);
        }
    }

    /**
     * Write to the save file.
     */
    void writeSaveFile(String saveFilename) throws IOException {
        Formatter formatter = new Formatter(new StringBuilder());
        BufferedWriter bw = null;
        bw = new BufferedWriter(new FileWriter(saveFilename));
        if (yScaling) {
            formatter.format("y\n");
            formatter.format("%.16g %.16g\n", yLower, yUpper);
            formatter.format("%.16g %.16g\n", yMin, yMax);
        }
        formatter.format("x\n");
        formatter.format("%.16g %.16g\n", lower, upper);
        for (int i=0; i<max_index; i++) {
            if (featureMin[i]!=featureMax[i]) formatter.format("%d %.16g %.16g\n", i, featureMin[i], featureMax[i]);
        }
        bw.write(formatter.toString());
        bw.close();
    }

    /**
     * Command-line version, loads data from files.
     */
    public static void main(String[] args) throws IOException {

        Options options = new Options();
        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd = null;

	Option svmFileOption = new Option("d", "svmfile", true, "Data file in SVM format");
	svmFileOption.setRequired(true);
	options.addOption(svmFileOption);
	//
        Option saveFileOption = new Option("s", "savefile", true, "file to save scaling parameters to (exclusive from -r)");
        saveFileOption.setRequired(false);
        options.addOption(saveFileOption);
	//
	Option restoreFileOption = new Option("r", "restorefile", true, "file to read scaling parameters from (exclusive from -s)");
	restoreFileOption.setRequired(false);
	options.addOption(restoreFileOption);
	//
        Option lowerOption = new Option("l", "xlowerlimit", true, "x scaling lower limit [-1]");
        lowerOption.setRequired(false);
        options.addOption(lowerOption);
	//
        Option upperOption = new Option("u", "xupperlimit", true, "x scaling upper limit [+1]");
        upperOption.setRequired(false);
        options.addOption(upperOption);
	//
        Option yScalingOption = new Option("y", "ylimits", true, "y scaling limits yLower,yUpper [no y scaling]");
        yScalingOption.setRequired(false);
        options.addOption(yScalingOption);

        if (args.length==0) {
            formatter.printHelp("Scaler [options]", options);
            System.exit(1);
        }
        
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println(e.getMessage());
            formatter.printHelp("Scaler [options]", options);
            System.exit(1);
        }

        // instantiate default instance from data file
        String svmFilename = cmd.getOptionValue("svmfile");
	// other files
        String saveFilename = null;
        if (cmd.hasOption("s")) {
            saveFilename = cmd.getOptionValue("s");
        }
	String restoreFilename = null;
	if (cmd.hasOption("r")) {
	    restoreFilename = cmd.getOptionValue("r");
	}

	// instantiate the Scaler and validate
        Scaler scaler = new Scaler(svmFilename);
	scaler.validateFiles(svmFilename, saveFilename, restoreFilename);
	scaler.restoreFilename = restoreFilename;

        // update scaler with options and validate
        if (cmd.hasOption("l")) {
            scaler.lower = Double.parseDouble(cmd.getOptionValue("l"));
        }
        if (cmd.hasOption("u")) {
            scaler.upper = Double.parseDouble(cmd.getOptionValue("u"));
        }
        if (cmd.hasOption("y")) {
            String[] parts = cmd.getOptionValue("y").split(",");
            scaler.yLower = Double.parseDouble(parts[0]);
            scaler.yUpper = Double.parseDouble(parts[1]);
            scaler.yScaling = true;
        }
        scaler.validateParameters();

        // run scaling and output
        scaler.findMaxIndex();
        scaler.findMinMaxValues();
        scaler.scaleToOutput(System.out);
	
        // write save file
        if (saveFilename!=null) {
            scaler.writeSaveFile(saveFilename);
        }
    }

    /**
     * Validate files.
     */
    void validateFiles(String svmFilename, String saveFilename, String restoreFilename) {
        if (svmFilename==null) {
            System.err.println("You have not provided a data file to scale.");
            System.exit(1);
	}
        if (saveFilename!=null && restoreFilename!=null) {
	    System.err.println("ERROR: you may not specify both a save file (-s) and a restore file (-r)");
	    System.exit(1);
	}
    }

    /**
     * Validate parameters
     */
    void validateParameters() {
        if (!(upper>lower)) {
            System.err.println("You have provided inconsistent lower/upper values.");
            System.exit(1);
        }
        if (yScaling && !(yUpper>yLower)) {
            System.err.println("You have provided inconsistent y-scaling values.");
            System.exit(1);
        }
    }
}
