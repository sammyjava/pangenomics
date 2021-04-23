package org.ncgr.weka;

import java.util.List;
import java.util.LinkedList;
import java.util.Random;

import weka.classifiers.Evaluation;
import weka.classifiers.rules.JRip;
import weka.core.Instances;
import weka.core.converters.ConverterUtils.DataSource;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
 
/**
 * Run Weka JRip on an ARFF file. Valid JRip options are:
 *
 *  -F <number of folds>
 *   Set number of folds for REP
 *   One fold is used as pruning set.
 *   (default 3)
 *
 *  -N <min. weights>
 *   Set the minimal weights of instances
 *   within a split.
 *   (default 2.0)
 *
 *  -O <number of runs>
 *   Set the number of runs of
 *   optimizations. (Default: 2)
 *
 *  -D
 *   Set whether turn on the
 *   debug mode (Default: false)
 *
 *  -S <seed>
 *   The seed of randomization
 *   (Default: 1)
 *
 *  -E
 *   Whether NOT check the error rate>=0.5
 *   in stopping criteria  (default: check)
 *
 *  -P
 *   Whether NOT use pruning
 *   (default: use pruning)
 */
public class WekaJRip {

    public static void main(String[] args) throws Exception {

	// JRip defaults
	int F = 3;         // number of folds for REP; one fold is used as pruning set.
	double N = 2.0;    // minimal weights of instances within a split.
	int O = 2;         // the number of runs of optimizations.
	boolean D = false; // true = debug mode
	int S = 1;         // seed of randomization
	boolean E = false; // true = do NOT check the error rate>=0.5 in stopping criteria
	boolean P = false; // true = do NOT use pruning

	// cross-validation defaults
	int kfold = 10;
	
	Options options = new Options();
        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd;

        Option arffFileOption = new Option("a", "arfffile", true, "ARFF input file (required)");
        arffFileOption.setRequired(true);
        options.addOption(arffFileOption);
	//
	Option kfoldOption = new Option("k", "kfold", true, "cross-validation k-fold ("+kfold+")");
	kfoldOption.setRequired(false);
	options.addOption(kfoldOption);
	//
	Option FOption = new Option("F", "folds", true, "number of folds for REP; one fold is used as pruning set ("+F+")");
	FOption.setRequired(false);
	options.addOption(FOption);
	//
	Option NOption = new Option("N", "minimalweights", true, "minimal weights of instances within a split ("+N+")");
	NOption.setRequired(false);
	options.addOption(NOption);
	//
	Option OOption = new Option("O", "nruns", true, "number of runs of optimizations ("+O+")");
	OOption.setRequired(false);
	options.addOption(OOption);
	//
	Option DOption = new Option("D", "debug", false, "enable debug option ("+D+")");
	DOption.setRequired(false);
	options.addOption(DOption);
	//
	Option SOption = new Option("S", "seed", true, "seed of randomization ("+S+")");
	SOption.setRequired(false);
	options.addOption(SOption);
	//
	Option EOption = new Option("E", false, "do NOT check the error rate>=0.5 in stopping criteria ("+E+")");
	EOption.setRequired(false);
	options.addOption(EOption);
	//
	Option POption = new Option("P", false, "do NOT use pruning ("+P+")");
	POption.setRequired(false);
	options.addOption(POption);
	
	try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println(e.getMessage());
            formatter.printHelp("WekaJRip", options);
            System.exit(1);
            return;
        }

        String arffFile = cmd.getOptionValue("arfffile");

	// cross-validation options
	if (cmd.hasOption("kfold")) kfold = Integer.parseInt(cmd.getOptionValue("kfold"));

	// JRip options
	if (cmd.hasOption("F")) F = Integer.parseInt(cmd.getOptionValue("F"));
	if (cmd.hasOption("N")) N = Double.parseDouble(cmd.getOptionValue("N"));
	if (cmd.hasOption("O")) O = Integer.parseInt(cmd.getOptionValue("O"));
	if (cmd.hasOption("S")) S = Integer.parseInt(cmd.getOptionValue("S"));
	D = cmd.hasOption("D");
	E = cmd.hasOption("E");
	P = cmd.hasOption("P");

	// form the JRip option string
	List<String> jripOptionsList = new LinkedList<>();
	jripOptionsList.add("-F"); jripOptionsList.add(String.valueOf(F));
	jripOptionsList.add("-N"); jripOptionsList.add(String.valueOf(N));
	jripOptionsList.add("-O"); jripOptionsList.add(String.valueOf(O));
	jripOptionsList.add("-S"); jripOptionsList.add(String.valueOf(S));
	if (D) jripOptionsList.add("-D");
	if (E) jripOptionsList.add("-E");
	if (P) jripOptionsList.add("-P");
	
	String[] jripOptions = new String[jripOptionsList.size()];
	int i = 0;
	for (String jripOption : jripOptionsList) {
	    jripOptions[i++] = jripOption;
	}

	// Read all the instances in the file (ARFF, CSV, XRFF, ...)
 	Instances data = Util.rearrange(new DataSource(arffFile).getDataSet());
        // remove the ID attribute
        data.deleteAttributeAt(0);
        // set the class attribute index
        data.setClassIndex(data.numAttributes() - 1);
	Evaluation evaluation = new Evaluation(data);
 
	// JRip
	JRip model = new JRip();
	model.setOptions(jripOptions);
	
	System.err.print("Starting JRip evaluation with: ");
	for (String option : model.getOptions()) {
	    System.err.print(option+" ");
	}
	System.err.println("");

	// run the cross-validation
	System.err.println("Running "+kfold+" cross-validations...");
	evaluation.crossValidateModel(model, data, kfold, new Random(1));

	// output
	System.out.println("--------------------------------------------------------------------------------------------------------------------");
	System.out.println(model.getClass().getName()+"\t"+evaluation.matthewsCorrelationCoefficient(0));
	System.out.println(evaluation.toSummaryString());
	System.out.println(evaluation.toClassDetailsString());
	System.out.println("--------------------------------------------------------------------------------------------------------------------");
    }
}
