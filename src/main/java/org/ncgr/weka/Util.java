package org.ncgr.weka;

import java.text.DecimalFormat;

import java.util.List;
import java.util.LinkedList;
import java.util.Random;
import java.util.TreeMap;

import java.util.concurrent.ThreadLocalRandom;

import weka.classifiers.Evaluation;
import weka.classifiers.IterativeClassifier;
import weka.core.Instance;
import weka.core.Instances;

import org.apache.commons.math3.stat.StatUtils;

/**
 * Static utility methods.
 */
public class Util {
    // number formats
    static final DecimalFormat sci = new DecimalFormat("0.00E00");
    static final DecimalFormat perc = new DecimalFormat("0.0%");
    static final DecimalFormat rate = new DecimalFormat("0.000");
    static final DecimalFormat round = new DecimalFormat("0");
    static final String pm = "\u00B1";

    /**
     * Print evaluation result numbers.
     */
    public static void printResults(double correct, double accuracy, double TPR, double FPR, double MCC) {
    	System.out.println("correct\taccur\tTPR\tFPR\tMCC");
	System.out.println(round.format(correct)+"\t"+perc.format(accuracy)+"\t"+rate.format(TPR)+"\t"+rate.format(FPR)+"\t"+rate.format(MCC));
    }

    /**
     * Print evaluation result numbers plus an extra two doubles, without heading.
     */
    public static void printResults(double correct, double accuracy, double TPR, double FPR, double MCC, double d1, double d2) {
	System.out.println(round.format(correct)+"\t"+perc.format(accuracy)+"\t"+rate.format(TPR)+"\t"+rate.format(FPR)+"\t"+rate.format(MCC)+"\t"+sci.format(d1)+"\t"+sci.format(d2));
    }
    
    /**
     * Print evaluation results.
     */
    public static void printResults(Evaluation evaluation) {
	double correct = evaluation.correct();
	double accuracy = evaluation.pctCorrect()/100.0;
	double TPR = evaluation.truePositiveRate(1);
	double FPR = evaluation.falsePositiveRate(1);
	double MCC = evaluation.matthewsCorrelationCoefficient(1);
	printResults(correct, accuracy, TPR, FPR, MCC);
    }

    /**
     * Print evaluation results along with two doubles (like C,gamma). Skip the header.
     */
    public static void printResults(Evaluation evaluation, double d1, double d2) {
	double correct = evaluation.correct();
	double accuracy = evaluation.pctCorrect()/100.0;
	double TPR = evaluation.truePositiveRate(1);
	double FPR = evaluation.falsePositiveRate(1);
	double MCC = evaluation.matthewsCorrelationCoefficient(1);
	printResults(correct, accuracy, TPR, FPR, MCC, d1, d2);
    }

    /**
     * Reduce the number of cases and controls from those in the given data.
     */
    public static Instances reduceData(Instances data, int nCases, int nControls) {
	if (nCases==0 || nControls==0) return data;
	System.err.println("Reducing data to "+nCases+" cases and "+nControls+" controls...");
	Instances reducedData = new Instances(data, (nCases+nControls));
	Random r = new Random(ThreadLocalRandom.current().nextInt());
	List<Integer> usedIndexes = new LinkedList<>();
	double caseValue = 1.0;
	double ctrlValue = 0.0;
	// load cases
	int caseCount = 0;
	while (caseCount<nCases) {
	    int i = r.nextInt(data.size());
	    if (usedIndexes.contains(i)) continue;
	    Instance instance = data.get(i);
	    if (instance.classValue()==caseValue) {
		caseCount++;
		usedIndexes.add(i);
		reducedData.add(instance);
	    }
	}
	// load controls
	int ctrlCount = 0;
	while (ctrlCount<nControls) {
	    int i = r.nextInt(data.size());
	    if (usedIndexes.contains(i)) continue;
	    Instance instance = data.get(i);
	    if (instance.classValue()==ctrlValue) {
		ctrlCount++;
		usedIndexes.add(i);
		reducedData.add(instance);
	    }
	}
	return reducedData;
    }

    /**
     * Rearrage the data, removing the first sample ID column and setting the class attribute (case/ctrl) to the last column.
     */
    static Instances rearrange(Instances data) {
	// remove the sample ID attribute
        data.deleteAttributeAt(0);
        // set the class attribute index at the end (0/1)
        data.setClassIndex(data.numAttributes() - 1);
	return data;
    }

    /**
     * Print out stats for a TreeMap of Evaluations.
     */
    static void printStats(TreeMap<Double,Evaluation> evaluations) {
    	// summary stats
	double corrects[] = new double[evaluations.size()];
	double accuracys[] = new double[evaluations.size()];
	double TPRs[] = new double[evaluations.size()];
	double FPRs[] = new double[evaluations.size()];
	double MCCs[] = new double[evaluations.size()];
	double precisions[] = new double[evaluations.size()];
	double recalls[] = new double[evaluations.size()];
	int i = 0;
	for (Evaluation evaluation : evaluations.values()) {
	    corrects[i] = evaluation.correct();
	    accuracys[i] = evaluation.pctCorrect()/100;
	    // assume case=1
	    TPRs[i] = evaluation.truePositiveRate(1);
	    FPRs[i] = evaluation.falsePositiveRate(1);
	    MCCs[i] = evaluation.matthewsCorrelationCoefficient(1);
	    precisions[i] = evaluation.precision(1);
	    recalls[i] = evaluation.recall(1);
	    i++;
	}
	// means with error bars
	System.out.println("# EVALUATION STATS");
	System.out.println("correct:\t"+round.format(StatUtils.mean(corrects))+pm+round.format(Math.sqrt(StatUtils.variance(corrects))));
	System.out.println("% correct:\t"+perc.format(StatUtils.mean(accuracys))+pm+perc.format(Math.sqrt(StatUtils.variance(accuracys))));
	System.out.println("TPR:\t\t"+rate.format(StatUtils.mean(TPRs))+pm+rate.format(Math.sqrt(StatUtils.variance(TPRs))));
	System.out.println("FPR:\t\t"+rate.format(StatUtils.mean(FPRs))+pm+rate.format(Math.sqrt(StatUtils.variance(FPRs))));
	System.out.println("MCC:\t\t"+rate.format(StatUtils.mean(MCCs))+pm+rate.format(Math.sqrt(StatUtils.variance(MCCs))));
	System.out.println("precision:\t"+rate.format(StatUtils.mean(precisions))+pm+rate.format(Math.sqrt(StatUtils.variance(precisions))));
	System.out.println("recall:\t\t"+rate.format(StatUtils.mean(recalls))+pm+rate.format(Math.sqrt(StatUtils.variance(recalls))));
	// worst result
	System.out.println("# WORST RESULT");
	printResults(evaluations.firstEntry().getValue());
	// mean result in printResults format
	System.out.println("# MEAN RESULT");
	printResults(StatUtils.mean(corrects), StatUtils.mean(accuracys), StatUtils.mean(TPRs), StatUtils.mean(FPRs), StatUtils.mean(MCCs));
	// best result
	System.out.println("# BEST RESULT");
	printResults(evaluations.lastEntry().getValue());
    }


    /**
     * Run boosting rounds on an iterative classifier.
     */
    public static int runBoosting(IterativeClassifier boost) throws Exception {
	int round = 0;
	while (boost.next()) {
	    round++;
	}
	return round;
    }
}
