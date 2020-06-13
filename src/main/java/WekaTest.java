import java.util.Random;
import java.util.TreeMap;
import java.util.concurrent.ConcurrentHashMap;

import java.text.DecimalFormat;

import weka.classifiers.Classifier;
import weka.classifiers.Evaluation;
import weka.classifiers.rules.*;
import weka.classifiers.trees.*;
import weka.core.Instances;
import weka.core.converters.ConverterUtils.DataSource;
 
public class WekaTest {

    public static int KFOLD = 10;
    public static DecimalFormat pf = new DecimalFormat("0.0%");

    public static void main(String[] args) throws Exception {
        if (args.length!=1) {
            System.err.println("Usage: WekaTest <arff filename>");
            System.exit(1);
        }

        String arffFile = args[0];

	// Read all the instances in the file (ARFF, CSV, XRFF, ...)
	DataSource source = new DataSource(arffFile);
 	Instances data = source.getDataSet();
        // remove the ID attribute
        data.deleteAttributeAt(0);
        // set the class attribute index
        data.setClassIndex(data.numAttributes() - 1);
 
	// load classifiers and their evaluations in a map
	ConcurrentHashMap<Classifier,Evaluation> modelEvaluations = new ConcurrentHashMap<>();
	// trees
	modelEvaluations.put(new DecisionStump(), new Evaluation(data)); // one-level decision tree	
	modelEvaluations.put(new HoeffdingTree(), new Evaluation(data)); //
	modelEvaluations.put(new J48(), new Evaluation(data));           // a decision tree
	modelEvaluations.put(new LMT(), new Evaluation(data));           //
	modelEvaluations.put(new RandomForest(), new Evaluation(data));  //
	modelEvaluations.put(new RandomTree(), new Evaluation(data));    // random tree
	modelEvaluations.put(new REPTree(), new Evaluation(data));       //
	// rules
	modelEvaluations.put(new DecisionTable(), new Evaluation(data)); // decision table majority classifier
	modelEvaluations.put(new JRip(), new Evaluation(data));          //
	modelEvaluations.put(new OneR(), new Evaluation(data));          //
	modelEvaluations.put(new PART(), new Evaluation(data));          // PART decision list
	modelEvaluations.put(new ZeroR(), new Evaluation(data));         // baseline: everything is the majority class

	// Run through the models in a parallel stream
	modelEvaluations.entrySet().parallelStream().forEach(entry -> {
		Classifier model = entry.getKey();
		Evaluation evaluation = entry.getValue();
		try {
		    evaluation.crossValidateModel(model, data, KFOLD, new Random(1));
		} catch (Exception e) {
		    System.err.println(e);
		    System.exit(1);
		}
		System.out.println("Finished "+model.getClass().getName());
	    });

	// spit out the results in order of number correct
	TreeMap<Double,Classifier> sortedModels = new TreeMap<>();
	for (Classifier model : modelEvaluations.keySet()) {
	    Evaluation evaluation = modelEvaluations.get(model);
	    sortedModels.put(evaluation.correct(), model);
	}
	for (Classifier model : sortedModels.values()) {
	    Evaluation evaluation = modelEvaluations.get(model);
	    System.out.println("--------------------------------------------------------------------------------------------------------------------");
	    System.out.println(model.getClass().getName());
	    System.out.println(evaluation.toSummaryString());
	    System.out.println(evaluation.toClassDetailsString());
	}
    }
}
