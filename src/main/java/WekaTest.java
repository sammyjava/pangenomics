import java.util.List;
import java.util.LinkedList;
import java.util.Random;

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
 
	// choose our classifiers
	List<Classifier> models = new LinkedList<>();
	// trees
	models.add(new DecisionStump()); // one-level decision tree	
	models.add(new HoeffdingTree()); //
	models.add(new J48());           // a decision tree
	models.add(new LMT());           //
	models.add(new RandomForest());  //
	models.add(new RandomTree());    // random tree
	models.add(new REPTree());       //
	// rules
	models.add(new DecisionTable()); // decision table majority classifier
	models.add(new JRip());          //
	models.add(new OneR());          //
	models.add(new PART());          // PART decision list
	models.add(new ZeroR());         // baseline: everything is the majority class
		
	// Run once for each model
	for (Classifier model :models) {
	    Evaluation evaluation = new Evaluation(data);
	    evaluation.crossValidateModel(model, data, KFOLD, new Random(1));
	    System.out.println("--------------------------------------------------------------------------------------------------------------------");
	    System.out.println(model.getClass().getName());
	    System.out.println(evaluation.toSummaryString());
	    System.out.println(evaluation.toClassDetailsString());
	}
    }
}
