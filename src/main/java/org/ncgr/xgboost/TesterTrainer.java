package org.ncgr.xgboost;

import java.util.Map;
import java.util.HashMap;

import ml.dmlc.xgboost4j.java.Booster;
import ml.dmlc.xgboost4j.java.DMatrix;
import ml.dmlc.xgboost4j.java.XGBoost;
import ml.dmlc.xgboost4j.java.XGBoostError;

/**
 * Test/train from a couple of SVM format files.
 *
 * @author Sam Hokin
 */
public class TesterTrainer {

    public static void main(String[] args) throws XGBoostError {

	DMatrix trainMat = new DMatrix("pathfrs.train.svm");
	DMatrix testMat = new DMatrix("pathfrs.test.svm");
	
	// JSON-ish parameters
	Map<String,Object> params = new HashMap<String,Object>() {
		{
		    put("eta", 1.0);
		    put("max_depth", 2);
		    put("objective", "binary:logistic");
		    put("eval_metric", "logloss");
		}
	    };
	
	// Specify a watch list to see model accuracy on data sets
	Map<String,DMatrix> watches = new HashMap<String, DMatrix>() {
		{
		    put("train", trainMat);
		    put("test", testMat);
		}
	    };
	int nround = 2;
	Booster booster = XGBoost.train(trainMat, params, nround, watches, null, null);
	
	// save a model
	// booster.saveModel("model.bin");
	
	// load a model
	// Booster booster = XGBoost.loadModel("model.bin");
	
	// predict
	float[][] predicts = booster.predict(testMat);
	// predict leaf
	float[][] leafPredicts = booster.predictLeaf(testMat, 0);

	
	// spit out the results
	int TP = 0;
	int FP = 0;
	int TN = 0;
	int FN = 0;
	for (int i=0; i<predicts.length; i++) {
	    boolean isCase = testMat.getLabel()[i]>0.5;
	    for (int j=0; j<predicts[i].length; j++) {
		double predict = predicts[i][j];
		if (predict>0.5) {
		    // positive call
		    if (isCase) {
			TP++;
		    } else {
			FP++;
		    }
		} else {
		    // negative call
		    if (isCase) {
			FN++;
		    } else {
			TN++;
		    }
		}
	    }
	}
	// stats
	int totCase = TP + FN;
	int totCtrl = FP + TN;
	int correct = TP + TN;
	int total = totCase + totCtrl;
	double accuracy = (double)correct / (double)total;
	double TPR = (double)TP / (double)totCase;
	double FPR = (double)FP / (double)totCtrl;
	double MCC = ((double)(TP*TN) - (double)(FP*FN)) / Math.sqrt((double)(TP+FP)*(double)(TP+FN)*(double)(TN+FP)*(double)(TN+FN));
	System.out.println("TP\tFP\tTN\tFN\tAcc\tTPR\tFPR\tMCC");
	System.out.println(TP+"\t"+FP+"\t"+TN+"\t"+FN+"\t"+accuracy+"\t"+TPR+"\t"+FPR+"\t"+MCC);
    }
}
