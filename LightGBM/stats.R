## pathfrs.train
## pathfrs.test
## LightGBM_model.txt
## LightGBM_predict_result.txt

test = read.table("pathfrs.test")
predict = read.table("LightGBM_predict_result.txt")

TP = nrow(test[test$V1==1 & predict$V1>0.5,])
FP = nrow(test[test$V1==0 & predict$V1>0.5,])

TN = nrow(test[test$V1==0 & predict$V1<0.5,])
FN = nrow(test[test$V1==1 & predict$V1<0.5,])

TPR = TP/(TP+FN)
FPR = FP/(TN+FP)

correct = (TP+TN)/(TP+TN+FP+FN)

