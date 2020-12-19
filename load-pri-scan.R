##
## load data from a priority scan
##

run1.graph = "HLAB"
run1.totFRs = 975
run1.minPri = c(1,461,429,121,73,44)
run1.numFRs = c(975,106,201,500,600,701)
run1.correct = c(.572,.515,.515,.548,.552,.551)
run1.TPR = c(.559,.982,.952,.672,.564,.553)
run1.FPR = c(.414,.952,.922,.575,.460,.451)
run1.MCC = c(.145,.083,.070,.099,.105,.102)

run1.sensitivity = run1.TPR
run1.specificity = 1.0 - run1.FPR
run1.fracFRs = run1.numFRs/run1.totFRs

run2.graph="SCZ6A+SCZ14C"
run2.totFRs = 1173
run2.minPri = c(1,673,662,647,622,595)
run2.numFRs = c(1173,100,209,503,805,900)
run2.correct = c(.595,.543,.542,.540,.539,.539)
run2.TPR = c(.570,.121,.113,.103,.100,.099)
run2.FPR = c(.380,.035,.028,.024,.022,.020)
run2.MCC = c(.190,.160,.164,.162,.163,.166)

run2.sensitivity = run2.TPR
run2.specificity = 1.0 - run2.FPR
run2.fracFRs = run2.numFRs/run2.totFRs
