### Matthews Correlation Coefficient ---
calcMCC <- function (CM){
	TP <- CM[1,1]
	FP <- CM[1,2]
	TN <- CM[2,2]
	FN <- CM[2,1]
	MCC <- ((TP*TN) - (FP*FN)) / sqrt((TP+FP)*(TP+FN)*(FP+TN)*(TN+FN))
	if (is.na(MCC)) {
		MCC <- 0
	}
	return(MCC)
}

### Accuracy ---	
calcAC <- function (CM) {
	TP <- CM[1,1]
	FP <- CM[1,2]
	TN <- CM[2,2]
	FN <- CM[2,1]
	AC <- (TP+TN) / (TP+FP+TN+FN)
	return(AC)
}

calculatePerformance <- function(data){
	SN <- sensitivity(data$PredClass, data$TrueClass)
	SP <- specificity(data$PredClass, data$TrueClass)
	PPV <- posPredValue(data$PredClass, data$TrueClass)
	NPV <- negPredValue(data$PredClass, data$TrueClass)
	LRP <- SN/(1-SP) # likelihood ratio positive
	LRN <- (1-SN)/SP # likelihood ratio negative	
	
	CM.table <- confusionMatrix(data$PredClass, data$TrueClass)$table
	TP <- CM.table[1,1]
	FP <- CM.table[1,2]
	TN <- CM.table[2,2]
	FN <- CM.table[2,1]	

	MCC <- calcMCC(CM.table)
	AC <- calcAC(CM.table)
	#F <- (2*SN*PPV) / (SN+PPV)
	AUC <- 0.5 * (SN + PPV)
					
	FPR <- FP / (FP+TN)
	FNR <- FN / (FN + TP)	

	performance <- data.frame(AC, AUC, MCC, SN, SP, PPV, NPV, FPR, FNR,  LRP, LRN, TP, FP, TN, FN)
	return(performance)
}
