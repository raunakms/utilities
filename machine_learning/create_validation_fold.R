########################################################### ------------
# createValidationFolds_Index(sample.class, file.output, nFolds, nReps) 
#
# Description:
# 	For 'i' repetitions of n-fold CV, generate the sample indices for the training dataset and write to file.
#
# Arguments:
#	file.class - vector of sample class labels, each sample label in a line
# 	file.output - full name of file that training indices will be written to
#	nFolds - # of folds in CV
#	nReps - # of times to repeat CV 
########################################################### ------------
createValidationFolds_Index <- function(sample.class, nFolds, nReps){
	require("caret")
	
	cat(paste(Sys.time()), "CREATING VALIDATION FOLDS ...", "\n", sep=" ")
	
	all.y <- as.factor(sample.class)

	folds <- createMultiFolds(all.y, k=nFolds, times=nReps)
	
	# Find max size of all training set. Must be done to create matrix of fixed size.
	maxLength <- 0
	for (n in (1:length(folds))) {
		if (length(unlist(folds[n])) > maxLength) {
			maxLength <- length(unlist(folds[n]))
		}
	}
	
	indices <- matrix(nrow=nReps*nFolds,ncol=maxLength)
	for (x in 1:nReps) {
		for (i in 1:nFolds) {
			index <- as.vector(t(unlist(folds[(x-1)*nFolds+i])))
			for (n in 1:length(index)) {
				indices[((x-1)*nFolds+i), n] <- index[n]
			}
		}
	}
	
	#write.table(indices, file.output, sep="\t", row.names=F, col.names=F, quote=F, eol="\n")

	#cat("FILE CREATED:", file.output, "\n", sep=" ")

	cat(paste(Sys.time()), "DONE ...", "\n", sep=" ")
	return(indices)
}

