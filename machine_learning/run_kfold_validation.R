### LOAD LIBRARIES ---
library("caret")
library("foreach")
library("doParallel")

### DEFINE PATH ----
dir.function <- file.path("../machine_learning")

### SOURCE FUNCTIONS ---
source(file.path(dir.function, "create_validation_fold.R"))
source(file.path(dir.function, "calculate_performance.R"))

### FUNCTION: RUN k-FOLD VALIDATION ----
run.kfold.validation <- function(univeral.matrix, cls, features, nFolds, nReps, batch.name, dir.output){
	cat("", "\n", sep=" ")
	cat(paste(Sys.time()), "START k-FOLD VALIDATION ...", "\n", sep=" ")

	# CREATE OUTPUT FOLDER ---
	dir.batch <- file.path(dir.output, batch.name)
	dir.create(dir.batch, showWarnings=FALSE)	

	cat(paste(Sys.time()), "OUTPUT PATH:", dir.batch, "\n", sep=" ")

	# CREATE FOLDS ---
	ind.matrix <- createValidationFolds_Index(sample.class=cls, nFolds, nReps) 

	cat(paste(Sys.time()), "TOTAL FEATURE =", length(features), "\n", sep=" ")
	cat(paste(Sys.time()), "START LOOP FOR EACH FEATURE ...", "\n", sep=" ")
	
	# LOOP FOR EACH FEATURE ---
	for(i in 1:length(features)){
		feature <- features[i]

		cat("", "\n", sep=" ")
		cat(paste(Sys.time()), "START FOR FEATURE =", feature, "\n", sep=" ")

		# CALL FUNCTION ---
		list.df <- get.kfold.validation(univeral.matrix, cls, ind.matrix, feature)

		# DEFINE OUTPUT FILE ---
		file.output1 <- file.path(dir.batch, paste("performance_all_", feature, ".tsv", sep=""))
		file.output2 <- file.path(dir.batch, paste("performance_summary_", feature, ".tsv", sep=""))

		# WRITE OUTPUT ---
		write.table(list.df[[1]], file.output1, sep="\t", row.names=F, col.names=T, quote=F)
		write.table(list.df[[2]], file.output2, sep="\t", row.names=F, col.names=T, quote=F)

		cat(paste(Sys.time()), "FEATURE PROCESSED:", feature, i, "OF", length(features), "\n", sep=" ")
	}
	
	cat(paste(Sys.time()), "ALL DONE ...", "\n", sep=" ")
}


### FUNCTION: GET k-FOLD VALIDATION ----
get.kfold.validation <- function(univeral.matrix, cls, ind.matrix, feature){
	cat(paste(Sys.time()), "CREATING TRAINING MATRIX ...", "\n", sep=" ")

	# Prepare Training Data ---
	train.matrix <- data.frame(SampleClass=cls,  Feature=as.numeric(univeral.matrix[feature,]))
	train.matrix$SampleClass <- as.factor(train.matrix$SampleClass)

	cat(paste(Sys.time()), "PARALLELIZING START ...", "\n", sep=" ")	

	# Declate Cluster 
	#no_cores <- detectCores() - 10
	no_cores <- 40
	cat(paste(Sys.time()), "No. OF CORES IN USE =", no_cores, "\n", sep=" ")	

	cl <- makeCluster(no_cores)
	registerDoParallel(cl)

	# PARALLELIZATION ---
	.func <- c("calculatePerformance", "calcMCC", "calcAC")
	performance <- foreach(ctr=1:nrow(ind.matrix), .combine='rbind', .packages="caret", .export=.func)  %dopar% {
						#set.seed(400)

						# Prepare Testing Data ---
						ind.test <- as.numeric(na.omit(ind.matrix[ctr,]))
						test.matrix <- data.frame(SampleClass=cls[ind.test], Feature=as.numeric(univeral.matrix[feature,ind.test]))
						del.ind <- which(is.na(test.matrix$Feature))
							
						if(length(del.ind) != 0){
							test.matrix <- test.matrix[-del.ind,]
						}

						test.matrix$SampleClass <- as.factor(test.matrix$SampleClass)
					
						# Train Classifier ---
						fitControl <- trainControl(method="LGOCV", repeats=10, p=0.75, search="grid") 
						knnFit <- train(SampleClass ~ ., data=train.matrix, 
										method="knn", trControl=fitControl, tuneGrid=expand.grid(.k=3),
										preProcess=NULL, weights = NULL, metric="Accuracy", na.action=na.omit)

						# Predict classifier in Test data ---
						knnPredict <- predict(knnFit, newdata=test.matrix)
						predictions <- data.frame(TrueClass=test.matrix$SampleClass, PredClass=knnPredict)

						# Compute Performance Metric---
						performance <- calculatePerformance(data=predictions)
					}
	stopCluster(cl)	

	cat(paste(Sys.time()), "PARALLELIZING END ...", "\n", sep=" ")	

	cat(paste(Sys.time()), "PREPARING RESULT SUMMARY ...", "\n", sep=" ")
	is.na(performance) <- do.call(cbind,lapply(performance, is.infinite))

	df.performance <- cbind(Feature=feature, Model="knn3", SubSampleIter=c(1:nrow(ind.matrix)), performance)
		
	### Find Confidence Interval
	p <- performance[,1:11]
	avg <- apply(p, 2, mean, na.rm=TRUE)
	sd <- apply(p, 2, sd, na.rm=TRUE)

	error <- qnorm(0.975)*sd/sqrt(nrow(p))
	left <- avg - error
	right <- avg + error

	df.metric.summary <- data.frame(Feature=feature, Model="knn3", Metric=colnames(p), Mean=avg, StdDev=sd, lower.CI=left, upper.CI=right)

	# COMPILE RESULTS ---
	#all.performance <- do.call(rbind.data.frame, df.performance)
	#all.metric.summary <- do.call(rbind.data.frame, df.metric.summary)

	list.output <- list(df.performance, df.metric.summary)

	cat(paste(Sys.time()), "DONE ...", "\n", sep=" ")

	return(list.output)
}
