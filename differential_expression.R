get.parametric.ttest <- function(dat, class1, class2, var.equal=FALSE){
	cat("COMPUTE DIFFERENTIAL EXPRESSION ...", "\n", sep="")
	#dat <- dat[-which(is.na(dat)),]
	
	#Get class indices
	index1 <- which(colnames(dat) %in% class1)
	index2 <- which(colnames(dat) %in% class2)
	
	# Check if the data points are constant ----
	exp <- dat[,c(index1, index2)]
	len <- apply(exp, 1, function(x) length(unique(as.numeric(x))))
	del.index <- which(len == 1)

	if(length(del.index) != 0){
		dat <- dat[-del.index,]
	}

	# Check if each group has at least 2 samples represented ---
	length_not_na.A <- as.numeric(apply(dat[,index1], 1, function(x) length(which(!is.na(as.numeric(x))))))
	length_not_na.B <- as.numeric(apply(dat[,index2], 1, function(x) length(which(!is.na(as.numeric(x))))))
	del.index_A <- which(length_not_na.A < 2)
	del.index_B <- which(length_not_na.B < 2)

	del.index_AB <- unique(c(del.index_A, del.index_B))
	
	if(length(del.index_AB) != 0){
		dat <- dat[-del.index_AB,]
	}
		
	# Compute Median and StdDev for the Sample Class
	median1 <- apply(as.matrix(dat)[,index1], 1, median, na.rm=TRUE)
	median2 <- apply(as.matrix(dat)[,index2], 1, median, na.rm=TRUE)
	cat("MEDIAN COMPUTED ...", "\n", sep="")
	
	stdev1 <- apply(as.matrix(dat)[,index1], 1, sd, na.rm=TRUE)
	stdev2 <- apply(as.matrix(dat)[,index2], 1, sd, na.rm=TRUE)
	cat("STANDARD DEVIATION COMPUTED ...", "\n", sep="")
	
	#Compute FoldChange
	foldchange <- ifelse(median2 > median1, 2^median2/2^median1, -2^median1/2^median2)
	cat("FOLD CHANGE COMPUTED ...", "\n", sep="")
	
	cat("PERFORMING STUDENT'S T-TEST ...", "\n", sep="")
	# Un-paired Student's t-test between class1 and class2 samples
	func.ttest <- function(x, index1, index2){
		ttestGene <- t.test(x[index1], x[index2], var.equal=var.equal, alternative = "two.sided", paired = FALSE, na.action=na.omit)
		return(ttestGene)
	}

	ttestGene.list <- apply(as.matrix(dat), 1, function(x) func.ttest(x, index1, index2))
	cat("T-TEST COMPUTED ...", "\n", sep="")
	
	pvalue <- do.call(c,lapply(ttestGene.list,function(x) x$p.value))

	index0 <- which(pvalue == 0)
	if(length(index0) != 0){
		pvalue[index0] <- .Machine$double.eps #smallest value
	}

	fdr <- p.adjust(pvalue, method = "BH", n = length(pvalue))
	cat("P-VALUE COMPUTED ...", "\n", sep="")

	dat.summary <- data.frame(Gene=rownames(dat), 
								MedianA=median1, StdevA=stdev1, 
								MedianB=median2, StdevB=stdev2, 
								FoldChange=foldchange, 
								pvalue=pvalue, fdr=fdr)
	rownames(dat.summary) <- c(1:nrow(dat))
	dat.summary$Gene <- as.character(dat.summary$Gene)
	dat.summary <- dat.summary[order(dat.summary$pvalue, decreasing=F),]
	cat("DATA SUMMARY COMPUTED ...", "\n", sep="")
		
	return(dat.summary)
}

get.diffexpr.SN <- function(expr){
	cat("COMPUTE DIFFERENTIAL EXPRESSION ...", "\n", sep="")
	
	classLabel <- as.character(expr[1,])
	expr.matrix <- data.matrix(expr[-1,])

	expr_subset <- expr.matrix

	#Compute Median and StdDev for the Sample Class
	median1 <- apply(expr_subset[,classLabel==unique(classLabel)[1]], 1, median)
	median2 <- apply(expr_subset[,classLabel==unique(classLabel)[2]], 1, median)
	cat("MEDIAN COMPUTED ...", "\n", sep="")
	
	sd1 <- apply(expr_subset[,classLabel==unique(classLabel)[1]], 1, sd)
	sd2 <- apply(expr_subset[,classLabel==unique(classLabel)[2]], 1, sd)
	cat("STANDARD DEVIATION COMPUTED ...", "\n", sep="")
	
	fc <- ifelse(median2 > median1, 2^median2/2^median1, -2^median1/2^median2)
	cat("FOLD CHANGE COMPUTED ...", "\n", sep="")
	
	cat("PERFORMING STUDENT'S T-TEST ...", "\n", sep="")
	ttestFun <- function(x, labels){
		result <- t.test(x~as.factor(labels), na.rm=T, paired=F, var.equal=F, , na.action=na.omit)
		return(result$p.value)
	}
	pVal <- apply(expr_subset, 1, ttestFun, classLabel)
	cat("P-VALUE COMPUTED ...", "\n", sep="")
	
	fdr <- p.adjust(pVal, method="BH")

	expr_summary <- data.frame(Gene=rownames(expr_subset), 
								MedianA=median1, StdevA=sd1, 
								MedianB=median2, StdevB=sd2, 
								FoldChange=fc, 
								pvalue=pVal, fdr=fdr)
	cat("DATA SUMMARY COMPUTED ...", "\n", sep="")
	
	return(expr_summary)
}


### SHAPIROTEST FOR NORMALITY ----
get.shapirotest <- function(dat){

	func.shapirotest <- function(x){
		st <- shapiro.test(x)
		return(st)
	}

	shapirotest.list <- apply(as.matrix(dat), 1, function(x) func.shapirotest(x))
	pvalue <- do.call(c,lapply(shapirotest.list, function(x) x$p.value))
	statistic <- do.call(c,lapply(shapirotest.list, function(x) x$statistic))
	fdr <- p.adjust(pvalue, method = "BH", n = length(pvalue))

	dat.summary <- data.frame(Gene=rownames(dat), 
								statistic=statistic,
								pvalue=pvalue, 
								fdr=fdr)	
	rownames(dat.summary) <- c(1:nrow(dat))
	dat.summary$Gene <- as.character(dat.summary$Gene)
	dat.summary <- dat.summary[order(dat.summary$pvalue, decreasing=F),]

	dat.summary$NormalityStatus <- ifelse(dat.summary$fdr > 0.1, 1, 0)	
	return(dat.summary)
}


### WILCOXN RANK SUM TEST ----
get.wilcox.rank.test <- function(dat, class1, class2){
	cat("COMPUTE WILCOX RANK TEST ...", "\n", sep="")
	#dat <- dat[-which(is.na(dat)),]
	
	#Get class indices
	index1 <- which(colnames(dat) %in% class1)
	index2 <- which(colnames(dat) %in% class2)
	
	# Check if the data points are constant ----
	exp <- dat[,c(index1, index2)]
	len <- apply(exp, 1, function(x) length(unique(as.numeric(x))))
	del.index <- which(len == 1)

	if(length(del.index) != 0){
		dat <- dat[-del.index,]
	}
	
	# Compute Median and StdDev for the Sample Class
	median1 <- apply(as.matrix(dat)[,index1], 1, median, na.rm=TRUE)
	median2 <- apply(as.matrix(dat)[,index2], 1, median, na.rm=TRUE)
	cat("MEDIAN COMPUTED ...", "\n", sep="")
	
	stdev1 <- apply(as.matrix(dat)[,index1], 1, sd, na.rm=TRUE)
	stdev2 <- apply(as.matrix(dat)[,index2], 1, sd, na.rm=TRUE)
	cat("STANDARD DEVIATION COMPUTED ...", "\n", sep="")
	
	#Compute FoldChange
	#foldchange <- ifelse(median2 > median1, 2^median2/2^median1, -2^median1/2^median2)
	#cat("FOLD CHANGE COMPUTED ...", "\n", sep="")
	
	cat("PERFORMING WILCOX RANK TEST ...", "\n", sep="")
	#Wilcox.rank.test
	func.wilcoxtest <- function(x, index1, index2){
		wilcoxtest <- wilcox.test(x[index1], x[index2], alternative="two.sided", na.action=na.omit)
		return(wilcoxtest)
	}

	wilcoxtest.list <- apply(as.matrix(dat), 1, function(x) func.wilcoxtest(x, index1, index2))
	cat("WILCOX RANK TEST COMPUTED ...", "\n", sep="")
	
	pvalue <- do.call(c,lapply(wilcoxtest.list,function(x) x$p.value))

	index0 <- which(pvalue == 0)
	if(length(index0) != 0){
		pvalue[index0] <- .Machine$double.eps #smallest value
	}

	fdr <- p.adjust(pvalue, method = "BH", n = length(pvalue))
	cat("P-VALUE COMPUTED ...", "\n", sep="")

	dat.summary <- data.frame(Gene=rownames(dat), 
								MedianA=median1, StdevA=stdev1, 
								MedianB=median2, StdevB=stdev2, 
								#FoldChange=foldchange, 
								pvalue=pvalue, fdr=fdr)
	rownames(dat.summary) <- c(1:nrow(dat))
	dat.summary$Gene <- as.character(dat.summary$Gene)
	dat.summary <- dat.summary[order(dat.summary$pvalue, decreasing=F),]
	cat("DATA SUMMARY COMPUTED ...", "\n", sep="")
		
	return(dat.summary)
}
