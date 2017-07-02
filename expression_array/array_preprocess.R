### GET UNIQUE GENE EXPRESSION ---
getunique.gene.expression <- function(expr){
	unique.genes <- sort(unique(as.character(expr[,1])), decreasing=F)
	unique.genes <- unique.genes[unique.genes != ""]
	expr <- subset(expr, expr[,1] != "")
	
	rown <- unique.genes
	coln <- colnames(expr[2:ncol(expr)])
	gene.expr <- matrix(nrow=length(unique.genes), ncol=ncol(expr)-1, dimnames=list(rown, coln)) 

	ids <- as.character(expr[,1])
	for(i in 1:length(unique.genes)){
		gene <- unique.genes[i]
		index <- which(ids == gene)
		expr.sub <- as.matrix(expr[index,-1])
		
		if(nrow(expr.sub) == 1){
			gene.expr[i,] <- as.numeric(expr.sub)
		} else {
			gene.expr[i,] <- apply(expr.sub, 2, mean, na.rm=TRUE)
		}	
		
		cat("GENES PROCESSED :", i, "OF", length(unique.genes), "\n", sep="\t")
	}
	return(gene.expr)
}


### QUANTILE NORMALIZE ---
library(aroma.light)
getQuantile.normalize <- function(dat){
	dQNorm <- data.frame(normalizeQuantile(as.matrix(dat)))
	return(dQNorm)
}	

### COMPUTE Z-SCORE ---
#Note: Row=Gene, Column=Sample
getZscore <- function(dat){
	z.dat <- matrix(nrow=nrow(dat), ncol=ncol(dat), dimnames=list(rownames(dat),colnames(dat)))

	dat.mean <- apply(dat, 1, mean)
	#cat("MEAN COMPUTED","\n", sep="\t")
	
	dat.stdev <- apply(dat, 1, sd)
	#cat("STANDARD DEVIATION COMPUTED","\n", sep="\t")

	#cat("COMPUTING Z-SCORE","\n", sep="\t")
	for(i in 1:nrow(dat)){
		x <- as.numeric(dat[i,])
		z.dat[i,] <- (x - dat.mean[i])/dat.stdev[i]
		#cat("Z-SCORE COMPUTED FOR ROW :", i, "OF", nrow(dat), "\n", sep="\t")
	}
	z.dat <- as.data.frame(z.dat)
	return(z.dat)
}

### CONVERT EXPRESSION MATRIX TO GCT FORMAT ----
get.gct <- function(dat, file.gct){
	names.row <- rownames(dat)
	num.row <- nrow(dat)
	num.col <- ncol(dat)
	description <- rep("", nrow(dat))
	
	names.header <- c("Name","Description", colnames(dat))
	
	dat1 <- cbind(names.row, description, dat)
	rownames(dat1) <- c(1:nrow(dat1))
	colnames(dat1) <- c(1:ncol(dat1))
	
	#dat.head <- matrix("", nrow=3, ncol=length(names.header))
	dat.head1 <- dat.head2 <- rep("", length(names.header))
	dat.head1[1] <- "#1.2"
	dat.head2[1] <- num.row
	dat.head2[2] <- num.col
	dat.head3 <- names.header
	#colnames(dat.head) <- c(1:ncol(dat1))
	
	file.temp1 <- tempfile(pattern = "file1", tmpdir = tempdir())
	file.temp2 <- tempfile(pattern = "file1", tmpdir = tempdir())
	
	dat.head <- rbind(dat.head1, dat.head2, dat.head3)
	write.table(dat.head, file.temp1, sep="\t", row.names=F, col.names=F, quote=F)
	write.table(dat1, file.temp2, sep="\t", row.names=F, col.names=F, quote=F)
	
	cmd <- paste("cat", file.temp1, file.temp2, ">", file.gct, sep=" ")
	system(cmd)
	cat("FILE GENERATED: ", file.gct, "\n", sep=" ")
}
