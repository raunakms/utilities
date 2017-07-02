#####################################################
# visualizeBatch
# - batch: numeric vector assigning each sample to a batch corresponding to a number (i.e. 1,2)
# - expr: matrix of sample exprs
# expr: expr matrix
# outcome: class labels
# dates: CEL run dates
# Note:
# - assumes that samples in outcome/batch vector are same as samples in expression matrix
#####################################################

### FUNCTION:  visualizeBatch -----
visualizeBatch <- function(expr, outcome, batch=NULL, batchIntervals=NULL, outDir) {
	# LOAD LIBRARIES 
	require("corpcor")
	require("affy")

	if (!is.matrix(expr))
		stop("Expression file needs to be a data matrix")
		
	# batch=as.numeric(cut(as.numeric(batch),batchIntervals)

	# Hierarchical Cluster
	##compute distance between and perform clustering 
	mydist <- dist(t(expr), method = "euclidean", diag = TRUE, upper = TRUE)

	hc <- hclust(mydist, method="ward.D2")

	##make cluster. show outcome in text, batch in color
	pdf(file.path(outDir, "PC1.pdf"))
		boxplot(split(s$v[,1],batch), main="PC1")
	dev.off()
	

	pdf(file.path(outDir, "PC2.pdf"))
		boxplot(split(s$v[,2],batch), main="PC2")
	dev.off()
	
	###and confounding between batch and outcome
	table(outcome,batch)
	
	# Hierarchical Clustering of PC1
	mydist <- dist(as.matrix(s$v[,1]))
	hc=hclust(mydist)
	pdf(file.path(outDir, "HClust_PC1.pdf"), width=30, height=6)
		myplclust(hc,lab=outcome,lab.col=batch, main="HClust_PC1")
	dev.off()
	
	# Hierarchical Clustering of PC2
	mydist <- dist(as.matrix(s$v[,2]))
	hc=hclust(mydist)

	pdf(file.path(outDir, "HClust_PC2.pdf"), width=10, height=6)
		myplclust(hc,lab=outcome,lab.col=batch, main="HClust_PC2")
	dev.off()
	
	# Hierarchical Clustering of PC3
	mydist <- dist(as.matrix(s$v[,3]))
	hc=hclust(mydist)	

	pdf(file.path(outDir, "HClust_PC3.pdf"), width=10, height=6)
		myplclust(hc,lab=outcome,lab.col=batch, main="HClust_PC3")
	dev.off()	
	
	# Hierarchical Clustering of PC4
	mydist <- dist(as.matrix(s$v[,4]))
	hc=hclust(mydist)	

	pdf(file.path(outDir, "HClust_PC4.pdf"), width=10, height=6)
		myplclust(hc,lab=outcome,lab.col=batch, main="HClust_PC4")
	dev.off()	
	
	# Hierarchical Clustering of Top 2 PC
	mydist <- dist(as.matrix(s$v[,1:2]))
	hc=hclust(mydist)	

	pdf(file.path(outDir, "HClust_Top2PC.pdf"), width=10, height=6)
		myplclust(hc,lab=outcome,lab.col=batch, main="HClust_Top2PC")
	dev.off()	
	
	# Hierarchical Clustering of Top 4 PC
	mydist <- dist(as.matrix(s$v[,1:10]))
	hc=hclust(mydist)	

	pdf(file.path(outDir, "HClust_Top10PC.pdf"), width=10, height=6)
		myplclust(hc,lab=outcome,lab.col=batch, main="HClust_Top10PC")
	dev.off()	
}

#####################################################
## myplclust
## modifiction of plclust for plotting hclust objects *in colour*!
## Arguments:
##    hclust:    hclust object
##    lab:        a character vector of labels of the leaves of the tree
##    lab.col:    colour for the labels; NA=default device foreground colour
##    hang:     as in hclust & plclust
## Side effect:
##    A display of hierarchical cluster with coloured leaf labels.
#####################################################

myplclust <- function(hclust, lab=hclust$labels, lab.col=rep(1,length(hclust$labels)), hang=0.1,...){
	# LOAD LIBRARIES ---
	require("RColorBrewer")

	# DEFINE COLOR PALETTE ---
	jColFun <- colorRampPalette(brewer.pal(n = 9, "Set1"))
	
	# GET COLORS ---
	print(lab.col)
	cols <- jColFun(length(unique(lab.col)))
	#cols <- rainbow(length(unique(lab.col)))
 
 	# GET COMPONENTS ---
	y <- rep(hclust$height,2)
	x <- as.numeric(hclust$merge)
	y <- y[which(x<0)]
	x <- x[which(x<0)]
	x <- abs(x)
	y <- y[order(x)]
	x <- x[order(x)]
	plot( hclust, labels=FALSE, hang=hang, ... )
	
	#text( x=x, y=y[hclust$order]-(max(hclust$height)*hang), labels=lab[hclust$order], col=lab.col[hclust$order], srt=90, adj=c(1,0.5), xpd=NA, ... )
	text(x=x, y=y[hclust$order]-(max(hclust$height)*hang), labels=lab[hclust$order], col=cols[lab.col[hclust$order]], srt=90, adj=c(1,0.5), xpd=NA, cex=0.5,...)
}

