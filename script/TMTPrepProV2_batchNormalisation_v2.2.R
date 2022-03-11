

###############
# TMTPrePro
# Developed around 2017 by Dana Pascovici and Jemma Wu 
# while at the Australian Proteome Analysis Facility, Macquarie University
# Not actively developed, provided "as is" with no warranty whatsoever
###############



setwd("C:/testTMTPrePro/results")

source("../script/functions.R")

	##############################################
	# main function for overall and targeted analysis
	##############################################

	library(heatmap3)
	library(openxlsx)
	library(scatterplot3d) 
	library(limma)
	library(lattice)
	library(tools)

	print(sessionInfo())
	
	############
	# parameters
	############

	
	zipfname = "../data/run12.zip"
	designfname = "../data/mehdiDesign12.xlsx"
	FCCutoff = 1.2
	PvalCutoff = 0.05
	CountsCutoff = 0
	Clean <- TRUE
	MasterFilter = 'IsMasterProtein' # 'Master Protein' #
	KeepREF = TRUE
	SampleLoadNorm = "total" # "median"
		
	
	args <- list(...)

	for(i in 1:length(args)) {
		flag <- substring(args[[i]], 0, 2)
		value <- substring(args[[i]], 3, nchar(args[[i]]))

		if(flag=='-f') zipfname <- value;
		if(flag=='-d') designfname <- value;
		if(flag=='-r') FCCutoff <- as.numeric(value);
		if(flag=='-c') CountsCutoff <- as.numeric(value);
		if(flag=='-p') PvalCutoff <- as.numeric(value);
		if(flag=='-l') Clean <- ifelse(tolower(value)=='yes', TRUE, FALSE)
		if(flag=='-m') MasterFilter <- value;
		if(flag=='-k') KeepREF <- value;
		if(flag=='-s') SampleLoadNorm <- value;

	} 

	dat.para = data.frame(Parameter=c('P value', 'Fold change', 'Clean', 'FileName', 'DesignFile', 'Date'), 
			Cutoff=c(PvalCutoff, FCCutoff, as.character(Clean), basename(zipfname), 
			basename(designfname), date() ))	
		

	# unzip files and check zip
	files <- unzip(zipfname, junkpaths=TRUE)
	if (length(files) < 1) stop("The zip file is empty")
		
	# Error checking, have two or more tabs etc?	

	designSheets = try(loadWorkbook(designfname) )
	if(inherits(designSheets, 'try-error')) stop("Error with loading design file.")

	if(length(names(designSheets)) < 2) stop("Design file needs to have at least 2 tab.")

	designAll = try( readWorkbook(designfname, 1) )
	if (inherits(designAll, "try-error")) stop("Could not extract the first tab from the design file.")


	 references = readWorkbook(designfname, 2)
	 if (inherits(references, "try-error")) stop("Could not extract the second tab from the design file.")

	
	files = readWorkbook(designfname, 2)[,1]
	refLabels = references[,2]
	
	if (length(unique(refLabels)) > 1) cat("Warning: the references belong to different groups; data interpretation may not make sense.");
	
	# Long format of design
	NRuns = length(grep("Group", colnames(designAll)))
	
	designLong = designAll[,c(1:3)]	
	designLong$Run = "R1"
	
	if(NRuns > 1) 
		for(ii in 2:NRuns) {
			designii = designAll[, c(1,ii*2,ii*2+1)]	
			designii$Run = paste0("R",ii)
			designLong = rbind(designLong, designii)
		}
	
	designLong$Label = gsub("\\_", "", designLong$Label)
  refLabels = gsub("\\_", "", refLabels)
	
	refInRuns = paste(paste0("R", 1:NRuns), refLabels, sep=".")
	
	
	# Read in files
	file.list = lapply(files, FUN=read.delim, as.is=T)

	file.list.clean = file.list
	
	# Cleaning
	
	file.list.clean = lapply(file.list, function(x) {dat1 = x; dat1 = dat1[!duplicated(dat1$Accession),];
		
		# data check & cleaning
		
		if(!'Master' %in% colnames(dat1) || !'Protein.FDR.Confidence' %in% colnames(dat1) ) 
			stop('Must have Master and Protein.FDR.Confidence columns');
		
		if(Clean==TRUE) dat1 = dat1[dat1$Master == MasterFilter & dat1$Protein.FDR.Confidence == "High",];
		dat1
		} )
	
	# Extract abundance, accession and count
	
	list.dat = lapply(1:length(file.list.clean), 
		function(ii) {x=file.list.clean[[ii]];
				col.abun = grepl("Abundance.", colnames(x)) & !grepl("Ratio", colnames(x)) & !grepl("Count", colnames(x));		
				abun1 = x[,col.abun]; 
				names(abun1) = gsub("(\\.)+", "\\.", names(abun1))
				names(abun1) = gsub("\\.$", "", names(abun1))

				# remove F? in the ratio names -- added 8 May 2018
				# names(abun1) = gsub("(F\\d+\\.)", "", names(abun1))
		
				colnames(abun1) = paste0("R",ii, colnames(abun1));
				
				counts.idx = grep("X..Peptides", colnames(x));		
				Counts1 = x[,counts.idx][,1]; 
		
				data.frame(x[,"Accession"], Counts1, abun1)} )
	
	# Merge abundance with counts and description	
	mg.dat = list.dat[[1]]
	if(length(list.dat) > 1)
	for(ii in 2:length(list.dat)) mg.dat = merge(mg.dat, list.dat[[ii]], by=1, all=T, sort=FALSE)
	
	# Fill in descriptions from files
	acc.descp.map = file.list.clean[[1]][,c("Accession","Description")]
	
	if(length(list.dat) > 1)
	for(ii in 2:length(file.list.clean)) acc.descp.map = rbind(acc.descp.map, file.list.clean[[ii]][,c("Accession","Description")])
	acc.descp.map = acc.descp.map[!duplicated(acc.descp.map$Accession),]
	
	mg.dat$Description = acc.descp.map[match(mg.dat[,1], acc.descp.map[,1]), 2]
	
	# Remove dots in col names
	colnames(mg.dat) = gsub("^x\\.+", "", colnames(mg.dat) )
	colnames(mg.dat) = gsub("\\.+$", "", colnames(mg.dat) )
	
	# Reorder columns
	colfirst = c("Accession", "Description", colnames(mg.dat)[grep("Counts", colnames(mg.dat))])
	mg.dat = mg.dat[,c(which(colnames(mg.dat)%in%colfirst), which(!colnames(mg.dat)%in%colfirst))]
				
	mg.mat = mg.dat[,-c(1:length(colfirst))]	
	rownames(mg.mat) = mg.dat[,"Accession"]

#
# Naming bit different each time
#
	# colnames(mg.mat) = paste(gsub("^(R\\d+).*", "\\1", colnames(mg.mat)),
				# gsub("^R.(.*)\\.(.*)","\\2", colnames(mg.mat)), sep="." )
				
	colnames(mg.mat) = paste(gsub("^(R\\d+).*", "\\1", colnames(mg.mat)),
				gsub("^R.(.*)\\.(.*).Sample","\\2", colnames(mg.mat)), sep="." )
				
			
	
	mg.dat = mg.dat[,colfirst]

	# if(KeepREF == FALSE) mg.mat = mg.mat[,!colnames(mg.mat) %in% refInRuns]
	
	# Remove rows with all NAs
	idx.allna = apply(mg.mat, 1, function(x) sum(is.na(x))==ncol(mg.mat))
	mg.mat = mg.mat[!idx.allna,]
	
	
	data_raw = mg.mat #data_no_na
		
	stopifnot(nrow(data_raw) > 0)

	Group = as.factor(designLong[match(gsub("\\.", "", colnames(data_raw)),
				paste0(designLong$Run, designLong$Label)),"Group"] )
				
				
		
	Replicate = designLong[match(gsub("\\.", "", colnames(data_raw)),
				paste0(designLong$Run, designLong$Label)),"Replicate"] 
	
	Replicate = as.factor(gsub("\\_\\d+", "", Replicate))
	
	keep.idx = Group != "X"
	
	
	Group = Group[keep.idx, drop=TRUE]
	mg.mat = mg.mat[,keep.idx]
	data_raw = data_raw[,keep.idx]
	Replicate = Replicate[keep.idx]
	
	
	
	Cols = rainbow(nlevels(Group))
	
	#################################
	# Sample loading normalistation
	#################################

	# Sample loading median or total normalisation
	if(tolower(SampleLoadNorm) == "total") {
		tot = apply(data_raw, 2, function(x) sum(na.omit(x)))
	} else {
		tot = apply(data_raw, 2, function(x) median(na.omit(x)))
	}

	
	data_sl = sweep(data_raw, 2, tot/max(tot), "/")

	format(round(colSums(na.omit(data_sl)), digits = 0), big.mark = ",")
	
	
	# see what the SL normalized data look like
	png(paste("Boxplot raw and", SampleLoadNorm, "norm.png"), 2000, 2000, res=300)
	layout(matrix(1:4, ncol=2))	
	boxplot(log2(data_raw), col = Cols[Group], 
		notch = TRUE, main = "Plain", las=2, cex.names=.8)
	plotDensities(log2(data_raw), col = Cols[Group], main = "Plain", legend=FALSE)
	
	boxplot(log2(data_sl), col = Cols[Group], 
			notch = TRUE, main = paste(SampleLoadNorm, "normalized data"), las=2, cex.names=.8)
	plotDensities(log2(data_sl), col = Cols[Group], main = paste(SampleLoadNorm, "normalization"), legend=FALSE)
	
	dev.off()

	#########
	# IRS
	########
	
	list.exp = list()
	for(ii in 1:length(list.dat)) list.exp[[ii]]=data_sl[,grep(paste0("R",ii), colnames(data_sl))]
	
	# Calculate the protein sum for each batch
	list.rowsum = lapply(list.exp, function(x) apply(x,1, function(y) (sum(na.omit(y)) ) ) )
	
	irs = as.data.frame(list.rowsum[[1]])
	if(length(list.rowsum) > 1)
	for(ii in 2:length(list.rowsum)) irs = cbind(irs, as.data.frame(list.rowsum[[ii]]) )
	
	# convert 0 to NA	
	irs[irs==0] = NA
	
	colnames(irs) <- paste0("sum", 1:ncol(irs))
	
	rowsum.average <- apply(irs, 1, function(x) exp(mean(log(na.omit(x)))))
		

	# compute the scaling factor vectors
	
	irs.fac = sweep(irs, 1, rowsum.average, "/")
	
	list.irs.scaled = lapply(1:length(list.exp), function(x) sweep(list.exp[[x]], 1, irs.fac[,x], "/") ) 
		
	# make new data frame with normalized data
	data_irs <- list.irs.scaled[[1]]
	if(length(list.irs.scaled) > 1)
	for(ii in 2:length(list.irs.scaled)) data_irs = cbind(data_irs, list.irs.scaled[[ii]])

	# in case of 0
	data_irs = data_irs+1
	
	format(round(colSums(na.omit(data_irs)), digits = 0), big.mark = ",")
	
	#########
	# END IRS
	########

	# remove reference ?
	data_irs_full = data_irs
	if(KeepREF == FALSE) data_irs = data_irs[,!colnames(data_irs) %in% refInRuns]
	
	Group = as.factor(designLong[match(gsub("\\.", "", colnames(data_irs)),
				paste0(designLong$Run, designLong$Label)),"Group"] )
				
	Replicate = designLong[match(gsub("\\.", "", colnames(data_irs)),
				paste0(designLong$Run, designLong$Label)),"Replicate"] 
	
	Replicate = as.factor(gsub("\\_\\d+", "", Replicate))	
	
	Cols = rainbow(nlevels(Group))
	
	####################################################
	# Some overall data look metrics
	# -- correlation
	# -- boxplots and density plots
	# -- within group correlation for each level of the group
	# -- PCA
	# -- heatmap
	# -- Anova, followed by clustering of DE proteins
	# -- Diff exp proteins _to the reference_, e.g. Mod/Control, ..., etc, via 1-sample t-test,
	# -- Also combine z-scores
	# -- Venn diagrams of overlap, and also barplots
	####################################################


	# Correlation heatmap
	png("Correlation heatmap IRS.png", 2000, 2000,res=300)
	heatmap3(cor(log(na.omit(data_irs+.5)), use="pairwise.complete.obs"), distfun=cordist,
		 col=colorRampPalette(c("green","black", "red"))(1024),
		 main="IRS correlation",
		ColSideColors=rainbow(nlevels(Group))[Group], margins=c(10,10))
	legend("topright", fill=Cols[1:nlevels(Group)], legend=levels(Group), xpd=TRUE, cex=.6)
	dev.off()


	png("BoxplotDensity.png", width=3500, height=1700,res=300)
	layout(matrix(1:2, nrow=1))
	par(mar=c(13,4,4,2)+.1)
	plotDensities(log(na.omit(data_irs+.5)), col=Cols[Group], legend=FALSE, main=paste(SampleLoadNorm, "and IRS normalised"))
	legend('topright', fill=Cols[1:nlevels(Group)], legend=levels(Group))
	
	# boxplots and density plots
	boxplot(log(data_irs[, order(Group)]+.5), las=2, col=Cols[Group[order(Group)]], 
		main=paste(SampleLoadNorm, "and IRS normalised"),
		cex.axis=0.6)	

	dev.off()


	cat('Begin correlations for all groups\n')
	


	# Correlations for all levels of the group
	for (lev in levels(Group)) {
	
		dd = na.omit(data_irs[,Group == lev, drop=FALSE])
		if(ncol(dd) > 1) {
		png(paste("Cor", lev, ".png", sep=""), 2000, 2000,res=300)
		pairs(log(dd+.5), lower.panel = panel.smooth, upper.panel = panel.cor, main=lev)
		dev.off()
		}
	}	

	cat('End correlations for all groups\n')
	
	############################################
	# unsupervised analysis: clustering and PCA
	############################################

		cat('Begin unsupervised\n')
		
	png("HeatmapAll.png", 2000, 2000,res=300)
	heatmap3(as.matrix(log(na.omit(data_irs+.5))), distfun=cordist,col=colorRampPalette(c("green","black","red"))(100),
		ColSideColors=rainbow(nlevels(Group))[Group], margins=c(10,10))
	legend("topright", fill=Cols[1:nlevels(Group)], legend=levels(Group),  xpd=TRUE, cex=.6 )
	dev.off()


	#pca.res <- PCA(log(t(na.omit(data_irs+.5))), Group, k=5, scaleC = FALSE)
	pca.res <- PCA(log(t(na.omit(data_irs+.5))), Group, k=5)
	z <- pca.res$componentScores

	ld = pca.res$componentLoadings

	# proportion of variance of the top 3 components
	props = round(100*pca.res$summary$importance[2,1:3], 1)
	
	png("PCA3dPlot.png", 2000, 2000, res=300)
	s3d <- scatterplot3d(z[,1:3], color = Cols[Group], col.axis=gray(0.85),col.grid="lightblue",
		box = T, angle = 26, pch=20 )
	s3d$points3d(z[,1:3], pch=21)
	legend("topright", fill=Cols[1:nlevels(Group)], legend=levels(Group), cex=.8)
	text(s3d$xyz.convert(3+z[,1], 3+z[,2], z[,3]), labels = colnames(data_irs), cex=0.4)
	dev.off()

	ord.list = list()

	png("PCA2DAll.png", 2000, 2000, res=300)
	layout(matrix(1:4, ncol=2))
	plot(z[,1], z[,2], col=Cols[Group], pch=20, xlab=paste0("PC1(", props[1],"%)"), 
		ylab=paste0("PC2(", props[2], "%)"))
	points(z[,1], z[,2], pch=21, cex=1.1, lwd=1.3)
	text(z[,1], z[,2], colnames(data_irs), pos=3, cex=.5)
	
	plot(z[,1], z[,3], col=Cols[Group], pch=20, xlab=paste0("PC1(", props[1],"%)"), 
		ylab=paste0("PC3(", props[3],"%)"))
	points(z[,1], z[,3], pch=21, cex=1.1, lwd=1.3)
	text(z[,1], z[,3], colnames(data_irs), pos=3, cex=.5)
	
	plot(z[,2], z[,3], col=Cols[Group], pch=20, xlab=paste0("PC2(", props[2], "%)"), 
		ylab=paste0("PC3(", props[3],"%)"))
	points(z[,2], z[,3], pch=21, cex=1.1, lwd=1.3)
	text(z[,2], z[,3], colnames(data_irs), pos=3, cex=.5)
	
	plot(z[,2], z[,3], col=Cols[Group], pch=20, xlab="", ylab="", axes=FALSE, type='n')
		
	legend("center", fill=Cols[1:nlevels(Group)], legend=levels(Group))	

	dev.off()
	
	png("PCATopLoadings.png", width=2000, height=700, res=300)
	par(oma=c(2,1,1,1))
	layout(matrix(1:3, nrow=1))

	for (ii in 1:3) {
	 ord = order(abs(ld[,ii]), decreasing=TRUE)[1:5]
	 barplot(sort(ld[ord, ii]), las=2, main=paste("Top loadings PC", ii))
	 ord.list[[ii]]=ord
	}
	dev.off()

	png("PCATopLoadingsProteinPatterns.png", width=2500, height=2500, res=300)
	par(mar=c(5,2,3,1))
	layout(matrix(1:15, nrow=3, byrow=T))
	for (ii in 1:3) {
		ord = ord.list[[ii]]
		for (xx in 1:5) {
			boxplot(as.vector(as.matrix(data_irs[match(rownames(ld)[ord[xx]], rownames(data_irs)),])) ~ Group, 
			boxwex=0.5, main=rownames(ld)[ord[xx]], col="gray", las=2)
		}
	}
	dev.off()


	cat('End unsupervised \n')
	############################################
	# End unsupervised analysis
	############################################


	################
	# ANOVA
	################

		cat('Begin anova\n')
	Anova = rep(NA, nrow(data_irs))

	# compute Group means (in log space, geometric)
	data.ag = aggregate(t(data_irs), by=list(Group=Group), 
			FUN=function(v){exp(mean(log(na.omit(v))))} )
	Means = t( data.ag[,-1])
	colnames(Means) = paste("Means",data.ag[,1])
	MaxFC = apply(Means, 1, FUN=function(v){max(v)/min(v)})


	for (i in 1:nrow(data_irs)) {
		
		v= t(data_irs[i,])
		nna.idx = !is.na(v)

		an.res =  try(anova(lm( log(v[nna.idx]+.5) ~ Group[nna.idx, drop=TRUE] ))[1,"Pr(>F)"])

		if (!inherits(an.res, "try-error")) Anova[i] = an.res;

	}



	Anova.adj  = p.adjust(Anova, method = "fdr")
	Anova.idx =  !is.na(MaxFC) & ( MaxFC > FCCutoff ) & !is.na(Anova) & (Anova < PvalCutoff)


	if(nrow(data_irs[Anova.idx,]) > 3) {
  	png("Heatmap - Anova DE.png", 2000, 2000, res=300)
  	hm1 <- heatmap3(as.matrix(na.omit(log(data_irs[Anova.idx,]+.5))),  margins=c(15,15), cexRow=1,
  		 col=colorRampPalette(c("green", "black", "red"))(120), 
  		ColSideColors=Cols[Group]  )
  	legend("topright", fill=Cols[1:nlevels(Group)], legend=levels(Group),
  	 xpd=TRUE,cex=.6 )
  	dev.off()
	}
		cat('End anova\n')

		cat('Begin clustering\n')
	#cluster.data = na.omit(log(data_irs[Anova.idx,]))
	cluster.data = na.omit(log(data_irs_full[Anova.idx,]))
	NotOmitted = setdiff(1:sum(Anova.idx), attr(cluster.data, "na.action"))

	gp = Group


	Cluster = rep(NA, sum(Anova.idx))
	
	res1 <- try(HClust(cluster.data, metric="pearsonCorrelation", method="complete", cutNumber=4))
	
	if(!inherits(res1, "try-error")) {
		clustID <- res1$clustID

		Cluster[NotOmitted] = clustID

		noClusters <- 4

		# scale cluster data for visulisation only	
		mat.ref = cluster.data[,colnames(cluster.data)%in%refInRuns, drop=FALSE]
			
		list.runs = list()
		for(ii in 1:ncol(mat.ref)) list.runs[[ii]] = cluster.data[,grep(paste0("R",ii), colnames(cluster.data))]
		
		list.scaledruns = lapply(1:length(list.runs), function(x) list.runs[[x]]/mat.ref[,x])
		
		scaled.cluster.data = list.scaledruns[[1]]
		if(length(list.scaledruns) > 1) 
		for(ii in 2:length(list.scaledruns)) scaled.cluster.data = cbind(scaled.cluster.data, list.scaledruns[[ii]])
		
		if(KeepREF == FALSE) scaled.cluster.data = 
			scaled.cluster.data[,!colnames(scaled.cluster.data)%in%refInRuns]
		
		plotClusterProfile(log(scaled.cluster.data), clustID, Group, k=noClusters, ylab="Average log ratio")

	}

	Clusters = rep(NA, nrow(data_irs))
	Clusters[Anova.idx] = Cluster

	# reorder abundance data
	data_irs_orderbygroup = data_irs[,order(Group)]
	colnames(data_irs_orderbygroup) = gsub("\\.", "\\.Abundance.", colnames(data_irs_orderbygroup))
	
	full.res = data.frame(Accession=rownames(data_irs_orderbygroup), data_irs_orderbygroup, 
		Means, MaxFC, Anova, Anova.adj, Clusters)
		
	full.res = merge(mg.dat, full.res, by=1, all.x=FALSE, all.y=TRUE, sort=FALSE)
	
	# Output all samples and corresponding groups in the new order
	write.csv(data.frame(Sample=colnames(data_irs_orderbygroup), Group=Group[order(Group)]), "samplegroup.csv")
	
	
	
	cat('end clustering\n')
	################
	# End ANOVA
	################

	
	# Output results
	
	wb <- createWorkbook()
	full.res[full.res=='NaN'] = NA

	# add parameter cutoffs
	addWorksheet(wb, sheet='Parameter')
	writeData(wb, 'Parameter', dat.para)
		
	
	
	ps <- try(printOpenxlsxStyle(full.res, ratios=grep("MaxFC", names(full.res)), 
					pvals=c(grep("Anova", names(full.res))), 
					wb = wb, 
					tabName = "AllData", hiCutoff=FCCutoff, lowCutoff=1/FCCutoff, pvalCutoff=PvalCutoff) )

	if(inherits(ps, 'try-error') ) warning('Error with print overall xlsx file')
					
	ps <- try(printOpenxlsxStyle(data.frame(rownames(pca.res$componentScores),pca.res$componentScores), ratios=NULL, pvals=NULL, wb = wb, tabName = "PCAScores"))

	if(inherits(ps, 'try-error') ) warning('Error with print overall PCA component scores tab')

	ps <- try(printOpenxlsxStyle(data.frame(rownames(ld), ld), ratios=NULL, pvals=NULL, wb = wb, tabName = "PCALoadings"))

	if(inherits(ps, 'try-error') ) warning('Error with print overall PCA loading tab')


	saveWorkbook(wb, file="ResultsOverall.xlsx", overwrite=TRUE)

	
	
	##############
	# Pairwise 
	##############

	if(sum(grepl("comp", tolower(names(designSheets)) )) == 1 ) {
	
	dat.comparisons = readWorkbook(designfname, 3)

	tarcompres.list = list()

	for(idx.comp in 1:nrow(dat.comparisons) ) {
		comp = dat.comparisons[idx.comp,-1]
		
		idx.group1 = which(Group %in% comp[1,1])

		idx.group2 = which(Group %in% comp[1,2])

		idx.mean1 = which(comp[1,1] == gsub('Means ', '', colnames(Means), fixed=TRUE) ) 
		  
		idx.mean2 = which(comp[1,2] == gsub('Means ', '', colnames(Means), fixed=TRUE) )

		FC = Means[,idx.mean1]/Means[,idx.mean2]
		
		TwoSplTTest = rep(NA, nrow(data_irs))
		
		for(ii in 1:nrow(data_irs)) {	
			v = t(data_irs[ii,])
			
			# TWo sample t test
			temp = try(t.test(log(na.omit(v[idx.group1])), 
				log(na.omit(v[idx.group2])), var.equal=TRUE) )
				
			if(!inherits(temp, "try-error")) TwoSplTTest[ii] = temp$p.value		
		
		}
		
		TwoSplTTestAdj = p.adjust(TwoSplTTest, 'fdr')
		
		idx.sig = (!is.na(TwoSplTTest)) & (!is.na(FC)) & (TwoSplTTest < PvalCutoff) & (abs(log(FC))) > log(FCCutoff)
	
		png(paste('Volcano plot for targeted', dat.comparisons[idx.comp,1], '.png', sep=''), 2000, 2000, res=300)
		plot(log(FC), -log(TwoSplTTest), xlab='log FC', ylab='-log p value', 
					main=paste('Protein volcano plot', dat.comparisons[idx.comp,1]))
		abline(h=-log(PvalCutoff), col="red")
		abline(v=log(FCCutoff), col="blue")
		abline(v=-log(FCCutoff), col="blue")
		
		if(sum(idx.sig) > 0) 			
			points(log(FC[idx.sig]), -log(TwoSplTTest[idx.sig]), col='red', pch=20)

		dev.off()	

		
		dat.annot = mg.dat[match(rownames(data_irs), mg.dat$Accession), c(grep("Description", colnames(mg.dat)), grep("Count", colnames(mg.dat)))]
		
		dat.comp = data.frame(Accession=rownames(data_irs), dat.annot, data_irs[,c(idx.group1, idx.group2)], Means[,c(idx.mean1, idx.mean2)],
					FC, TwoSplTTest, TwoSplTTestAdj, Significant = idx.sig)

		tarcompres.list[[idx.comp]] = list(Result=dat.comp, ComparisonName=dat.comparisons[idx.comp,1] )

	}
	
	wb <- createWorkbook()
	
	# add parameter cutoffs
	addWorksheet(wb, sheet='Parameter')
	writeData(wb, 'Parameter', dat.para)
	
	for(idx.comp in 1:length(tarcompres.list)) {
	
		dd = tarcompres.list[[idx.comp]][[1]]
		dd[dd=='NaN'] = NA
		
		colnames(dd)[1] = 'Accession'
		
		# sort by significant
		dd = dd[order(dd$Significant, decreasing=TRUE), ]
		
		printOpenxlsxStyle(dd, ratios=grep('FC', names(tarcompres.list[[idx.comp]][[1]])),
				pvals=grep('TwoSplTTest', names(tarcompres.list[[idx.comp]][[1]])), wb=wb,
				tabName = tarcompres.list[[idx.comp]]$ComparisonName, hiCutoff=FCCutoff, lowCutoff=1/FCCutoff, pvalCutoff=PvalCutoff)
	}
	
		
	# images
	addWorksheet(wb, sheet='images')	
	startCol = 1
	for(idx.comp in 1:length(tarcompres.list)) {
		if(file.exists(paste('Volcano plot for targeted', dat.comparisons[idx.comp,1], '.png', sep=''))) {			
			insertImage(wb, sheet='images', 
				file=paste('Volcano plot for targeted', dat.comparisons[idx.comp,1], '.png', sep=''), 
				width=12, height=15, startRow=2, startCol=startCol, units='cm')
			
			startCol = startCol + 10
		}
	}
	
	# comparisons
	addWorksheet(wb, sheet='comparisons')
	writeData(wb, sheet='comparisons', dat.comparisons)
	
	saveWorkbook(wb, file="ResultsTargetedPaiwise.xlsx", overwrite=TRUE)
	}
	
	#output design
	saveWorkbook(designSheets, file='Design.xlsx', overwrite=TRUE)
	
	# generate MD5 checksums
		
	flist = list.files()

	md5list = lapply(flist, FUN=md5sum)

	# format •	MD5 checksum
	#        •	two spaces ('  ')
	#        •	the filename

	write.table(data.frame(unlist(md5list), names(unlist(md5list))), sep="  ", quote=F, row.names=F, file = "checksums.txt")
	
		
}	
	


###############
# Functions
###############

#### test cases

if(FALSE) {
  rm(list=ls())
  setwd("\\\\APAF-HPV-FILE\\BioInfo_Project\\Projects\\TMT\\TMTPrepPro_V2.1_Batch_Normalisation\\example 18 - s50008 Philip")
  source("..\\TMTPrepProV2_batchNormalisation_v2.2.R")
  TMTPrepProV2.1_IRS_Norm_GP ("-frawData.zip", "-dDesign_s5008_chm.xlsx", "-r1.2", "-mMaster Protein", "-p0.05", "-lYes", "-kTRUE", "-stotal")
  
}

