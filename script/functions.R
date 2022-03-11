
###############
# Functions
# Developed over time by Dana Pascovici and Jemma Wu 
# while at the Australian Proteome Analysis Facility, Macquarie University
###############

##############################
# function printOpenxlsxStyle
##############################

printOpenxlsxStyle  <- function (dat, ratios, pvals, wb, tabName = "results", hiCutoff = 1.5, lowCutoff=0.67, pvalCutoff=0.05) 
{
  
  addWorksheet(wb, sheet=tabName)
  
  upReg <- createStyle(fgFill = "violet")
  downReg <- createStyle(fgFill = "forestgreen")
  sigStyle <- createStyle(fgFill = "khaki1")
  
  writeData(wb, tabName, dat, keepNA=FALSE)
  
  for (rat in ratios) {
    up.idx <- which(!is.na(dat[, rat]) & (dat[, rat] > hiCutoff))
    if (length(up.idx) > 0) 
      addStyle(wb, tabName, style=upReg, rows = 1 + up.idx, cols = rat)
    
    down.idx <- which(!is.na(dat[, rat]) & (dat[, rat] < 
                                              lowCutoff))
    if (length(down.idx) > 0) 
      addStyle(wb, tabName, style=downReg, rows = 1 + down.idx, cols = rat)
  }
  
  for (pval in pvals) {
    sig.idx <- which(!is.na(dat[, pval]) & (dat[, pval] < 
                                              pvalCutoff))
    if (length(sig.idx) > 0) 
      addStyle(wb, tabName, style=sigStyle, rows = 1 + sig.idx, cols = pval)
  }
}


# Added on 22/08/18 to fix the overlapped labels
plotErrorBarsLines <- function (v, barSizes, lines, labels = NULL, col = "blue", ylim = c(min(lines), 
                                                                                          max(lines)), ...) 
{
  barSizes[is.na(barSizes)] <- 0
  topBars <- v + 0.5 * barSizes
  bottomBars <- v - 0.5 * barSizes
  N <- length(v)
  if (is.null(labels)) 
    labels <- 1:N
  ylims <- c(min(bottomBars, ylim[1], min(lines)), max(topBars, 
                                                       ylim[2], max(lines)))
  par(pch = 19, xaxt = "n")
  plot(as.numeric(labels), v, ylim = ylims, col = col, type = "b", 
       lwd = 3, ...)
  par(xaxt = "s")
  
  for (i in 1:N) {
    lines(c(i, i), c(topBars[i], bottomBars[i]))
  }
  for (i in 1:ncol(lines)) {
    lines(as.numeric(labels), lines[, i], lwd = 0.5, lty = "dotted", 
          col = "gray")
  }
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

make_CVs = function(df, replicate) {
  
  cvs = matrix(NA, nrow=nrow(df), ncol=nlevels(as.factor(replicate) ))
  colnames(cvs) = levels(as.factor(replicate))
  rownames(cvs) = rownames(df)
  
  for(ii in 1:nrow(df))
    cvs[ii,] = aggregate(t(df[ii,]), by=list(replicate), function(x) sd(x)/mean(x))[,2]
  
  return(cvs)
  
}

plotClusterProfile <- function(cluster.data, clustID, group, k=4, ylab="Abundance") {
	
	gp = group
	noClusters <- k

	r.temp <- aggregate(t(cluster.data), by=list(gp=gp), FUN=mean)
	ag.sample <- r.temp[,-1]
	rownames(ag.sample) <- r.temp[,1]
	ag.genes <- aggregate(t(ag.sample), by=list(Cluster=clustID), FUN=mean)
	ag.sd <- aggregate(t(ag.sample), by=list(Cluster=clustID), FUN=sd)
	ag.matrix <- as.matrix(ag.genes[,-1])
	ag.counts <- summary(as.factor(clustID))
	ag.bars <- as.matrix(ag.sd[,-1])
	
	png("ClusterPatterns.png", 2000, 2000, res=300)
	par(bg=gray(.95), fg=gray(0.3), oma= c(5, 2, 2, 1) + 0.1, col.main="black", col.sub="black", col.lab="black", col.axis="black")
	layout(matrix(1:4, ncol=2, byrow=TRUE))
	NSig <- noClusters
	for(i in 1:NSig) {
		cols <- rainbow(4) 
		# cols <- rep("gray", 6)
		gname <-  paste("Cluster", i, "(", ag.counts[i], "proteins )")
		lines <- ag.sample[, clustID==i, drop=FALSE]
		plotErrorBarsLines(ag.matrix[i,], 2*ag.bars[i,], lines, 
			labels=1:ncol(ag.matrix), 
			col=cols[i],  main=gname, # bgcol="gray", split=split,
			ylab=ylab, xlab="",
			ylim=c(min(ag.matrix), max(ag.matrix)))
		axis(1,at=1:ncol(ag.matrix), las=2, labels=colnames(ag.matrix), col="black")
		abline(h=0, lty="dotted")
	}
	
	dev.off()

}




plotDensityNA <- function (data, group = rownames(data), xlab = "Abundance",
		main="Sample Densities", legend=TRUE) 
{
    group <- as.factor(group)
    colours <- rainbow(length(levels(group)))
    col <- colours[group]
    # par(mfrow = c(1, 1))
   

        x <- as.matrix(data)
	  y <- x[!is.na(x)]

        try(d1 <- density(y))
        if (inherits(d1, "try-error")) 
            Error("Failed to generate the density plots")
        ymx <- max(d1$y)
        plot(d1, type = "n", xlab = xlab, ylab = "PDF", main = main)
         #   ylim = c(0, 2 * ymx), yaxp = c(0, 2 * ymx, 5))
        for (i in 1:ncol(x)) {
            try(d1 <- density(x[!is.na(x[,i]), i]))
            if (inherits(d1, "try-error")) 
                Error(paste("Failed to generate the density plot for sample", 
                  i))
            # lines(d1, lty = i, col = col[i])
		lines(d1, col = col[i])

        }
	
	  if (legend) {
	    par(cex=0.5)
        legend("topright", legend = levels(group), col = colours, 
            lty = 1)
         par(cex=1)
		 }

    
}




# Jemma Wu - IRS; needs raw data, not log transformed
IrsNormalisation = function (Mat, Run, method="total") {

	Run = as.factor(Run)
	
	list.exp = list()
	
	for(ii in 1:nlevels(Run)) list.exp[[ii]]=Mat[,Run==levels(Run)[ii] ]
	
	# Calculate the protein sum for each batch
	if(tolower(method) =="total" ) {
		list.rowsum = lapply(list.exp, function(x) apply(x,1, function(y) (sum(na.omit(y)) ) ) )
	} else {
	
		list.rowsum = lapply(list.exp, function(x) apply(x,1, function(y) (median(na.omit(y)) ) ) )
	}
	
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

	data_irs
}



# function printOpenxlsxStyle

printOpenxlsxStyle  <- function (dat, ratios, pvals, wb, tabName = "results", hiCutoff = 1.5, lowCutoff=0.67, pvalCutoff=0.05) 
{
  
  addWorksheet(wb, sheet=tabName)
  
  upReg <- createStyle(fgFill = "violet")
  downReg <- createStyle(fgFill = "forestgreen")
  sigStyle <- createStyle(fgFill = "khaki1")
  
  writeData(wb, tabName, dat, keepNA=FALSE)
  
  for (rat in ratios) {
    up.idx <- which(!is.na(dat[, rat]) & (dat[, rat] > hiCutoff))
    if (length(up.idx) > 0) 
      addStyle(wb, tabName, style=upReg, rows = 1 + up.idx, cols = rat)
    
    down.idx <- which(!is.na(dat[, rat]) & (dat[, rat] < 
                                              lowCutoff))
    if (length(down.idx) > 0) 
      addStyle(wb, tabName, style=downReg, rows = 1 + down.idx, cols = rat)
  }
  
  for (pval in pvals) {
    sig.idx <- which(!is.na(dat[, pval]) & (dat[, pval] < 
                                              pvalCutoff))
    if (length(sig.idx) > 0) 
      addStyle(wb, tabName, style=sigStyle, rows = 1 + sig.idx, cols = pval)
  }
}


# Added on 22/08/18 to fix the overlapped labels
plotErrorBarsLines <- function (v, barSizes, lines, labels = NULL, col = "blue", ylim = c(min(lines), 
                                                                                          max(lines)), ...) 
{
  barSizes[is.na(barSizes)] <- 0
  topBars <- v + 0.5 * barSizes
  bottomBars <- v - 0.5 * barSizes
  N <- length(v)
  if (is.null(labels)) 
    labels <- 1:N
  ylims <- c(min(bottomBars, ylim[1], min(lines)), max(topBars, 
                                                       ylim[2], max(lines)))
  par(pch = 19, xaxt = "n")
  plot(as.numeric(labels), v, ylim = ylims, col = col, type = "b", 
       lwd = 3, ...)
  par(xaxt = "s")
  
  for (i in 1:N) {
    lines(c(i, i), c(topBars[i], bottomBars[i]))
  }
  for (i in 1:ncol(lines)) {
    lines(as.numeric(labels), lines[, i], lwd = 0.5, lty = "dotted", 
          col = "gray")
  }
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

make_CVs = function(df, replicate) {
  
  cvs = matrix(NA, nrow=nrow(df), ncol=nlevels(as.factor(replicate) ))
  colnames(cvs) = levels(as.factor(replicate))
  rownames(cvs) = rownames(df)
  
  for(ii in 1:nrow(df))
    cvs[ii,] = aggregate(t(df[ii,]), by=list(replicate), function(x) sd(x)/mean(x))[,2]
  
  return(cvs)
  
}

plotClusterProfile <- function(cluster.data, clustID, group, k=4, ylab="Abundance") {
	
	gp = group
	noClusters <- k

	r.temp <- aggregate(t(cluster.data), by=list(gp=gp), FUN=mean)
	ag.sample <- r.temp[,-1]
	rownames(ag.sample) <- r.temp[,1]
	ag.genes <- aggregate(t(ag.sample), by=list(Cluster=clustID), FUN=mean)
	ag.sd <- aggregate(t(ag.sample), by=list(Cluster=clustID), FUN=sd)
	ag.matrix <- as.matrix(ag.genes[,-1])
	ag.counts <- summary(as.factor(clustID))
	ag.bars <- as.matrix(ag.sd[,-1])
	
	png("ClusterPatterns.png", 2000, 2000, res=300)
	par(bg=gray(.95), fg=gray(0.3), oma= c(5, 2, 2, 1) + 0.1, col.main="black", col.sub="black", col.lab="black", col.axis="black")
	layout(matrix(1:4, ncol=2, byrow=TRUE))
	NSig <- noClusters
	for(i in 1:NSig) {
		cols <- rainbow(4) 
		# cols <- rep("gray", 6)
		gname <-  paste("Cluster", i, "(", ag.counts[i], "proteins )")
		lines <- ag.sample[, clustID==i, drop=FALSE]
		plotErrorBarsLines(ag.matrix[i,], 2*ag.bars[i,], lines, 
			labels=1:ncol(ag.matrix), 
			col=cols[i],  main=gname, # bgcol="gray", split=split,
			ylab=ylab, xlab="",
			ylim=c(min(ag.matrix), max(ag.matrix)))
		axis(1,at=1:ncol(ag.matrix), las=2, labels=colnames(ag.matrix), col="black")
		abline(h=0, lty="dotted")
	}
	
	dev.off()

}


PCA <-
function (data, labelValue, scaleR = FALSE, scaleC = TRUE, k = min(dim(data)) - 
    1) 
{
    if (k > min(dim(data) - 1)) 
        warning("The number of components was too large compared to the data and was adjusted accordingly")
    k <- min(k, min(dim(data)) - 1)
    if (scaleR) {
        row.nrm <- apply(data, 1, sd)
        row.nrm <- pmax(row.nrm, 1e-04)
        data <- sweep(data, 1, row.nrm, FUN = "/")
    }
    result <- try(prcomp(data, retx = TRUE, scale = scaleC), 
        silent = TRUE)
    if (inherits(result, "try-error")) 
        stop("Failed to Calculate Principal Components")
    componentVariances <- result$sdev^2
    componentLoadings <- result$rotation[, 1:k]
    componentScores <- result$x[, 1:k]
    totalVariance <- sum(componentVariances)
    componentVariances <- componentVariances[1:k]
    z <- componentScores
    plot(cloud(z[, 1] ~ z[, 3] + z[, 2], groups = as.factor(labelValue), 
        auto.key = list(points = TRUE, pch = 19, space = "right"), 
        xlab = "PC 3", ylab = "PC 2", zlab = "PC 1", distance = 0.1, 
        main = "Projection in the space of the first 3 princial components"))
    value <- list(componentVariances = componentVariances, componentScores = componentScores, 
        componentLoadings = componentLoadings, summary = summary(result))
    value
}


cordist <- function(x) { as.dist((1-cor(t(x)))/2) }




