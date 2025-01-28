

PCA_calc <- function(Data){
	pca <- dudi.pca(Data, center=T, scale=T, scannf=F, nf=5)
	propvar <- 100 * pca$eig/sum(pca$eig)
	co.df <- as.data.frame(pca$co)
	li.df <- as.data.frame(pca$li)
	return(list("pca"=pca, "propvar"=propvar, "co.df"=co.df, "li.df"=li.df))
}

plot_PCA <- function(main, vec1, lab1, vec2, lab2, colors, samples, legcolors, leglabels, xlim, ylim){
	vec1 <- as.numeric(vec1)
	vec2 <- as.numeric(vec2)
	yflank <- (ylim[2]-ylim[1])*0.1
	xflank <- (xlim[2]-xlim[1])*0.1
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(ylim[1]-yflank,ylim[2]+yflank), xlim=c(xlim[1]-xflank,xlim[2]+xflank), col=NA)
	mtext(lab1, side = 1, line = 4, cex=2)
	mtext(lab2, side = 2, line = 4, cex=2)
	mtext(main, side = 3, line = 1, cex=2.5)
	points(vec1, vec2, pch=21, bg=modif_alpha(colors,0.3), col=colors, cex=4)
	text(vec1, vec2, labels=samples, pos=3, font=1, cex=1.5)
	if(length(leglabels)>1 & !is.na(leglabels[1])){
		legend("bottomri", as.character(leglabels), pch=21, text.col="black", col=legcolors, pt.bg=modif_alpha(legcolors,0.1), bty = "n", cex=2, xjust = 0, yjust = 0)
	}
	#axis(1, at = seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/5), lwd.ticks=1, las=2, cex.axis=1.5)
	#axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1.5)
	box()
}

modif_alpha <- function(col, alpha=.5){
  if(missing(col))
  stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))  
}

plot_x4_PCA <- function(pca, main, metad, pcol, Pop, scol, Sex, lim1, lim2, lim3, lim4){
	plot_PCA(main, 
		pca$co.df$Comp1, paste("PC1", round(pca$propvar[1], digits = 2), "%"), 
		pca$co.df$Comp2, paste("PC2", round(pca$propvar[2], digits = 2), "%"), 
		metad$ColorP, metad$Sample, pcol, Pop, lim1, lim2)
	plot_PCA(main, 
		pca$co.df$Comp1, paste("PC1", round(pca$propvar[1], digits = 2), "%"), 
		pca$co.df$Comp2, paste("PC2", round(pca$propvar[2], digits = 2), "%"), 
		metad$ColorS, metad$Sample, scol, Sex, lim1, lim2)
	plot_PCA(main, 
		pca$co.df$Comp3, paste("PC3", round(pca$propvar[3], digits = 2), "%"), 
		pca$co.df$Comp4, paste("PC4", round(pca$propvar[4], digits = 2), "%"), 
		metad$ColorP, metad$Sample, pcol, Pop, lim3, lim4)
	plot_PCA(main, 
		pca$co.df$Comp3, paste("PC3", round(pca$propvar[3], digits = 2), "%"), 
		pca$co.df$Comp4, paste("PC4", round(pca$propvar[4], digits = 2), "%"), 
		metad$ColorS, metad$Sample, scol, Sex, lim3, lim4)
}

plot_Het_PCA <- function(pca, main, metad, hcol){
	plot_PCA(main, 
		pca$co.df$Comp1, paste("PC1", round(pca$propvar[1], digits = 2), "%"), 
		pca$co.df$Comp2, paste("PC2", round(pca$propvar[2], digits = 2), "%"), 
		metad$ColorHet, metad$Sample, NA, NA, range(pca$co.df$Comp1), range(pca$co.df$Comp2))
	plot_PCA(main, 
		pca$co.df$Comp3, paste("PC3", round(pca$propvar[3], digits = 2), "%"), 
		pca$co.df$Comp4, paste("PC4", round(pca$propvar[4], digits = 2), "%"), 
		metad$ColorHet, metad$Sample, NA, NA, range(pca$co.df$Comp3), range(pca$co.df$Comp4))
}

calc_chr_means <- function(pca, DF, cs){
	DF <- cbind(DF, pca$li.df)
	chrData <- as.data.frame(matrix(ncol = length(colnames(pca$li.df))+1, nrow = 0))
	colnames(chrData) <- c("chr", colnames(pca$li.df))
	for(c in c(1:length(cs))){
		chr.mean.axis <- c()
		for(axis in colnames(pca$li.df)){
			chr.mean.axis <- c(chr.mean.axis,mean(DF[which(DF$chr==cs[c]),axis]))
		}
		chrData[c,] <- c(cs[c],chr.mean.axis)
	}
	return(chrData)
}

plot_Chrs_PCA <- function(chrData, main, cs, lim1, lim2, lim3, lim4){
	plot_PCA(main, 
		chrData$Axis1, paste("PC1", round(pca$propvar[1], digits = 2), "%"), 
		chrData$Axis2, paste("PC2", round(pca$propvar[2], digits = 2), "%"), 
		"black", cs, NA, NA, lim1, lim2)
	plot_PCA(main, 
		chrData$Axis3, paste("PC3", round(pca$propvar[3], digits = 2), "%"), 
		chrData$Axis4, paste("PC4", round(pca$propvar[4], digits = 2), "%"), 
		"black", cs, NA, NA, lim3, lim4)
}


