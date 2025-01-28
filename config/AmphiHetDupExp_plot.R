

ppal.color <- "slateblue3"
ppal.color2 <- "orangered"

types.of.features <- c("intergenic", "promoter", "exon", "intron")
types.of.CodingSites <- c("S", "N", "D")

types.of.features.color <- c("#B2533E", "#FCE09B", "#B5CB99", "#186F65")
types.of.CodingSites.color <- c("#B5CB99", "#B5CB99", "#B5CB99")

pvalthresh <- 0.001
population.colors <- c("orange", "purple") 
sex.colors <- c("black", "firebrick", "gold")
chr.colors <- c("gold","orange","darkred")



modif_color <- function(col, change){
	apply(sapply(col, col2rgb)/255, 2,
	function(x)
	return(rgb(max(0,min(x[1]+change,1)), max(0,min(x[2]+change,1)), max(0,min(x[3]+change,1)))))
}

modif_alpha <- function(col, alpha=.5){
  if(missing(col))
  stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))  
}
# Plots for Real data vs. Simulations comparison
plot_ObsSim <- function(main, ylab, mvec, avec, sdf, column, colNe, ylim){
	w <- 0.5
	adj <- 2
	uNemu <- unique(paste(sdf$Ne, sdf$mu, sep=" "))
	Nevalues <- unique(sdf$Ne)[order(unique(sdf$Ne))]
	muvalues <- unique(sdf$mu)[order(unique(sdf$mu))]
	print(c(0,length(Nevalues)+2))
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(0,length(Nevalues)+3), col=NA)
	mtext("", side = 1, line = 4, cex=1.2)
	mtext(ylab, side = 2, line = 4, cex=1.2)
	mtext(main, side = 3, line = 1, cex=1.2)
	#points(jitter(rep(2,length(avec)), amount=0.25), avec, pch=19, cex=1, col=population.colors[1])
	d <- density(avec, adjust = adj)
	ynorm <- d$y/max(d$y)*w/2
	polygon(c(1+ynorm, rev(1-ynorm)), c(d$x,rev(d$x)), col=modif_alpha(population.colors[1]), border=population.colors[1], lwd=1)
	points(1, median(avec), pch=19, cex=1, col=population.colors[1])
	#points(jitter(rep(1,length(mvec)), amount=0.25), mvec, pch=19, cex=1, col=population.colors[2])
	d <- density(mvec, adjust = adj)
	ynorm <- d$y/max(d$y)*w/2
	polygon(c(2+ynorm, rev(2-ynorm)), c(d$x,rev(d$x)), col=modif_alpha(population.colors[2]), border=population.colors[2], lwd=1)
	points(2, median(mvec), pch=19, cex=1, col=population.colors[2])
	for(nm in c(1:length(uNemu))){
		vec <- c()
		for(rep in unique(sdf$Rep)){
			points(jitter(2+nm, amount=0.25), sdf[which(paste(sdf$Ne, sdf$mu, sep=" ")==uNemu[nm] & sdf$Rep==rep), column], pch=19, cex=1, col=colNe[nm])
			vec <- c(vec, sdf[which(paste(sdf$Ne, sdf$mu, sep=" ")==uNemu[nm] & sdf$Rep==rep), column])
		}
		d <- density(vec, adjust = adj)
		ynorm <- d$y/max(d$y)*w/2
		polygon(c(2+nm+ynorm, rev(2+nm-ynorm)), c(d$x,rev(d$x)), col=modif_alpha(colNe[nm]), border=colNe[nm], lwd=1)
		#points(2+nm, median(vec), pch=19, cex=1, col=colNe[nm])
	}
	#axis(1, at = seq(1,length(Nevalues)+1), labels = NA, lwd.ticks=1, las=1, cex.axis=1)
	#axis(1, at = seq(1,length(Nevalues)+1), labels = c("Real", gsub(" ", "\n", uNemu)), lwd.ticks=NA, lwd=NA, line=1.5, las=1, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	box()
}

plot_Ne_mu <- function(sdf, colNe, leglabs){
	uNemu <- unique(paste(sdf$Ne, sdf$mu, sep=" "))
	Nevalues <- unique(sdf$Ne)[order(unique(sdf$Ne))]
	muvalues <- unique(sdf$mu)[order(unique(sdf$mu))]
	rNe <- c(min(log2(sdf$Ne)),max(log2(sdf$Ne)))
	print(rNe)
	rmu <- c(min(log2(sdf$mu)),max(log2(sdf$mu)))
	print(rmu)
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(-5,105), xlim=c(0,length(Nevalues)+3), col=NA)
	mtext("Ne", side = 2, line = 4, cex=1.2)
	mtext("mu", side = 4, line = 4, cex=1.2, col="darkred")
	Neline <- c()
	muline <- c()
	for(nm in c(1:length(uNemu))){
		vecNemu <- as.numeric(as.character(unlist(strsplit(uNemu[nm], " "))))
		points(2+nm, (log2(vecNemu[1])-rNe[1])/(rNe[2]-rNe[1])*100, pch=19, cex=1, col="black")
		points(2+nm, (log2(vecNemu[2])-rmu[1])/(rmu[2]-rmu[1])*100, pch=19, cex=1, col="darkred")
		Neline <- c(Neline, vecNemu[1])
		muline <- c(muline, vecNemu[2])
	}
	lines(c(1:length(uNemu))+2, (log2(Neline)-rNe[1])/(rNe[2]-rNe[1])*100, col="black")
	lines(c(1:length(uNemu))+2, (log2(muline)-rmu[1])/(rmu[2]-rmu[1])*100, col="darkred")
	axis(1, at = seq(1,length(Nevalues)+2), labels = NA, lwd.ticks=1, las=1, cex.axis=.5)
	axis(1, at = seq(1,length(Nevalues)+2), labels = c("Atlantic", "Mediterranean", gsub(" ", "\n", uNemu)), lwd.ticks=NA, lwd=NA, line=1, las=1, cex.axis=.5)
	axis(2, at = (log2(Nevalues)-rNe[1])/(rNe[2]-rNe[1])*100, labels=Nevalues, lwd.ticks=1, las=1, cex.axis=.5)
	axis(4, at = (log2(muvalues)-rmu[1])/(rmu[2]-rmu[1])*100, labels=muvalues, lwd.ticks=1, las=1, cex.axis=.5, col="darkred", col.axis="darkred")
	box()
}
