


plot_scatter <- function(main, vec1, lab1, vec2, lab2, colors, xlim, ylim, cexp=1, text=NULL){
	plot(c(1:10), c(1:10), xlab="", ylab="", ylim=ylim, xlim=xlim, col=NA)
	mtext(lab1, side = 1, line = 3, cex=1)
	mtext(lab2, side = 2, line = 3, cex=1)
	mtext(main, side = 3, line = 1, cex=1.5)
	points(vec1, vec2, pch=21, bg=modif_alpha(colors), col=colors, cex=cexp)
	text(vec1, vec2, labels=text, font=1, cex=.5)
	box()
}


modif_alpha <- function(col, alpha=.5){
  if(missing(col))
  stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))  
}


barplot_samples <- function(ss, df, column, ylab, ylim, col=rep("black",length(ss)), hline=NULL){
	w <- 0.8
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(1,length(ss)), col=NA)
	mtext(ylab, side = 2, line = 4, cex=1.2)
	for(s in c(1:length(ss))){
		polygon(c(s-w/2,s+w/2,s+w/2,s-w/2), c(0,0,df[s,column],df[s,column]), col=col[s], border=NA)
	}
	abline(h=hline, lty=2, lwd=2)
	axis(1, at = c(1:length(ss)), labels=ss, lwd.ticks=1, las=2, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	box()
}

barplot_samples_chrs <- function(ss, cs, df, column, ylab, ylim, col=rep("black",length(ss)), hline=NULL){
	ns <- length(ss)
	twc <- 0.9
	w <- twc/ns
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(1,length(cs)), col=NA)
	mtext(ylab, side = 2, line = 4, cex=1.2)
	for(c in c(1:length(cs))){
		cdf <- df[which(df$chr==cs[c]),]
		stp <- c-twc/2
		for(s in c(1:length(ss))){
			polygon(c(stp+s*w-w,stp+s*w,stp+s*w,stp+s*w-w), c(0,0,cdf[s,column],cdf[s,column]), col=col[s], border=NA)
		}
	}
	abline(h=hline, lty=2, lwd=2)
	axis(1, at = c(1:length(cs)), labels=cs, lwd.ticks=1, las=2, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	box()
}

points_samples_chrs <- function(ss, cs, df, column, ylab, ylim, cols=c("gold", "forestgreen")){
	colfunc <- colorRampPalette(cols)
	ccols <- colfunc(length(cs))
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(1,length(ss)), col=NA)
	mtext(ylab, side = 2, line = 4, cex=1.2)
	for(s in c(1:length(ss))){
		for(c in c(1:length(cs))){
			val <- df[which(df$chr==cs[c] & df$Sample==ss[s]),column]
			points(s, val, col=ccols[c], pch=19)
		}
	}
	axis(1, at = c(1:length(ss)), labels=ss, lwd.ticks=1, las=2, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	#legend("topright", cs, pch=19, text.col="black", col=ccols, bty = "n", cex=1, xjust = 0, yjust = 0)
	box()
}

points_samples_chrs_legend <- function(cs, cols=c("gold", "forestgreen")){
	colfunc <- colorRampPalette(cols)
	ccols <- colfunc(length(cs))
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", col=NA)
	legend("topright", cs, pch=19, text.col="black", col=ccols, bty = "n", cex=1, xjust = 0, yjust = 0)
}


pi_perChr <- function(cs, cpi, pi, ylab, ylim){
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(1,length(cs)), col=NA)
	mtext(ylab, side = 2, line = 4, cex=1.2)
	points(c(1:length(cs)), cpi, col="black", pch=19)
	abline(h=pi, lty=2, lwd=2)
	axis(1, at = c(1:length(cs)), labels=cs, lwd.ticks=1, las=2, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	box()
}
