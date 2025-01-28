

draw_chr_row <- function(h, p, starts, ends, color){
	for(r in c(1:length(starts))){
		polygon(c(starts[r],ends[r],ends[r],starts[r]), c(p-h/2,p-h/2,p+h/2,p+h/2), col=color, border=NA)
	}
}

plot_all_chr_rows_len_bins <- function(lenDF, regDF, ss, es, cs){
	height <- 0.6
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,length(lenDF$Chr)), xlim=c(0,10), col=NA)
	for(c in c(1:length(lenDF$Chr))){
		print(lenDF$Chr[c])
		pos <- length(lenDF$Chr)-c
		for(b in c(1:length(ss))){
			sts <- regDF$st[which(regDF$chr==lenDF$Chr[c] & regDF$len>=ss[b] & regDF$len<es[b])]/max(lenDF$Length)*10
			ends <- regDF$end[which(regDF$chr==lenDF$Chr[c] & regDF$len>=ss[b] & regDF$len<es[b])]/max(lenDF$Length)*10
			draw_chr_row(height, pos, sts, ends, cs[b])
		}
		polygon(c(0,lenDF$Length[c]/max(lenDF$Length)*10,lenDF$Length[c]/max(lenDF$Length)*10,0), c(pos-height/2,pos-height/2,pos+height/2,pos+height/2), col=NA, border="black")
	}
	axis(2, at = c(1:length(lenDF$Chr))-1, labels=rev(lenDF$Chr), lwd=NA, lwd.ticks=NA, las=1, cex.axis=1)
	legend("bottomright", paste(ss, "<= x <", es), pch=19, text.col="black", col=cs, bty = "n", cex=1, xjust = 0, yjust = 0)
}


barplot_chromosomes <- function(chrs, df, type, ss, es, cs, ylab, ylim){
	print(paste("barplot_chromosomes",type))
	w <- 0.8
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(1-w/2,length(chrs)+w/2), col=NA)
	mtext(type, side = 3, line = 3, cex=1.5)
	mtext(ylab, side = 2, line = 3, cex=1.5)
	mtext("Chromosomes", side = 1, line = 3, cex=1.5)
	for(c in c(1:length(chrs))){
		print(chrs[c])
		cumst <- 0
		for(b in c(1:length(ss))){
			print(paste(b, df[c,paste0("Perc", type, b)], cs[b]))
			polygon(c(c-w/2,c+w/2,c+w/2,c-w/2), c(cumst,cumst,cumst+df[c,paste0("Perc", type, b)],cumst+df[c,paste0("Perc", type, b)]), col=cs[b], border=NA)
			cumst <- cumst+df[c,paste0("Perc", type, b)]
		}
	}
	axis(1, at = c(1:length(chrs)), labels=chrs, lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	legend("topright", paste(ss, "<= x <", es), pch=19, text.col="black", col=cs, bty = "n", cex=1, xjust = 0, yjust = 0)
	box()
}

plot_density <- function(vec, lab, color, xlim, ylim, thresh=NULL, thresh2=NULL, printlab=TRUE){
	print(lab)
	d <- density(na.omit(vec), adjust = 3)
	print(c(min(d$x), max(d$x)))
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=xlim, col=NA, xaxs="i", yaxs="i")	
	if(printlab==TRUE){
		mtext(lab, side = 1, line = 3, cex=1.5)
		mtext("Density", side = 2, line = 3, cex=1.5)
	}	
	polygon(c(min(d$x)-10,d$x,max(d$x)+10), c(-10,d$y,10), col=color, border=NA)
	lines(d$x, d$y, lwd=2)
	abline(v=thresh, lty=2, lwd=2)
	abline(v=thresh2, lty=2, lwd=2)
	axis(1, at = seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	box()
}












