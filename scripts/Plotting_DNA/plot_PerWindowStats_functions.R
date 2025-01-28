


AlongSequence_AllChr <- function(DF, values, laby, ylim, cs, colors, niDF){
	DF$value <- values
	DF$midp <- DF$st + (DF$end-DF$st)/2
	xlim <- c(0, max(DF$end))
	layout(matrix(seq(1,length(cs),1),nrow=length(cs),ncol=1,byrow=T), widths=c(1), heights=rep(1/length(cs),length(cs)), TRUE)
	for(c in c(1:length(cs))){
		AlongSequence_OneChr(DF[which(DF$chr==cs[c]),], ylim, xlim, colors[c], niDF[which(niDF$chr==cs[c]),])
	}
	text(x=-xlim[2]/10, y=ylim[2]*length(cs)/2, laby, xpd=NA, srt = 90, cex=2)
}

AlongSequence_OneChr <- function(cDF, ylim, xlim, color, cniDF){
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=xlim, col=NA)
	lines(cDF$midp, cDF$value, lwd=2, col=color)
	for(i in c(1:length(cniDF[,1]))){
		polygon(c(cniDF$st[i], cniDF$end[i], cniDF$end[i], cniDF$st[i]), c(ylim[1],ylim[1],ylim[2],ylim[2]), col="white", border=NA, lwd=3)		
	}
	mtext(cDF$chr[1], side = 4, line = 1, cex=1)
	axis(1, at = seq(0,max(cDF$end), 10000000), labels=paste0(seq(0,max(cDF$end)/1000000,10),"Mb"), lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/2), lwd.ticks=1, las=1, cex.axis=1)
	box()
}

AlongSequence_AllChr_AllSamples <- function(DF, denominator, laby, ylim, cs, ss, scolbasename, scolors){
	DF$denominator <- denominator
	DF$midp <- DF$st + (DF$end-DF$st)/2
	for(s in ss){
		DF[,s] <- DF[,paste0(scolbasename,s)]/DF$denominator
	}
	xlim <- c(0, max(DF$end))
	layout(matrix(seq(1,length(cs),1),nrow=length(cs),ncol=1,byrow=T), widths=c(1), heights=rep(1/length(cs),length(cs)), TRUE)
	for(c in c(1:length(cs))){
		AlongSequence_OneChr_AllSamples(DF[which(DF$chr==cs[c]),], ylim, xlim, ss, scolors)
	}
	text(x=-xlim[2]/10, y=ylim[2]*length(cs)/2, laby, xpd=NA, srt = 90, cex=2)
}

AlongSequence_OneChr_AllSamples <- function(cDF, ylim, xlim, ss, colors){
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=xlim, col=NA)
	for(s in c(1:length(ss))){
		lines(cDF$midp, cDF[,ss[s]], lwd=2, col=colors[s])
	}
	mtext(cDF$chr[1], side = 4, line = 1, cex=1)
	axis(1, at = seq(0,max(cDF$end), 10000000), labels=paste0(seq(0,max(cDF$end)/1000000,10),"Mb"), lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/2), lwd.ticks=1, las=1, cex.axis=1)
	box()
}


plot_scatter <- function(main, vec1, lab1, vec2, lab2, colors, xlim, ylim, cexp=1){
	plot(c(1:10), c(1:10), xlab="", ylab="", ylim=ylim, xlim=xlim, col=NA)
	mtext(lab1, side = 1, line = 4, cex=1.2)
	mtext(lab2, side = 2, line = 4, cex=1.2)
	mtext(main, side = 3, line = 1, cex=1.2)
	points(vec1, vec2, pch=21, bg=modif_alpha(colors,0.01), col=NA, cex=cexp)
	box()
}


plot_Histogram <- function(values, xlab, ylab, ylim, xlim, step, vabline=NULL){
	h <- hist(values, plot=F, breaks=seq(0,100,step))
	w <- h$mids[2]-h$mids[1]
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=xlim, col=NA)
	mtext(xlab, side = 1, line = 3, cex=1.2)
	mtext(ylab, side = 2, line = 3, cex=1.2)
	h$counts <- h$counts/sum(h$counts)
	for(i in c(1:length(h$mids))){
		polygon(c(h$mids[i]-w/2,h$mids[i]+w/2,h$mids[i]+w/2,h$mids[i]-w/2), c(0,0,h$counts[i],h$counts[i]), col="black", border=NA)
	}
	abline(v=vabline)
	axis(1, at = seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	box()
}



