
Extract_AbsMajorFreq <- function(x){
	vfreq <- unlist(strsplit(x, ";|:"))
	return(as.numeric(vfreq[length(vfreq)]))
}

Extract_MajorAllele <- function(x){
	vfreq <- unlist(strsplit(x, ";|:"))
	return(vfreq[length(vfreq)-1])
}

Extract_SingletonSample <- function(x){
	if(x[grep("AbsMajorFreq", names(x))] == 71){
		genovalues <- x[grep("geno", names(x))]
		genosamples <- names(genovalues[grep(genovalues[1], genovalues, invert=TRUE)])
		if(length(genosamples)==1){
			return(substring(genosamples,5))
		}else{
			return(substring(names(genovalues)[1],5))
		}
	}else{return(NA)}
}

plot_scatter <- function(main, vec1, lab1, vec2, lab2, colors, xlim, ylim, cexp=1){
	plot(c(1:10), c(1:10), xlab="", ylab="", ylim=ylim, xlim=xlim, col=NA)
	mtext(lab1, side = 1, line = 4, cex=1.2)
	mtext(lab2, side = 2, line = 4, cex=1.2)
	mtext(main, side = 3, line = 1, cex=1.2)
	points(vec1, vec2, pch=21, bg=modif_alpha(colors,0.01), col=NA, cex=cexp)
	box()
}

plot_scatter_heatmap <- function(main, vec1, lab1, vec2, lab2, colors, xlim, ylim, by1, by2, cexp=1){
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=xlim, col=NA)
	mtext(lab1, side = 1, line = 4, cex=1.2)
	mtext(lab2, side = 2, line = 4, cex=1.2)
	mtext(main, side = 3, line = 1, cex=1.2)
	heatmapFromVectors(vec1, vec2, seq(xlim[1],xlim[2],by1), seq(ylim[1],ylim[2],by2), colors)
	axis(1, at = seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	box()
}

plot_RelativeFreqDiscreteValues_histogram <- function(values, total, mids, xlab, ylab, main, ylim, col){
	w <- mids[2]-mids[1]
	h <- hist(values, plot=F, breaks=seq(mids[1]-w/2, mids[length(mids)]+w/2, w))
	h$counts <- h$counts/total
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(mids[1]-w/2, mids[length(mids)]+w/2), col=NA)
	mtext(xlab, side = 1, line = 3, cex=1.2)
	mtext(ylab, side = 2, line = 3, cex=1.2)
	mtext(main, side = 3, line = 1, cex=1.2)
	for(i in c(1:length(h$mids))){
		polygon(c(h$mids[i]-w/2,h$mids[i]+w/2,h$mids[i]+w/2,h$mids[i]-w/2), c(0,0,h$counts[i],h$counts[i]), col=col, border=NA)
	}
	axis(1, at = h$mids, labels=, lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	box()
}

plot_histogram_absoluteSFS <- function(values, mids, xlab, main, ylim, col){
	w <- mids[2]-mids[1]
	h <- hist(values, plot=F, breaks=seq(mids[1]-w/2, mids[length(mids)]+w/2, w))
	print(main)
	print(h$counts)
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(mids[1]-w/2, mids[length(mids)]+w/2), col=NA)
	mtext(xlab, side = 1, line = 3, cex=1.2)
	mtext("Proportion of sites", side = 2, line = 3, cex=1.2)
	mtext(main, side = 3, line = 1, cex=1.2)
	h$counts <- h$counts/sum(h$counts)
	# Draw expectation of SFS
	eh <- h$counts[1]
	for(i in c(2:length(h$counts))){
		eh <- c(eh, h$counts[1]/i)
	}
	for(i in c(1:length(h$mids))){
		polygon(c(h$mids[i]-w/2,h$mids[i]+w/2,h$mids[i]+w/2,h$mids[i]-w/2), c(0,0,eh[i],eh[i]), col=modif_color(col, .3), border=NA)
	}
	# Draw real values
	for(i in c(1:length(h$mids))){
		polygon(c(h$mids[i]-w/2,h$mids[i]+w/2,h$mids[i]+w/2,h$mids[i]-w/2), c(0,0,h$counts[i],h$counts[i]), col=col, border=NA)
	}
	abline(v=32)
	axis(1, at = seq(h$mids[1],h$mids[length(h$mids)],5), lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	box()
	#return(h)
}



plot_joined_histogram_SFS <- function(values, feat, types, breaks, xlab, main, ylim, colors){
	w <- breaks[2]-breaks[1]
	df <- as.data.frame(cbind(values, feat))
	print(head(df))
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(min(breaks),max(breaks)), col=NA)
	mtext(xlab, side = 1, line = 3, cex=1.2)
	mtext("Proportion of sites", side = 2, line = 3, cex=1.2)
	mtext(main, side = 3, line = 1, cex=1.2)
	## Draw expectation of SFS
	#eh <- h$counts[1]
	#for(i in c(2:length(h$counts))){
	#	eh <- c(eh, h$counts[1]/i)
	#}
	#for(i in c(1:length(h$mids))){
	#	polygon(c(h$mids[i]-w/2,h$mids[i]+w/2,h$mids[i]+w/2,h$mids[i]-w/2), c(0,0,eh[i],eh[i]), col=modif_color(col, .3), border=NA)
	#}
	# Draw real values
	for(i in c(1:length(types))){
		print(i)
		h <- hist(as.numeric(na.omit(df$values[which(df$feat == types[i])])), plot=F, breaks=breaks)
		h$counts <- h$counts/sum(h$counts)
		print(breaks)
		print(h$counts)
		stepsx <- c(rbind(breaks,breaks))
		stepsy <- c(0,rbind(h$counts,h$counts)[1:length(c(rbind(h$mids,h$mids)))],0)
		print(stepsx)
		print(stepsy)
		lines(stepsx, stepsy, lwd=2, col=colors[i])
	}
	axis(1, at = seq(min(breaks),max(breaks),(max(breaks)-min(breaks))/5), lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	legend("topright", types, pch=19, text.col="black", col=colors, bty = "n", cex=2, xjust = 0, yjust = 0)
	box()
	#return(h)
}

plot_propVar_propLength_PerFeat <- function(DF, types, vcolumn, lcolumn, main, ylim, col=c("forestgreen","gold")){
	w <- .3
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(.5,length(types)+.5), col=NA)
	mtext(main, side = 3, line = 1, cex=1.5)
	mtext("%", side = 2, line = 3, cex=1.2)
	for(i in c(1:length(types))){
		propLen <- DF[which(DF$feature==types[i]), lcolumn]
		polygon(c(i-w,i,i,i-w), c(0,0,propLen,propLen), col=col[1], lwd=3)
		text(i-w/2, propLen+5,  labels =paste0(format(round(propLen, 1), nsmall = 1),"%"), pos=3, srt = 90, cex=1.5)
		propVar <- DF[which(DF$feature==types[i]), vcolumn]
		polygon(c(i,i+w,i+w,i), c(0,0,propVar,propVar), col=col[2], lwd=3)	
		text(i+w/2, propVar+5,  labels =paste0(format(round(propVar, 1), nsmall = 1),"%"), pos=3, srt = 90, cex=1.5)
	}
	axis(1, at = seq(1,length(types),1), labels=types, lwd.ticks=1, las=2, cex.axis=1.5)
	axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=1)
	legend("topleft", c("% of callable length", "% of variants"), pch=19, text.col="black", col=col, bty = "n", cex=2, xjust = 0, yjust = 0)
	box()
}

plot_value_PerFeat <- function(DF, types, vcolumn, ylab, main, ylim, col){
	w <- .3
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(.5,length(types)+.5), col=NA)
	mtext(main, side = 3, line = 1, cex=1.5)
	mtext(ylab, side = 2, line = 3, cex=1.2)
	for(i in c(1:length(types))){
		value <- DF[which(DF$feature==types[i]), vcolumn]
		polygon(c(i-w/2,i+w/2,i+w/2,i-w/2), c(0,0,value,value), col=col[i], lwd=3)
	}
	axis(1, at = seq(1,length(types),1), labels=types, lwd.ticks=1, las=2, cex.axis=1.5)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	box()	
}

plot_distribution_PerFeat <- function(values, features, types, xlab, main, col, xlim, ptresh, ofolder){
	w <- .6
	df <- as.data.frame(cbind(values, features))
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,10), xlim=xlim, col=NA)
	mtext("Density", side = 2, line = 3, cex=1.2)
	mtext(xlab, side = 1, line = 3, cex=1.2)
	mtext(main, side = 3, line = 1, cex=1.5)
	intergenic <- as.numeric(df$values[which(df$features=="intergenic")])
	write("Comparison\tp.value", file = paste0(ofolder, "/T.test_", gsub(" ", "_", xlab),".txt"))
	for(i in c(1:length(types))){
		print(types[i])
		v <- as.numeric(df$values[grep(types[i], df$features)])
		t<-t.test(v, intergenic)
		if(t$p.value<=ptresh){
			print(t)
			#text(i, xlim[2],  labels =t$p.value, pos=1, cex=1)
			write(paste0("intergenic-", types[i], "\t" ,t$p.value), file = paste0(ofolder, "/T.test_", gsub(" ", "_", xlab),".txt"), append=TRUE)
		}
		d <- density(v, adjust = 7)
		lines(d$x, d$y, col=col[i], lwd=3)
	}
	axis(1, at = seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/5), lwd.ticks=1, las=1, cex.axis=1.5)
	axis(2, at = seq(0,10,2), lwd.ticks=1, las=1, cex.axis=1.5)
	legend("topright", types, pch=19, text.col="black", col=col, bty = "n", cex=2, xjust = 0, yjust = 0)
	box()
}

plot_boxplot_PerFeat <- function(values, features, types, ylab, main, col, ylim, ptresh, ofolder){
	w <- .6
	df <- as.data.frame(cbind(values, features))
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(.5,length(types)+.5), col=NA)
	mtext(ylab, side = 2, line = 3, cex=1.2)
	mtext(main, side = 3, line = 1, cex=1.5)
	intergenic <- as.numeric(df$values[which(df$features=="intergenic")])
	write("Comparison\tp.value", file = paste0(ofolder, "/T.test_", gsub(" ", "_", ylab),".txt"))
	for(i in c(1:length(types))){
		print(types[i])
		v <- as.numeric(df$values[grep(types[i], df$features)])
		t<-t.test(v, intergenic)
		if(t$p.value<=ptresh){
			print(t)
			text(i, ylim[2],  labels =t$p.value, pos=1, cex=1)
			write(paste0("intergenic-", types[i], "\t" ,t$p.value), file = paste0(ofolder, "/T.test_", gsub(" ", "_", ylab),".txt"), append=TRUE)
		}
		DrawBox(v, i, col[i], .6)
		#d <- density(v, adjust = 2)
		#polygon(c(i+d$y/max(d$y)*w/2, rev(i-d$y/max(d$y)*w/2)), c(d$x,rev(d$x)), col=NA, border=col, lwd=2)
	}
	abline(h=median(intergenic), lty=2, lwd=2)
	axis(1, at = seq(1,length(types),1), labels=types, lwd.ticks=1, las=2, cex.axis=1.5)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	box()
}

plot_boxplot_PerFeat_SND <- function(values, features, SND, types, typescod, ylab, main, col, ylim, ptresh, ofolder){
	w <- .6
	df <- as.data.frame(cbind(values, features, SND))
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(.5,(length(types)-1)*1+length(typescod)*.5+1), col=NA)
	mtext(ylab, side = 2, line = 3, cex=1.2)
	mtext(main, side = 3, line = 1, cex=1.5)
	intergenic <- as.numeric(df$values[which(df$features=="intergenic")])
	write("Comparison\tp.value", file = paste0(ofolder, "/T.test_", gsub(" ", "_", ylab),".txt"))
	pos <- 0
	axislabels <- c()
	axispos <- c()
	for(i in c(1:length(types))){
		print(types[i])
		if(types[i] != "exon"){
			pos <- pos+1
			v <- as.numeric(df$values[grep(types[i], df$features)])
			n <- as.numeric(length(grep(types[i], df$features)))
			t<-t.test(v, intergenic)
			if(t$p.value<=ptresh){
				print(t)
				text(pos, ylim[2],  labels =t$p.value, pos=1, cex=1)
				write(paste0("intergenic-", types[i], "\t" ,t$p.value), file = paste0(ofolder, "/T.test_", gsub(" ", "_", ylab),".txt"), append=TRUE)
			}
			DrawBox(v, pos, col[i], w)
			text(pos, ylim[1],  labels =n, pos=1, cex=1)
			axislabels <- c(axislabels, types[i])
			axispos <- c(axispos, pos)
		}else{
			pos <- pos+.5
			for(j in c(1:length(typescod))){
				pos <- pos+.5
				v <- as.numeric(df$values[intersect(grep(types[i], df$features), grep(typescod[j], df$SND))])
				n <- as.numeric(length(intersect(grep(types[i], df$features), grep(typescod[j], df$SND))))
				t<-t.test(v, intergenic)
				if(t$p.value<=ptresh){
					print(t)
					text(pos, ylim[2],  labels =t$p.value, pos=1, cex=1)
					write(paste0("intergenic-", typescod[j], "\t" ,t$p.value), file = paste0(ofolder, "/T.test_", gsub(" ", "_", ylab),".txt"), append=TRUE)
				}
				DrawBox(v, pos, col[i], w/2)
				text(pos, ylim[1],  labels =n, pos=1, cex=1)
				axislabels <- c(axislabels, typescod[j])
				axispos <- c(axispos, pos)
			}
		}
	}
	abline(h=median(intergenic), lty=2, lwd=2)
	axis(1, at = axispos, labels=axislabels, lwd.ticks=1, las=2, cex.axis=1.5)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	box()
}

DrawBox <- function(values, pos, col, w=.8, den=NULL, text=FALSE, cextext=1){
	s <- boxplot(values, plot=FALSE)
	#points(rep(pos,length(s$out)), s$out, col=modif_alpha("black", 0.3), pch=16, cex=1.5, xpd = NA)
	arrows(x0=pos, y0=s$stats[1], x1=pos, y1=s$stats[5], angle=90, code=3, length=w/10, lwd=2, xpd = NA)
	#lines(c(pos, pos),c(s$stats[1], s$stats[5]), lwd=2)
	if(!is.null(den)){
		polygon(c(pos-w/2, pos+w/2, pos+w/2, pos-w/2), c(s$stats[2],s$stats[2],s$stats[4],s$stats[4]), col=col, border=col, lwd=3)		
		polygon(c(pos-w/2, pos+w/2, pos+w/2, pos-w/2), c(s$stats[2],s$stats[2],s$stats[4],s$stats[4]), col="black", border="black", lwd=3, density=den)
	}else{
		polygon(c(pos-w/2, pos+w/2, pos+w/2, pos-w/2), c(s$stats[2],s$stats[2],s$stats[4],s$stats[4]), col=col, border="black", lwd=3, density=den)		
	}
	lines(c(pos-w/2, pos+w/2),c(s$stats[3], s$stats[3]), lwd=2)
	polygon(c(pos-w/2, pos+w/2, pos+w/2, pos-w/2), c(s$stats[2],s$stats[2],s$stats[4],s$stats[4]), col=NA, border="black", lwd=3)		
	if(text){
		par(xpd=TRUE) 
		text(pos, 0,  labels =length(values), pos=1, cex=cextext)
		par(xpd=FALSE) 		
	}
}

heatmapFromVectors <- function(vec1, vec2, breaks1, breaks2, cols){
	ncol <- 200
	colfunc <- colorRampPalette(cols)
	hcol <- colfunc(ncol)
	mat <- cbind(vec1[!is.na(vec1) & !is.na(vec2)], vec2[!is.na(vec1) & !is.na(vec2)])
	bined <- t(bin2(mat, ab=rbind(c(min(breaks1),max(breaks1)),c(min(breaks2),max(breaks2))), nbin=c(length(breaks1)-1,length(breaks2)-1))$nc)
	for(c in c(1:(length(breaks1)-1))){
		for(r in c(1:(length(breaks2)-1))){
			if(bined[r,c]>0){
				cellcol <- hcol[unlist(cut(bined[r,c], breaks=seq(min(bined[which(bined !=0)]),max(bined),(max(bined)-min(bined[which(bined !=0)]))/ncol), labels=c(1:ncol)))]
				polygon(c(breaks1[c], breaks1[c+1], breaks1[c+1], breaks1[c]), c(breaks2[r], breaks2[r], breaks2[r+1], breaks2[r+1]), col=cellcol, border=NA)
			}
		}
	}
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
	#polygon(c(min(d$x)-10,d$x,max(d$x)+10), c(-10,d$y,10), col=color, border=NA)
	lines(d$x, d$y, lwd=2)
	abline(v=thresh, lty=2, lwd=2)
	abline(v=thresh2, lty=2, lwd=2)
	axis(1, at = seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	box()
}

plot_histogram <- function(vec, ylab, xlab, xlim, ylim, breaks, color, log=FALSE){
	w <- breaks[2]-breaks[1]
	h <- hist(na.omit(vec), breaks=breaks, plot=FALSE)
	if(log){
		counts = log(h$counts)
	}else{
		counts = h$counts
	}
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=xlim, col=NA, xaxs="i", yaxs="i")	
	mtext(xlab, side = 1, line = 3, cex=1.2)
	mtext(ylab, side = 2, line = 3, cex=1.2)
	for(i in c(1:length(h$mids))){
		if(counts[i]>0){
			polygon(c(h$mids[i]-w/2, h$mids[i]+w/2, h$mids[i]+w/2, h$mids[i]-w/2), c(0,0,counts[i],counts[i]), col=color, border=NA, lwd=3)		
		}
	}
	axis(1, at = seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/4), lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/4), lwd.ticks=1, las=1, cex.axis=1)
	box()
}


Check.EqGenotypes <- function(x){
	len <- length(unlist(unique(as.vector(x))))
	if(len==1){
		return(1)
	}else{
		return(0)
	}
}


AlongSequence_AllChr <- function(DF, values, laby, ylim, cs, colors){
	DF$value <- values
	xlim <- c(0, max(DF$end))
	layout(matrix(seq(1,length(cs),1),nrow=length(cs),ncol=1,byrow=T), widths=c(1), heights=rep(1/length(cs),length(cs)), TRUE)
	for(c in c(1:length(cs))){
		AlongSequence_OneChr(DF[which(DF$chr==cs[c]),], ylim, xlim, colors[c])
	}
	text(x=-xlim[2]/10, y=ylim[2]*length(cs)/2, laby, xpd=NA, srt = 90, cex=2)
}

AlongSequence_OneChr <- function(cDF, ylim, xlim, color){
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=xlim, col=NA)
	lines(cDF$st, cDF$value, lwd=2, col=color)
	mtext(cDF$chr[1], side = 4, line = 1, cex=1)
	axis(1, at = seq(0,max(cDF$end), 10000000), labels=paste0(seq(0,max(cDF$end)/1000000,10),"Mb"), lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/2), lwd.ticks=1, las=1, cex.axis=1)
	box()
}


