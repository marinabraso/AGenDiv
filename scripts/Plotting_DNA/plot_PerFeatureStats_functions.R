
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



plot_perGenomicFeature_boxplot <- function(values, features, types, ylab, col, ylim, ptresh, ofolder){
	w <- .6
	df <- as.data.frame(cbind(values, features))
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(.5,length(types)+.5), col=NA)
	mtext(ylab, side = 2, line = 4, cex=1.5)
	intergenic <- as.numeric(df$values[which(df$features=="intergenic")])
	write("Comparison\tp.value", file = paste0(ofolder, "/T.test_", gsub(" ", "_", gsub("/", "_", ylab)),".txt"))
	for(i in c(1:length(types))){
		print(types[i])
		v <- as.numeric(df$values[which(df$features==types[i])])
		DrawBox(v, i, col[i], w)
		t<-t.test(v, intergenic)
		stars <- 0
		if(t$p.value<=ptresh[1]){
			stars <- stars+1
			write(paste0("intergenic-", types[i], "\t" ,t$p.value), file = paste0(ofolder, "/T.test_", gsub(" ", "_", gsub("/", "_", ylab)),".txt"), append=TRUE)
		}
		if(t$p.value<=ptresh[2]){
			stars <- stars+1
		}
		if(t$p.value<=ptresh[3]){
			stars <- stars+1
		}
		print(t)
		if(stars > 0){
			text(i, ylim[2],  labels =paste0(rep("*", stars), collapse=""), pos=1, cex=2)
		}
		#d <- density(v, adjust = 2)
		#polygon(c(i+d$y/max(d$y)*w/2, rev(i-d$y/max(d$y)*w/2)), c(d$x,rev(d$x)), col=NA, border=col, lwd=2)
	}
	abline(v=3.5, lty=2, lwd=.5)
	axis(1, at = seq(1,length(types),1), labels=types, lwd.ticks=1, las=1, cex.axis=1.5)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	box()
}


plot_perGenomicFeature_violin <- function(values, features, types, ylab, col, ylim, ptresh, ofolder){
	w <- .6
	df <- as.data.frame(cbind(values, features))
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(.5,length(types)+.5), col=NA)
	mtext(ylab, side = 2, line = 4, cex=1.5)
	intergenic <- as.numeric(df$values[which(df$features=="intergenic")])
	write("Comparison\tp.value", file = paste0(ofolder, "/T.test_", gsub(" ", "_", gsub("/", "_", ylab)),".txt"))
	for(i in c(1:length(types))){
		print(types[i])
		v <- as.numeric(df$values[which(df$features==types[i])])
		d <- density(v, adjust = 1.5)
		ynorm <- d$y/max(d$y)*w/2
		polygon(c(i-ynorm, rev(i+ynorm)), c(d$x,rev(d$x)), col=modif_alpha(col[i],0.4), border=col[i], lwd=2)
		#points(jitter(rep(i, length(featvalA)), amount=w/2), featvalA, pch=21, bg=modif_alpha(col[i],0.2), col=modif_alpha(col[i]))
		t<-t.test(v, intergenic)
		stars <- 0
		if(t$p.value<=ptresh[1]){
			stars <- stars+1
			write(paste0("intergenic-", types[i], "\t" ,t$p.value), file = paste0(ofolder, "/T.test_", gsub(" ", "_", gsub("/", "_", ylab)),".txt"), append=TRUE)
		}
		if(t$p.value<=ptresh[2]){
			stars <- stars+1
		}
		if(t$p.value<=ptresh[3]){
			stars <- stars+1
		}
		print(t)
		if(stars > 0){
			text(i, ylim[2],  labels =paste0(rep("*", stars), collapse=""), pos=1, cex=2)
		}
		#d <- density(v, adjust = 2)
		#polygon(c(i+d$y/max(d$y)*w/2, rev(i-d$y/max(d$y)*w/2)), c(d$x,rev(d$x)), col=NA, border=col, lwd=2)
	}
	#abline(h=median(intergenic), lty=2, lwd=2)
	axis(1, at = seq(1,length(types),1), labels=types, lwd.ticks=1, las=2, cex.axis=1.5)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	box()
}

plot_propSND_boxplot <- function(values, SNDvalues, features, types, SNDtypes, ylab, col, ylim, ptresh, ofolder){
	w <- .6
	df <- as.data.frame(cbind(values, SNDvalues, features))
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(.5,(length(types)-1)*1+length(SNDtypes)*.5+1), col=NA)
	mtext(ylab, side = 2, line = 3, cex=1.5)
	intergenic <- as.numeric(df$values[which(df$features=="intergenic")])
	write("Comparison\tp.value", file = paste0(ofolder, "/T.test_", gsub(" ", "_", gsub("/", "_", ylab)),".txt"))
	pos <- 0
	axislabels <- c()
	axispos <- c()
	for(i in c(1:length(types))){
		print(types[i])
		if(types[i]!="exon"){
			pos <- pos+1
			v <- as.numeric(df$values[which(df$features==types[i])])
			t<-t.test(v, intergenic)
			stars <- 0
			if(t$p.value<=ptresh[1]){
				stars <- stars+1
				write(paste0("intergenic-", types[i], "\t" ,t$p.value), file = paste0(ofolder, "/T.test_", gsub(" ", "_", gsub("/", "_", ylab)),".txt"), append=TRUE)
			}
			if(t$p.value<=ptresh[2]){
				stars <- stars+1
			}
			if(t$p.value<=ptresh[3]){
				stars <- stars+1
			}
			print(t)
			if(stars > 0){
				text(pos, ylim[2],  labels =paste0(rep("*", stars)), pos=1, cex=1)
			}
			DrawBox(v, pos, col[i], w)
			#d <- density(v, adjust = 2)
			#polygon(c(i+d$y/max(d$y)*w/2, rev(i-d$y/max(d$y)*w/2)), c(d$x,rev(d$x)), col=NA, border=col, lwd=2)
			axislabels <- c(axislabels, types[i])
			axispos <- c(axispos, pos)
		}else{
			pos <- pos+.5
			for(j in c(1:length(SNDtypes))){
				print(SNDtypes[j])
				pos <- pos+.5
				v <- as.numeric(as.character(df[which(df$features==types[i]),paste0("prop", SNDtypes[j])]))
				print(head(v))
				stars <- 0
				if(t$p.value<=ptresh[1]){
					stars <- stars+1
					write(paste0("intergenic-", types[i], "\t" ,t$p.value), file = paste0(ofolder, "/T.test_", gsub(" ", "_", gsub("/", "_", ylab)),".txt"), append=TRUE)
				}
				if(t$p.value<=ptresh[2]){
					stars <- stars+1
				}
				if(t$p.value<=ptresh[3]){
					stars <- stars+1
				}
				print(t)
				if(stars > 0){
					text(pos, ylim[2],  labels =paste0(rep("*", stars)), pos=1, cex=1)
				}
				DrawBox(v, pos, col[i], w/2)
				axislabels <- c(axislabels, SNDtypes[j])
				axispos <- c(axispos, pos)
			}
		}
	}
	#abline(h=median(intergenic), lty=2, lwd=2)
	axis(1, at = axispos, labels=axislabels, lwd.ticks=1, las=2, cex.axis=1.5)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	box()
}

plot_perGenomicFeature_perSample <- function(DF, scolumns, fcolumn, totalcolumn, types, ylab, pop, pcol, ylim){
	w <- .3
	print(pop)
	DF[,scolumns] <- DF[,scolumns]/DF[,totalcolumn]
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(.5,length(types)+.5), col=NA)
	mtext(ylab, side = 2, line = 3, cex=1.5)
	for(i in c(1:length(types))){
		print(types[i])
		featval <- colMeans(DF[grep(types[i], DF[,fcolumn]),scolumns])
		featvalM <- featval[grep("Banyuls", pop)]
		d <- density(featvalM, adjust = 2)
		ynorm <- d$y/max(d$y)*w/2
		polygon(c(i-w*.6-ynorm, rev(i-w*.6+ynorm)), c(d$x,rev(d$x)), col=modif_alpha(pcol[grep("Banyuls", pop)],0.2), border=pcol[grep("Banyuls", pop)], lwd=1)
		points(jitter(rep(i-w*.6, length(featvalM)), amount=w/2), featvalM, pch=21, bg=modif_alpha(pcol[grep("Banyuls", pop)],0.2), col=modif_alpha(pcol[grep("Banyuls", pop)]))

		featvalA <- featval[grep("Roscoff", pop)]
		d <- density(featvalA, adjust = 2)
		ynorm <- d$y/max(d$y)*w/2
		polygon(c(i+w*.6-ynorm, rev(i+w*.6+ynorm)), c(d$x,rev(d$x)), col=modif_alpha(pcol[grep("Roscoff", pop)],0.2), border=pcol[grep("Roscoff", pop)], lwd=1)
		points(jitter(rep(i+w*.6, length(featvalA)), amount=w/2), featvalA, pch=21, bg=modif_alpha(pcol[grep("Roscoff", pop)],0.2), col=modif_alpha(pcol[grep("Roscoff", pop)]))
		
		#points(jitter(rep(i, length(scolumns)),amount=.3), featval, pch=19, col=pcol)
	}
	axis(1, at = seq(1,length(types),1), labels=types, lwd.ticks=1, las=2, cex.axis=1.5)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	box()
}

DrawBox <- function(values, pos, col, w=.8, den=NULL, text=FALSE, cextext=1){
	s <- boxplot(values, plot=FALSE)
	#points(rep(pos,length(s$out)), s$out, col=modif_alpha("black", 0.3), pch=16, cex=1.5, xpd = NA)
	arrows(x0=pos, y0=s$stats[1], x1=pos, y1=s$stats[5], angle=90, code=3, length=w/10, lwd=2, xpd = NA)
	#lines(c(pos, pos),c(s$stats[1], s$stats[5]), lwd=2)
	polygon(c(pos-w/2, pos+w/2, pos+w/2, pos-w/2), c(s$stats[2],s$stats[2],s$stats[4],s$stats[4]), col="white", border="white", lwd=2, density=den)		
	polygon(c(pos-w/2, pos+w/2, pos+w/2, pos-w/2), c(s$stats[2],s$stats[2],s$stats[4],s$stats[4]), col=modif_alpha(col), border=col, lwd=2, density=den)		
	lines(c(pos-w/2, pos+w/2),c(s$stats[3], s$stats[3]), lwd=2, col=col)
	polygon(c(pos-w/2, pos+w/2, pos+w/2, pos-w/2), c(s$stats[2],s$stats[2],s$stats[4],s$stats[4]), col=NA, border=col, lwd=2)		
	if(text){
		par(xpd=TRUE) 
		text(pos, 0,  labels =length(values), pos=1, cex=cextext)
		par(xpd=FALSE) 		
	}
}

