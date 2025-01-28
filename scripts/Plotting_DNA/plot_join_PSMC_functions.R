

modif_alpha <- function(col, alpha=.5){
  if(missing(col))
  stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))  
}


plot_PSMC_samples_wo_mu <- function(df, samples, murange, g, colors, vabline, sldf){
	Xcolumn1 <- "T_k_lowmu"
	Xcolumn2 <- "T_k_highmu"
	Ycolumn1 <- "N_k_lowmu"
	Ycolumn2 <- "N_k_highmu"
	yscaling <- 100000
	ylim <- c(0, max(df[,Ycolumn1])/yscaling)
	xlim <- log(c(round(min(df[which(df[,Xcolumn1]!=0),Xcolumn1])/1000, digits=0)*1000, round(max(df[,Xcolumn1])/100000, digits=0)*100000))
	ylimright <- c(-100, 100)
	sldf$OnAxisValues <- (sldf$SeaLevel_m-ylimright[1])/(ylimright[2]-ylimright[1])*(ylim[2]-ylim[1])+ylim[1]
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=xlim, col=NA)		
	mtext(paste0("Years ago (g=", g, ", ", expression(mu), "=", murange[1], "-", murange[2],")"), side = 1, line = 3, cex=1)
	mtext(paste0("Effective population size (x",yscaling,")"), side = 2, line = 3, cex=1)
	#arrows(log(sldf$Age), sldf$OnAxisValues, x1 = log(sldf$Age*murange[1]/murange[2]), y1 = sldf$OnAxisValues, length = 0)
	#lines(log(sldf$Age), sldf$OnAxisValues, lwd=2, col=modif_alpha("black"))
	#lines(log(sldf$Age*murange[1]/murange[2]), sldf$OnAxisValues, lwd=2, col=modif_alpha("black"))
	polygon(log(c(vabline*murange[2]/murange[1],vabline,vabline,vabline*murange[2]/murange[1])), c(ylim[1],ylim[1],ylim[2],ylim[2]), col=modif_alpha("forestgreen", .1), border=NA)
	for(s in c(1:length(samples))){
		sdf <- df[which(df$sample==samples[s]),]
		x <- sdf[,Xcolumn1]
		y <- sdf[,Ycolumn1]
		stepsx <- c(rbind(x,x))[2:length(c(rbind(x,x)))]
		stepsy <- c(rbind(y,y))[1:length(c(rbind(x,x)))-1]
		lines(log(stepsx), stepsy/yscaling, lwd=2, col=colors[s])
	}
	xaxlab1 <- c(10^3,10^4,10^5,10^6,10^7,10^8,10^9)
	xaxlab2 <- xaxlab1*murange[1]/murange[2]
	yaxlab1 <- seq(0,12,2)
	yaxlab2 <- yaxlab1*murange[1]/murange[2]

	axis(1, at = log(xaxlab1), labels=paste0(xaxlab1,"-",xaxlab2), lwd.ticks=1, las=1, cex.axis=1.5)
	axis(2, at = yaxlab1, labels=paste0(yaxlab1,"-",yaxlab2), lwd.ticks=1, las=1, cex.axis=1.5)
	#axis(4, at = seq(ylimright[1], ylimright[2], 50), lwd.ticks=1, las=1, cex.axis=1.5)
	box()
}

plot_PSMC_samples_w_mu <- function(df, samples, type, mu, g, colors){
	print(type)
	ylim <- c(0, round(max(df[,paste0("N_k_", type)])/1000000, digits=0)*100)
	xlim <- log(c(round(min(df[which(df[,paste0("T_k_", type)]!=0),paste0("T_k_", type)])/1000, digits=0)*1000, round(max(df[,paste0("T_k_", type)])/100000, digits=0)*100000))

	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=xlim, col=NA)		
	mtext(paste0("Years ago (g=", g, ", mu=", mu, ")"), side = 1, line = 3, cex=1)
	mtext("Effective population size (x10^4)", side = 2, line = 3, cex=1)
	for(s in c(1:length(samples))){
		sdf <- df[which(df$sample==samples[s]),]
		x <- sdf[,paste0("T_k_", type)]
		y <- sdf[,paste0("N_k_", type)]
		stepsx <- c(rbind(x,x))[2:length(c(rbind(x,x)))]
		stepsy <- c(rbind(y,y))[1:length(c(rbind(x,x)))-1]
		lines(log(stepsx), stepsy/10000, lwd=2, col=colors[s])
	}
	axis(1, at = log(c(10^3,10^4,10^5,10^6,10^7,10^8,10^9)), labels=c(10^3,10^4,10^5,10^6,10^7,10^8,10^9), lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(0,ylim[2],round(ylim[2]/100,digits=0)*10), lwd.ticks=1, las=1, cex.axis=1)
	box()
}
