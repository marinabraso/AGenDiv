

ppal.color <- "slateblue3"
ppal.color2 <- "orangered"

types.of.features <- c("Exon", "Intron", "Promoter", "Intergenic")
types.of.CodingSites <- c("S", "N", "D")

types.of.features.color <- c("#B5CB99", "#186F65", "#FCE09B", "#B2533E")
types.of.CodingSites.color <- c("#B5CB99", "#B5CB99", "#B5CB99")

pvalthresh <- 0.001
populations <- c("Atlantic", "Mediterranean")
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

violin <- function(values, pos, width, col, border, ad, lwd=1){
	d <- density(as.numeric(values), adjust = ad)
	ynorm <- d$y/max(d$y)*width/2
	polygon(c(pos+ynorm, rev(pos-ynorm)), c(d$x,rev(d$x)), col=col, border=border, lwd=lwd)
}

writePlotLabel <- function(txt){
    old.par <- par(xpd=NA)
    di <- dev.size("in")
    x <- grconvertX(c(0, di[1]), from="in", to="user")
    y <- grconvertY(c(0, di[2]), from="in", to="user")
    fig <- par("fig")
    x <- x[1] + (x[2] - x[1]) * fig[1:2]
    y <- y[1] + (y[2] - y[1]) * fig[3:4]
    x <- x[1] + strwidth(txt, cex=3) / 2
    y <- y[2] - strheight(txt, cex=3) / 2
    text(x, y, txt, cex=3)
    par(old.par)
}

