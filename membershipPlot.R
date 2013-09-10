membershipPlot <- function(tab, 
                           method=c("bar","pie"), 
                           fill.method=c("bw","color"),
                           distinguish.multiton=FALSE,
                           postscript.file=NULL, 
                           pdf.file=NULL, 
                           file.dim=c(5.25,2.00),
                           xlab="Seed cache",
                           ylab="Seed source proportion",
                           label.rotation=1,
                           x.mar.width=ifelse(label.rotation >= 2, 4, 3),  # method="bar"
                           l1=NULL,  # the first line of annotation to appear above the bars
                                     # assumed to be real (e.g. alpha values)
                           l2=NULL,  # the second line of annotation to appear above the bars
                                     # assumed to be integer (e.g., total counts)
                           ...)
{

  .eps <- function(...) {
      postscript(onefile=FALSE, horizontal=FALSE, paper="special", ...)
  }
  .pdf <- function(...) {
      pdf(onefile=FALSE, paper="special", ...)
  }

  method <- match.arg(method)
  fill.method <- match.arg(fill.method)
  stopifnot(label.rotation >= 0, label.rotation <= 3)
  stopifnot(is.null(postscript.file) || is.null(pdf.file))
  to.file = ! is.null(postscript.file) || ! is.null(pdf.file)

  # data comes in "backwards" to what we expect, correct the order
  tab <- tab[rev(1:nrow(tab)), ]
  tab <- t(tab)

  # Generate colors: white for singleton types, grayscale for
  # types that appear in a single site, color for multi-site 
  # types

  fill.colors <- rep("", nrow(tab))
  names(fill.colors) <- rownames(tab)
  n.types <- apply(tab, 1, sum)
  names(n.types) <- rownames(tab)

  # sort the types in descending order of number of appearances of type
  n.types <- rev(sort(n.types))
  tab <- tab[names(n.types), ]

  # number of sites in which each type can be found
  n.sites <- apply(tab, 1, function(x) sum(x > 0))
  singleton.idx <- n.types == 1
  multiton.idx <- n.sites == 1 & n.types > 1
  multisite.idx <- n.sites > 1

  if (fill.method == "color") {

    # White for singleton types; grayscale for types that appear in a
    # single site; color for multi-site types.

    fill.colors[singleton.idx] <- "white"

    fill.colors[multiton.idx] <- if (distinguish.multiton)
      gray(seq(from=0.5, to=0.8, length=sum(multiton.idx)))
    else
      "white"

    if (sum(multisite.idx) <= 8) 
      #fill.colors[multisite.idx] <- brewer.pal(sum(multisite.idx), "Set1")
      fill.colors[multisite.idx] <- brewer.pal(sum(multisite.idx), "Set1")
    else 
      fill.colors[multisite.idx] <- rainbow(sum(multisite.idx))

    fill.args <- list(col=fill.colors)

    # colors = c("red", "orange", "green", "blue" "indigo", "violet")
    n.cols = min(sum(multisite.idx), 8)
    colors = brewer.pal(n.cols, "Dark2")
    ang <- dens <- numeric(nrow(tab))
    ang[! multisite.idx] <- 0
    dens[! multisite.idx] <- -1
    fill.colors[! multisite.idx] <- "white"
    lo <- sum(multisite.idx) - n.cols 
    # ang[multisite.idx] <- c(rep(90, 6), rep(c(90,45,135), length.out=lo))
    # dens[multisite.idx] <- c(rep(-1, 6), rep(c(40,25), each=3, length.out=lo))
    ang[multisite.idx] <- c(rep(90, n.cols), rep(c(45, 75, 105, 135), each=n.cols/2, length.out=lo))
    dens[multisite.idx] <- c(rep(-1, n.cols), rep(c(40,25), each=n.cols/2, length.out=lo))
    fill.colors[multisite.idx] <- c(colors,
      rep(colors, length.out=lo))
    fill.args <- list(angle=ang, density=dens, col=fill.colors)

  } else if (fill.method == "bw") {

    # White for singleton types and for types that appear in a single site;
    # b+w patterns for multi-site types, following black, gray, very
    # slightly off vertical, 45 degrees right, very slightly off
    # horizontal, 45 degrees left.

    ang <- dens <- numeric(nrow(tab))
    ang[! multisite.idx] <- 0
    dens[! multisite.idx] <- -1
    fill.colors[! multisite.idx] <- "white"
    lo <- sum(multisite.idx) - 2 
    ang[multisite.idx] <- c(90, 90, rep(c(90,45,135), length.out=lo))
    dens[multisite.idx] <- c(-1, -1, rep(c(40,25), each=3, length.out=lo))
    fill.colors[multisite.idx] <- c("black", gray(0.5), 
      rep(c(gray(0.5), "black"), each=3, length.out=lo))
    fill.args <- list(angle=ang, density=dens, col=fill.colors)

  }

  par.cex = 0.7
  header.cex = 0.9
  names.cex = 0.85

  if (method == "bar") {

    # adjust to proportions so they sum to 1
    tab <- scale(tab, scale=apply(tab, 2, sum), center=FALSE)
    if (any(abs(apply(tab, 2, sum) - 1) > (.Machine$double.eps ^ 0.5))) 
      stop("didn't scale tab correctly")

    if (to.file) {
      if (! is.null(postscript.file))
        .eps(file=postscript.file, width=file.dim[1], height=file.dim[2])
      else 
        .pdf(file=pdf.file, width=file.dim[1], height=file.dim[2])
      par(mar=c(x.mar.width,3,3,0), mgp=c(1.9,0.4,0), 
          bty="n", xpd=NA, las=label.rotation, tcl=-0.25, cex=par.cex)
    } else {
      par(mar=c(x.mar.width,7,3,1), bty="n", xpd=NA, las=label.rotation)
    }

    tab <- tab[, rev(1:ncol(tab))]
    do.call("barplot", c(list(height=tab, space=0.3, horiz=FALSE),
                          cex.names=names.cex, fill.args, ...))
    title(xlab=, ylab=ylab)
    xl <- seq(0.8, by=1.3, length.out=length(l1))
    if (! is.null(l1)) {
      text(-0.6, 1.250, labels=expression(italic(alpha[g])==""), cex=header.cex)
      text(xl, 1.250, labels=sapply(l1, function(x) sprintf("%.1f",x)), cex=header.cex)
    }
    if (! is.null(l2)) {
      text(-0.6, 1.125, labels=expression(italic(n[g])==""), cex=header.cex)
      text(xl, 1.125, labels=l2, cex=header.cex)
    }
    if (! is.null(postscript.file) || ! is.null(pdf.file)) {
      dev.off()
    }

  } else if (method == "pie") {

    n.pie <- ncol(tab)
    n.pie.col <- 5
    n.pie.row <- ceiling(n.pie/n.pie.col)
    if (! is.null(postscript.file)) {
      .eps(file=postscript.file, width=file.dim[1], 
          height=file.dim[2])
    }
    par(mfrow=c(n.pie.row, n.pie.col), mar=c(1,1,1,1), bty="n", xpd=NA)
    for (p in 1:n.pie) {
      ttl <- paste(sep="",colnames(tab)[p]," (n=",sum(tab[,p]),")")
      do.call("pie", c(list(x=tab[,p], radius=0.95, main=ttl, labels=""),
                          fill.args, ...))
    }
    if (! is.null(postscript.file)) {
      dev.off()
    }

  }
}
