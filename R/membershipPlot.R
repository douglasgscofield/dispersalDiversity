#' Plot a series of bars or pies indicating membership in groups
#'
#' Singleton groups, those with just one member, appear as white regions
#' within their respective sites.  Multiton groups, with more than one
#' member but all members appearing within a single site, also appear as
#' white regions unless \code{fill.method = "color"} with
#' \code{distinguish.multiton = TRUE}, then they are grey.  Groups with 
#' members in more than one site are given a unique colour/hatching 
#' combination.  An example if its use is in Figure 2A-C of
#' Scofield et al. (2012).
#'
#' @param tab Matrix of counts or proportions for observations of 
#' group membership (columns) within sites (rows).  The table is 
#' standardized by site so it doesn't matter whether the
#' table values are counts or proportions.
#'
#' @param method \code{"bar"} (the default) or \code{"pie"}, for the form of
#' plot produced.  Bar plots have received considerably more attention than 
#' pie plots.
#'
#' @param fill.method \code{"color"} (the default, with synonym 
#' \code{"colour"}, or \code{"bw"}, for the method of colouring different 
#' groups within the plots.  For \code{"color"}, no more than eight separate
#' colours are chosen; if there are more than eight groups, the
#' colours for lower-frequency groups are modified by hatching lines of
#' varying angles and densities.  For \code{"bw"}, fewer grades of grey are
#' used, together with black-and-white patterns and hatching.
#'
#' For \code{"color"}, if the package \code{\link{RColorBrewer}} is 
#' available, colours will be chosen from its \code{"Dark2"} palette.
#' This may be modified with the \code{fill.pallete} option.  The aesthetics
#' of color choice are important for distinguishing among groups in 
#' membership plots.  If \code{\link{RColorBrewer}} is not available, colours
#' are chosen using \code{\link{rainbow}}.
#'
#' @param fill.palette For \code{fill.method = "color"}, if the 
#' package \code{\link{RColorBrewer}} is available, use this palette to
#' choose the base colours to use
#'
#' @param distinguish.multiton For \code{fill.method = "color"}, whether 
#' to distinguish multiton groups (see Description) with grey rather than 
#' the default white colour.
#'
#' @param pdf.file,postscript.file Filename for output of PDF file, using
#' \code{\link{pdf}}, or PostScript file, using \code{\link{postscript}}.
#' file.  Only one can be set.  Both are produced using the options
#' \code{onefile = FALSE, paper = "special"}.
#'
#' @param file.dim Dimensions of the plot, in inches
#'
#' @param xlab,ylab,las,cex Plot parameters (see \code{\link{par}}
#'
#' @param x.mar.width  Width of the X-margin, in a bar plot
#'
#' @param l1,l2 First and second lines of annotation to appear above
#' bars in a bar plot
#'
#' @param header.cex Character expansion factor for plot header
#'
#' @param cex.names Character expansion factor for the names plotted below
#' the bars in a bar plot
#'
#' @return No value is returned
#
# @references
#
# Scofield, D. G., Smouse, P. E., Karubian, J. and Sork, V. L. (2012)
# Use of alpha, beta and gamma diversity measures to characterize seed
# dispersal by animals.  \emph{American Naturalist} 180:719-732.
#
# @example
#
#
#
#
#
#
#
#'
#' @export membershipPlot
#'
membershipPlot <- function(tab, method = c("bar", "pie"), 
    fill.method = c("color", "bw", "colour"), fill.palette = "Dark2",
    distinguish.multiton = FALSE, 
    pdf.file = NULL, postscript.file = NULL, file.dim = c(5.25, 2), 
    xlab = "Site", ylab = "Membership", las = 1, 
    x.mar.width = ifelse(las >= 2, 4, 3), 
    cex = 0.7, header.cex = 0.9, cex.names = 0.85, 
    l1 = NULL, l2 = NULL, ...) {
    
    .eps <- function(...) {
        postscript(onefile = FALSE, horizontal = FALSE, paper = "special", ...)
    }
    .pdf <- function(...) {
        pdf(onefile = FALSE, paper = "special", ...)
    }
    
    method <- match.arg(method)
    fill.method <- match.arg(fill.method)
    if (fill.method == "colour") fill.method <- "color"
    stopifnot(is.null(pdf.file) || is.null(postscript.file))
    to.file <- !is.null(pdf.file) || !is.null(postscript.file)
    
    # data comes in 'backwards' to what we expect, correct the order
    tab <- tab[rev(1:nrow(tab)), ]
    tab <- t(tab)
    
    # Generate colors: white for singleton types, grayscale for types that
    # appear in a single site, color for multi-site types
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
        
        # White for singleton types
        # Grayscale for types that appear in a single site,
        #     if distinguish.multiton, else white
        # Colour for multi-site types, use RColorBrewer if available

        fill.colors[singleton.idx] <- "white"
        fill.colors[multiton.idx] <- if (distinguish.multiton) 
            gray(seq(from = 0.5, to = 0.8, length = sum(multiton.idx)))
        else "white"
        
        n.cols <- min(sum(multisite.idx), 8)
        # if RColorBrewer is available, use it for colours, much nicer,
        colors <- if (requireNamespace("RColorBrewer", quietly = TRUE))
            colors <- RColorBrewer::brewer.pal(n.cols, fill.palette)
        else rainbow(n.cols)
        ang <- dens <- numeric(nrow(tab))
        ang[!multisite.idx] <- 0
        dens[!multisite.idx] <- -1
        fill.colors[!multisite.idx] <- "white"
        lo <- sum(multisite.idx) - n.cols
        # calculate hatching angles and densities
        ang[multisite.idx] <- c(rep(90, n.cols), 
                                rep(c(45, 75, 105, 135), each = n.cols/2,
                                    length.out = lo))
        dens[multisite.idx] <- c(rep(-1, n.cols), 
                                 rep(c(40, 25), each = n.cols/2, 
                                     length.out = lo))
        fill.colors[multisite.idx] <- c(colors, rep(colors, length.out = lo))
        fill.args <- list(angle = ang, density = dens, col = fill.colors)
        
    } else if (fill.method == "bw") {

        # White for singleton types and for types that appear in a single site;
        # b+w patterns for multi-site types, following black, gray, very
        # slightly off vertical, 45 degrees right, very slightly off
        # horizontal, 45 degrees left.

        ang <- dens <- numeric(nrow(tab))
        ang[!multisite.idx] <- 0
        dens[!multisite.idx] <- -1
        fill.colors[!multisite.idx] <- "white"
        lo <- sum(multisite.idx) - 2
        ang[multisite.idx] <- c(90, 90, rep(c(90, 45, 135), length.out = lo))
        dens[multisite.idx] <- c(-1, -1, rep(c(40, 25), each = 3, 
                                             length.out = lo))
        fill.colors[multisite.idx] <- c("black", gray(0.5), rep(c(gray(0.5), 
            "black"), each = 3, length.out = lo))
        fill.args <- list(angle = ang, density = dens, col = fill.colors)
        
    }
    
    if (method == "bar") {
        
        # adjust to proportions so they sum to 1
        tab <- scale(tab, scale = apply(tab, 2, sum), center = FALSE)
        if (any(abs(apply(tab, 2, sum) - 1) > (.Machine$double.eps^0.5))) 
            stop("didn't scale tab correctly")
        
        if (to.file) {
            if (!is.null(pdf.file)) 
                .pdf(file = pdf.file, 
                     width = file.dim[1], height = file.dim[2])
            else .eps(file = postscript.file, 
                      width = file.dim[1], height = file.dim[2])
            par(mar = c(x.mar.width, 3, 3, 0), mgp = c(1.9, 0.4, 0), 
                bty = "n", xpd = NA, las = las, tcl = -0.25, cex = cex)
        } else par(mar = c(x.mar.width, 7, 3, 1), bty = "n", xpd = NA, 
                   las = las)
        
        tab <- tab[, rev(1:ncol(tab))]
        do.call("barplot", c(list(height = tab, space = 0.3, horiz = FALSE), 
            cex.names = cex.names, fill.args, ...))
        title(xlab = , ylab = ylab)
        xl <- seq(0.8, by = 1.3, length.out = max(length(l1), length(l2)))
        if (!is.null(l1)) 
            text(xl, 1.25, labels = l1, cex = header.cex)
        if (!is.null(l2)) 
            text(xl, 1.125, labels = l2, cex = header.cex)

        if (to.file)
            dev.off()
        
    } else if (method == "pie") {
        
        n.pie <- ncol(tab)
        n.pie.col <- 5
        n.pie.row <- ceiling(n.pie/n.pie.col)
        if (to.file) {
            if (!is.null(pdf.file)) 
                .pdf(file = pdf.file, 
                     width = file.dim[1], height = file.dim[2])
            else .eps(file = postscript.file, 
                      width = file.dim[1], height = file.dim[2])
        }
        par(mfrow = c(n.pie.row, n.pie.col), mar = c(1, 1, 1, 1), bty = "n", 
            xpd = NA, cex = cex)
        for (p in 1:n.pie) {
            ttl <- paste(sep = "", colnames(tab)[p], " (n=", sum(tab[, p]), 
                         ")")
            do.call("pie", c(list(x = tab[, p], radius = 0.95, main = ttl, 
                                  labels = ""), 
                             fill.args, ...))
        }
        if (to.file)
            dev.off()
        
    }
    
}
 
