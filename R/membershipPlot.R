#' Plot a series of stacked bars indicating membership in groups
#'
#' Singleton groups, those with just one member, appear as white regions
#' within their respective sites.  Multiton groups, with more than one
#' member but all members appearing within a single site, also appear as
#' white regions unless \code{fill.method = "colour"} with
#' \code{distinguish.multiton = TRUE}, then they are grey.  Groups with 
#' members in more than one site are given a unique colour/hatching 
#' combination.  An example if its use is in Figure 2A-C of
#' Scofield et al. (2012).
#'
#' @param tab   Object of class \code{\link{divtable}} of counts of sites
#' (rows) by group membership (columns), or any similar object
#' (\code{table}, \code{matrix}, \code{data.frame}, \code{xtabs}) that
#' can be converted with \code{\link{as.divtable}}.  This argument is
#' standardized by site so values may be counts or proportions.
#'
#' @param fill.method   \code{"colour"} (the default, with synonym 
#' \code{"color"}, or \code{"bw"}, for the method of colouring different 
#' groups within the plots.  For \code{"colour"}, no more than eight separate
#' colours are chosen; if there are more than eight groups, the
#' colours for lower-frequency groups are modified by hatching lines of
#' varying angles and densities.  For \code{"bw"}, fewer grades of grey are
#' used, together with black-and-white patterns and hatching.
#' For \code{"colour"}, if the package \code{\link{RColorBrewer}} is 
#' available, colours will be chosen from its \code{"Dark2"} palette.
#' This may be modified with the \code{fill.pallete} option.  The aesthetics
#' of colour choice are important for distinguishing among groups in 
#' membership plots.  If \code{\link{RColorBrewer}} is not available, colours
#' are chosen using \code{\link{rainbow}}.
#'
#' @param fill.palette   For \code{fill.method = "colour"}, if the 
#' package \code{\link{RColorBrewer}} is available, use this palette to
#' choose the base colours to use
#'
#' @param distinguish.multiton   For \code{fill.method = "colour"}, whether 
#' to distinguish multiton groups (see Description) with grey rather than 
#' the default white colour
#'
#' @param file   If specified, produce PDF output to \code{file}, using
#' \code{\link{pdf}}.  Options \code{onefile = FALSE, paper = "special"}
#' are used.
#'
#' @param file.dim   Dimensions of the PDF plot if \code{file} is specified,
#' in inches
#'
#' @param xlab,ylab,las,cex   Plot parameters (see \code{\link{par}}
#'
#' @param x.mar.width  Width of the X-margin, in a bar plot
#'
#' @param l1,l2   First and second lines of annotation to appear above
#' bars in a bar plot
#'
#' @param header.cex Character expansion factor for plot header
#'
#' @param cex.names Character expansion factor for the names plotted below
#' the bars in a bar plot
#'
#' @return No value is returned
#'
#' @references
#'
#' Scofield, D. G., Smouse, P. E., Karubian, J. and Sork, V. L. (2012)
#' Use of alpha, beta and gamma diversity measures to characterize seed
#' dispersal by animals.  \emph{American Naturalist} 180:719-732.
#'
#' @examples
#'
#' ## create table of random membership data
#' n.sites <- 5
#' n.sources <- 10
#' n.samples <- 160
#' ## data frame of site-source pairs
#' set.seed(75333)
#' t <- data.frame(site = sample(n.sites, n.samples, replace = TRUE),
#'                 source = round(runif(n.samples) * n.sources + 0.5))
#' ## create site-by-source table
#' tab <- as.divtable(table(t))
#' membershipPlot(tab, distinguish.multiton = TRUE)
#'
#' @export
#'
#' @name membershipPlot
#'
NULL

membershipPlot <- function(tab, ...) UseMethod("membershipPlot")



#' @rdname membershipPlot
#'
#' @export
#'
membershipPlot.default <- function(tab, ...)
{
    membershipPlot.divtable(as.divtable(tab), ...)
}



#' @rdname membershipPlot
#'
#' @export
#'
membershipPlot.divtable <- function(tab,
    fill.method = c("colour", "bw", "color"), fill.palette = "Dark2",
    distinguish.multiton = FALSE, xlab = "Site", ylab = "Membership",
    las = 1, x.mar.width = ifelse(las >= 2, 4, 3), cex = 0.7,
    header.cex = 0.9, cex.names = 0.85, l1 = NULL, l2 = NULL,
    file = NULL, file.dim = c(5.25, 2),
    ...)
{
    fill.method <- match.arg(fill.method)
    if (fill.method == "color") fill.method <- "colour"
    
    # data comes in 'backwards' to what we expect when organising for the plot
    tab <- tab[rev(1:nrow(tab)), , drop = FALSE]
    tab <- t(tab)
    # so now types are along the rows, groups along the columns
    
    # Generate colours: white for singleton types, grayscale for types that
    # appear in a single site, colour for multi-site types.  Sort the types in
    # descending order of number of appearances of type
    fill.colours <- rep("", nrow(tab))
    names(fill.colours) <- rownames(tab)
    n.types <- apply(tab, 1, sum)
    n.types <- rev(sort(n.types))
    tab <- tab[names(n.types), , drop = FALSE]
    # Number of sites in which each type can be found
    n.sites <- apply(tab, 1, function(x) sum(x > 0))
    # Indices into type categories
    singleton.idx <- n.types == 1
    multiton.idx <- n.sites == 1 & n.types > 1
    multisite.idx <- n.sites > 1
    
    if (fill.method == "colour") {
        
        fill.colours[singleton.idx] <- "white"
        fill.colours[multiton.idx] <- if (distinguish.multiton) 
            gray(seq(from = 0.5, to = 0.8, length = sum(multiton.idx)))
        else "white"
        ang <- dens <- numeric(nrow(tab))
        ang[! multisite.idx] <- 0
        dens[! multisite.idx] <- -1
        fill.colours[! multisite.idx] <- "white"
        n.multisite <- sum(multisite.idx)
        if (n.multisite) {  # we need colours
            n.cols <- min(sum(multisite.idx), 8)
            # if RColorBrewer is available, use it for colours, much nicer,
            colours <- if (requireNamespace("RColorBrewer", quietly = TRUE))
                RColorBrewer::brewer.pal(max(n.cols, 3), fill.palette)
            else rainbow(n.cols)
            lo <- n.multisite - n.cols
            # calculate hatching angles and densities
            ang[multisite.idx] <- c(rep(90, n.cols), 
                                    rep(c(45, 75, 105, 135), each = n.cols / 2,
                                        length.out = lo))
            dens[multisite.idx] <- c(rep(-1, n.cols), 
                                     rep(c(40, 25), each = n.cols / 2, 
                                         length.out = lo))
            fill.colours[multisite.idx] <- c(colours, rep(colours, length.out = lo))
        }
        fill.args <- list(angle = ang, density = dens, col = fill.colours)
        
    } else if (fill.method == "bw") {

        # White for singleton types and for types that appear in a single site;
        # b+w patterns for multi-site types, following black, gray, very
        # slightly off vertical, 45 degrees right, very slightly off
        # horizontal, 45 degrees left.
        ang <- dens <- numeric(nrow(tab))
        ang[!multisite.idx] <- 0
        dens[!multisite.idx] <- -1
        fill.colours[!multisite.idx] <- "white"
        n.multisite <- sum(multisite.idx)
        if (n.multisite) {
            lo <- sum(multisite.idx) - 2
            ang[multisite.idx] <- c(90, 90, rep(c(90, 45, 135), length.out = lo))
            dens[multisite.idx] <- c(-1, -1, rep(c(40, 25), each = 3, 
                                                 length.out = lo))
            fill.colours[multisite.idx] <- c("black", gray(0.5), rep(c(gray(0.5), 
                "black"), each = 3, length.out = lo))
        }
        fill.args <- list(angle = ang, density = dens, col = fill.colours)
        
    }
    
    # adjust to proportions so they sum to 1
    tab <- scale(tab, scale = apply(tab, 2, sum), center = FALSE)
    if (any(abs(apply(tab, 2, sum) - 1) > (.Machine$double.eps^0.5))) 
        stop("didn't scale tab correctly")
    
    if (! is.null(file)) {
        pdf(file = file, width = file.dim[1], height = file.dim[2],
            onefile = FALSE, paper = "special")
        par(mar = c(x.mar.width, 3, 3, 0), mgp = c(1.9, 0.4, 0), 
            bty = "n", xpd = NA, las = las, tcl = -0.25, cex = cex)
    } else {
        par(mar = c(x.mar.width, 7, 3, 1), bty = "n", xpd = NA, las = las)
    }
    
    tab <- tab[, rev(1:ncol(tab)), drop = FALSE]
    do.call("barplot", c(list(height = tab, space = 0.3, horiz = FALSE), 
        cex.names = cex.names, fill.args, ...))
    title(xlab = , ylab = ylab)
    xl <- seq(0.8, by = 1.3, length.out = max(length(l1), length(l2)))
    if (!is.null(l1)) 
        mtext(l1, at = xl, line = 1, cex = header.cex, xpd = NA)
    if (!is.null(l2)) 
        mtext(l2, at = xl, line = 0.5, cex = header.cex, xpd = NA)

    if (! is.null(file)) dev.off()
}
 
