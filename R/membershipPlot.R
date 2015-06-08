#' Plot a series of stacked bars indicating membership
#'
#' Plot a series of stacked bars indicating membership of species across
#' sites, items across groups, or however the contents of \code{tab} should
#' be interpreted.  In the descriptions below we will treat \code{tab} as
#' describing species counts across sites.
#'
#' In the plot, species are ordered bottom to top in order of decreasing total
#' abundance.  Species which are present in more than one site are given a
#' unique colour/hatching combination, with hatching used only if the number of
#' species is greater than the number of readily distinguishable colours,
#' currently set to 8.  Multiton species, with more than one individual but all
#' individual appearing within a single site, are coloured grey when
#' \code{fill.method = "colour"}, and are coloured white when \code{fill.method
#' = "bw"} or \code{distinguish.multiton = FALSE}.  Singleton groups, those with
#' just one member, are coloured white within their respective sites.
#'
#' Examples of membership plots in publications can be found in
#' Scofield \emph{et al}. (2010), Scofield \emph{et al}. (2011) and Scofield
#' \emph{et al}. (2012).
#'
#' @param tab   Object of class \code{\link{divtable}} containing counts of
#' sites (rows) by species membership (columns), or any similar object
#' (\code{table}, \code{matrix}, \code{data.frame}, \code{xtabs}) that
#' can be converted with \code{\link{as.divtable}}.  This argument is
#' normalized to sum to 1 within each site so values may be counts or
#' proportions.
#'
#' @param fill.method   \code{"colour"} (the default, with synonym
#' \code{"color"}), or \code{"bw"}, for the method of colouring different
#' species within the plots.  For \code{"colour"}, no more than eight separate
#' colours are chosen; if there are more than eight species, the
#' colours for lower-frequency groups are modified by hatching lines of
#' varying angles and densities.  For \code{"bw"}, fewer grades of grey are
#' used, together with black-and-white patterns and hatching.
#' For \code{"colour"}, if the package \code{\link{RColorBrewer}} is
#' available, colours will be chosen from the pallete designated by the
#' \code{fill.pallete} option.
#' The aesthetics of colour choice are important for distinguishing among
#' species in membership plots.  If \code{\link{RColorBrewer}} is not
#' available, colours are chosen using \code{\link{rainbow}}.
#'
#' @param fill.palette   For \code{fill.method = "colour"}, if the
#' package \code{\link{RColorBrewer}} is available, use this palette to
#' choose the base colours to use.  Default is \code{"Dark2"}.
#'
#' @param distinguish.multiton   For \code{fill.method = "colour"}, whether
#' to distinguish multiton groups (see Description) with grey rather than
#' the default white colour.  Default is \code{TRUE}.
#'
#' @param file   If specified, produce PDF output to \code{file}, using
#' \code{\link{pdf}}.  Options \code{onefile = FALSE, paper = "special"}
#' are used.
#'
#' @param file.dim   If \code{file} is specified, dimensions of the PDF
#' plot in inches
#'
#' @param xlab,ylab,las,cex   Plot parameters (see \code{\link{par}}
#'
#' @param x.mar.width  Width of the X-margin
#'
#' @param l1,l2,annotate.cex,annotate   First and second lines of annotation
#' to appear above bars in a bar plot, the character expansion factor to use,
#' and the function to place them there.  Alternate annotations may be
#' enabled by redefining \code{annotate}.
#'
#' @param cex.names Character expansion factor for the names plotted below
#' the bars
#'
#' @return No value is returned
#'
#' @references
#'
#' Scofield, D. G.,  Sork, V. L. and Smouse, P. E. (2010) Influence of acorn
#' woodpecker social behaviour on transport of coast live oak
#' (\emph{Quercus agrifolia}) acorns in a southern California oak savanna.
#' \emph{Journal of Ecology} 98:561-572.
#'
#' Scofield, D. G., Alfaro, V. R., Sork, V. L., Grivet, D., Martinez, E.,
#' Papp, J., Pluess, A. R., Koenig, W. D. and Smouse, P. E. (2011) Foraging
#' patterns of acorn woodpeckers (\emph{Melanerpes formicivorus}) on valley
#' oak (\emph{Quercus lobata} N\'{e}e) in two California oak savanna-woodlands.
#' Oecologia 166:187-196.
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
#' membershipPlot(tab)
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
    distinguish.multiton = TRUE, xlab = "Site", ylab = "Membership",
    las = 1, x.mar.width = ifelse(las >= 2, 4, 3), cex = 0.7, cex.names = 0.85,
    l1 = NULL, l2 = NULL, annotate.cex = 0.9,
    annotate = function() {
        xl <- seq(0.8, by = 1.3, length.out = max(length(l1), length(l2)))
        if (!is.null(l1)) mtext(l1, at = xl, line = 0.8, cex = annotate.cex, xpd = NA)
        if (!is.null(l2)) mtext(l2, at = xl, line = 0.4, cex = annotate.cex, xpd = NA)
    },
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
        par(mar = c(x.mar.width, 4, 3, 1), bty = "n", xpd = NA, las = las)
    }

    tab <- tab[, rev(1:ncol(tab)), drop = FALSE]
    do.call("barplot", c(list(height = tab, space = 0.3, horiz = FALSE),
        cex.names = cex.names, fill.args, ...))
    title(xlab = , ylab = ylab)

    #
    if (! is.null(annotate))
        annotate()

    if (! is.null(file)) dev.off()
}

