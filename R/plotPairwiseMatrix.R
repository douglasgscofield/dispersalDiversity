#' Plot pairwise divergence or overlap as calculated by diversity
#'
#' Plot lower-triangular matrix of pairwise divergence or overlap values
#' using a colour scale and \code{levelplot} from the \code{lattice} package.
#' The mean value may be included.
#'
#' The default colour scheme used is grayscale, with white and black used
#' for values very near the extremes (within 0.01 of 0 and 1, respectively).
#' Options are available for modifying the colour scheme.
#'
#' @param x    Pairwise divergence or overlap matrix as found at,
#' e.g., \code{diversity()$q$diversity.mat}.  Currently, prior to plotting
#' this matrix has its diagonal and upper triangle zeroed, and is then
#' rotated prior to passing to \code{\link{lattice::levelplot}}.
#'
#' @param statistic    \code{"divergence"} or \code{"overlap"}, indicating
#' the statistic being supplied for plotting
#'
#' @param pairwise.mean    Mean pairwise divergence or overlap as found at,
#' e.g., \code{diversity()$q$divergence}.  If provided, this is added to
#' the plot in the upper triangle, rounded to three digits, with
#' positions specified by \code{mean.pos}.  The value is plotted
#' together with '\eqn{\bar{\delta} = }' if \code{statistic} is
#' \code{"divergence"}, and '\eqn{\bar{\omega} = }' if \code{statistic} is
#' \code{"overlap"}.
#'
#' @param mean.pos    If \code{pairwise.mean} is given, the relative
#' X and Y positions within the panel at which the value is plotted, in
#' a two-element vector.  \code{adj = c(0, 0)} is used when plotting the
#' value.
#'
#' @param axis.label    Label used for the X and Y axes, the X and Y axis
#' labels can be changed with \code{xlab} and \code{ylab}
#'
#' @param bty,aspect,col.regions,colorkey,at,scales,xlab,ylab    Additional
#' plot options passed to \code{\link{lattice::levelplot}}
#'
#' @param digits        Number of digits to use when printing mean value,
#' taken from the \code{"digits"} option if not supplied
#'
#' @return The class \code{\link{lattice}} plot object is returned invisibly
#'
#' @references
#'
#' Scofield, D. G., Smouse, P. E., Karubian, J. and Sork, V. L. (2012)
#' Use of alpha, beta and gamma diversity measures to characterize seed
#' dispersal by animals.  \emph{American Naturalist} 180:719-732.
#'
#' @examples
#'
#' data(granaries_2002_Qlob)
#' d <- diversity(granaries_2002_Qlob)
#' plotPairwiseMatrix(x = d$r$divergence.mat,
#'                    pairwise.mean = d$r$divergence,
#'                    statistic = "divergence",
#'                    axis.label = "Granary")
#'
#' @seealso \code{\link{diversity}}, \code{\link{lattice::levelplot}}
#'
#' @export plotPairwiseMatrix
#'
plotPairwiseMatrix <- function(x, statistic = c("divergence", "overlap"),
    pairwise.mean = NULL, mean.position = c(0.45, 0.3),
    bty = "L", aspect = "iso",
    col.regions = function(x) gray(c(1, seq(.9, .6, length.out=(x - 2)), 0)),
    colorkey = list(labels = list(cex = 1.0)),
    at = c(0, 0.01, seq(0.2, 0.8, 0.2), 0.99, 1.0),
    scales = list(draw = FALSE, tck = 0, cex = 0.7, x = list(rot = 90)),
    axis.label = "Group",
    xlab = list(axis.label, cex=1.0), ylab = list(axis.label, cex=1.0),
    digits = getOption("digits"),
    ...)
{
    stopifnot(is.matrix(x))
    statistic <- match.arg(statistic)
    # Diagonal is the identity value
    diag(x) <- if (statistic == "divergence") 0 else 1
    # Zeroeing for plotting
    x[upper.tri(x)] <- 0
    rotateMatrix <- function(.x) t(.x[nrow(.x):1, , drop = FALSE])
    x = rotateMatrix(x)
    opa <- par(mar = c(0, 0, 0, 5), ps = 10, xpd = NA)
    lp <- lattice::levelplot(x, bty = bty, aspect = aspect,
        regions = TRUE, col.regions = col.regions, colorkey = colorkey,
        at = at, scales = scales, xlab = xlab, ylab = ylab, ...)
    plot(lp, ...)
    if (! is.null(pairwise.mean)) {
        lims <- lapply(lattice::current.panel.limits(), round)
        rr <- abs(sapply(lims, diff))
        val <- signif(pairwise.mean, digits)
        xpr <- if (statistic == "divergence")
            substitute(bar(delta) == OBS, list(OBS = val))
        else substitute(bar(omega) == OBS, list(OBS = val))
        lattice::panel.text(xpr, x = mean.position[1] * rr[1],
                            y = mean.position[2] * rr[2], adj = c(0, 0),
                            cex = 1.0)
    }
    par(opa)
    invisible(lp)
}

