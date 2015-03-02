#' Create a visual plot of pairwise divergence or overlap as calculated by \code{pmiDiversity}
#'
#' Plot pairwise values using \code{levelplot} from the \code{lattice} 
#' package.  See Scofield et al. Am. Nat, Figure 4A-C to see examples 
#' of this function in use.
#'
#' @param pairwise.mat Pairwise divergence or overlap matrix as found at
#' \code{pmiDiversity()$r.diversity.mat}.  Currently, prior to plotting
#' this matrix has its diagonal and upper triangle zeroed, and is then
#' rotated prior to passing to \code{lattice::levelplot}
#'
#' @param statistic \code{"divergence"} or \code{"overlap"}, for which
#' statistic is being presented
#'
#' @param pairwise.mean Mean pairwise divergence or overlap as found at
#' \code{pmiDiversity()$r.divergence}.  If provided, this is added to
#' the plot in the upper triangle, with positions specified by
#' \code{mean.positions}.  The value is plotted together with
#' '$\bar{\delta} =$' if \code{statistic} is \code{"divergence"}, and
#' '$\bar{\omega} =$' if \code{statistic} is \code{"overlap"}.
#'
#' @param mean.positions If \code{pairwise.mean} is given, the relative
#' X and Y positions within the panel at which the value is plotted, in 
#' a two-element vector.  \code{adj = c(0, 0)} is used when plotting the
#' value.
#'
#' @param axis.label Label used for the X and Y axes, the X and Y axis
#' labels can be changed with \code{xlab} and \code{ylab}
#'
#' @param bty,aspect,col.regions,colorkey,at,scales,xlab,ylab Additional
#' plot options passed to \code{lattice::levelplot}
#'
#' @return The lattice plot is returned invisibly
#'
#' @examples
#' pmiD = pmiDiversity(tab)
#' plotPairwiseMatrix(pairwise.mat = pmiD$r.divergence.mat, 
#'                    pairwise.mean = pmiD$r.divergence, 
#'                    statistic = "divergence", 
#'                    axis.label = "Seed Pool")
#'
#' @seealso \code{\link{lattice::levelplot}}
#'
#' @export plotPairwiseMatrix
#'
# library(lattice)
#
plotPairwiseMatrix <- function(pairwise.mat, 
    statistic = c("divergence", "overlap"), 
    pairwise.mean = NULL, mean.position = c(0.45, 0.7),
    bty = "L", aspect = "iso", 
    col.regions = function(.x) gray(c(1, seq(.9, .6, length.out=(.x - 2)), 0)),
    colorkey = list(labels = list(cex = 1.0)),
    at = c(0, 0.01, seq(0.2, 0.8, 0.2), 0.99, 1.0),
    scales = list(draw = FALSE, tck = 0, cex = 0.7, x = list(rot = 90)),
    axis.label = "Species Pool",
    xlab = list(axis.label, cex=1.0), ylab = list(axis.label, cex=1.0),
                               ...) {
    # the matrix to be passed in is found at pmiDiversity()$r.diversity.mat 
    statistic = match.arg(statistic)
    # should I be doing this zeroeing?
    diag(pairwise.mat) <- 0
    pairwise.mat[upper.tri(pairwise.mat)] <- 0
    rotateMatrix = function(mat) t(mat[nrow(mat):1, , drop=FALSE])
    pairwise.mat = rotateMatrix(pairwise.mat)
    opa <- par(mar = c(0, 0, 0, 5), ps = 10, xpd = NA)
    lp <- lattice::levelplot(pairwise.mat, 
              bty = bty, aspect = aspect, 
              regions = TRUE, col.regions = col.regions, 
              colorkey = colorkey,
              at = at, 
              scales = scales,
              xlab = xlab, ylab = ylab,
              ...)
    ## this needs a lattice::, right?
    lattice::plot(lp, ...)
    if (! is.null(pairwise.mean)) {
        lims <- lapply(lattice::current.panel.limits(), round)
        # is it sapply that i should use instead of unlist lapply?
        rr <- abs(unlist(lapply(lims, diff)))
        xpr <- if (statistic == "divergence")
            substitute(bar(delta) == OBS, list(OBS = round(pairwise.mean, 3)))
        else substitute(bar(omega) == OBS, list(OBS = round(pairwise.mean, 3)))

        lattice::panel.text(xpr, 
                            x = mean.position[1] * rr[1], 
                            y = mean.position[2] * rr[2], 
                            adj = c(0, 0), cex = 1.0)
        #} else {
        #    lattice::panel.text(x=0.45*rr[1], y=0.7*rr[2], #expression(omega[italic(gh)]),
        #            substitute(bar(omega)==OBS, list(OBS=round(pairwise.mean,3))),
        #            adj=c(0,0), cex=1.0)
        #}
    }
    par(opa)
    invisible(lp)
}
