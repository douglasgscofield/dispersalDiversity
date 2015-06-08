#' @include diversity-divtable.R
NULL



#' Plot the gamma accumulation result from \code{gammaAccum}
#'
#' For an example of its use see Figure 4D-F of Scofield \emph{et al}. (2012).
#'
#' @param g  Object of class \code{'gamma_accum'} returned by
#' \code{\link{gammaAccum}}
#'
#' @param xmax   Maximum extent of the X-axis in the plot
#'
#' @param ymax   Maximum extent of the Y-axis in the plot
#'
#' @param obs.omega   If provided, the observed omega value for the
#' data underlying the gamma accumulation curve; this value is added to
#' the plot
#'
#' @param omega.pos  If \code{obs.omega} is given, the location at
#' which it is plotted, relative to \code{xmax} and \code{ymax},
#' in a two-element vector.  1 is added to each dimension.
#'
#' @param pch  Plotting character used for mean gamma values accumulated,
#' see \code{?pch}
#'
#' @param pch.gamma  Plotting character used for total gamma value,
#' \code{?pch}
#'
#' @param xlab,ylab,col,bty,lwd,bg  Options passed to \code{\link{plot}},
#' \code{\link{lines}} and/or \code{\link{points}}
#'
#' @param \dots  Additional options passed to \code{\link{plot}},
#' \code{\link{lines}} and/or \code{\link{points}}
#'
#' @return Nothing returned
#'
#' @references
#'
#' Scofield, D. G., Smouse, P. E., Karubian, J. and Sork, V. L. (2012)
#' Use of alpha, beta and gamma diversity measures to characterize seed
#' dispersal by animals.  \emph{American Naturalist} 180:719-732.
#'
#' @examples
#'
#' ## generate random divtable
#' n.sites <- 5
#' n.sources <- 10
#' n.samples <- 160
#' set.seed(75333)
#' t <- data.frame(site = sample(n.sites, n.samples, replace = TRUE),
#'                 source = round(runif(n.samples) * n.sources + 0.5))
#' tab <- as.divtable(table(t))
#' ## generate gamma accumulation results
#' rga.result = gammaAccum(tab)
#' ## plot gamma accumulation curve
#' plot(rga.result)
#' ## generate random allele_divtables
#' ##
#' ## plot gamma accumulation curve for allele data
#' ##
#'
#'
#' @seealso \code{\link{gammaAccum}}
#'
#' @export
#'
plot.gamma_accum <- function(g, xmax = length(g$simple.results$mns),
    ymax = max(g$simple.results$mns + g$simple.results$SE),
    obs.omega = NULL, omega.pos = c(0.15, 0.97),
    xlab = "Number of groups",
    ylab = expression("Accumulated  " * gamma * "  diversity"),
    pch = 19, pch.gamma = 24, col = "black", bty = "L", lwd = 1.5,
    bg = "white", lty = 2, ...)
{
    gtr <- g$simple.results
    obs.gamma <- g$obs.gamma
    axis.xmax <- length(gtr$mns)
    x <- 1:axis.xmax
    #lim <- c(1, max(axis.xmax + 3, gtr$mns + gtr$SE))
    #xlim <- ylim <- lim
    xlim <- c(1, xmax)
    ylim <- c(1, ymax)

    opa <- par(mar = c(2.7, 2.7, 0.5, 0.5), mgp = c(1.7, 0.4, 0),
               las = 1, ps = 10, tcl = -0.4, xpd = NA)
    plot.separate.gamma <- FALSE
    if (abs(gtr$mns[axis.xmax] - obs.gamma) < 0.1) {
        # include gamma in means sequence
        pch <- c(rep(pch, axis.xmax - 1), pch.gamma)
    } else {
        plot.separate.gamma <- TRUE
    }
    plot(x, gtr$mns, pch = pch, col = col, bty = bty, xlim = xlim,
         ylim = ylim, lwd = lwd, xlab = xlab, ylab = ylab, bg = bg, ...)
    lines(x, gtr$mns + gtr$SE, lty = lty, ...)
    lines(x, gtr$mns - gtr$SE, lty = lty, ...)
    if (plot.separate.gamma)
        points(axis.xmax, obs.gamma, pch = pch.gamma, lwd = lwd, bg = bg,
               cex = 1.2, ...)
    if (! is.null(obs.omega)) {
        text(x = (omega.pos[1] * xlim[2]) + 1,
             y = (omega.pos[2] * ylim[2]) + 1,
             substitute(bar(omega) == OMEGA, list(OMEGA = obs.omega)),
             cex = 1.2)
    }
    par(opa)
}



#' Perform gamma diversity accumulation on \code{divtable} or \code{allele_divtables} objects
#'
#' Perform a gamma diversity accumulation on site-by-source data in \code{tab},
#' an object of class \code{\link{divtable}}, or a set of site-by-alleles count
#' data for several loci, an object of classs \code{\link{allele_divtables}}.
#' Several arguments control the method of accumulation and value of gamma
#' calculated.  Only the defaults have been tested; the others were developed
#' while exploring the data and must be considered experimental.  The result is
#' returned in a list, which may be passed to \code{plotGammaAccum} to plot the
#' result.
#'
#' @param tab    Site-by-source table of class \code{\link{divtable}}
#'
#' @param adt    Allele diversity dataset of class
#' \code{\link{allele_divtables}}
#'
#' @param resample.method  \code{"permute"} (default) or \code{"bootstrap"}.
#' Whether to resample sites without replacement (\code{"permute"}) or with
#' replacement (\code{"bootstrap"}).
#'
#' @param n.resample Number of resamples for accumulation curve
#'
#' @param gamma.method Calculate gamma using \code{"q"} (default for class
#' \code{\link{divtable}}), \code{"r"} (default for class
#' \code{\link{allele_divtables}}) or \code{"q.nielsen"} method
#'
#' @param \dots Additional parameters, currently unused
#'
#' @note Only the defaults have been tested; other options were developed while
#' exploring the data and must be considered experimental.
#'
#' @references
#'
#' Scofield, D. G., Smouse, P. E., Karubian, J. and Sork, V. L. (2012)
#' Use of alpha, beta and gamma diversity measures to characterize seed
#' dispersal by animals.  \emph{American Naturalist} 180:719-732.
#'
#' Sork, V. L., Smouse, P. E., Grivet, D. and Scofield, D. G. (In press)
#' Impact of asymmetric male and female gamete dispersal on allelic
#' diversity and spatial genetic structure in valley oak
#' (\emph{Quercus lobata} N\'{e}e).  \emph{Evolutionary Ecology}.
#'
#' @seealso \code{\link{plot.gamma_accum}}, \code{\link{gammaAccumSimple}}
#'
#' @examples
#'
#' ## generate random divtable
#' n.sites <- 5
#' n.sources <- 10
#' n.samples <- 160
#' set.seed(75333)
#' t <- data.frame(site = sample(n.sites, n.samples, replace = TRUE),
#'                 source = round(runif(n.samples) * n.sources + 0.5))
#' tab <- as.divtable(table(t))
#' ## generate gamma accumulation results
#' rga.result = gammaAccum(tab)
#'
#' ## generate random allele_divtables
#' ##
#' ##
#'
#' @export
#'
#' @name gammaAccum
#'
NULL

gammaAccum <- function(x, ...) UseMethod("gammaAccum")



#' @rdname gammaAccum
#'
#' @export
#'
gammaAccum.divtable <- function(tab,
    resample.method = c("permute", "bootstrap"), n.resample = 1000,
    gamma.method = c("q", "r", "q.nielsen"), ...)
{
    resample.method <- match.arg(resample.method)
    gamma.method <- match.arg(gamma.method)
    d <- diversity(tab)
    ans <- list()
    #TODO: account for new diversity return value?
    ans$obs.gamma <- d[[gamma.method]][["d.gamma"]]
    ans$obs.omega.mean <- d[[gamma.method]][["overlap"]]
    ans$obs.delta.mean <- d[[gamma.method]][["divergence"]]
    ans$simple.results <- gammaAccumSimple.divtable(tab,
        resample.method = resample.method, gamma.method = gamma.method)
    structure(ans, class = c('gamma_accum', 'list'))
}



gammaAccumSimple <- function(x, ...) UseMethod("gammaAccumSimple")


# .allele_divtables method in gammaAccum-allele_divtables.R
#
gammaAccumSimple.divtable <- function(tab, ...)
{
    return(gammaAccumStats(gammaAccumWorker.divtable(tab, ...)))
}



gammaAccumStats <- function(ga)
{
    f = function(x, FUN) {
        recip = 1/x
        rmax = max(recip[recip != Inf])
        recip[recip == Inf] = rmax
        FUN(recip)
    }
    mns <- unlist(lapply(ga, f, mean))
    vars <- unlist(lapply(ga, f, var))
    SE <- sqrt(vars)
    quants <- lapply(ga, function(x) quantile(1/x, c(0.05, 0.5, 0.95)))
    ####
    return(list(mns = mns, vars = vars, SE = SE, quants = quants))
}



gammaAccumWorker <- function(x, ...) UseMethod("gammaAccumWorker")

# .allele_divtables method in gammaAccum-allele_divtables.R
#
gammaAccumWorker.divtable <- function(tab, n.sites=dim(tab)[1],
    n.resample=1000, resample.method=c("permute", "bootstrap"),
    gamma.method=c("q", "r", "q.nielsen"), ...)
{
    resample.method <- match.arg(resample.method)
    gamma.method <- match.arg(gamma.method)
    pool.names <- dimnames(tab)[[1]]
    G <- dim(tab)[1]
    K <- dim(tab)[2]
    N <- sum(tab)
    ans <- lapply(1:n.sites, function(x) numeric(0))
    for (i in 1:n.resample) {
        row.order <- sample(1:G,
                            size = n.sites,
                            replace = ifelse(resample.method == "bootstrap",
                                             TRUE, FALSE))
        for (g in 1:n.sites) {
            this.gamma <- .calculateGammaAccum(apply(tab[row.order[1:g], , drop=FALSE], 2, sum),
                                              gamma.method)
            ans[[g]] <- c(ans[[g]], this.gamma)
        }
    }
    ans
}


#--------------------------------------------- local functions


.calculateGammaAccum = function(this.accum, gamma.method)
{
    sum.accum <- sum(this.accum)
    this.prop <- this.accum / sum.accum
    sum.prop <- sum(this.prop * this.prop)
    switch(gamma.method,
           q         = sum.prop,
           q.nielsen = nielsenTransform(sum.prop, sum.accum),
           r         = (((sum.accum * sum.prop) - 1) / (sum.accum - 1)))
}


