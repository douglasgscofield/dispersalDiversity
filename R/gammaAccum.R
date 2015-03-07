#' @include diversity.R
# for collation order
NULL



#' Plot the gamma accumulation result from \code{runGammaAccum}
#'
#' For an example of its use see Figure 4D-F of Scofield et al. (2012).
#'
#' @param gamma.accum Result from \code{\link{runGammaAccum}}
#' 
#' @param plot.xmax   Maximum extent of the X-axis in the plot
#' 
#' @param plot.ymax   Maximum extent of the Y-axis in the plot
#' 
#' @param obs.omega   If provided, the observed omega value for the
#' data underlying the gamma accumulation curve; this value is added to 
#' the plot
#' 
#' @param omega.positions  If \code{obs.omega} is given, the location at
#' which it is plotted, relative to \code{plot.xmax} and \code{plot.ymax},
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
#
# @references
#
# Scofield, D. G., Smouse, P. E., Karubian, J. and Sork, V. L. (2012)
# Use of alpha, beta and gamma diversity measures to characterize seed
# dispersal by animals.  \emph{American Naturalist} 180:719-732.
#
# @examples
#
# rga.result = runGammaAccum(tab)
# plotGammaAccum(rga.result)
#
#' @seealso \code{\link{runGammaAccum}}
#'
#' @export plotGammaAccum
#'
plotGammaAccum <- function(gamma.accum, 
                           plot.xmax = length(gamma.accum$simple.results$mns), 
                           plot.ymax = max(gamma.accum$simple.results$mns + 
                                           gamma.accum$simple.results$SE), 
                           obs.omega = NULL,
                           omega.positions = c(0.15, 0.97),
                           xlab = "Number of seed pools",
                           ylab = expression("Accumulated  " * gamma * 
                                             "  diversity"),
                           pch = 19, pch.gamma = 24,
                           col = "black", bty = "L", lwd = 1.5, bg = "white",
                           lty = 2,
                           ...) {
    gtr <- gamma.accum$simple.results
    obs.gamma <- gamma.accum$obs.gamma
    xmax <- length(gtr$mns)
    x <- 1:xmax
    #lim <- c(1, max(xmax + 3, gtr$mns + gtr$SE))
    #xlim <- ylim <- lim
    xlim <- c(1, plot.xmax)
    ylim <- c(1, plot.ymax)

    opa <- par(mar = c(2.7, 2.7, 0.5, 0.5), mgp = c(1.7, 0.4, 0), 
               las = 1, ps = 10, tcl = -0.4, xpd = NA)
    plot.separate.gamma <- FALSE
    if (abs(gtr$mns[xmax] - obs.gamma) < 0.1) {
        # include gamma in means sequence
        pch <- c(rep(pch, xmax - 1), pch.gamma)
    } else {
        plot.separate.gamma <- TRUE
    }
    plot(x, gtr$mns, pch = pch, col = col, bty = bty, xlim = xlim, 
         ylim = ylim, lwd = lwd, xlab = xlab, ylab = ylab, bg = bg, ...)
    lines(x, gtr$mns + gtr$SE, lty = lty, ...)
    lines(x, gtr$mns - gtr$SE, lty = lty, ...)
    if (plot.separate.gamma)
        points(xmax, obs.gamma, pch = pch.gamma, lwd = lwd, bg = bg, 
               cex = 1.2, ...)
    if (! is.null(obs.omega)) {
        text(x = (omega.positions[1] * xlim[2]) + 1, 
             y = (omega.positions[2] * ylim[2]) + 1,
             substitute(bar(omega) == OMEGA, list(OMEGA = obs.omega)),
             cex = 1.2)
    }
    par(opa)
}



#' Perform gamma diversity accumulation on site-by-source data
#'
#' Perform a gamma diversity accumulation on the site-by-source data in tab.
#' Several arguments control the method of accumulation and value of gamma
#' calculated.  Only the defaults have been tested; the others were developed
#' while exploring the data and must be considered experimental.  The result is
#' returned in a list, which may be passed to \code{plotGammaAccum} to plot the
#' result.
#'
#' @param tab           Site-by-source table, as passed to diversity()
#'
#' @param accum.method  \code{"random"} (default) or \code{"proximity"}.  
#' Method for accumulating sites.  If \code{"random"}, then the next site is
#' chosen at random. If \code{"proximity"}, then the next site chosen is the
#' next closest site, and \code{distance.file} must be supplied.
#'
#' @param resample.method  \code{"permute"} (default) or \code{"bootstrap"}.
#' Whether to resample sites without replacement (\code{"permute"}) or with
#' replacement (\code{"bootstrap"}).
#'
#' @param distance.file  A file or data.frame containing three columns of data,
#' with the header/column names being \code{"pool"}, \code{"X"}, and \code{"Y"},
#' containing the spatial locations of the seed pools named in the row names of
#' tab.  Only used if the \code{accum.method=="proximity"}.
#'
#' @param gamma.method Calculate gamma using \code{"r"} (default),
#' \code{"q"} or \code{"q.nielsen"} method (see paper).
#'
#' @param \dots Additional parameters passed to \code{runGammaAccumSimple}
#
# TODO does runGammaAccumSimple need to be exposed?  does it need more options?
#
#'
#' @seealso \code{\link{plotGammaAccum}}, \code{\link{runGammaAccumSimple}}
#'
# @examples
#'
#' @export runGammaAccum
#'
runGammaAccum <- function(tab, 
                          accum.method = c("random", "proximity"),
                          resample.method = c("permute", "bootstrap"),
                          gamma.method = c("r", "q.nielsen", "q"),
                          distance.file = NULL,
                          ...)
{
    accum.method <- match.arg(accum.method)
    resample.method <- match.arg(resample.method)
    gamma.method <- match.arg(gamma.method)
    pmiD <- diversity(tab)
    ans <- list()
    ans$obs.gamma <- pmiD[[paste(sep="", "d.gamma.", gamma.method)]]
    ans$obs.omega.mean <- pmiD[[paste(sep="", gamma.method, ".overlap")]]
    ans$obs.delta.mean <- pmiD[[paste(sep="", gamma.method, ".divergence")]]
    ans$simple.results <- runGammaAccumSimple(tab,
        accum.method = accum.method, resample.method = resample.method,
        gamma.method = gamma.method, distance.file = distance.file, ...)
    return(ans)
}



#
# gammaAccum()            : workhorse function for gamma accumulation
# gammaAccumStats()       : extracts stats from result of gammaAccum()
# runGammaAccumSimple()   : wrapper that runs and then return stats from gammaAccum()
#
#
#
# CHANGELOG
#
# 0.1.3: Fix typo and simplify gammaAccum() by adding local functions
# 0.1.2: Add calculation of gamma via Nielsen et al. transform of q_gg
# 0.1.1: Minor bugfix for runGammaAccum arguments
# 0.1: First release
#
#
# TODO
#
# * Turn this into an actual R package



#---------------------------------------------


runGammaAccumSimple <- function(tab, ...)
{
  return(gammaAccumStats(gammaAccum(tab, ...)))
}


#---------------------------------------------


gammaAccumStats <- function(gammaAccum.results)
{
  f = function(x, FUN) {
    recip = 1/x; rmax = max(recip[recip != Inf]); recip[recip == Inf] = rmax; FUN(recip)
  }
  mns <- unlist(lapply(gammaAccum.results, f, mean))
  vars <- unlist(lapply(gammaAccum.results, f, var))
  SE <- sqrt(vars)
  quants <- lapply(gammaAccum.results, 
                   function(x) quantile(1/x, c(0.05, 0.5, 0.95)))
  ####
  return(list(mns=mns, vars=vars, SE=SE, quants=quants))
}


#---------------------------------------------


gammaAccum <- function(tab, 
                       n.sites=dim(tab)[1],
                       n.resample=1000,
                       accum.method=c("random", "proximity"),
                       resample.method=c("permute", "bootstrap"),
                       distance.file=NULL,
                       gamma.method=c("r", "q.nielsen", "q"))
{
  # If used, the distance.file has three columns, with names: pool, X, Y
  # It is either a filename to read, or a data.frame
  accum.method <- match.arg(accum.method)
  resample.method <- match.arg(resample.method)
  gamma.method <- match.arg(gamma.method)
  pool.names <- dimnames(tab)[[1]]
  G <- dim(tab)[1]
  K <- dim(tab)[2]
  N <- sum(tab)
  if (accum.method == "proximity") {
      if (! is.null(distance.file))
          gxy = .loadDistanceFile(distance.file, pool.names)
      else stop("distance.file must be provided")
  }

  ans = lapply(1:n.sites, function(x) numeric(0))
  for (i in 1:n.resample) {
    row.order <- sample(1:G,
                        size=n.sites,
                        replace=ifelse(resample.method=="bootstrap", TRUE, FALSE))
    if (accum.method == "proximity")
      row.order = .reorderByProximity(row.order, gxy, pool.names)
    for (g in 1:n.sites) {
      this.gamma = .calculateGammaAccum(apply(tab[row.order[1:g], , drop=FALSE], 2, sum),
                                        gamma.method)
      ans[[g]] = c(ans[[g]], this.gamma)
    }
  }
  return(ans)
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


.reorderByProximity = function(row.order, gxy, pool.names)
{
  new.row.order <- c()
  pools <- gxy[pool.names[row.order], ]
  dm <- .createSpatialDistmat(pools$X, pools$Y, names.1=pools$pool)
  diag(dm) <- +Inf  # distance to same pool is always farthest
  next.row.i <- 1
  new.row.order <- c(new.row.order, row.order[next.row.i])
  dm[, next.row.i] <- Inf  # distance to same pool is always farthest
  for (r in 2:n.sites) {
    next.row.i <- which(dm[next.row.i, ] == min(dm[next.row.i, ]))[1]
    dm[, next.row.i] <- Inf
    new.row.order <- c(new.row.order, row.order[next.row.i])
  }
  new.row.order
}


.createSpatialDistmat <- function(x.1, y.1, x.2, y.2, names.1, names.2)
{
    if (! missing(names.1)) names(x.1) = names.1
    if (missing(x.2)) x.2 <- x.1
    if (! missing(names.2)) names(x.2) = names.2
    if (missing(y.2)) y.2 <- y.1
    xx <- outer(x.1, x.2, "-")
    xx <- xx * xx
    yy <- outer(y.1, y.2, "-")
    yy <- yy * yy
    ans <- sqrt(xx + yy)
    dimnames(ans) <- list(names(x.1), names(x.2))
    ans
}


.loadDistanceFile = function(distance.file, pool.names)
{
  if (is.character(distance.file))
    gxy <- read.delim(file=distance.file)
  else if ("data.frame" %in% class(distance.file))
    gxy <- distance.file
  else
    stop("unknown type of distance information for", 
        deparse(substitute(distance.file)), " with class ", 
        paste(sep=", ", class(distance.file)))
  subset(gxy, pool %in% pool.names)
}

