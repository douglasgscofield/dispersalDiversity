#' @include pmiDiversity.R
# for collation order
NULL



# Provide functions for calculating gamma accumulation across sites, and
# plotting the results.  Used during data analysis for Scofield et al 2012
# American Naturalist 180(6) 719-732,
# http://www.jstor.org/stable/10.1086/668202.  Requires as input one or more
# tables of counts in sites (rows) X sources (columns) format.
#
# A typical workflow might look like
#
#    rga.result = runGammaAccum(tab)  # where tab is site-by-source
#    plotGammaAccum(rga.result)


# FUNCTIONS PROVIDED 
# 
# See Scofield et al. Am. Nat, Figure 4D-F to see figures derived from using
# these functions.
# 

#' Plot the gamma accumulation result from \code{runGammaAccum}
#'
#' See Scofield et al. Am. Nat, Figure 4D-F to see figures derived from using
#' these functions.
#'
#' @param gamma.accum Result from \code{\link{runGammaAccum}}
#' 
#' @param plot.xmax   If provided, the maximum extent of the X-axis
#' in the plot.  If missing, calculated from the data.
#' 
#' @param plot.ymax   If provided, the maximum extent of the Y-axis
#' in the plot.  If missing, calculated from the data.
#' 
#' @param obs.omega   If provided, the observed omega value for the
#' data underlying the gamma accumulation curve is added to the plot
#' 
#' @param xlab,ylab,col,bty,lwd,bg  Options passed to \code{\link{plot}},
#' \code{\link{lines}} and/or \code{\link{points}}
#'
#' @param \dots  Additional options passed to \code{\link{plot}},
#' \code{\link{lines}} and/or \code{\link{points}}
#'
#' @return Nothing returned
#'
#' @seealso \code{\link{runGammaAccum}}
#'
#' @export plotGammaAccum
#'
plotGammaAccum <- function(gamma.accum, plot.xmax, plot.ymax, obs.omega,
                           xlab = "Number of seed pools",
                           ylab = expression("Accumulated  " * gamma * 
                                             "  diversity"),
                           col = "black", bty = "L", lwd = 1.5, bg = "white",
                           lty = 2,
                           ...) {
  gtr <- gamma.accum$simple.results
  obs.gamma <- gamma.accum$obs.gamma
  xmax <- length(gtr$mns)
  x <- 1:xmax
  lim <- c(1, max(xmax + 3, gtr$mns + gtr$SE))
  xlim <- ylim <- lim
  xlim <- c(1, ifelse(missing(plot.xmax), xmax, plot.xmax))
  ylim <- c(1, ifelse(missing(plot.ymax), max(gtr$mns + gtr$SE), plot.ymax))

  #par(mar = c(3.5, 3.5, 0.5, 0.5), mgp = c(2.0, 0.5, 0), las = 1, xpd = NA)
  par(mar = c(2.7, 2.7, 0.5, 0.5), mgp = c(1.7, 0.4, 0), las = 1, ps = 10, tcl = -0.4, xpd = NA)

  pch.gamma <- 24
  pch.means <- 19
  plot.separate.gamma <- FALSE
  if (abs(gtr$mns[xmax] - obs.gamma) < 0.1) { # gamma in means sequence
    pch <- c(rep(pch.means, xmax - 1), pch.gamma)
  } else {
    pch <- pch.means
    plot.separate.gamma <- TRUE
  }
  plot(x, gtr$mns,
       pch = pch, col = col, bty = bty, xlim = xlim, ylim = ylim,
       lwd = lwd, xlab = xlab, ylab = ylab, bg = bg, ...)
  lines(x, gtr$mns + gtr$SE, lty = lty, ...)
  lines(x, gtr$mns - gtr$SE, lty = lty, ...)
  if (plot.separate.gamma) {
    points(xmax, obs.gamma, pch = pch.gamma, lwd = lwd, bg = bg, cex = 1.2)
  }
  if (! missing(obs.gamma)) {
    text(x = (0.15 * xlim[2]) + 1, y = (0.97 * ylim[2]) + 1,
         substitute(bar(omega) == OMEGA, list(OMEGA = obs.omega)),
         cex = 1.2)
  }
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
#' @param tab           Site-by-source table, as passed to pmiDiversity
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
#
# plotGammaAccum(rga.result)
#                         : plot the gamma accumulation result from runGammaAccum()
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



runGammaAccum <- function(tab, 
                          accum.method=c("random", "proximity"),
                          resample.method=c("permute", "bootstrap"),
                          gamma.method=c("r", "q.nielsen", "q"),
                          distance.file=NA,
                          ...)
{
  accum.method = match.arg(accum.method)
  resample.method = match.arg(resample.method)
  gamma.method = match.arg(gamma.method)
  pmiD = pmiDiversity(tab)
  ans = list()
  ans$obs.gamma <- pmiD[[paste(sep="", "d.gamma.", gamma.method)]]
  ans$obs.omega.mean <- pmiD[[paste(sep="", gamma.method, ".overlap")]]
  ans$obs.delta.mean <- pmiD[[paste(sep="", gamma.method, ".divergence")]]
  ans$simple.results <- runGammaAccumSimple(tab,
                                            accum.method=accum.method,
                                            resample.method=resample.method,
                                            gamma.method=gamma.method,
                                            distance.file=distance.file,
                                            ...)
  ####          
  ans
}


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
                       distance.file=NA,
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
  if (accum.method == "proximity")
    gxy = .loadDistanceFile(distance.file, pool.names)

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

