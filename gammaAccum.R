# gammaAccum.R

# Provide functions for calculating gamma accumulation across sites, and
# plotting the results.  Used during data analysis for Scofield et al 2012
# American Naturalist 180(6) 719-723,
# http://www.jstor.org/stable/10.1086/668202.  Requires as input one or more
# tables of counts in sites (rows) X sources (columns) format.

# A typical workflow might look like
#
#    rga.result = runGammaAccum(tab)  # where tab is site-by-source
#    plotGammaAccum(rga.result)

.gammaAccum.Version = "0.1.2"

# Copyright (c) 2012 Douglas G. Scofield, Umeå Plant Science Centre, Umeå, Sweden
#
# douglas.scofield@plantphys.umu.se
# douglasgscofield@gmail.com
#
# These statistical tools were developed in collaboration with Peter Smouse
# (Rutgers University) and Victoria Sork (UCLA) and were funded by U.S. National
# Science Foundation awards NSF-DEB-0514956 and NSF-DEB-0516529.
#
# Use as you see fit.  No warranty regarding this code is implied nor should be
# assumed.  Send bug reports etc. to one of the above email addresses.
#
#
# FUNCTIONS PROVIDED 
# 
# See Scofield et al. Am. Nat, Figure 4D-F to see figures derived from using
# these functions.
# 
#
# runGammaAccum(tab)      : Perform a gamma diversity accumulation on the
#                           site-by-source data in tab.  Several arguments
#                           control the method of accumulation and value of
#                           gamma calculated.  Only the defaults have been
#                           tested; the others were developed while exploring
#                           the data and must be considered experimental.
#                           The result is returned in a list, which may be
#                           passed to plotGammaAccum() to plot the result.
#
#       Arguments:
#           tab 
#                 Site-by-source table, as passed to pmiDiversity()
#           accum.method
#                 "random" (default) or "proximity".  If "proximity" is used,
#                 then distance.file must be supplied (see below)
#           resample.method
#                 "permute" (default) or "bootstrap"; whether to resample sites
#                 without ("permute") or with ("bootstrap") replacement
#           distance.file
#                 A file or data.frame containing three columns of data, with
#                 the header/column names being "pool", "X", and "Y",
#                 containing the spatial locations of the seed pools named in
#                 the row names of tab.  Only used if the
#                 accum.method=="proximity"
#           gamma.method
#                 Calculate gamma using "r" (default), "q" or "q.nielsen" method
#                 (see paper).
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
# 0.1.2: Add calculation of gamma via Nielsen et al. transform of q_gg
# 0.1.1: Minor bugfix for runGammaAccum arguments
# 0.1: First release
#
#
# TODO
#
# * Turn this into an actual R package


source("pmiDiversity.R")


#---------------------------------------------


plotGammaAccum <- function(gamma.accum, 
                           plot.xmax=NA, 
                           plot.ymax=NA,
                           add.regression=FALSE, 
                           obs.omega=NA, 
                           ...)
{
  gtr = gamma.accum$simple.results
  obs.gamma = gamma.accum$obs.gamma
  xmax <- length(gtr$mns)
  x <- 1:xmax
  lim <- c(1, max(xmax + 3, gtr$mns + gtr$SE))
  xlim <- ylim <- lim
  xlim <- c(1, ifelse(is.na(plot.xmax), xmax, plot.xmax))
  ylim <- c(1, ifelse(is.na(plot.ymax), max(gtr$mns + gtr$SE), plot.ymax))

  #par(mar=c(3.5, 3.5, 0.5, 0.5), mgp=c(2.0, 0.5, 0), las=1, xpd=NA)
  par(mar=c(2.7, 2.7, 0.5, 0.5), mgp=c(1.7, 0.4, 0), las=1, ps=10, tcl=-0.4, xpd=NA)

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
       pch=pch, col="black", bty="L",
       xlim=xlim,
       ylim=ylim,
       #asp=1,
       lwd=1.5,
       tcl=-0.3,
       xlab="Number of Seed Pools",
       ylab=expression("Accumulated  "*gamma*"  diversity"),
       bg="white",
       ...
       )
  lines(x, gtr$mns + gtr$SE, lty=2, ...)
  lines(x, gtr$mns - gtr$SE, lty=2, ...)
  if (plot.separate.gamma) {
    points(xmax, obs.gamma, pch=pch.gamma, lwd=1.5, bg="white", cex=1.2)
  }
  if (! is.na(obs.omega)) {
    text(x=(0.15*xlim[2])+1, y=(0.97*ylim[2])+1,
         substitute(bar(omega)==OMEGA, list(OMEGA=obs.omega)),
         cex=1.2)
  }
}


#---------------------------------------------


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
  mns <- unlist(lapply(gammaAccum.results, 
                       function(x) { recip = 1/x; rmax = max(recip[recip != Inf]); recip[recip == Inf] = rmax; ; mean(recip) } ))
  vars <- unlist(lapply(gammaAccum.results, 
                        function(x) { recip = 1/x; rmax = max(recip[recip != Inf]); recip[recip == Inf] = rmax; ; var(recip) } ))
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
  # If used, the distance.file has three columns, with the header line pool, X, Y
  accum.method <- match.arg(accum.method)
  resample.method <- match.arg(resample.method)
  gamma.method <- match.arg(gamma.method)
  pool.names <- dimnames(tab)[[1]]
  G <- dim(tab)[1]
  K <- dim(tab)[2]
  N <- sum(tab)
  if (accum.method == "proximity") {
    if (is.character(distance.file))
      gxy <- read.delim(file=distance.file)
    else if ("data.frame" %in% class(distance.file))
      gxy <- distance.file
    else
      stop("unknown type of distance information for", 
          deparse(substitute(distance.file)), " with class ", 
          paste(sep=", ", class(distance.file)))
    gxy <- subset(gxy, granary %in% pool.names)
  }

  ans <- list()
  for (g in 1:n.sites)
    ans[[g]] <- numeric(0)
  for (i in 1:n.resample) {
    row.order <- sample(1:G,
                        size=n.sites,
                        replace=ifelse(resample.method=="bootstrap", TRUE, FALSE))
    if (accum.method == "proximity") {
      new.row.order <- c()
      pools <- gxy[pool.names[row.order], ]
      dm <- CreateSpatialDistmat(pools$X, pools$Y, names.1=pools$pool)
      diag(dm) <- +Inf  # distance to same pool is always farthest
      next.row.i <- 1
      new.row.order <- c(new.row.order, row.order[next.row.i])
      dm[, next.row.i] <- Inf  # distance to same pool is always farthest
      for (r in 2:n.sites) {
        next.row.i <- which(dm[next.row.i, ] == min(dm[next.row.i, ]))[1]
        dm[, next.row.i] <- Inf
        new.row.order <- c(new.row.order, row.order[next.row.i])
      }
      row.order <- new.row.order
    }
    for (g in 1:n.sites) {
      this.accum <- apply(tab[row.order[1:g], , drop=FALSE], 2, sum)
      sum.accum <- sum(this.accum)
      this.prop <- this.accum / sum.accum
      sum.prop <- sum(this.prop * this.prop)
      ans[[g]] <- c(ans[[g]], switch(gamma.method,
        q         = sum.prop,
        q.nielsen = nielsenTransform(sum.prop, sum.accum),
        r         = (((sum.accum * sum.prop) - 1) / (sum.accum - 1))))
    }
  }
  return(ans)
}


#---------------------------------------------


