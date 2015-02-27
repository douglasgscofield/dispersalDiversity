# alleleGammaAccum.R
#
# Copyright (c) 2015 Douglas G. Scofield, Uppsala University
#
# douglas.scofield@ebc.uu.se
# douglasgscofield@gmail.com
#
# These statistical tools were developed in collaboration with Peter Smouse
# (Rutgers University) and Victoria Sork (UCLA) and were funded by U.S. National
# Science Foundation awards NSF-DEB-0514956 and NSF-DEB-0516529.
#
# Use as you see fit.  No warranty regarding this code is implied nor should be
# assumed.  Send bug reports etc. to one of the above email addresses.

# Provide functions for calculating gamma accumulation across sites for 
# multiple loci.  Based on concepts developed in Scofield et al 2012 American
# Naturalist 180(6) 719-732, http://www.jstor.org/stable/10.1086/668202.
# Requires as input one list of tables produced by allele.createTableList().

# A typical workflow might look like
#
#    dat = readGenalex("genotypes.txt")
#    lst = allele.createTableList(dat)
#    allele.rga.result = allele.runGammaAccum(lst)
#    plotGammaAccum(allele.rga.result)

.alleleGammaAccum.Version = "0.1"

# FUNCTIONS PROVIDED 
# 
# allele.runGammaAccum(lst): Perform a gamma diversity accumulation on the
#                           site-by-source data in tab.  Several arguments
#                           control the method of accumulation and value of
#                           gamma calculated.  Only the defaults have been
#                           tested; the others were developed while exploring
#                           the data and must be considered experimental.
#                           The result is returned in a list, which may be
#                           passed to plotGammaAccum() to plot the result.
#
#       Arguments:
#           lst 
#                 List of site-by-source tables for loci, as produced by
#                 allele.createTableList()
#           accum.method
#                 "random" (default) or "proximity".  If "proximity" is used,
#                 then distance.file must be supplied (see below)
#           resample.method
#                 "permute" (default) or "bootstrap"; whether to resample sites
#                 without ("permute") or with ("bootstrap") replacement
#           n.resample
#                 Number of resamples for accumulation curve (default 1000)
#           distance.file
#                 A file or data.frame containing three columns of data, with
#                 the header/column names being "pool", "X", and "Y",
#                 containing the spatial locations of the seed pools named in
#                 the row names of tab.  Only used if the
#                 accum.method=="proximity"
#           gamma.method
#                 Calculate gamma using "r" (default), "q" or "q.nielsen" method
#                 (see Scofield et al. 2012).  The allele.* functions are all based
#                 on "r" versions and undefined results may happen otherwise.
#
# allele.gammaAccum()                : workhorse function for gamma accumulation
# allele.runAlleleGammaAccumSimple() : wrapper that runs and then return stats from allele.gammaAccum()
#
#
#
# CHANGELOG
#
# 0.1: First release
#
#
# TODO
#
# * Turn this into an actual R package


# uses functions from both these files
source("allelePmiDiversity.R")
source("gammaAccum.R")


#---------------------------------------------


allele.runGammaAccum <- function(lst, 
                          accum.method=c("random", "proximity"),
                          resample.method=c("permute", "bootstrap"),
                          n.resample=1000,
                          gamma.method=c("r", "q.nielsen", "q"),
                          distance.file=NA,
                          ...)
{
  accum.method = match.arg(accum.method)
  resample.method = match.arg(resample.method)
  gamma.method = match.arg(gamma.method)
  if (gamma.method != "r")
    warning("allele.runGammaAccum() should be used with gamma.method=='r'", immediate.=TRUE)
  pmiD = allele.pmiDiversity(lst)
  ans = list()
  ans$obs.gamma <- pmiD$gamma
  ans$obs.omega.mean <- pmiD$overlap
  ans$obs.delta.mean <- pmiD$divergence
  ans$simple.results <- allele.runGammaAccumSimple(lst,
                                            accum.method=accum.method,
                                            resample.method=resample.method,
                                            n.resample=n.resample,
                                            gamma.method=gamma.method,
                                            distance.file=distance.file,
                                            ...)
  ####          
  ans
}


#---------------------------------------------


allele.runGammaAccumSimple <- function(lst, ...)
{
  return(gammaAccumStats(allele.gammaAccum(lst, ...)))
}


#---------------------------------------------


allele.gammaAccum <- function(lst, 
                       n.sites=dim(lst[[1]])[1],
                       accum.method=c("random", "proximity"),
                       resample.method=c("permute", "bootstrap"),
                       n.resample=1000,
                       distance.file=NA,
                       gamma.method=c("r", "q.nielsen", "q"))
{
  # If used, the distance.file has three columns, with names: pool, X, Y
  # It is either a filename to read, or a data.frame
  accum.method <- match.arg(accum.method)
  resample.method <- match.arg(resample.method)
  gamma.method <- match.arg(gamma.method)
  pool.names <- dimnames(lst[[1]])[[1]]
  L = names(lst)
  G <- dim(lst[[1]])[1]
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
      gamma.all.loci = lapply(lst, 
                              function(x) {
                                .calculateGammaAccum(apply(x[row.order[1:g], , drop=FALSE], 2, sum),
                                                     gamma.method)
                              })
      ans[[g]] = c(ans[[g]], mean(unlist(gamma.all.loci)))
    }
  }
  return(ans)
}


#---------------------------------------------

