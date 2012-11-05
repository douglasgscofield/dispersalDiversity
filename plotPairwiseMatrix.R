# plotPairwiseMatrix.R

# Provide function for plotting pairwise diversity matrices as returned by the
# pmiDiversity() function.  Used during data analysis for Scofield et al
# American Naturalist, http://www.jstor.org/stable/10.1086/668202.  Requires as
# input one or more tables of counts in sites (rows) X sources (columns)
# format.

.plotPairwiseMatrix.Version = "0.1"

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
# 
# See Scofield et al. Am. Nat, Figure 4A-C to see examples of this function in use.
# 
#
# plotPairwiseMatrix() : Create a visual plot of pairwise divergence or overlap
#                        values as calculated by pmiDiversity().  For example, with
#                        tab defined as a sites-by-sources matrix as for pmiDiversity(),
#                        plot the divergence matrix based on r_gg calculations, 
#                        labelling the axes "Seed Pool", using the following code:
#
#                        pmiD = pmiDiversity(tab)
#                        plotPairwiseMatrix(pairwise.mat=pmiD$r.divergence.mat, 
#                                           pairwise.mean=pmiD$r.divergence, 
#                                           statistic="divergence", 
#                                           axis.label="Seed Pool")
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


library(lattice)


#---------------------------------------------


plotPairwiseMatrix <- function(pairwise.mat, 
                               pairwise.mean=NA, 
                               statistic=c("divergence","overlap"), 
                               axis.label="Seed Pool",
                               ...)
{
  statistic = match.arg(statistic)
  # the matrix to be passed in is found at pmiDiversity()$r.diversity.mat 
  diag(pairwise.mat) <- 0
  pairwise.mat[upper.tri(pairwise.mat)] <- 0
  rotateMatrix = function(mat) t(mat[nrow(mat):1,,drop=FALSE])
  pairwise.mat = rotateMatrix(pairwise.mat)
  par(mar=c(0,0,0,5), ps=10, xpd=NA)
  lp <- levelplot(pairwise.mat, 
                  at=c(0,0.01,seq(0.2,0.8,0.2),0.99,1.0),
                  ### Uncomment below for the 5X figure containing the horizontal scale and delta_gh
                  ### colorkey=list(space="top", labels=list(cex=1.0)),
                  colorkey=list(labels=list(cex=1.0)),
                  aspect="iso", 
                  bty="L",
                  regions=TRUE,
                  col.regions=function(.x) gray(c(1,seq(.9,.6,length.out=(.x-2)),0)),
                  xlab=list(axis.label, cex=1.0),
                  ylab=list(axis.label, cex=1.0),
                  scales=list(draw=FALSE, tck=0, cex=0.7, x=list(rot=90))
                  )
  plot(lp, ...)
  if (! is.na(pairwise.mean)) {
    lims <- lapply(current.panel.limits(), round)
    rr <- abs(unlist(lapply(lims, diff)))
    if (statistic == "divergence") {
      panel.text(x=0.45*rr[1], y=0.7*rr[2],
                substitute(bar(delta)==OBS, list(OBS=round(pairwise.mean,3))),
                adj=c(0,0), cex=1.0)
      ### Uncomment below for the 5X figure containing the scale and delta_gh
      ### panel.text(x=0.45*rr[1], y=0.7*rr[2], expression(delta[italic(gh)]),
      ###           adj=c(0,0), cex=1.0)
    } else {
      panel.text(x=0.45*rr[1], y=0.7*rr[2], #expression(omega[italic(gh)]),
                substitute(bar(omega)==OMEGA, list(OMEGA=round(pairwise.mean,3))),
                adj=c(0,0), cex=1.0)
    }
  }
}


#---------------------------------------------

