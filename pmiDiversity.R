# pmiDiversity.R

# Calculate PMI statistics (Grivet et al. 2005, Scofield et al.  2010, 2011) as
# well as alpha, beta, gamma based on both q_qq and r_gg (Scofield et al
# American Naturalist, http://www.jstor.org/stable/10.1086/668202), as well as
# q_gg adjusted following Nielsen et al 2003.  Used during data analysis for
# Scofield et al Am Nat; earlier versions (pre-github) were used for Scofield
# et al 2010 J Ecol and for Scofield et al 2011 Oecologia.  The single argument
# is a table of counts in sites (rows) X sources (columns) format.  If the
# argument is not a matrix, it is converted to one, and rows are ordered
# numerically by rowname if the rownames are numeric.

.pmiDiversityVersion = "0.3"

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
# CHANGELOG
#
# 0.3: Pull Nielsen et al. out into separate function for use in diversity tests,
#      add collaborators/funding statement
# 0.2: Finalize Nielsen et al. calculations of diversity measures
# 0.1: First release


nielsenTransform = function(q.gg, n.g)
{  
    # better bias correction, Nielsen et al 2003 Mol Ecol
    (((q.gg*(n.g+1)*(n.g-2))+(3-n.g))/((n.g-1)^2))
}


pmiDiversity <- function(tab)
{
  .thisis = paste("pmiDiversity v", .pmiDiversityVersion, "Douglas G. Scofield, douglasgscofield@gmail.com")
  tab = as.matrix(tab)
  #if (! any(is.na(suppressWarnings(as.numeric(rownames(tab))))))
  #  tab <- tab[order(as.numeric(rownames(tab))), ]
  G <- dim(tab)[1]
  K <- dim(tab)[2]
  N <- sum(tab)
  n.g <- apply(tab, 1, sum)
  # PMI statistics (Grivet et al 2005 Mol Ecol, Scofield et al 2010 J Ecol)

  r.gg <- numeric(G)
  for (g in 1:G) {
    tmp.numer <- (tab[g,]*tab[g,]) - tab[g,]
    tmp.denom <- (n.g[g]*n.g[g]) - n.g[g]
    r.gg[g] <- sum(tmp.numer / tmp.denom)
  }
  # for q.gh (same as r.gh), create table of relative genotype frequencies
  reltab <- sweep(tab, 1, n.g, FUN="/")
  Q.mat <- reltab %*% t(reltab)
  q.gg <- diag(Q.mat) # note biased estimator is diagonal of this matrix
  names(r.gg) <- names(q.gg)
  q.nielsen.gg <- nielsenTransform(q.gg, n.g)
  q.bar.0 <- sum(n.g * n.g * q.gg) / sum(n.g * n.g)
  q.nielsen.bar.0 <- sum(n.g * n.g * q.nielsen.gg) / sum(n.g * n.g)
  r.bar.0 <- sum((n.g*n.g*r.gg) - (n.g*r.gg)) / sum((n.g*n.g) - n.g)
  # more pairwise PMI stuff
  q.gh <- lower.tri(Q.mat)
  # pooled PMI (Scofield et al 2010 J Ecol)
  q.0.gh <- prop.y.0.gh <- matrix(0,G,G)
  dimnames(q.0.gh) <- dimnames(prop.y.0.gh) <- dimnames(Q.mat)
  for (g in 1:(G-1)) {
    for (h in (g+1):G) {
      shared.k <- tab[g,] & tab[h,]
      y.gh <- prop.y.0.gh[g,h] <- sum(tab[g,shared.k])
      y.hg <- prop.y.0.gh[h,g] <- sum(tab[h,shared.k])
      prop.y.0.gh[g,h] <- y.gh / n.g[g]
      prop.y.0.gh[h,g] <- y.hg / n.g[h]
      q.0.gh[g,h] <- (y.gh * y.hg)/(n.g[g]*n.g[h])
    }
  }

  # Diversity measures (Scofield et al Am Nat http://www.jstor.org/stable/10.1086/668202)
  # alpha
  alpha.q <- 1/q.gg
  alpha.q.unweighted.mean <- (1/G) * sum(alpha.q)
  q.gg.unweighted.mean <- mean(q.gg)
  d.alpha.q <- 1/q.gg.unweighted.mean

  alpha.q.nielsen <- 1/q.nielsen.gg
  alpha.q.nielsen.unweighted.mean <- (1/G) * sum(alpha.q.nielsen)
  q.nielsen.gg.unweighted.mean <- mean(q.nielsen.gg)
  d.alpha.q.nielsen <- 1/q.nielsen.gg.unweighted.mean

  alpha.r <- 1/r.gg
  alpha.r.unweighted.mean <- (1/G) * sum(alpha.r)
  r.gg.unweighted.mean <- mean(r.gg)
  d.alpha.r <- 1/r.gg.unweighted.mean

  # gamma
  n.k <- apply(tab, 2, sum)
  Q.k <- n.k / N
  d.gamma.q <- 1/sum(Q.k * Q.k)
  d.gamma.q.nielsen <- 1/nielsenTransform(sum(Q.k * Q.k), N)
  d.gamma.r <- 1/sum((n.k*(n.k - 1)) / (N*(N-1)))

  # form pairwise omega matrix, replace its diagonal with alpha
  q.diversity.mat <- (2 * Q.mat) / outer(q.gg, q.gg, FUN="+")
  diag(q.diversity.mat) <- alpha.q
  q.divergence.mat <- 1 - (2 * Q.mat) / outer(q.gg, q.gg, FUN="+")
  diag(q.divergence.mat) <- 0
  q.overlap.mat <- 1 - q.divergence.mat

  q.nielsen.diversity.mat <- (2 * Q.mat) / outer(q.nielsen.gg, q.nielsen.gg, FUN="+")
  diag(q.nielsen.diversity.mat) <- alpha.q.nielsen
  q.nielsen.divergence.mat <- 1 - (2 * Q.mat) / outer(q.nielsen.gg, q.nielsen.gg, FUN="+")
  diag(q.nielsen.divergence.mat) <- 0
  q.nielsen.overlap.mat <- 1 - q.nielsen.divergence.mat

  r.diversity.mat <- (2 * Q.mat) / outer(r.gg, r.gg, FUN="+")
  diag(r.diversity.mat) <- alpha.r
  r.divergence.mat <- 1 - (2 * Q.mat) / outer(r.gg, r.gg, FUN="+")
  diag(r.divergence.mat) <- 0
  r.overlap.mat <- 1 - r.divergence.mat

  ## mean overlap and divergence
  C <- Q.mat
  diag(C) <- 0

  q.overlap <- sum(C)/((G - 1)*sum(q.gg))
  q.divergence <- 1 - q.overlap

  q.nielsen.overlap <- sum(C)/((G - 1)*sum(q.nielsen.gg))
  q.nielsen.divergence <- 1 - q.nielsen.overlap

  r.overlap <- sum(C)/((G - 1)*sum(r.gg))
  r.divergence <- 1 - r.overlap

  # return value
  list(produced.by=.thisis,
     table=tab,   # table passed in, as a matrix
     num.groups=G,  # number of rows (sites)
     num.sources=K, # number of columns (sources)
     num.samples=N, 
     num.samples.group=n.g,
     num.sources.group=apply(tab, 1, function(x)sum(x > 0)),
     num.samples.source=n.k,
     num.groups.source=apply(tab, 2, function(x)sum(x > 0)),

     # PMI measures by group/site
     q.gg=q.gg, 
     q.nielsen.gg=q.nielsen.gg,
     r.gg=r.gg, 

     # Matrix of q_gg and q_gh values
     Q.mat=Q.mat,  # this is q.gh

     q.0.gh=q.0.gh,

     # Pooled values (Scofield et al 2010 J Ecol)
     y.gh=y.gh,
     prop.y.0.gh=prop.y.0.gh,  # 

     # Weighted grand means
     q.bar.0=q.bar.0, 
     q.nielsen.bar.0=q.nielsen.bar.0, 
     r.bar.0=r.bar.0,

     # Unweighted means
     q.gg.unweighted.mean=q.gg.unweighted.mean,
     q.nielsen.gg.unweighted.mean=q.nielsen.gg.unweighted.mean,
     r.gg.unweighted.mean=r.gg.unweighted.mean,

     # Alpha diversities at each site
     alpha.q=alpha.q,
     alpha.q.nielsen=alpha.q.nielsen,
     alpha.r=alpha.r,

     # Mean alpha diversity as reciprocal of unweighted mean site PMI
     d.alpha.q=d.alpha.q,
     d.alpha.q.nielsen=d.alpha.q.nielsen,
     d.alpha.r=d.alpha.r,

     # Mean alpha diversity as mean of individual site alpha diversities
     alpha.q.unweighted.mean=alpha.q.unweighted.mean,
     alpha.q.nielsen.unweighted.mean=alpha.q.nielsen.unweighted.mean,
     alpha.r.unweighted.mean=alpha.r.unweighted.mean,

     # Gamma diversity
     d.gamma.q=d.gamma.q,
     d.gamma.q.nielsen=d.gamma.q.nielsen,
     d.gamma.r=d.gamma.r,

     # Beta diversity
     d.beta.q=d.gamma.q/d.alpha.q,
     d.beta.q.nielsen=d.gamma.q.nielsen/d.alpha.q.nielsen,
     d.beta.r=d.gamma.r/d.alpha.r,

     # Diversity matrix: q_gg on diagonal, q_gh off-diagonal
     # Divergence matrix: 0 diagonal, delta_gh off-diagonal
     # You can construct the overlap matrix (1 diagonal,
     #   omega_gh off-diagonal) via 1 - Divergence matrix
     q.diversity.mat=q.diversity.mat,
     q.divergence.mat=q.divergence.mat,
     # Mean divergence and overlap
     q.divergence=q.divergence,
     q.overlap=q.overlap,

     q.nielsen.diversity.mat=q.nielsen.diversity.mat,
     q.nielsen.divergence.mat=q.nielsen.divergence.mat,
     q.nielsen.divergence=q.nielsen.divergence,
     q.nielsen.overlap=q.nielsen.overlap,

     r.diversity.mat=r.diversity.mat,
     r.divergence.mat=r.divergence.mat,
     r.divergence=r.divergence,
     r.overlap=r.overlap
    )
}


