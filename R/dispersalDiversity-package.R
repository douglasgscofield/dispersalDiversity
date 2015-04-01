#' Analyse individual and genetic diversity of dispersal processes in biological populations
#'
#' A collection of R functions to analyse individual and genetic diversity in
#' biological populations, primarily in the context of studying dispersal
#' processes.  The basic data format is a types x sites matrix.  Functions are
#' provided for calculating alpha, beta, delta/omega and gamma diversity using a
#' variety of measures; calculating alpha, beta, delta/omega and gamma diversity
#' of allelic variation; and comparing diversity measures between samples.
#'
#' A similar set of functions provides analyses of two similar dispersal data
#' sets.  The first is ecological sites-by-counts data, and in this context the
#' alpha, beta and gamma quantities produced should be similar.  The second is
#' genetic data such as allele identities, which is analysed to partition
#' diversity into a structure similar to alpha, beta etc.  See Sork et al.
#'
#' The function \code{createAlleleTables} accepts an argument of class
#' \code{'genalex'} and creates an object of class \code{'allele_tables'}, a
#' list of site-by-allele count tables, one per locus.  Further analyses may be
#' applied to objects of class \code{'allele_tables'}, including structural
#' tests and comparisons.
#'
#' These statistical tools were developed in collaboration with Peter E.
#' Smouse (Rutgers University) and Victoria L. Sork (UCLA) and were funded by
#' U.S. National Science Foundation awards NSF-DEB-0514956 and
#' NSF-DEB-0516529.
#'
#' @references
#'
#' Grivet, D., Smouse, P. E. and Sork, V. L. (2005) A novel approach to an old
#' problem: tracking dispersed seeds.  \emph{Molecular Ecology} 14:3585-3595.
#'
#' Nielsen, R., Tarpy, D. R. and Reeve, H. K. (2003) Estimating effective
#' paternity number in social insects and the effective number of alleles in
#' a population.  \emph{Molecular Ecology} 12:3157-3164.
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
#' Sork, V. L., Smouse, P. E., Grivet, D. and Scofield, D. G.  Impact of
#' asymmetric male and female gamete dispersal on allelic diversity and
#' spatial genetic structure in valley oak (\emph{Quercus lobata} N\'{e}e).
#' \emph{In revision}.
#'
#' \url{https://github.com/douglasgscofield/dispersalDiversity}
#'
#' @docType package
#'
#' @name dispersalDiversity-package
#'
NULL

