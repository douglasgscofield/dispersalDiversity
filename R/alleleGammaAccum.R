#' @include allelePmiDiversity.R
#' @include gammaAccum.R
# for collation
NULL



#' Perform gamma diversity accumulation on allele counts
#'
#' Several arguments control the method of accumulation and value of gamma
#' calculated.  Only the defaults have been tested; the others were developed
#' while exploring the data and must be considered experimental.  The result is
#' returned in a list, which may be passed to \code{\link{plotGammaAccum}} to
#' plot the result.
#'
#' @param lst List of site-by-source tables for loci, as produced by 
#' \code{allele.createTableList}
#'
#' @param accum.method \code{"random"} (default) or \code{"proximity"}.  If 
#' \code{"proximity"} is used, then \code{distance.file} must be supplied
#' (see below).
#'
#' @param resample.method \code{"permute"} (default) or \code{"bootstrap"}, 
#' whether to resample sites without (\code{"permute}") or with 
#' (\code{"bootstrap"}) replacement
#'
#' @param n.resample Number of resamples for accumulation curve
#'
#' @param distance.file A file or data frame containing three columns of data,
#' with the header/column names being \code{pool}, \code{X}, and \code{Y},
#' containing the spatial locations of the seed pools named in the row names 
#' of tab.  Only used if the \code{accum.method == "proximity"}.
#'
#' @param gamma.method Calculate gamma using the \code{"r"} (default), 
#' \code{"q"} or \code{"q.nielsen"} method, see Scofield et al. (2012).
#' The allele.* functions are all based on \code{"r"} versions and undefined 
#' results may happen otherwise.
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
#' @examples
#'
#' dat <- readGenalex("genotypes.txt")
#' lst <- allele.createTableList(dat)
#' allele.rga.result <- allele.runGammaAccum(lst)
#' plotGammaAccum(allele.rga.result)
#'
#' @export allele.runGammaAccum
#'
allele.runGammaAccum <- function(lst, 
    accum.method = c("random", "proximity"),
    resample.method = c("permute", "bootstrap"), n.resample = 1000,
    gamma.method = c("r", "q.nielsen", "q"), distance.file=NULL,
    ...) {
    accum.method <- match.arg(accum.method)
    resample.method <- match.arg(resample.method)
    gamma.method <- match.arg(gamma.method)
    if (gamma.method != "r")
        warning("allele.runGammaAccum() should be used with gamma.method == 'r'", 
                immediate. = TRUE)
    pmiD <- allele.diversity(lst)
    ans <- list()
    ans$obs.gamma <- pmiD$gamma
    ans$obs.omega.mean <- pmiD$overlap
    ans$obs.delta.mean <- pmiD$divergence
    ans$simple.results <- allele.runGammaAccumSimple(lst, 
                            accum.method = accum.method,
                            resample.method = resample.method,
                            n.resample = n.resample,
                            gamma.method = gamma.method,
                            distance.file = distance.file,
                            ...)
    return(ans)
}



# Wrapper that runs and then return stats from allele.gammaAccum()
#
allele.runGammaAccumSimple <- function(lst, ...) {
    return(gammaAccumStats(allele.gammaAccum(lst, ...)))
}



# Actually perform the gamma accumulation
#
allele.gammaAccum <- function(lst, n.sites = dim(lst[[1]])[1],
    accum.method = c("random", "proximity"),
    resample.method = c("permute", "bootstrap"), n.resample = 1000,
    distance.file = NULL, gamma.method = c("r", "q.nielsen", "q")) {

    # If used, the distance.file has three columns, with names: pool, X, Y
    # It is either a filename to read, or a data.frame
    accum.method <- match.arg(accum.method)
    resample.method <- match.arg(resample.method)
    gamma.method <- match.arg(gamma.method)
    pool.names <- dimnames(lst[[1]])[[1]]
    L <- names(lst)
    G <- dim(lst[[1]])[1]
    if (accum.method == "proximity")
        gxy <- .loadDistanceFile(distance.file, pool.names)

    ans <- lapply(1:n.sites, function(x) numeric(0))
    for (i in 1:n.resample) {
        row.order <- sample(1:G,
                            size = n.sites,
                            replace = ifelse(resample.method == "bootstrap", 
                                             TRUE, FALSE))
        if (accum.method == "proximity")
            row.order <- .reorderByProximity(row.order, gxy, pool.names)
        for (g in 1:n.sites) {
            gamma.all.loci <- sapply(lst, function(x) {
                .calculateGammaAccum(apply(x[row.order[1:g], , drop=FALSE], 
                                           2, sum),
                                     gamma.method) })
            ans[[g]] <- c(ans[[g]], mean(gamma.all.loci))
        }
    }
    return(ans)
}

