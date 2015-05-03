#' @include diversity-allele_divtables.R
#' @include gammaAccum-divtable.R
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
#' @param resample.method \code{"permute"} (default) or \code{"bootstrap"}, 
#' whether to resample sites without (\code{"permute"}) or with 
#' (\code{"bootstrap"}) replacement
#'
#' @param n.resample Number of resamples for accumulation curve
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
#' lst <- createAlleleTables(dat)
#' allele.rga.result <- gammaAccum(lst)
#' plot(allele.rga.result)
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
gammaAccum.allele_divtables <- function(lst, 
    resample.method = c("permute", "bootstrap"), n.resample = 1000,
    gamma.method = c("r", "q.nielsen", "q"), ...)
{
    accum.method <- match.arg(accum.method)
    resample.method <- match.arg(resample.method)
    gamma.method <- match.arg(gamma.method)
    d <- diversityMultilocus(lst)
    ans <- list()
    ans$obs.gamma <- d$gamma
    ans$obs.omega.mean <- d$overlap
    ans$obs.delta.mean <- d$divergence
    ans$simple.results <- 
        gammaAccumSimple(lst,
                         resample.method = resample.method,
                         n.resample = n.resample,
                         gamma.method = gamma.method,
                         ...)
    structure(ans, class = c('gamma_accum', 'list'))
}



gammaAccumSimple <- function(x, ...) UseMethod("gammaAccumSimple")

# Wrapper that runs and then return stats from allele.gammaAccum()
#
gammaAccumSimple.allele_divtables <- function(lst, ...)
{
    gammaAccumStats(gammaAccumWorker.allele_divtables(lst, ...))
}


gammaAccumWorker <- function(x, ...) UseMethod("gammaAccumWorker")

# Actually perform the gamma accumulation
#
gammaAccumWorker.allele_divtables <- function(lst, n.sites = dim(lst[[1]])[1],
    resample.method = c("permute", "bootstrap"), n.resample = 1000,
    gamma.method = c("r", "q.nielsen", "q"), ...)
{
    # If used, the distance.file has three columns, with names: pool, X, Y
    # It is either a filename to read, or a data.frame
    accum.method <- match.arg(accum.method)
    resample.method <- match.arg(resample.method)
    gamma.method <- match.arg(gamma.method)
    pool.names <- dimnames(lst[[1]])[[1]]
    L <- names(lst)
    G <- dim(lst[[1]])[1]
    ans <- lapply(1:n.sites, function(x) numeric(0))
    for (i in 1:n.resample) {
        row.order <- sample(1:G,
                            size = n.sites,
                            replace = ifelse(resample.method == "bootstrap", 
                                             TRUE, FALSE))
        for (g in 1:n.sites) {
            gamma.all.loci <- sapply(lst, function(x) {
                .calculateGammaAccum(apply(x[row.order[1:g], , drop=FALSE], 
                                           2, sum),
                                     gamma.method) })
            ans[[g]] <- c(ans[[g]], mean(gamma.all.loci))
        }
    }
    ans
}

