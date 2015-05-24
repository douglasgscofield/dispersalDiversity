#' @include diversity-allele_divtables.R
#' @include gammaAccum-divtable.R
# for collation
NULL



#' @rdname gammaAccum
#'
#' @export
#'
gammaAccum.allele_divtables <- function(adt, 
    resample.method = c("permute", "bootstrap"), n.resample = 1000,
    gamma.method = c("q", "r", "q.nielsen"), ...)
{
    accum.method <- match.arg(accum.method)
    resample.method <- match.arg(resample.method)
    gamma.method <- match.arg(gamma.method)
    d <- diversityMultilocus(adt)
    ans <- list()
    ans$obs.gamma <- d$gamma
    ans$obs.omega.mean <- d$overlap
    ans$obs.delta.mean <- d$divergence
    ans$simple.results <- 
        gammaAccumSimple.allele_divtables(adt,
                         resample.method = resample.method,
                         n.resample = n.resample,
                         gamma.method = gamma.method,
                         ...)
    structure(ans, class = c('gamma_accum', 'list'))
}



gammaAccumSimple.allele_divtables <- function(adt, ...)
{
    gammaAccumStats(gammaAccumWorker.allele_divtables(adt, ...))
}


gammaAccumWorker.allele_divtables <- function(adt, n.sites = dim(adt[[1]])[1],
    resample.method = c("permute", "bootstrap"), n.resample = 1000,
    gamma.method = c("q", "r", "q.nielsen"), ...)
{
    accum.method <- match.arg(accum.method)
    resample.method <- match.arg(resample.method)
    gamma.method <- match.arg(gamma.method)
    pool.names <- dimnames(adt[[1]])[[1]]
    L <- names(adt)
    G <- dim(adt[[1]])[1]
    ans <- lapply(1:n.sites, function(x) numeric(0))
    for (i in 1:n.resample) {
        row.order <- sample(1:G,
                            size = n.sites,
                            replace = ifelse(resample.method == "bootstrap", 
                                             TRUE, FALSE))
        for (g in 1:n.sites) {
            gamma.all.loci <- sapply(adt, function(x) {
                .calculateGammaAccum(apply(x[row.order[1:g], , drop=FALSE], 
                                           2, sum),
                                     gamma.method) })
            ans[[g]] <- c(ans[[g]], mean(gamma.all.loci))
        }
    }
    ans
}

