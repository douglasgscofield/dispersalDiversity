#' @include diversity-divtable.R
NULL



#' Calculate diversity for \code{allele_divtables} object using \code{diversityMultilocus}
#'
#' @param x Object of class \code{\link{allele_divtables}}
#'
#' @param \dots Additional arguments passed to
#' \code{\link{diversityMultilocus}}
#'
#' @return The value from \code{\link{diversityMultilocus}}
#'
#' @seealso \code{\link{diversityMultilocus}}
#'
#' @rdname diversity
#'
#' @export
#'
diversity.allele_divtables <- function(x, ...)
{
    diversityMultilocus.allele_divtables(x, ...)
}



#' Calculate diversity across a set of allele tables
#'
#' Calculate diversity for a set of loci encoded in an object of class
#' \code{'allele_tables'} as created with \code{\link{createAlleleTables}},
#' and diversity for each locus is calculated separately
#' using \code{\link{diversitySingleLocus}}.
#'
#' @param x Allele diversity data contained in an object of class
#' \code{\link{allele_divtables}}
#'
#' @param ploidy Ploidy of underlying genotypes
#'
#' @param method Method to use for calculating diversity statistics.
#' See \code{\link{diversity}} and Scofield \emph{et al}. (2012),
#' Sort \emph{et al}. (in press).
#'
#' @return A list of class \code{multilocus_diversity}.  Elements of the
#' list are calculated from values returned by
#' \code{\link{diversitySingleLocus}} and include:
#' \itemize{
#'     \item \code{q.gg}, matrix of \code{allele.alpha.g} values across sites
#'           (rows), for each locus (columns)
#'     \item \code{alpha.gk}, matrix of \code{allele.alpha.g} values across
#'           sites (rows), for each locus (columns) (\emph{are these simply
#'           reciprocals of q.gg?})
#'     \item \code{q.k}, vector of mean \code{q.gg} values across loci, the
#'           mean values of the rows of \code{q.gg}
#'     \item \code{alpha.k}, vector of mean \code{allele.alpha.g} values
#'           across loci, the mean values of the rows of \code{alpha.gk}
#'     \item \code{q.bar.weighted.a}, vector of \code{allele.q.bar.weighted}
#'           for each locus
#'     \item \emph{\code{alpha.bar.weighted.a}, vector of
#'           \code{allele.alpha.bar.weighted} for each locus ????}
#'     \item \emph{\code{overlap.a}, vector of \code{allele.overlap} for each locus ????}
#'     \item \emph{\code{divergence.a}, vector of \code{allele.divergence} for each locus ????}
#'     \item \code{Q.0.a}, vector of \code{allele.Q.0} values for each locus
#'     \item \code{q.bar.weighted}, mean of \code{q.bar.weighted.a}
#'     \item \code{alpha.bar.weighted}, reciprocal of \code{q.bar.weighted}
#'     \item \emph{\code{Q.0}, mean of \code{Q.0.a} ????}
#'     \item \code{gamma}, reciprocal of mean of \code{Q.0.a}, \emph{Q.0 ????}
#'     \item \code{beta}, \code{gamma} divided by \code{alpha.bar.weighted}
#'     \item \code{G}, number of loci
#'     \item \code{N}, allele count for in each locus
#'     \item \code{K}, number of sites
#'     \item \code{alpha.scaled}, \code{alpha.bar.weighted} scaled to fall
#'           in the interval [0, 1]
#'     \item \code{gamma.scaled}, \code{gamma} scaled to fall in the interval
#'           [0, 1]
#'     \item \code{beta.scaled}, \code{beta} scaled to fall in the interval
#'           [0, 1]
#'     \item \code{q.gh}, site-by-site matrix of mean of \code{q.gh} values
#'           across all loci
#'     \item \code{overlap}, sum of \code{q.gh} divided by \code{K} - 1 times
#'           the sum of \code{q.k}
#'     \item \code{divergence}, 1 - \code{overlap}
#'     \item \code{overlap.mean}, mean of \code{allele.divergence} values
#'           for each locus (\emph{\code{overlap.a} ???})
#'     \item \code{divergence.mean}, mean of \code{allele.overlap} values
#'           for each locus (\emph{\code{divergence.a} ???})
#' }
#'
#' These items were left out for no apparent reason or are confusing:
#' \itemize{
#'     \item \code{Q.0}
#'     \item \code{alpha.bar.weighted.a}
#'     \item \code{G}
#'     \item \code{N}
#'     \item \code{K}
#'     \item \code{overlap.a}
#'     \item \code{divergence.a}
#' }
#'
#'   ans <- list(q.gg               = q.gg,
#'               alpha.gk           = alpha.gk,
#'               q.k                = q.k,
#'               alpha.k            = alpha.k,
#'               q.bar.weighted.a   = q.bar.weighted.a,
#'               Q.0.a              = Q.0.a,
#'               q.bar.weighted     = q.bar.weighted,
#'               alpha.bar.weighted = alpha.bar.weighted,
#'               gamma              = gamma,
#'               beta               = beta,
#'               G                  = G,
#'               N                  = N,
#'               K                  = K,
#'               alpha.scaled       = alpha.scaled,
#'               gamma.scaled       = gamma.scaled,
#'               beta.scaled        = beta.scaled,
#'               q.gh               = q.gh,
#'               overlap            = overlap,
#'               divergence         = divergence,
#'               overlap.mean       = overlap.mean,
#'               divergence.mean    = divergence.mean)
#'
#' @references
#'
#' Scofield, D. G., Smouse, P. E., Karubian, J. and Sork, V. L. (2012)
#' Use of alpha, beta and gamma diversity measures to characterize seed
#' dispersal by animals.  \emph{American Naturalist} 180:719-732.
#'
#' Sork, V. L., Smouse, P. E., Grivet, D. and Scofield, D. G. (In press)
#' Impact of asymmetric male and female gamete dispersal on allelic
#' diversity and spatial genetic structure in valley oak
#' (\emph{Quercus lobata} N\'{e}e).  \emph{Evolutionary Ecology}.
#'
#' @examples
#'
#' ## Workflow to calculate basic allelic diversity statistics:
#' library(readGenalex)
#' data(Qagr_pericarp_genotypes)
#' gt <- createAlleleTables(Qagr_pericarp_genotypes)
#' div <- diversityMultilocus(gt)
#'
#' @export
#'
#' @name diversityMultilocus
#'
NULL

diversityMultilocus <- function(x, ...) UseMethod("diversityMultilocus")



diversityMultilocus.default <- function(x, ...)
{
    stop("intended to operate on an object of class 'allele_divtables', ",
         "perhaps you need to use createAlleleTables() first?")
}



# TODO: now with method=, the returned values are not correctly named
#
#' @rdname diversityMultilocus
#'
#' @export
#'
diversityMultilocus.allele_divtables <- function(x, ploidy = 2,
    method = c("r", "q.nielsen", "q"), ...)
{
    # calculates diversity statistics for a collection of loci, the argument
    # is produced by createAlleleTables()

    method <- match.arg(method)

    # create utility zero vectors and matrices
    pop.locus.0 <- matrix(0, nrow = nrow(lst[[1]]), ncol = length(names(lst)),
                          dimnames = list(Pop = rownames(lst[[1]]),
                                          Locus = names(lst)))
    pop.pop.0 <- matrix(0, nrow = nrow(lst[[1]]), ncol = nrow(lst[[1]]),
                        dimnames = list(Pop = rownames(lst[[1]]),
                                        Pop = rownames(lst[[1]])))
    pop.0 <- pop.locus.0[, 1]
    locus.0 <- pop.locus.0[1, ]
    # accumulator vectors and matrices
    alpha.gk <- q.gg <- pop.locus.0
    q.gh <- pop.pop.0
    overlap.a <- divergence.a <- Q.0.a <- q.bar.weighted.a <-
        alpha.bar.weighted.a <- locus.0
    q.k <- alpha.k <- pop.0
    p <- list()
    # Go through loci, storing diversity values
    for (a in names(lst)) {
        p[[a]] <- diversitySingleLocus(lst[[a]], method)
        alpha.gk[, a] <- p[[a]]$allele.alpha.g
        q.gg[, a] <- p[[a]]$allele.q.gg
        Q.0.a[a] <- p[[a]]$allele.Q.0
        q.bar.weighted.a[a] <- p[[a]]$allele.q.bar.weighted
        alpha.bar.weighted.a[a] <- p[[a]]$allele.alpha.bar.weighted
        q.gh <- q.gh + p[[a]]$allele.q.gh
        overlap.a[a] <- p[[a]]$allele.overlap
        divergence.a[a] <- p[[a]]$allele.divergence
    }
    # summary variables
    G <- length(locus.0)
    N <- p[[1]]$allele.N  # number of alleles, so for diploid this will already be 2*N
    K <- p[[1]]$allele.G
    # derived variables
    alpha.k[] <- apply(alpha.gk, 1, mean)
    q.k[] <- apply(q.gg, 1, mean)
    q.bar.weighted <- mean(q.bar.weighted.a)
    alpha.bar.weighted <- 1 / q.bar.weighted
    Q.0 <- mean(Q.0.a)
    gamma <- 1 / Q.0
    beta <- gamma / alpha.bar.weighted
    # derived overlaps
    # mean of r.gh elements across loci, note matrix is symmetrical
    # containing both r_gh and r_hg so no '2*' required
    q.gh <- q.gh / G
    overlap <- sum(q.gh) / ((K - 1) * sum(q.k))
    divergence <- 1 - overlap
    overlap.mean <- mean(overlap.a)
    divergence.mean <- mean(divergence.a)
    # scaled values, note that N here is already 2*N when diploid
    alpha.max <- ((ploidy * N) - G) / G
    gamma.max <- ((ploidy * N) - G)
    beta.max <- G
    alpha.scaled <- ((alpha.bar.weighted - 1) / alpha.bar.weighted) *
                    (alpha.max / (alpha.max - 1))
    gamma.scaled <- ((gamma - 1) / gamma) * (gamma.max / (gamma.max - 1))
    beta.scaled <- ((beta - 1) / beta) * (beta.max / (beta.max - 1))
    # return value
    ans <- list(q.gg               = q.gg,
                alpha.gk           = alpha.gk,
                q.k                = q.k,
                alpha.k            = alpha.k,
                q.bar.weighted.a   = q.bar.weighted.a,
                Q.0.a              = Q.0.a,
                q.bar.weighted     = q.bar.weighted,
                alpha.bar.weighted = alpha.bar.weighted,
                gamma              = gamma,
                beta               = beta,
                G                  = G,
                N                  = N,
                K                  = K,
                alpha.scaled       = alpha.scaled,
                gamma.scaled       = gamma.scaled,
                beta.scaled        = beta.scaled,
                q.gh               = q.gh,
                overlap            = overlap,
                divergence         = divergence,
                overlap.mean       = overlap.mean,
                divergence.mean    = divergence.mean)
    structure(ans, class = c('multilocus_diversity', 'list'))
}


#' Calculate allelic diversity for a single locus
#'
#' The single argument is, for a single locus, a table of site-by-allele
#' counts, with row names being the site names, and column names being the
#' names given to the individual alleles.
#'
#' For specific formulas and discussion of each statistic, see
#' Scofield \emph{et al}. (2015).
#'
#' @param tab    Table of site-by-allele counts, with row names being the site
#' names and column names being the names given to the individual alleles
#'
#' @param method   Method to use when calculating diversity within the locus.
#' See \code{\link{diversity}}.
#'
#' @return A list of class \code{singlelocus_diversity}, containing values
#' returned by \code{\link{diversity}} or derived from them:
#' \itemize{
#'     \item \code{allele.table}, class \code{\link{divtable}} table of sites
#'           by alleles
#'     \item \code{allele.N}, the total allele count
#'     \item \code{allele.G}, the number of unique alleles
#'     \item \code{allele.q.gg}, vector of length number of sites, containing
#'           squared frequencies of alleles in each site
#'     \item \code{allele.n.g}, vector of the number of alleles in each site
#'     \item \code{allele.n.k}, vector of allele counts for each unique allele
#'     \item \code{allele.alpha.g}, \code{alpha.g} values for each site (EXPAND from diversity)
#'     \item \code{allele.q.bar.weighted}, weighted mean allele frequency
#'           calculated across sites; must this always be an "r" type mean?
#'     \item \code{allele.alpha.bar.weighted}, reciprocal of
#'           \code{allele.q.bar.weighted}
#'     \item \code{allele.Q.0}, reciprocal of \code{allele.gamma}
#'     \item \code{allele.q.gh}, site-by-site matrix of \code{q.gh} values
#'           across all loci
#'     \item \code{allele.overlap}, mean site-by-site overlap as calculated by
#'           \code{diversity}
#'     \item \code{allele.divergence}, mean site-by-site divergence as
#'           calculated by \code{diversity}
#'     \item \code{allele.gamma}, gamma diversity across sites as calculated
#'           by \code{diversity}
#'     \item \code{allele.beta}, \code{allele.gamma} divided by
#'           \code{allele.alpha.bar.weighted}
#'     \item \code{allele.func.scale.alpha}, function for scaling
#'           \code{allele.alpha.g} and \code{allele.alpha.bar.weighted} with
#'           given \code{allele.N} and \code{allele.G} and arbitrary ploidy
#'     \item \code{allele.func.scale.gamma}, function for scaling
#'           \code{allele.gamma} with
#'           given \code{allele.N} and \code{allele.G} and arbitrary ploidy
#'     \item \code{allele.func.scale.beta}, function for scaling
#'           \code{allele.beta} with given \code{allele.G} and arbitrary ploidy
#'     \item \code{allele.scaled.ploidy1}, a list containing values for
#'           \code{allele.alpha.g}, \code{allele.alpha.bar.weighted},
#'           \code{allele.gamma} and \code{allele.beta} scaled using the above
#'           functions and ploidy of 1
#'     \item \code{allele.scaled.ploidy2}, a list containing values for
#'           \code{allele.alpha.g}, \code{allele.alpha.bar.weighted},
#'           \code{allele.gamma} and \code{allele.beta} scaled using the above
#'           functions and ploidy of 2
#' }
#'
#' \emph{Do we want to generalise allele.q.bar.weighted to compute whatever
#' form of mean is appropriate for \code{method}???}
#'
#' @export
#'
#' @name diversitySingleLocus
#'
NULL

diversitySingleLocus <- function(tab, ...) UseMethod("diversitySingleLocus")



diversitySingleLocus.default <- function(tab, ...)
{
    stop("intended to operate on an object of class allele_divtable describing a single locus")
}



#' @rdname diversitySingleLocus
#'
#' @export
#'
diversitySingleLocus.divtable <- function(tab,
    method = c("r", "q.nielsen", "q"), ...)
{
    method <- match.arg(method)

    d <- diversity(tab)

    allele.N <- d$num.samples
    allele.G <- d$num.groups
    allele.n.g <- d$num.samples.group
    allele.n.k <- d$num.samples.source
    # recreate C matrix not returned from diversity()
    allele.q.gh <- d$Q.mat; diag(allele.q.gh) <- 0 # C matrix here

    p <- d[[method]]

    allele.q.gg <- p$q.gg
    allele.alpha.g <- p$alpha.g

    # the allelic version in eq.2 is a bit different from the ecological version
    #allele.q.bar.weighted <- p$q.bar.0
    allele.q.bar.weighted <- sum((allele.n.g - 1) * allele.q.gg) /
                             (allele.N - allele.G)
    allele.alpha.bar.weighted <- 1 / allele.q.bar.weighted
    #allele.R.0 <- sum(allele.n.k * (allele.n.k - 1)) /
    #              (allele.N * (allele.N - 1))
    allele.gamma <- p$d.gamma
    allele.Q.0 <- 1 / allele.gamma
    allele.beta <- allele.gamma / allele.alpha.bar.weighted

    allele.overlap <- p$overlap
    allele.divergence <- p$divergence

    # functions to scale alpha, beta, gamma, see Sork et al
    # the use of body() replaces N and G in the parse tree with current values
    allele.func.scale.alpha <- function(a, ploidy = 2) {
        N <- allele.N
        G <- allele.G
        alpha.max <- ((ploidy * N) - G) / G
        ((a - 1) / a) * (alpha.max / (alpha.max - 1))
    }
    body(allele.func.scale.alpha)[[2]][[3]] <- allele.N
    body(allele.func.scale.alpha)[[3]][[3]] <- allele.G
    allele.func.scale.gamma <- function(g, ploidy=2) {
        N <- allele.N
        G <- allele.G
        gamma.max <- (ploidy * N) - G
        ((g - 1) / g) * (gamma.max / (gamma.max - 1))
    }
    body(allele.func.scale.gamma)[[2]][[3]] <- allele.N
    body(allele.func.scale.gamma)[[3]][[3]] <- allele.G
    allele.func.scale.beta <- function(b, ploidy = 2) {
        G <- allele.G
        beta.max <- G
        ((b - 1) / b) * (beta.max / (beta.max - 1))
    }
    body(allele.func.scale.beta)[[2]][[3]] <- allele.G
    # scaled to [0, 1] assuming ploidy of 1 and 2
    allele.scaled.ploidy1 <-
        list(allele.alpha.g            = allele.func.scale.alpha(allele.alpha.g, 1),
             allele.alpha.bar.weighted = allele.func.scale.alpha(allele.alpha.bar.weighted, 1),
             allele.gamma              = allele.func.scale.gamma(allele.gamma, 1),
             allele.beta               = allele.func.scale.beta(allele.beta, 1))
    allele.scaled.ploidy2 <-
        list(allele.alpha.g            = allele.func.scale.alpha(allele.alpha.g, 2),
             allele.alpha.bar.weighted = allele.func.scale.alpha(allele.alpha.bar.weighted, 2),
             allele.gamma              = allele.func.scale.gamma(allele.gamma, 2),
             allele.beta               = allele.func.scale.beta(allele.beta, 2))

    ans <- list(allele.table              = tab,
                allele.N                  = allele.N,
                allele.G                  = allele.G,
                allele.q.gg               = allele.q.gg,
                allele.n.g                = allele.n.g,
                allele.n.k                = allele.n.k,
                allele.alpha.g            = allele.alpha.g,
                allele.q.bar.weighted     = allele.q.bar.weighted,
                allele.alpha.bar.weighted = allele.alpha.bar.weighted,
                allele.Q.0                = allele.Q.0,
                allele.q.gh               = allele.q.gh,
                allele.overlap            = allele.overlap,
                allele.divergence         = allele.divergence,
                allele.gamma              = allele.gamma,
                allele.beta               = allele.beta,
                allele.func.scale.alpha   = allele.func.scale.alpha,
                allele.func.scale.gamma   = allele.func.scale.gamma,
                allele.func.scale.beta    = allele.func.scale.beta,
                allele.scaled.ploidy1     = allele.scaled.ploidy1,
                allele.scaled.ploidy2     = allele.scaled.ploidy2)
    structure(ans, class = c('singlelocus_diversity', 'list'))
}

