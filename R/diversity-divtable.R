#' @include nielsenTransform.R
# for collation
NULL

#' Calculate a collection of diversity statistics from a site x source matrix
#'
#' Given a site-by-source matrix, calculate a collection of diversity
#' statistics.  Statistics include PMI (probability of maternal inheritance)
#' from Grivet \emph{et al}. (2005) and Scofield \emph{et al}. (2010, 2011) 
#' as well as alpha, beta, gamma and delta/omega from classical ecology and
#' applied to a dispersal context by Scofield \emph{et al}. (2012).
#'
#' @param tab     \code{\link{divtable}} of sites-by-sources counts.  An
#' attempt is make to convert it to class \code{\link{divtable}} if it is not
#' already.  Rows are not reordered.
#'
#' @param \dots   Additional arguments, currently ignored
#'
#' @return A list of class \code{diversity} containing three forms of
#' diversity statistics, each returned as lists under separate named
#' elements:
#' \itemize{
#'     \item \code{q}, based on squared frequencies, known to be biased at 
#'           smaller sample sizes
#'     \item \code{r}, unbiased versions of \code{q} (Simpson, 1949), but
#'           singletons are dropped and values are nonsensical if a site
#'           contains only singletons
#'     \item \code{q.nielsen}, a better unbiased transformation of \code{q}
#'           following Nielsen \emph{et al}. (2003) that does not have the
#'           singleton properties of \code{r}
#' }
#' Each such list \code{q}, \code{r} and \code{q.nielsen} contains several 
#' elements representing different diversity statistics.  Several 
#' quantities are included:
#' \itemize{
#'     \item \code{q.gg}, vector of length number of sites, containing
#'           squared frequencies of groups in each site.  In form \code{q},
#'           these are untransformed, while in form \code{r} and
#'           \code{q.nielsen}, these are transformed accordingly
#'     \item \code{q.bar.0}, scalar, weighted mean of \code{q.gg} 
#'     \item \code{q.variance.0}, scalar, weighted sample variance of
#'           \code{q.bar.0} 
#'     \item \code{q.unweighted.mean}, mean of \code{q.gg}
#'     \item \code{alpha.g}, vector of length number of sites, containing
#'           reciprocals of \code{q.gg}
#'     \item \code{alpha.unweighted.mean}, scalar, mean of \code{alpha.g}
#'     \item \code{d.alpha}, scalar, reciprocal of \code{q.unweighted.mean}
#'     \item \code{d.gamma}, scalar, reciprocal of sum of grand squared
#'           frequencies
#'     \item \code{d.beta}, scalar, \code{d.gamma / d.alpha}
#'     \item \code{diversity.mat}, sites-by-sites matrix containing Chao
#'           pairwise divergences off-diagonal, \code{alpha.g} on-diagonal
#'     \item \code{divergence.mat}, sites-by-sites matrix matching
#'           \code{diversity.mat} but with with \code{0} on-diagonal
#'     \item \code{overlap.mat}, sites-by-sites matrix containing
#'           \code{1 - divergence.mat}
#'     \item \code{overlap}, scalar, mean overlap,
#'           \eqn{\sum{C} / ((G - 1) \sum{q_{gg}})}
#'     \item \code{divergence}, scalar, mean divergence, \code{1 - overlap}
#' }
#' The list also contains the original \code{\link{divtable}} passed to
#' \code{diversity}, together with several descriptive quantities:
#' \itemize{
#'     \item \code{table}, the class \code{\link{divtable}} object passed in
#'           as \code{tab}
#'     \item \code{num.groups}, the number of rows in \code{tab}
#'     \item \code{num.sources}, the number of columns in \code{tab}
#'     \item \code{num.samples}, the sum of all elements in \code{tab}
#'     \item \code{num.samples.group}, a vector of length \code{num.groups},
#'           the sum of elements in each row of \code{tab}
#'     \item \code{num.sources.group}, a vector of length \code{num.groups},
#'           the number of distinct groups observed in each row of \code{tab}
#'     \item \code{num.samples.source}, a vector of length \code{num.sources},
#'           the sum of elements in each column of \code{tab}
#'     \item \code{num.groups.source}, a vector of length \code{num.groups},
#'           the number of distinct groups observed in each row of \code{tab}
#' }
#' Finally, some additional statistics representing pooled PMI measures
#' (Scofield \emph{et al}. 2010) and weighted overlaps (Scofield \emph{et al}.
#' 2011) are returned:
#' \itemize{
#'     \item \code{y.gh}, a pairwise matrix of counts of items in each group
#'           from any source that are shared with each other group
#'     \item \code{prop.y.0.gh}, a matrix of proportions of items in each
#'           group from any source that are shared with each other group,
#'           derived from \code{y.gh}
#'     \item \code{q.0.gh}, a matrix of probabilities that two items drawn
#'           from group g and h are from any of the sources shared between
#'           g and h
#'     \item \code{q.bar.gh}, weighted mean overlap (\eqn{q_{gh}}), 
#'           probabilities that two items, one drawn from each of groups
#'           \emph{g} and \emph{h}, are from the same source, weighted by
#'           sample sizes
#'     \item \code{q.variance.gh}, weighted sample variance of \code{q.bar.gh}
#' }
#'
#' @references
#'
#' Simpson, E. D. (1949) Measurement of diversity. \emph{Nature} 163:688.
#'
#' Chao, A., Jost, L., Chiang, S. C., Jiang, Y. H., and Chazdon, R. L. (2008)
#' A two-stage probabilistic approach to multiple-community similarity 
#' indices.  \emph{Biometrics} 64:1178-1186.
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
#' @examples
#'
#' ## create table of random memberships from data frame of site-source pairs
#' n.sites <- 5
#' n.sources <- 10
#' n.samples <- 160
#' set.seed(75333)
#' t <- data.frame(site = sample(n.sites, n.samples, replace = TRUE),
#'                 source = round(runif(n.samples) * n.sources + 0.5))
#' ## diversity accepts table via default method, which applies as.divtable
#' div <- diversity(table(t))
#' div$q
#' 
#' @export
#'
#' @name diversity
#'
NULL

diversity <- function(tab, ...) UseMethod("diversity")



#' @rdname diversity
#'
#' @export
#'
diversity.default <- function(tab, ...)
{
    diversity.divtable(as.divtable(tab), ...)
}



#' @rdname diversity
#'
#' @export
#'
diversity.divtable <- function(tab, ...)
{
    orig.tab <- tab
    if (storage.mode(tab) == "integer") storage.mode(tab) <- "double"
    G <- dim(tab)[1]
    K <- dim(tab)[2]
    N <- sum(tab)
    n.g <- apply(tab, 1, sum)
    reltab <- sweep(tab, 1, n.g, FUN="/")
    Q.mat <- reltab %*% t(reltab)
    C <- Q.mat
    diag(C) <- 0
    n.k <- apply(tab, 2, sum)
    Q.k <- n.k / N

    # pooled PMI (Scofield et al 2010 J Ecol)
    q.0.gh <- prop.y.0.gh <- matrix(0, G, G)
    dimnames(q.0.gh) <- dimnames(prop.y.0.gh) <- dimnames(Q.mat)
    for (g in 1:(G - 1)) {
        for (h in (g + 1):G) {
            shared.k <- tab[g, ] & tab[h, ]
            y.gh <- prop.y.0.gh[g, h] <- sum(tab[g, shared.k])
            y.hg <- prop.y.0.gh[h, g] <- sum(tab[h, shared.k])
            prop.y.0.gh[g, h] <- y.gh / n.g[g]
            prop.y.0.gh[h, g] <- y.hg / n.g[h]
            q.0.gh[g, h] <- (y.gh * y.hg)/(n.g[g] * n.g[h])
        }
    }

    # weighted q.gh (Scofield et al 2011 Oecologia)
    n.gh <- outer(n.g, n.g, '*'); diag(n.gh) <- 0
    q.bar.gh <- sum(n.gh * C) / sum(n.gh)
    # note n.gh is 0 along diagonal so diag(C - q.bar.0) is ignored
    q.variance.gh <- sum(n.gh) * sum(n.gh * (C - q.bar.gh)) /
        (sum(n.gh)^2 - sum(n.gh * n.gh))

    # Separate lists for standard calculations:
    #     q
    # classical bias correction:
    #     r
    # and Nielsen et al. 2003 bias correction:
    #     q.nielsen
    #
    # Each list has the same elements
    #
    #    q.gg                   based on squared frequencies
    #    q.bar.0                weighted mean of q.gg
    #    q.variance.0           weighted sample variance of q.gg
    #    q.unweighted.mean      unwieghted mean of q.gg
    #    alpha.g                reciprocal of individual q.gg
    #    alpha.unweighted.mean  mean of alpha.g
    #    d.alpha                reciprocal of q.unweighted.mean  (TODO RENAME?)
    #    d.gamma                reciprocal of sum of grand squared frequencies  (TODO RENAME?)
    #    d.beta                 d.gamma / d.alpha  (TODO RENAME?)
    #    diversity.mat          Chao pairwise off-diagonal, alpha.g on diagonal  (TODO RENAME?)
    #    divergence.mat         diversity.mat with 0 diagonal  (TODO RENAME?)
    #    overlap.mat            1 - divergence.mat  (TODO RENAME?)
    #    overlap                mean overlap  (TODO RENAME?)
    #    divergence             mean divergence  (TODO RENAME?)

    q = list()

    q$q.gg <- diag(Q.mat)
    q$q.bar.0 <- sum(n.g * n.g * q$q.gg) / sum(n.g * n.g)

    q$q.variance.0 <- (sum(n.g * n.g) * sum(n.g * n.g * (q$q.gg - q$q.bar.0))) /
        (sum(n.g * n.g)^2 - sum(n.g * n.g * n.g * n.g))

    q$q.unweighted.mean <- mean(q$q.gg)
    q$alpha.g <- 1 / q$q.gg
    q$alpha.unweighted.mean <- sum(q$alpha.g) / G
    q$d.alpha <- 1 / q$q.unweighted.mean
    q$d.gamma <- 1 / sum(Q.k * Q.k)
    q$d.beta <- q$d.gamma / q$d.alpha
    q$diversity.mat <- (2 * Q.mat) / 
                       outer(q$q.gg, q$q.gg, FUN="+")
    diag(q$diversity.mat) <- q$alpha.g
    q$divergence.mat <- 1 - q$diversity.mat
    diag(q$divergence.mat) <- 0
    q$overlap.mat <- 1 - q$divergence.mat
    q$overlap <- sum(C) / ((G - 1) * sum(q$q.gg))
    q$divergence <- 1 - q$overlap

    q.nielsen = list()

    q.nielsen$q.gg <- nielsenTransform(diag(Q.mat), n.g)
    q.nielsen$q.bar.0 <- sum(n.g * n.g * q.nielsen$q.gg) / sum(n.g * n.g)

    q.nielsen$q.variance.0 <- (sum(n.g * n.g) * sum(n.g * n.g *
                                   (q.nielsen$q.gg - q.nielsen$q.bar.0))) /
        (sum(n.g * n.g)^2 - sum(n.g * n.g * n.g * n.g))

    q.nielsen$q.unweighted.mean <- mean(q.nielsen$q.gg)
    q.nielsen$alpha.g <- 1 / q.nielsen$q.gg
    q.nielsen$alpha.unweighted.mean <- sum(q.nielsen$alpha.g) / G
    q.nielsen$d.alpha <- 1 / q.nielsen$q.unweighted.mean
    q.nielsen$d.gamma <- 1 / nielsenTransform(sum(Q.k * Q.k), N)
    q.nielsen$d.beta <- q.nielsen$d.gamma / q.nielsen$d.alpha
    q.nielsen$diversity.mat <- (2 * Q.mat) / 
                               outer(q.nielsen$q.gg, q.nielsen$q.gg, FUN="+")
    diag(q.nielsen$diversity.mat) <- q.nielsen$alpha.g
    q.nielsen$divergence.mat <- 1 - q.nielsen$diversity.mat
    diag(q.nielsen$divergence.mat) <- 0
    q.nielsen$overlap.mat <- 1 - q.nielsen$divergence.mat
    q.nielsen$overlap <- sum(C) / ((G - 1) * sum(q.nielsen$q.gg))
    q.nielsen$divergence <- 1 - q.nielsen$overlap

    r = list()

    r.gg <- numeric(G)
    for (g in 1:G) {
        tmp.numer <- (tab[g, ] * tab[g, ]) - tab[g, ]
        tmp.denom <- (n.g[g] * n.g[g]) - n.g[g]
        r.gg[g] <- sum(tmp.numer / tmp.denom)
    }
    names(r.gg) <- names(q$q.gg)
    r$q.gg <- r.gg
    r$q.bar.0 <- sum((n.g * n.g * r$q.gg) - (n.g * r$q.gg)) / 
                 sum((n.g * n.g) - n.g)

    M.g <- n.g * (n.g - 1)
    r$q.variance.0 <- (sum(M.g) * sum(M.g * (r$q.gg - r$q.bar.0))) /
        (sum(M.g)^2 - sum(M.g * M.g))

    r$q.unweighted.mean <- mean(r$q.gg)
    r$alpha.g <- 1 / r$q.gg
    r$alpha.unweighted.mean <- sum(r$alpha.g) / G
    r$d.alpha <- 1 / r$q.unweighted.mean
    r$d.gamma <- 1/sum((n.k * (n.k - 1)) / (N * (N - 1)))
    r$d.beta <- r$d.gamma / r$d.alpha
    r$diversity.mat <- (2 * Q.mat) / 
                       outer(r$q.gg, r$q.gg, FUN="+")
    diag(r$diversity.mat) <- r$alpha.g
    r$divergence.mat <- 1 - r$diversity.mat
    diag(r$divergence.mat) <- 0
    r$overlap.mat <- 1 - r$divergence.mat
    r$overlap <- sum(C) / ((G - 1) * sum(r$q.gg))
    r$divergence <- 1 - r$overlap

    # return value
    ans <- list(table       = orig.tab, # table passed in
                num.groups  = G,   # number of rows (sites)
                num.sources = K,   # number of columns (sources)
                num.samples = N, 
                num.samples.group  = n.g,
                num.sources.group  = apply(tab, 1, function(x) sum(x > 0)),
                num.samples.source = n.k,
                num.groups.source  = apply(tab, 2, function(x) sum(x > 0)),

                # Pooled PMI (Scofield et al 2010 J Ecol)
                y.gh        = y.gh,
                prop.y.0.gh = prop.y.0.gh,
                q.0.gh      = q.0.gh,

                # Weighted q.gh (Scofield et al 2010 Oecologia)
                q.bar.gh      = q.bar.gh,
                q.variance.gh = q.variance.gh,

                # Q matrix, used to reconstruct C matrix in (at least) allele code
                Q.mat     = Q.mat,

                # diversity calculations
                q         = q,
                r         = r, 
                q.nielsen = q.nielsen)
    structure(ans, class = c('diversity', 'list'))
}

