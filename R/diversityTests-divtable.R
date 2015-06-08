#' @include diversity-divtable.R
NULL



#' Print the result of a diversity test (alpha or gamma)
#'
#' Print an object of class \code{diversity_test}, the result of
#' \code{\link{alphaDiversityTest}}, \code{\link{alphaContrastTest}},
#' \code{\link{alphaContrastTest3}}, \code{\link{gammaContrastTest}} and
#' \code{\link{gammaContrastTest3}} when given an argument of class
#' \code{\link{divtable}}.
#'
#' @param x       Object of class \code{diversity_test}, returned by
#' \code{\link{alphaDiversityTest}}, \code{\link{alphaContrastTest}},
#' \code{\link{alphaContrastTest3}}, \code{\link{gammaContrastTest}} and
#' \code{\link{gammaContrastTest3}} when given an argument of class
#' \code{\link{allele_divtables}}.
#'
#' @param digits  Number of significant digits to use when printing
#' numeric values
#'
#' @return \code{x}, returned invisibly
#'
#' @seealso \code{\link{alphaDiversityTest}}, \code{\link{alphaContrastTest}}, \code{\link{alphaContrastTest3}}, \code{\link{gammaContrastTest}}, \code{\link{gammaContrastTest3}}
#'
#' @export
#'
print.diversity_test <- function(x, digits = getOption("digits"), ...)
{
    cat(x$method, "\n\n", sep = "")
    cat("data:  ", x$data.name, "\n", sep = "")
    # Sample sizes: those included depend on what is present in the object
    out <- c()
    .o.x.var <- function(o, v) {
        if (! is.null(x[[v]])) o <- c(o, paste(v, "=", x[[v]])) else o
    }
    out <- .o.x.var(out, "N")
    out <- .o.x.var(out, "N.a")
    out <- .o.x.var(out, "N.b")
    out <- .o.x.var(out, "N.c")
    out <- .o.x.var(out, "N.groups")
    cat("Samples", paste(out, collapse = ", "), "\n")
    cat("Observed log-likelihood ratio = ",
        format(signif(x$observed.ln.LR, max(1L, digits - 2L))), "\n\n",
        sep = "")
    cat("Test against empirical X^2 distribution:\n")
    cat("Iterations = ", x$n.resample,
        ", P = ", format.pval(x$P, digits = max(1L, digits - 3L)),
        "\n", sep = "")
    cat("\nQuantiles of the empirical distribution:\n")
    print(x$quantiles, digits = digits, ...)
    invisible(x)
}



#' Plot the result of a diversity test (alpha or gamma)
#'
#' Plot an object of class \code{diversity_test}, the result of
#' \code{\link{alphaDiversityTest}}, \code{\link{alphaContrastTest}},
#' \code{\link{alphaContrastTest3}}, \code{\link{gammaContrastTest}} and
#' \code{\link{gammaContrastTest3}}.
#'
#' @note The method used to compress the x-axis when observed values greatly
#' exceed the empirical distribution has not been well thought through.
#'
#' @param x  Object of class \code{diversity_test}, returned by
#' \code{\link{alphaDiversityTest}}, \code{\link{alphaContrastTest}},
#' \code{\link{alphaContrastTest3}}, \code{\link{gammaContrastTest}} and
#' \code{\link{gammaContrastTest3}}
#'
#' @param breaks  Number of breaks to use when plotting the histogram of
#' the empirical distribution, passed to \code{\link{hist(..., plot = FALSE)}}
#'
#' @param compress.x  Logical, if \code{TRUE} and the observed value is more
#' than \code{compress.ratio} times the maximum value of the empirical
#' distribution, the x-axis is compressed to include the observed value
#'
#' @param compress.ration See \code{compress.x}
#'
#' @param xlab,ylab,main  Labels added to the plot
#'
#' @param legend.text Text to use when printing legend containing the
#' observed value in the upper right of the plot. Set to \code{NULL} to
#' suppress printing the legend.
#'
#' @param digits  Number of significant digits to use when printing numeric
#' values
#'
#' @seealso \code{\link{alphaDiversityTest}}, \code{\link{alphaContrastTest}}, \code{\link{alphaContrastTest3}}, \code{\link{gammaContrastTest}}, \code{\link{gammaContrastTest3}}
#'
#' @export
#'
plot.diversity_test <- function(x, breaks = 50, compress.x = TRUE,
    compress.ratio = 1.2, xlab = "ln(LR) value", ylab = "Frequency",
    main = x$method, legend.text = "Observed ln(LR) = ",
    digits = getOption("digits"), ...)
{
    par(mar = c(3.5, 3.5, 1.5, 1), mgp = c(2, 0.5, 0))
    plotted.observed.ln.LR <- x$observed.ln.LR
    empdist.range <- diff(range(x$empdist[-x$n.resample]))
    full.range <- diff(range(x$empdist))
    if (compress.x && full.range > empdist.range * compress.ratio) {
        xlim <- range(x$empdist[-x$n.resample])
        h <- hist(x$empdist[-x$n.resample],
                  breaks = seq(xlim[1], xlim[2], length.out = breaks),
                  plot = FALSE)
        xlim <- xlim * c(1, compress.ratio)
        plotted.observed.ln.LR <- xlim[2]
    } else {
        h <- hist(x$empdist[-x$n.resample], breaks = breaks, plot = FALSE)
        xlim <- range(c(h$breaks, x$empdist[x$n.resample]))
    }
    ylim <- range(c(0, h$count))
    plot(h, xlim = xlim, ylim = ylim, freq = TRUE,
         xlab = ylab, ylab = ylab, main = main, ...)
    lines(rep(plotted.observed.ln.LR, 2), c(ylim[1], ylim[2]*0.5),
          col = "darkgray", lty = 1, lwd = 2)
    if (! is.null(legend.text))
        legend("topright", bty = "n", legend = paste0(legend.text,
               format(signif(x$observed.ln.LR, max(1L, digits - 2L)))))
}



#' Test for differences in alpha diversity among sites within a single data set
#'
#' The null hypothesis for this tests is that there is no difference in alpha
#' diversity between the sites represented in \code{tab} or \code{adt}.  The
#' initial (class \code{\link{divtable}}) version of this was described in
#' Scofield et al. (2012), while the allelic (class
#' \code{link{allele_divtables}} extension was described in Sort et al. (In
#' press).
#'
#' @param tab    Site-by-source table of class \code{\link{divtable}}
#'
#' @param adt    Allele diversity dataset of class
#' \code{\link{allele_divtables}}
#'
#' @param zero.div.adjust Logical, if \code{TRUE} (the default), then groups
#' with 0 within-group diversity are assigned a minimum diversity which is
#' half the empirical diversity possible given the group size
#'
#' @param n.resample Number of iterations for creation of the null distribution
#'
#' @param method \code{"bootstrap"} or \code{"permute"}, whether to create null
#' distribution iterations with (\code{"bootstrap"}) or without
#' (\code{"permute"}) replacement
#'
#' @param \dots  Additional parameters
#'
#' @return An \code{diversity_test} object with the result of the test
#
#' @seealso \code{\link{alphaContrastTest}}, \code{\link{alphaContrastTest3}}, \code{\link{print.diversity_test}}
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
#' ## Add example with class divtable
#' ##
#' ## Using allele diversity dataset of class allele_divtables.  Compare
#' ## allele diversity between sites in the same sample:
#' ##
#' ## library(readGenalex)  # already loaded as a prerequisite
#' data(Qagr_pericarp_genotypes)  # from readGenalex
#' gt <- createAlleleTables(Qagr_pericarp_genotypes)
#' alpha.test <- alphaDiversityTest(gt)
#'
#' @export
#'
#' @name alphaDiversityTest
#'
NULL

alphaDiversityTest <- function(tab, ...) UseMethod("alphaDiversityTest")



#' @rdname alphaDiversityTest
#'
#' @export
#'
alphaDiversityTest.divtable <- function(tab, zero.div.adjust = TRUE,
    n.resample = 10000, method = c("bootstrap", "permute"),
    test.quantiles = c(0.001, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 0.999),
    ...)
{
    method <- match.arg(method)
    .checkRowSums(tab)
    ans <- list(method = "Alpha diversity test, contrast among sites in single data set")
    ans$data.name <- deparse(substitute(tab))
    g.vardist <- .diversityTest.directGowerDiag(tab)
    n.g <- sapply(g.vardist, length)
    N <- sum(n.g)
    G <- length(n.g)
    ans$N.samples <- N
    ans$N.groups <- G
    terms = .diversityTest.CalcTerms(n.g, g.vardist, zero.div.adjust)
    ans$observed.ln.LR <- terms$ln.LR
    # Empirical distribution
    nulldist <- .diversityTest.NullDist(obs = terms$ln.LR, n.g = n.g,
        g.vardist = g.vardist, zero.div.adjust = zero.div.adjust,
        method = method, n.resample=n.resample)
    ans$n.resample <- n.resample
    ans$resample.method <- method
    ans$quantiles <- quantile(nulldist, test.quantiles)
    ans$P <- sum(terms$ln.LR <= nulldist) / n.resample
    ans$empdist <- nulldist
    structure(ans, class = c('diversity_test', 'list'))
}



#' @rdname alphaDiversityTest
#'
#' @export
#'
alphaDiversityTest.default <- function(tab, ...)
{
    stop(deparse(substitute(tab)),
         " must be of class divtable or class allele_divtables")
}



#' Test for differences in alpha diversity between two data sets
#'
#' The null hypothesis for this tests is that there is no difference in alpha
#' diversity between the two sets of datasets sites represented in \code{tab.a}
#' and \code{tab.b}, or \code{adt.a} and \code{adt.b}.  The method for the
#' initial class \code{\link{divtable}} version of this was described in
#' Scofield et al. (2012), while that for the allelic extension as class
#' \code{link{allele_divtables}} was described in Sork \emph{et al}.
#' (In press).
#'
#' @note \code{adt.a} and \code{adt.b} must contain the same loci
#'
#' @param tab.a Site-by-source table of class \code{\link{divtable}}
#'
#' @param tab.b Site-by-source table of class \code{\link{divtable}}
#'
#' @param adt.a Allelic diversity dataset of class
#' \code{\link{allele_divtables}}
#'
#' @param adt.b Allelic diversity dataset of class
#' \code{\link{allele_divtables}}
#'
#' @param zero.div.adjust Logical, if \code{TRUE} (the default), then groups
#' with 0 within-group diversity are assigned a minimum diversity which is
#' half the empirical diversity possible given the group size
#'
#' @param n.resample Number of iterations for creation of the null distribution
#'
#' @param method \code{"bootstrap"} or \code{"permute"}, whether to create null
#' distribution iterations with (\code{"bootstrap"}) or without
#' (\code{"permute"}) replacement
#'
#' @return An \code{diversity_test} object with the result of the test
#'
#' @seealso \code{\link{alphaDiversityTest}}, \code{\link{alphaContrastTest3}}
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
#' ## Comparing alpha diversity between two different sites:
#' ##
#' data(granaries_2002_Qlob)
#' data(granaries_2004_Qlob)
#' par(mfcol = c(2, 1))
#' membershipPlot(granaries_2002_Qlob)
#' membershipPlot(granaries_2004_Qlob)
#' alphaContrastTest(granaries_2002_Qlob, granaries_2004_Qlob)
#'
#' ## Comparing allele diversity between two different samples:
#' ##
#' # dat1 <- readGenalex("file-of-genotypes-sample-1.txt")
#' # dat2 <- readGenalex("file-of-genotypes-sample-2.txt")
#' # gt1 <- createAlleleTables(dat1)
#' # gt2 <- createAlleleTables(dat2)
#' # alpha.contrast <- alphaContrastTest(gt1, gt2)
#'
#' @export
#'
#' @name alphaContrastTest
#'
NULL

alphaContrastTest <- function(a, b, ...) UseMethod("alphaContrastTest")



#' @rdname alphaContrastTest
#'
#' @export
#"
alphaContrastTest.divtable <- function(tab.a, tab.b, zero.div.adjust = TRUE,
    n.resample = 10000, method = c("bootstrap", "permute"),
    test.quantiles = c(0.001, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 0.999))
{
    stopifnot(inherits(tab.b, 'divtable'))
    .checkRowSums(tab.a)
    .checkRowSums(tab.b)
    method <- match.arg(method)
    ans <- list(method = "Alpha diversity test, contrast between 2 datasets")
    ans$data.name <- paste(sep = ", ", deparse(substitute(tab.a)),
        deparse(substitute(tab.b)))
    a.vardist <- .diversityTest.directGowerDiag(tab.a)
    n.a <- sapply(a.vardist, length)
    N.a <- sum(n.a)
    G.a <- length(n.a)
    ans$N.a <- N.a
    ans$G.a <- G.a
    terms.a = .diversityTest.CalcTerms(n.a, a.vardist, zero.div.adjust)
    # V.a.p = terms.a$V.p

    b.vardist <- .diversityTest.directGowerDiag(tab.b)
    n.b <- sapply(b.vardist, length)
    N.b <- sum(n.b)
    G.b <- length(n.b)
    ans$N.b <- N.b
    ans$G.b <- G.b
    terms.b = .diversityTest.CalcTerms(n.b, b.vardist, zero.div.adjust)
    # V.b.p = terms.b$V.p
    V.a.b.p <- (((N.a - G.a) * terms.a$V.p) + ((N.b - G.b) * terms.b$V.p)) /
              (N.a + N.b - G.a - G.b)
    observed.ln.LR.a.b <- ((N.a + N.b - G.a - G.b) * log(V.a.b.p)) -
                          ((N.a - G.a) * log(terms.a$V.p)) -
                          ((N.b - G.b) * log(terms.b$V.p))

    # Combine A and B into strata for comparison
    n.a.b <- c(a = N.a, b = N.b)
    a.b.vardist <- list(a = unlist(a.vardist, use.names=FALSE),
                        b = unlist(b.vardist, use.names=FALSE))
    N <- sum(n.a.b)
    G <- length(n.a.b)
    DF <- G -1
    ans$N.samples <- N
    ans$N.groups <- G
    ans$observed.ln.LR <- observed.ln.LR.a.b
    # Empirical distribution
    nulldist <- .diversityTest.NullDist(obs=observed.ln.LR.a.b,
                                        n.g = n.a.b,
                                        g.vardist = a.b.vardist,
                                        zero.div.adjust, method,
                                        n.resample)
    ans$n.resample <- n.resample
    ans$resample.method <- method
    ans$quantiles <- quantile(nulldist, test.quantiles)
    ans$P <- sum(observed.ln.LR.a.b <= nulldist) / n.resample
    ans$empdist <- nulldist
    ans$a.b.vardist = a.b.vardist
    structure(ans, class = c('diversity_test', 'list'))
}



#' @rdname alphaContrastTest
#'
#' @export
#'
alphaContrastTest.default <- function(tab.a, tab.b, ...)
{
    args <- c()
    if (! inherits(tab.a, "divtable"))
        args <- c(args, deparse(substitute(tab.a)))
    if (! inherits(tab.b, "divtable"))
        args <- c(args, deparse(substitute(tab.b)))
    stop(paste(collapse = " and ", args),
         " must be of class divtable or class allele_divtables")
}



#' Test for differences in alpha diversity between three data sets
#'
#' The null hypothesis for this tests is that there is no difference in alpha
#' diversity between the three sets of datasets sites represented in \code{tab}
#' and \code{tab}, or or \code{adt.a} and \code{adt.b}.  The method for the
#' initial class \code{\link{divtable}} version of this was described in
#' Scofield et al. (2012).
#' There currently is no support for contrasting allelic data.
#'
#' The null hypothesis for this tests is that there is no difference in alpha
#' diversity between the three sets of datasets sites represented in
#' \code{tab.a}, \code{tab.b} and \code{tab.c}.  The method was described in
#' Scofield et al. (2012).
#'
#' @param tab.a First site-by-source table, of class \code{\link{divtable}}
#'
#' @param tab.b Second site-by-source table, of class \code{\link{divtable}}
#'
#' @param tab.c Third site-by-source table, of class \code{\link{divtable}}
#'
#' @param zero.div.adjust Logical, if \code{TRUE} (the default), then groups
#' with 0 within-group diversity are assigned a minimum diversity which is
#' half the empirical diversity possible given the group size
#'
#' @param n.resample Number of iterations for creation of the null distribution
#'
#' @param method \code{"bootstrap"} or \code{"permute"}, whether to create null
#' distribution iterations with (\code{"bootstrap"}) or without
#' (\code{"permute"}) replacement
#'
#' @return An \code{diversity_test} object with the result of the test
#'
#' @seealso \code{\link{alphaDiversityTest}}, \code{\link{alphaContrastTest}}
#'
#' @references
#'
#' Scofield, D. G., Smouse, P. E., Karubian, J. and Sork, V. L. (2012)
#' Use of alpha, beta and gamma diversity measures to characterize seed
#' dispersal by animals.  \emph{American Naturalist} 180:719-732.
#'
# @examples
#
# # Suitable examples would require three datasets
#'
#' @export
#'
#' @name alphaContrastTest3
#'
NULL

alphaContrastTest3 <- function(a, b, c, ...)
    UseMethod("alphaContrastTest3")



#' @rdname alphaContrastTest3
#'
#' @export
#'
alphaContrastTest3.divtable <- function(tab.a, tab.b, tab.c,
    zero.div.adjust = TRUE, n.resample = 10000,
    method = c("bootstrap", "permute"),
    test.quantiles = c(0.001, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 0.999),
    ...)
{
    stopifnot(inherits(tab.b, 'divtable'))
    stopifnot(inherits(tab.c, 'divtable'))
    .checkRowSums(tab.a)
    .checkRowSums(tab.b)
    .checkRowSums(tab.c)
    method <- match.arg(method)
    ans <- list(method = "Alpha diversity test, contrast between 3 datasets")
    ans$data.name <- paste(sep = ", ", deparse(substitute(tab.a)),
        deparse(substitute(tab.b)), deparse(substitute(tab.b)))
    a.vardist <- .diversityTest.directGowerDiag(tab.a)
    n.a <- sapply(a.vardist, length)
    N.a <- sum(n.a)
    G.a <- length(n.a)
    ans$N.a <- N.a
    ans$G.a <- G.a
    terms.a <- .diversityTest.CalcTerms(n.a, a.vardist, zero.div.adjust)
    # V.a.p <- terms.a$V.p

    b.vardist <- .diversityTest.directGowerDiag(tab.b)
    n.b <- sapply(b.vardist, length)
    N.b <- sum(n.b)
    G.b <- length(n.b)
    ans$N.b <- N.b
    ans$G.b <- G.b
    terms.b <- .diversityTest.CalcTerms(n.b, b.vardist, zero.div.adjust)
    # V.b.p <- terms.b$V.p

    c.vardist <- .diversityTest.directGowerDiag(tab.c)
    n.c <- sapply(c.vardist, length)
    N.c <- sum(n.c)
    G.c <- length(n.c)
    ans$N.c <- N.c
    ans$G.c <- G.c
    terms.c <- .diversityTest.CalcTerms(n.c, c.vardist, zero.div.adjust)
    # V.c.p <- terms.c$V.p

    V.a.b.c.p <- (((N.a - G.a) * terms.a$V.p) +
                  ((N.b - G.b) * terms.b$V.p) +
                  ((N.c - G.c) * terms.c$V.p)) /
                 (N.a + N.b + N.c - G.a - G.b - G.c)
    observed.ln.LR.a.b.c <-
        ((N.a + N.b + N.c - G.a - G.b - G.c) * log(V.a.b.c.p)) -
        ((N.a - G.a) * log(terms.a$V.p)) -
        ((N.b - G.b) * log(terms.b$V.p)) -
        ((N.c - G.c) * log(terms.c$V.p))


    # Combine A B C into strata for comparison
    n.a.b.c <- c(a = N.a, b = N.b, c = N.c)
    a.b.c.vardist <- list(a = unlist(a.vardist, use.names = FALSE),
                          b = unlist(b.vardist, use.names = FALSE),
                          c = unlist(c.vardist, use.names = FALSE))
    N <- sum(n.a.b.c)
    G <- length(n.a.b.c)
    DF <- G - 1
    ans$N.samples <- N
    ans$N.groups <- G
    ans$observed.ln.LR <- observed.ln.LR.a.b.c
    # Empirical distribution
    nulldist <- .diversityTest.NullDist(obs = observed.ln.LR.a.b.c,
                                        n.g = n.a.b.c,
                                        g.vardist = a.b.c.vardist,
                                        zero.div.adjust, method,
                                        n.resample)
    ans$n.resample <- n.resample
    ans$resample.method <- method
    ans$quantiles <- quantile(nulldist, test.quantiles)
    ans$P <- sum(observed.ln.LR.a.b.c <= nulldist) / n.resample
    ans$empdist <- nulldist
    ans$a.b.c.vardist <- a.b.c.vardist
    structure(ans, class = c('diversity_test', 'list'))
}



#' @rdname alphaContrastTest3
#'
#' @export
#'
alphaContrastTest3.default <- function(tab.a, tab.b, tab.c, ...)
{
    args <- c()
    if (! inherits(tab.a, "divtable"))
        args <- c(args, deparse(substitute(tab.a)))
    if (! inherits(tab.b, "divtable"))
        args <- c(args, deparse(substitute(tab.b)))
    if (! inherits(tab.c, "divtable"))
        args <- c(args, deparse(substitute(tab.c)))
    stop(paste(collapse = " and ", args),
         " must be of class divtable or class allele_divtables")
}



#' Test for structure in pairwise divergence between sites
#'
#' @export
#'
#' @name pairwiseMeanTest
#'
NULL

pairwiseMeanTest <- function(a, ...) UseMethod("pairwiseMeanTest")



#' @rdname pairwiseMeanTest
#'
#' @export
#'
pairwiseMeanTest.divtable <- function(tab, n.iter = 10000,
    method = c("r", "q", "q.nielsen"), statistic = c("divergence", "overlap"),
    test.quantiles = c(0.001, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 0.999),
    ...)
{
    method <- match.arg(method)
    statistic <- match.arg(statistic)
    d <- diversity(tab)
    gammafreq <- d$num.samples.source
    K <- d$num.sources
    n.g <- d$num.samples.group
    names.g <- rownames(tab)
    G <- length(n.g)
    obs <- NULL
    p <- d[[method]]
    obs <- p[[statistic]]

    nulldist <- numeric(n.iter)
    for (i in 1:(n.iter - 1)) {
        mat <- c()
        for (g in 1:G)
            mat <- rbind(mat, t(rmultinom(1, n.g[g], gammafreq)))
        reltab <- sweep(mat, 1, n.g, FUN="/")
        Q.mat <- reltab %*% t(reltab)
        cache.gg <- diag(Q.mat) # q.gg
        diag(Q.mat) <- 0
        cache.gg <- switch(method,
                           q = cache.gg,
                           r = (((n.g * cache.gg) - 1) / (n.g - 1)),
                           q.nielsen = nielsenTransform(cache.gg, n.g))
        stat <- sum(Q.mat) / ((G - 1)*sum(cache.gg))
        nulldist[i] <- if (statistic == "divergence") (1 - stat) else stat
    }
    nulldist[n.iter] <- obs
    P.lower <- sum(obs >= nulldist) / n.iter
    P.upper <- sum(obs <= nulldist) / n.iter
    q2 <- quantiles(nulldist, test.quantiles)
    ans <- list(n.cache = G, n.source = K, n.seed = sum(n.g),
                obs = obs, n.iter = n.iter, nulldist = nulldist,
                P.lower = P.lower, P.upper = P.upper,
                quantiles = q1, method = method, statistic = statistic)
    structure(ans, class = c('pairwise_mean_test', 'list'))
}



#' @rdname pairwiseMeanTest
#'
#' @export
#'
pairwiseMeanTest.default <- function(tab, ...)
{
    pairwiseMeanTest.divtable(as.divtable(tab), ...)
}



#' Plot the result of \code{\link{pairwiseMeanTest}}
#'
#' @seealso \code{\link{pairwiseMeanTest}}
#'
#' @export
#'
plot.pairwise_mean_test <- function(result, ...)
{
    if (! inherits(result, 'pairwise_mean_test'))
        stop(deparse(substitute(result)), "not a result of pairwiseMeanTest()")
    par(mar = c(2.8, 2.8, 0.5, 0), mgp = c(1.6, 0.4, 0), tcl = -0.25, cex = 0.9)
    xlim <- range(c(0, 1, result$nulldist))
    breaks <- seq(xlim[1], xlim[2], length.out=50)
    h <- hist(result$nulldist[-result$n.iter], breaks = breaks, plot = FALSE)
    ylim <- range(c(0, h$counts))
    xlab <- switch(result$statistic,
                overlap = switch(result$method,
                    q = expression(bar(omega)[q]*" value"),
                    q.nielsen = expression(bar(omega)[q^'*']*" value"),
                    r = expression(bar(omega)[r]*" value")),
                divergence = switch(result$method,
                    q = expression(bar(delta)[q]*" value"),
                    q.nielsen = expression(bar(delta)[q^'*']*" value"),
                    r = expression(bar(delta)[r]*" value")))
    plot(h, xlim = xlim, ylim = ylim, freq = TRUE, xlab = xlab,
         ylab = "Frequency", main = "", ...)
    lines(rep(result$obs, 2), c(ylim[1], ylim[2]*0.5), col="black", lty=2, lwd=4)
    OBS <- round(result$obs, 3)
    leg <- switch(result$statistic,
               overlap = switch(result$method,
                   q = bquote("Observed "*bar(omega)[q] == .(OBS)),
                   q.nielsen = bquote("Observed "*bar(omega)[q^'*'] == .(OBS)),
                   r = bquote("Observed "*bar(omega)[r] == .(OBS))),
               divergence = switch(result$method,
                   q = bquote("Observed "*bar(delta)[q] == .(OBS)),
                   q.nielsen = bquote("Observed "*bar(delta)[q^'*'] == .(OBS)),
                   r = bquote("Observed "*bar(delta)[r] == .(OBS))))
    legend("topleft", bty = "n", legend = leg)
}



#' Test for difference in gamma diversity between two datasets
#'
#' The null hypothesis for this tests is that there is no difference in gamma
#' diversity between the two sets of datasets sites represented in \code{tab.a}
#' and \code{tab.b}, or \code{adt.a} and \code{adt.b}.  The method for the
#' initial class \code{\link{divtable}} version of this was described in
#' Scofield et al. (2012), while that for the allelic extension as class
#' \code{link{allele_divtables}} was described in Sork \emph{et al}.
#' (in press).
#'
#' @note \code{adt.a} and \code{adt.b} must contain the same loci
#'
#' @param tab.a Site-by-source table of class \code{\link{divtable}}
#'
#' @param tab.b Site-by-source table of class \code{\link{divtable}}
#'
#' @param adt.a Allelic diversity dataset of class
#' \code{\link{allele_divtables}}
#'
#' @param adt.b Allelic diversity dataset of class
#' \code{\link{allele_divtables}}
#'
#' @param zero.div.adjust Logical, if \code{TRUE} (the default), then groups
#' with 0 within-group diversity are assigned a minimum diversity which is
#' half the empirical diversity possible given the group size
#'
#' @param n.resample Number of iterations for creation of the null distribution
#'
#' @param method \code{"bootstrap"} or \code{"permute"}, whether to create null
#' distribution iterations with (\code{"bootstrap"}) or without
#' (\code{"permute"}) replacement
#'
#' @return  Class \code{diversity_test} object with the result of the test
#'
#' @seealso \code{\link{gammaContrastTest3}}, \code{\link{alphaContrastTest}}
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
#' ## Compare gamma diversity between two diversity datasets
#' ##
#' data(granaries_2002_Qlob)
#' data(granaries_2004_Qlob)
#' par(mfrow = c(1, 2))
#' plot(gammaAccum(granaries_2002_Qlob))
#' plot(gammaAccum(granaries_2004_Qlob))
#' gammaContrastTest(granaries_2002_Qlob, granaries_2004_Qlob)
#'
#' ## Compare gamma diversity between two allele datasets
#' ##
#' # dat1 <- readGenalex("file-of-genotypes-sample-1.txt")
#' # dat2 <- readGenalex("file-of-genotypes-sample-2.txt")
#' # gt1 <- createAlleleTables(dat1)
#' # gt2 <- createAlleleTables(dat2)
#' # gamma.contrast <- gammaContrastTest(gt1, gt2)
#'
#' @export
#'
#' @name gammaContrastTest
#'
NULL

gammaContrastTest <- function(a, b, ...) UseMethod("gammaContrastTest")



#' @rdname gammaContrastTest
#'
#' @export
#"
gammaContrastTest.divtable <- function(tab.a, tab.b, zero.div.adjust = TRUE,
    n.resample = 10000, method = c("bootstrap", "permute"),
    test.quantiles = c(0.001, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 0.999),
    ...)
{
    stopifnot(inherits(tab.b, 'divtable'))
    .checkRowSums(tab.a)
    .checkRowSums(tab.b)
    ans <- list(method = "Gamma diversity test, contrast between two data sets")
    ans$data.name <- paste(sep = ", ", deparse(substitute(tab.a)),
        deparse(substitute(tab.b)))
    X.a.k <- apply(tab.a, 2, sum)
    X.b.k <- apply(tab.b, 2, sum)
    N.a <- sum(X.a.k)
    N.b <- sum(X.b.k)
    R.a.0 <- sum((X.a.k * (X.a.k - 1)) / (N.a * (N.a - 1)))
    R.b.0 <- sum((X.b.k * (X.b.k - 1)) / (N.b * (N.b - 1)))
    V.a.tot <- 1 - R.a.0
    V.b.tot <- 1 - R.b.0
    V.a.b.tot <- ((N.a - 1) * V.a.tot + (N.b - 1) * V.b.tot) / (N.a + N.b - 2)
    observed.ln.LR.a.b <- ((N.a + N.b - 2) * log(V.a.b.tot)) -
                          ((N.a - 1) * log(V.a.tot)) -
                          ((N.b - 1) * log(V.b.tot))
    #a.vardist <- list(b = diag(.diversityTest.gower(.diversityTest.distmat(X.a.k))))
    a.vardist <- .diversityTest.directGowerDiag(X.a.k)
    n.a <- sapply(a.vardist, length)
    stopifnot(sum(n.a) == N.a)
    G.a <- length(n.a)
    ans$N.a <- N.a
    ans$G.a <- G.a
    terms.a = .diversityTest.CalcTerms(n.a, a.vardist, zero.div.adjust)
    #cat(sprintf("terms.a$V.p = %f  V.a.tot = %f\n", terms.a$V.p, V.a.tot))

    #b.vardist <- list(b = diag(.diversityTest.gower(.diversityTest.distmat(X.b.k))))
    b.vardist <- .diversityTest.directGowerDiag(X.b.k)
    n.b <- lapply(b.vardist, length)
    stopifnot(sum(n.b) == N.b)
    G.b <- length(n.b)
    ans$N.b <- N.b
    ans$G.b <- G.b
    terms.b = .diversityTest.CalcTerms(n.b, b.vardist, zero.div.adjust)
    #cat(sprintf("terms.b$V.p = %f  V.b.tot = %f\n", terms.b$V.p, V.b.tot))

    # Combine A and B into stratta for comparison
    n.a.b <- c(a = N.a, b = N.b)
    a.b.vardist <- list(a = unlist(a.vardist, use.names = FALSE),
                        b = unlist(b.vardist, use.names = FALSE))
    N <- sum(n.a.b)
    G <- length(n.a.b)
    DF <- G - 1
    ans$N.samples <- N
    ans$N.groups <- G
    ans$observed.ln.LR <- observed.ln.LR.a.b
    # Empirical distribution
    nulldist <- .diversityTest.NullDist(obs = observed.ln.LR.a.b,
                                        n.g = n.a.b, g.vardist = a.b.vardist,
                                        zero.div.adjust, method, n.resample)
    ans$n.resample <- n.resample
    ans$resample.method <- method
    ans$quantiles <- quantile(nulldist, test.quantiles)
    ans$P <- sum(observed.ln.LR.a.b <= nulldist) / n.resample
    ans$empdist <- nulldist
    ans$a.b.vardist = a.b.vardist
    structure(ans, class = c('diversity_test', 'list'))
}



#' @rdname gammaContrastTest
#'
#' @export
#'
gammaContrastTest.default <- function(tab.a, tab.b, ...)
{
    args <- c()
    if (! inherits(tab.a, "divtable"))
        args <- c(args, deparse(substitute(tab.a)))
    if (! inherits(tab.b, "divtable"))
        args <- c(args, deparse(substitute(tab.b)))
    stop(paste(collapse = " and ", args),
         " must be of class divtable or class allele_divtables")
}



#' Test for difference in gamma diversity among three datasets
#'
#' The null hypothesis for this tests is that there is no difference in gamma
#' diversity between the three datasets in \code{tab.a}, \code{tab.b} and
#' and \code{tab.c}.  The method was described in Scofield et al. (2012).
#' There currently is no support for contrasting allelic data.
#'
#' @param tab.a Site-by-source table of class \code{\link{divtable}}
#'
#' @param tab.b Site-by-source table of class \code{\link{divtable}}
#'
#' @param tab.c Site-by-source table of class \code{\link{divtable}}
#'
#' @param zero.div.adjust Logical, if \code{TRUE} (the default), then groups
#' with 0 within-group diversity are assigned a minimum diversity which is
#' half the empirical diversity possible given the group size
#'
#' @param n.resample Number of iterations for creation of the null distribution
#'
#' @param method \code{"bootstrap"} or \code{"permute"}, whether to create null
#' distribution iterations with (\code{"bootstrap"}) or without
#' (\code{"permute"}) replacement
#'
#' @return  Class \code{diversity_test} object with the result of the test
#'
#' @seealso \code{\link{gammaContrastTest}}, \code{\link{alphaContrastTest3}}
#'
#' @references
#'
#' Scofield, D. G., Smouse, P. E., Karubian, J. and Sork, V. L. (2012)
#' Use of alpha, beta and gamma diversity measures to characterize seed
#' dispersal by animals.  \emph{American Naturalist} 180:719-732.
#'
# @examples
#
#' @export
#'
#' @name gammaContrastTest3
#'
NULL

gammaContrastTest3 <- function(tab.a, tab.b, ...)
    UseMethod("gammaContrastTest3")



#' @rdname gammaContrastTest3
#'
#' @export
#"
gammaContrastTest3.divtable <- function(tab.a, tab.b, tab.c,
    zero.div.adjust = TRUE, n.resample = 10000,
    method = c("bootstrap", "permute"),
    test.quantiles = c(0.001, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 0.999),
    ...)
{
    stopifnot(inherits(tab.b, 'divtable'))
    stopifnot(inherits(tab.c, 'divtable'))
    .checkRowSums(tab.a)
    .checkRowSums(tab.b)
    .checkRowSums(tab.c)
    method <- match.arg(method)
    ans <- list(method = "Gamma diversity test, contrast between three data sets")
    ans$data.name <- paste(sep = ", ", deparse(substitute(tab.a)),
        deparse(substitute(tab.b)), deparse(substitute(tab.c)))
    X.a.k <- apply(tab.a, 2, sum)
    X.b.k <- apply(tab.b, 2, sum)
    X.c.k <- apply(tab.c, 2, sum)
    N.a <- sum(X.a.k)
    N.b <- sum(X.b.k)
    N.c <- sum(X.c.k)
    R.a.0 <- sum((X.a.k * (X.a.k - 1)) / (N.a * (N.a - 1)))
    R.b.0 <- sum((X.b.k * (X.b.k - 1)) / (N.b * (N.b - 1)))
    R.c.0 <- sum((X.c.k * (X.c.k - 1)) / (N.c * (N.c - 1)))
    V.a.tot <- 1 - R.a.0
    V.b.tot <- 1 - R.b.0
    V.c.tot <- 1 - R.c.0
    V.a.b.c.tot <- ((N.a - 1) * V.a.tot +
                    (N.b - 1) * V.b.tot +
                    (N.c - 1) * V.c.tot) / (N.a + N.b + N.c - 3)
    observed.ln.LR.a.b.c <- ((N.a + N.b + N.c - 3) * log(V.a.b.c.tot)) -
                          ((N.a - 1) * log(V.a.tot)) -
                          ((N.b - 1) * log(V.b.tot)) -
                          ((N.c - 1) * log(V.c.tot))
    # distances: distance matrix, then the diagonal of a Gower matrix
    #a.distmat <- .diversityTest.distmat(X.a.k)
    #a.vardist <- list(a = diag(.diversityTest.gower(a.distmat)))
    a.vardist <- .diversityTest.directGowerDiag(X.a.k)
    n.a <- sapply(a.vardist, length)
    stopifnot(sum(n.a) == N.a)
    N.a <- sum(n.a)
    G.a <- length(n.a)
    ans$N.a <- N.a
    ans$G.a <- G.a
    terms.a <- .diversityTest.CalcTerms(n.a, a.vardist, zero.div.adjust)
    #cat(sprintf("terms.a$V.p = %f  V.a.tot = %f\n", terms.a$V.p, V.a.tot))

    #b.distmat <- .diversityTest.distmat(X.b.k)
    #b.vardist <- list(b = diag(.diversityTest.gower(b.distmat)))
    b.vardist <- .diversityTest.directGowerDiag(X.b.k)
    n.b <- sapply(b.vardist, length)
    stopifnot(sum(n.b) == N.b)
    N.b <- sum(n.b)
    G.b <- length(n.b)
    ans$N.b <- N.b
    ans$G.b <- G.b
    terms.b <- .diversityTest.CalcTerms(n.b, b.vardist, zero.div.adjust)
    #cat(sprintf("terms.b$V.p = %f  V.b.tot = %f\n", terms.b$V.p, V.b.tot))

    #c.distmat <- .diversityTest.distmat(X.c.k)
    #c.vardist <- list(c = diag(.diversityTest.gower(c.distmat)))
    c.vardist <- .diversityTest.directGowerDiag(X.c.k)
    n.c <- sapply(c.vardist, length)
    stopifnot(sum(n.c) == N.c)
    N.c <- sum(n.c)
    G.c <- length(n.c)
    ans$N.c <- N.c
    ans$G.c <- G.c
    terms.c <- .diversityTest.CalcTerms(n.c, c.vardist, zero.div.adjust)
    #cat(sprintf("terms.c$V.p = %f  V.c.tot = %f\n", terms.c$V.p, V.c.tot))

    # Combine A B C into stratta for comparison
    n.a.b.c <- c(a = N.a, b = N.b, c = N.b)
    a.b.c.vardist <- list(a = unlist(a.vardist, use.names = FALSE),
                          b = unlist(b.vardist, use.names = FALSE),
                          c = unlist(c.vardist, use.names = FALSE))
    N <- sum(n.a.b.c)
    G <- length(n.a.b.c)
    DF <- G -1
    ans$N.samples <- N
    ans$N.groups <- G
    ans$observed.ln.LR <- observed.ln.LR.a.b.c
    # Empirical distribution
    nulldist <- .diversityTest.NullDist(obs = observed.ln.LR.a.b.c,
                                        n.g = n.a.b.c,
                                        g.vardist = a.b.c.vardist,
                                        zero.div.adjust, method, n.resample)
    ans$n.resample <- n.resample
    ans$resample.method <- method
    ans$quantiles <- quantile(nulldist, test.quantiles)
    ans$P <- sum(observed.ln.LR.a.b.c <= nulldist) / n.resample
    ans$empdist <- nulldist
    ans$a.b.c.vardist <- a.b.c.vardist
    structure(ans, class = c('diversity_test', 'list'))
}



#' @rdname gammaContrastTest3
#'
#' @export
#'
gammaContrastTest3.default <- function(tab.a, tab.b, tab.c, ...)
{
    args <- c()
    if (! inherits(tab.a, "divtable"))
        args <- c(args, deparse(substitute(tab.a)))
    if (! inherits(tab.b, "divtable"))
        args <- c(args, deparse(substitute(tab.b)))
    if (! inherits(tab.c, "divtable"))
        args <- c(args, deparse(substitute(tab.c)))
    stop(paste(collapse = " and ", args),
         " must be of class divtable or class allele_divtables")
}



# ------- Internal functions and constants



# Maximum tolerance for row and column sums for Gower matrix
.diversityTest.epsilon <- 1.0e-12



# Return diagonal of matrix of centroid distances for variances based on Gower (1966)
# directly from a site-by-source table.
#
# Gower JC. 1966. Some distance properties of latent root and vector
# methods used in multivariate analysis.  Biometrika 53:325-338.
#
.diversityTest.directGowerDiag <- function(tab, group = dimnames(tab)[[1]], drop = TRUE)
{
    if (is.null(dim(tab))) {
        dim(tab) <- c(1, length(tab))
        dimnames(tab) <- list(Site = "onedim", Group = names(tab))
    }
    if (dim(tab)[1] > 1 && is.null(group))
        stop("must supply group(s), all groups not supported")
    else if (missing(group) && dim(tab)[1] == 1)
        group <- 1
    N.G <- rowSums(tab)
    D <- list()
    for (g in group) {
        this.N.G <- N.G[g]  # total N for site
        n.K <- unname(tab[g, ][tab[g, ] > 0])  # nonzero sources for site
        # protect against integer overflow
        storage.mode(this.N.G) <- storage.mode(n.K) <- "double"
        # Calculate total number of 1s that would be found in the full distance
        # matrix for this site.  This is the total elements in the matrix,
        # minus the number of elements in each of the 'self' 0-matrices along
        # the diagonal.  Then calculate the mean using this count.  Note some
        # algebraic simplification was used here, the mean is -0.5*(N^2 -
        # sum(n.k^2)) / N^2
        mean.d <- -0.5 + (sum(n.K * n.K) / (2 * this.N.G * this.N.G))
        # Calculate the mean of each row in what would be the full distance
        # matrix for this site.  We calculate the total number of 0s that would
        # be found in the row, for each row, and then from this calculate the
        # number of 1s and then the distance mean.
        row.means <- -0.5 * (this.N.G - rep(n.K, times = n.K)) / this.N.G
        D[[as.character(g)]] <- -2 * row.means + mean.d
    }
    if (length(D) == 1 && drop)
        D[[1]]
    else D
}



# If there is no diversity within a group (ss.g == 0), assign a
# minimum diversity.  The minimum distance prior to dividing by
# (n.g - 1) is 1/(2 * n.g * n.g), so replace 0 diversity with
# half this quantity, divided by (n.g - 1):
#
#        1 / ((4 * n.g * n.g) * (n.g - 1))
#
.diversityTest.ZeroVarAdjust <- function(ss.g, n.g)
{
    nn <- names(ss.g[ss.g == 0 & ! is.na(ss.g)]) # names of ss.g==0 elements
    if (length(nn)) # index by names
        ss.g[nn] <- 1 / (4 * n.g[nn] * n.g[nn] * (n.g[nn] - 1))
    ss.g
}



# Calculate terms of the variance, log-likelihood and degrees of freedom
.diversityTest.CalcTerms <- function(n.g, g.vardist, zero.div.adjust = TRUE)
{
    N <- sum(n.g)
    G <- length(n.g)
    V.g <- sapply(g.vardist, sum) / (n.g - 1)
    if (zero.div.adjust)
        V.g <- .diversityTest.ZeroVarAdjust(V.g, n.g)
    # ss.pooled
    V.p <- sum((n.g - 1) * V.g) / (N - G)
    term.V.g <- sum((n.g - 1) * log(V.g))
    term.V.p <- (N - G) * log(V.p)
    term.denom <- 1 + ((1 / (3 * (G - 1))) *
                       (sum(1 / (n.g - 1)) - (1 / (N - G))))
    ln.LR <- (term.V.p - term.V.g) / term.denom
    DF <- G - 1
    list(V.g = V.g, V.p = V.p, ln.LR = ln.LR, DF = DF)
}



# Construct null distribution of the variance
.diversityTest.NullDist <- function(obs, n.g, g.vardist,
    zero.div.adjust = TRUE, method = c("bootstrap", "permute"),
    n.resample = 10000)
{
    method <- match.arg(method)
    N <- sum(n.g)
    G <- length(n.g)
    cum.n.g <- cumsum(n.g)
    all.g.vardist <- unlist(g.vardist, use.names=FALSE)
    nulldist <- obs  # observed
    for (i in 2:n.resample) {
        p <- switch(method,
                    "permute" = sample(all.g.vardist),
                    "bootstrap" = sample(all.g.vardist, replace=TRUE))
        for (n in names(g.vardist)) {
            # peel off a slice of the distance permutation for each group
            slice <- (cum.n.g[n] - n.g[n] + 1):cum.n.g[n]
            g.vardist[[n]] <- p[slice]
        }
        terms <- .diversityTest.CalcTerms(n.g, g.vardist, zero.div.adjust)
        nulldist <- c(nulldist, terms$ln.LR)
    }
    sort(nulldist)
}



# After attempting to use these functions to contrast diversities of OTUs,
# where counts/site were on the order of 10^5s, it was found to be very slow
# and require a lot of memory for the distance matrix intermediate, as well as
# the Gower intermediate.  The first thought was to cut out the generation of
# the entire Gower matrix and just generate the diagonal.  This didn't remove
# the need for the full distance matrix, there was still a space issue, but it
# definitely sped up the calculations.
#
# Addressing the space issue, we can avoid creating the distance matrix
# completely.  We still need to calculate a per-site quantity, but never create
# a complete row, instead we count the 1s and 0s in the row and matrix, which
# is all that is required to calculate the means.
#
# The code below is legacy, containing both the original way of doing things,
# the diagonal shortcut, and some wrappers used while doing a bit of
# benchmarking.
#
# In the tests, the original method was
#
#   g.distmat <- .diversityTest.distmat(tab)
#   g.vardist <- lapply(g.distmat, function(x) diag(.diversityTest.gower(x)))
#
# Then the diag trick gave us
#
#   g.distmat <- .diversityTest.distmat(tab)
#   g.vardist <- lapply(g.distmat, .diversityTest.gowerDiag(x))
#
# Then the final optimisation saving lots of time and space became
#
#   g.vardist <- .diversityTest.directGowerDiag(tab)
#
#
# Benchmarking these against each other for a small example dataset:
#
# m = matrix(c(2,1,0,0,0,0,0,1,2,3,0,0,0,1,0,1,7,9),3,6,TRUE)
# dimnames(m) = list(site=c("s1","s2","s3"), source=c("A","B","C","D","E","F"))
#
# alt.gower1 is the original method, alt.gower2 is the diagonal trick, and
# alt.gowerDiag is the new optimisation.
#
# > microbenchmark(alt.gower1(m), alt.gower2(m), alt.gowerDiag(m), times = 100000)
# Unit: microseconds
#              expr     min       lq     mean   median       uq      max neval
#     alt.gower1(m) 464.887 496.0135 538.8139 510.5865 530.2535 148723.3 1e+05
#     alt.gower2(m) 245.627 264.0980 288.4119 271.9580 283.1540 157597.7 1e+05
#  alt.gowerDiag(m)  85.680  95.5160 105.2933  98.9100 103.8760 132955.0 1e+05
#
# alt.gower1 <- function(tab)
# {
#     g.distmat <- .diversityTest.distmat(tab)
#     lapply(g.distmat, function(x) diag(.diversityTest.gower(x)))
# }
#
# alt.gower2 <- function(tab)
# {
#     g.distmat <- .diversityTest.distmat(tab)
#     lapply(g.distmat, .diversityTest.gowerDiag)
# }
#
# For space, with a total of N items belonging to K sites and G groups, the
# first method was something like O(2*N^2), while the diagonal method was
# O(N+N^2) and the new method is O(K).

# Construct a 0-1 distance matrix from the site x group matrix, each
# entry is 1 of the group is identical and 0 if it is not.
#
.diversityTest.distmat <- function(tab, group = dimnames(tab)[[1]],
                                   drop = TRUE)
{
    if (is.null(dim(tab))) {
        dim(tab) <- c(1, length(tab))
        dimnames(tab) <- list(Site = "onedim", Genotype = names(tab))
    }
    if (dim(tab)[1] > 1 && is.null(group))
        stop("must supply group(s), all groups not supported")
    else if (missing(group) && dim(tab)[1] == 1)
        group <- 1
    G <- dim(tab)[1]
    K <- dim(tab)[2]
    N <- sum(tab)
    N.G <- apply(tab, 1, sum)
    D <- list()
    for (g in group) {
        Dmat <- matrix(1, N.G[g], N.G[g])
        n.K <- tab[g, ][tab[g, ] > 0]
        cum.n.K <- cumsum(n.K)
        for (src in 1:length(n.K)) {
            # which rows/cols to 0
            slice <- (cum.n.K[src] - n.K[src] + 1):cum.n.K[src]
            Dmat[slice, slice] <- 0
        }
        D[[as.character(g)]] <- Dmat
    }
    if (length(D) == 1 && drop)
        D[[1]]
    else D
}



# Create centroid distances for variances based on Gower (1966)
#
# Gower JC. 1966. Some distance properties of latent root and vector
# methods used in multivariate analysis.  Biometrika 53:325-338.
#
.diversityTest.gower <- function(dmat)
{
    if (is.null(dim(dmat))) {
        dim(dmat) <- c(1, length(dmat))
        dimnames(dmat) <- list(Site = "onedim", Genotype = names(dmat))
    }
    if (! all(dmat == t(dmat)))
        stop("dmat not symmetric")
    d <- -0.5 * dmat
    rd <- apply(d, 1, mean)
    if (! all(rd == apply(d, 2, mean)))
        stop("dmat not symmetric")
    gower.mat <- d + outer(-rd, -rd, "+") + mean(d)
    if (! all(abs(rowSums(gower.mat)) <= .diversityTest.epsilon))
        stop("abs(rowSums(gower.mat)) > .diversityTest.epsilon")
    if (! all(abs(colSums(gower.mat)) <= .diversityTest.epsilon))
        stop("abs(colSums(gower.mat)) > .diversityTest.epsilon")
    gower.mat
}



# Return diagonal of matrix of centroid distances for variances based on Gower (1966)
#
# Gower JC. 1966. Some distance properties of latent root and vector
# methods used in multivariate analysis.  Biometrika 53:325-338.
#
.diversityTest.gowerDiag <- function(dmat)
{
    if (is.null(dim(dmat))) {
        dim(dmat) <- c(1, length(dmat))
        dimnames(dmat) <- list(Site = "onedim", Group = names(dmat))
    }
    if (! all(dmat == t(dmat)))
        stop("dmat not symmetric")
    d <- -0.5 * dmat
    rd <- rowMeans(d)
    # diag(d) is always 0
    # diag(d) + (-rd + -rd) + mean(d)
    (-rd + -rd) + mean(d)
}



# Stop if any row contains just one item
#
.checkRowSums <- function(tab)
{
    if (any(rowSums(tab) == 1)) {
        stop("row(s) ", paste(names(which(rowSums(tab) == 1))),
             " contain 1 or fewer items and cannot be included in this analysis",
             call. = FALSE)
    }
}

