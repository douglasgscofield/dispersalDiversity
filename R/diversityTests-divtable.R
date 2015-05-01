#' @include diversity-divtable.R
# for collation
NULL



#' Print the result of a diversity test (alpha or gamma)
#'
#' Print an object of class \code{diversity_test}, the result of
#' \code{\link{alphaDiversityTest}}, \code{\link{alphaContrastTest}},
#' \code{\link{alphaContrastTest3}}, \code{\link{gammaContrastTest}} and
#' \code{\link{gammaContrastTest3}}.
#'
#' @param x       Object of class \code{diversity_test}, returned by
#' \code{\link{alphaDiversityTest}}, \code{\link{alphaContrastTest}},
#' \code{\link{alphaContrastTest3}}, \code{\link{gammaContrastTest}} and
#' \code{\link{gammaContrastTest3}}
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
    #cat("Samples N.a = ", x$N.a, ", N.b = ", x$N.b,
    #    ", groups = ", x$N.groups, "\n")
    cat("Observed log-likelihood ratio = ",
        format(signif(x$observed.ln.lR, max(1L, digits - 2L))), "\n\n",
        sep = "")
    cat("Test against analytic X-2 distribution (usually not appropriate):\n")
    cat("Degrees of freedom = (N.groups - 1) = ", x$df.X2,
        ", P = ", format.pval(x$P.analytic, digits = max(1L, digits - 3L)),
        "\n", sep = "")
    cat("Test against empirical X-2 distribution (much better test):\n")
    cat("Iterations = ", x$n.resample,
        ", P = ", format.pval(x$P.empirical, digits = max(1L, digits - 3L)),
        "\n", sep = "")
    cat("Quantiles of the empirical distribution:\n")
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
#' @param add.analytic  Logical, if \code{TRUE} a curve of the \eqn{\Chi^2}
#' distribution with the appropriate degrees of freedom is added to the plot
#' @param breaks  Number of breaks to use when plotting the histogram of
#' the empirical distribution, passed to \code{\link{hist(..., plot = FALSE)}}
#' @param compress.x  Logical, if \code{TRUE} and the observed value is more
#' than \code{compress.ratio} times the maximum value of the empirical
#' distribution, the x-axis is compressed to include the observed value
#' @param compress.ration See \code{compress.x}
#' @param xlab,ylab,main  Labels added to the plot
#' @param legend.text Text to use when printing legend containing the
#' observed value in the upper right of the plot. Set to \code{NULL} to
#' suppress printing the legend.
#' @param digits  Number of significant digits to use when printing numeric
#' values
#'
#' @seealso \code{\link{alphaDiversityTest}}, \code{\link{alphaContrastTest}}, \code{\link{alphaContrastTest3}}, \code{\link{gammaContrastTest}}, \code{\link{gammaContrastTest3}}
#'
#' @export
#'
plot.diversity_test <- function(x, add.analytic = FALSE, 
    breaks = 50, compress.x = TRUE, compress.ratio = 1.2,
    xlab = "ln(LR) value", ylab = "Frequency", main = x$method, 
    legend.text = "Observed ln(LR) = ", digits = getOption("digits"), ...)
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
    if (add.analytic) {
        analytic.x <- seq(xlim[1], xlim[2], length.out = 200)
        analytic.X2 <- dchisq(x = analytic.x, df = x$df.ln.LR)
        ylim <- range(c(0, h$count, analytic.X2))
    } else {
        ylim <- range(c(0, h$count))
    }
    plot(h, xlim = xlim, ylim = ylim, freq = TRUE, 
         xlab = ylab, ylab = ylab, main = main, ...)
    lines(rep(plotted.observed.ln.LR, 2), c(ylim[1], ylim[2]*0.5),
          col = "darkgray", lty = 1, lwd = 2)
    if (add.analytic)
        lines(analytic.x, analytic.X2, lty=3, lwd=1)
    if (! is.null(legend.text))
        legend("topright", bty = "n", legend = paste0(legend.text,
               format(signif(x$observed.ln.LR, max(1L, digits - 2L)))))
}



#' Test for differences in alpha diversity among sites within a single data set
#'
#' @param tab    Site-by-source table of class \code{\link{divtable}}
#' @param zero.var.adjust Logical, whether to adjust zero-variance groups
#' with minimum value, see \code{\link{BLAHBLAH}}
#' @param n.resample Number of iterations for creation of the null distribution
#' @param method \code{"bootstrap"} or \code{"permute"}, whether to create null
#' distribution iterations with (\code{"bootstrap"}) or without
#' (\code{"permute"}) replacement
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
# @examples
#
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
alphaDiversityTest.divtable <- function(tab, zero.var.adjust = TRUE,
    n.resample = 10000, method = c("bootstrap", "permute"),
    test.quantiles = c(0.001, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 0.999),
    ...)
{
    method <- match.arg(method)
    ans <- list(method = "Alpha diversity test, contrast among sites in single data set")
    ans$data.name <- deparse(substitute(tab))
    g.distmat <- .diversityTest.distmat(tab)
    g.vardist <- lapply(g.distmat, function(x) diag(.diversityTest.gower(x)))
    n.g <- unlist(lapply(g.vardist, length))
    N <- sum(n.g)
    G <- length(n.g)
    ans$N.samples <- N
    ans$N.groups <- G
    terms = .diversityTest.CalcTerms(n.g, g.vardist, zero.var.adjust)
    ans$observed.ln.LR <- terms$ln.LR
    # Analytic distribution
    ans$df.X2 <- terms$DF
    ans$P.analytic <- pchisq(terms$ln.LR, df = terms$DF, TRUE)
    # Empirical distribution
    nulldist <- .diversityTest.NullDist(obs = terms$ln.LR, n.g = n.g, 
        g.vardist = g.vardist, zero.var.adjust = zero.var.adjust,
        method = method, n.resample=n.resample)
    ans$n.resample <- n.resample
    ans$resample.method <- method
    ans$quantiles <- quantile(nulldist, test.quantiles)
    ans$P.empirical <- sum(terms$ln.LR <= nulldist) / n.resample
    ans$empdist <- nulldist
    structure(ans, class = c('diversity_test', 'list'))
}



#' Test for differences in alpha diversity between two data sets
#'
#' @param tab.a First ite x soure table
#' @param tab.b Second site x soure table
#' @param zero.var.adjust Logical, whether to adjust zero-variance groups
#' with minimum value, see \code{\link{BLAHBLAH}}
#' @param n.resample Number of iterations for creation of the null distribution
#' @param method \code{"bootstrap"} or \code{"permute"}, whether to create null
#' distribution iterations with (\code{"bootstrap"}) or without
#' (\code{"permute"}) replacement
#'
#' @return An \code{'htest_boot'} object with the result of the test. The
#' test also prints the result, so perhaps I should modify it to only
#' return the htest_boot object?  Or should it return a different object?
#'
#' @seealso \code{\link{alphaDiversityTest}}, \code{\link{alphaContrastTest3}}
#'
#' @references
#'
#' Scofield, D. G., Smouse, P. E., Karubian, J. and Sork, V. L. (2012)
#' Use of alpha, beta and gamma diversity measures to characterize seed
#' dispersal by animals.  \emph{American Naturalist} 180:719-732.
#'
# @examples
#'
#' @export
#'
#' @name alphaContrastTest
#'
NULL

alphaContrastTest <- function(tab.a, tab.b, ...) UseMethod("alphaContrastTest")



#' @rdname alphaContrastTest
#'
#' @export
#"
alphaContrastTest.divtable <- function(tab.a, tab.b, zero.var.adjust = TRUE,
    n.resample = 10000, method = c("bootstrap", "permute"),
    test.quantiles = c(0.001, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 0.999))
{
    stopifnot(inherits(tab.b, 'divtable'))
    method <- match.arg(method)
    #.RT = .diversityTest.ReverseTerms
    #.diversityTest.ReverseTerms = FALSE
    ans <- list(method = "Alpha diversity test, contrast between 2 datasets")
    ans$data.name <- paste(sep = ", ", deparse(substitute(tab.a)),
        deparse(substitute(tab.b)))
    a.distmat <- .diversityTest.distmat(tab.a)
    a.vardist <- lapply(a.distmat, function(x) diag(.diversityTest.gower(x)))
    n.a <- sapply(a.vardist, length)
    N.a <- sum(n.a)
    G.a <- length(n.a)
    ans$N.a <- N.a
    ans$G.a <- G.a
    terms.a = .diversityTest.CalcTerms(n.a, a.vardist, zero.var.adjust)
    # V.a.p = terms.a$V.p

    b.distmat <- .diversityTest.distmat(tab.b)
    b.vardist <- lapply(b.distmat, function(x) diag(.diversityTest.gower(x)))
    n.b <- sapply(b.vardist, length)
    N.b <- sum(n.b)
    G.b <- length(n.b)
    ans$N.b <- N.b
    ans$G.b <- G.b
    terms.b = .diversityTest.CalcTerms(n.b, b.vardist, zero.var.adjust)
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
    # Analytic distribution
    ans$df.X2 <- G - 1
    ans$P.analytic <- pchisq(observed.ln.LR.a.b, df = DF, lower.tail = TRUE)
    # Empirical distribution
    nulldist <- .diversityTest.NullDist(obs=observed.ln.LR.a.b,
                                        n.g = n.a.b,
                                        g.vardist = a.b.vardist,
                                        zero.var.adjust, method, 
                                        n.resample)
    ans$n.resample <- n.resample
    ans$resample.method <- method
    ans$quantiles <- quantile(nulldist, test.quantiles)
    ans$P.empirical <- sum(observed.ln.LR.a.b <= nulldist) / n.resample
    ans$empdist <- nulldist
    ans$a.b.vardist = a.b.vardist
    structure(ans, class = c('diversity_test', 'list'))
}



#' Test for differences in alpha diversity between three data sets
#'
#' @param tab.a First site-by-source table
#' @param tab.b Second site-by-source table
#' @param tab.c Third site-by-source table
#' @param zero.var.adjust Logical, whether to adjust zero-variance groups
#' with minimum value, see \code{\link{BLAHBLAH}}
#' @param n.resample Number of iterations for creation of the null distribution
#' @param method \code{"bootstrap"} or \code{"permute"}, whether to create null
#' distribution iterations with (\code{"bootstrap"}) or without
#' (\code{"permute"}) replacement
#'
#' @return An \code{'htest_boot'} object with the result of the test. The
#' test also prints the result, so perhaps I should modify it to only
#' return the htest_boot object?  Or should it return a different object?
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
#' @export
#'
#' @name alphaContrastTest3
#'
NULL

alphaContrastTest3 <- function(tab.a, tab.b, tab.c, ...)
    UseMethod("alphaContrastTest3")



#' @rdname alphaContrastTest3
#'
#' @export
#'
alphaContrastTest3.divtable <- function(tab.a, tab.b, tab.c,
    zero.var.adjust = TRUE, n.resample = 10000,
    method = c("bootstrap", "permute"),
    test.quantiles = c(0.001, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 0.999),
    ...)
{
    stopifnot(inherits(tab.b, 'divtable'))
    stopifnot(inherits(tab.c, 'divtable'))
    method <- match.arg(method)
    #.RT <- .diversityTest.ReverseTerms
    #.diversityTest.ReverseTerms <- FALSE
    ans <- list(method = "Alpha diversity test, contrast between 3 datasets")
    ans$data.name <- paste(sep = ", ", deparse(substitute(tab.a)),
        deparse(substitute(tab.b)), deparse(substitute(tab.b)))
    a.distmat <- .diversityTest.distmat(tab.a)
    a.vardist <- lapply(a.distmat, function(x) diag(.diversityTest.gower(x)))
    n.a <- sapply(a.vardist, length)
    N.a <- sum(n.a)
    G.a <- length(n.a)
    ans$N.a <- N.a
    ans$G.a <- G.a
    terms.a <- .diversityTest.CalcTerms(n.a, a.vardist, zero.var.adjust)
    # V.a.p <- terms.a$V.p

    b.distmat <- .diversityTest.distmat(tab.b)
    b.vardist <- lapply(b.distmat, function(x) diag(.diversityTest.gower(x)))
    n.b <- unlist(lapply(b.vardist, length))
    N.b <- sum(n.b)
    G.b <- length(n.b)
    ans$N.b <- N.b
    ans$G.b <- G.b
    terms.b <- .diversityTest.CalcTerms(n.b, b.vardist, zero.var.adjust)
    # V.b.p <- terms.b$V.p

    c.distmat <- .diversityTest.distmat(tab.c)
    c.vardist <- lapply(c.distmat, function(x) diag(.diversityTest.gower(x)))
    n.c <- unlist(lapply(c.vardist, length))
    N.c <- sum(n.c)
    G.c <- length(n.c)
    ans$N.c <- N.c
    ans$G.c <- G.c
    terms.c <- .diversityTest.CalcTerms(n.c, c.vardist, zero.var.adjust)
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
    # Analytic distribution
    ans$df.X2 <- G - 1
    ans$P.analytic <- pchisq(observed.ln.LR.a.b.c, df = DF, lower.tail = TRUE)
    # Empirical distribution
    nulldist <- .diversityTest.NullDist(obs = observed.ln.LR.a.b.c,
                                        n.g = n.a.b.c,
                                        g.vardist = a.b.c.vardist,
                                        zero.var.adjust, method, 
                                        n.resample)
    ans$n.resample <- n.resample
    ans$resample.method <- method
    ans$quantiles <- quantile(nulldist, test.quantiles)
    ans$P.empirical <- sum(observed.ln.LR.a.b.c <= nulldist) / n.resample
    ans$empdist <- nulldist
    ans$a.b.c.vardist <- a.b.c.vardist
    structure(ans, class = c('diversity_test', 'list'))
}




#' Test for structure in pairwise divergence between sites
#'
#' @export
#'
#' @name pairwiseMeanTest
#'
NULL

pairwiseMeanTest <- function(tab, ...) UseMethod("pairwiseMeanTest")



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



#' Plot the result of \code{'pairwiseMeanTest'}
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
#' @param tab.a First site-by-source table
#' @param tab.b Second site-by-source table
#' @param zero.var.adjust Logical, whether to adjust zero-variance groups
#' with minimum value, see \code{\link{BLAHBLAH}}
#' @param n.resample Number of iterations for creation of the null distribution
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
# @examples
#
#' @export
#'
#' @name gammaContrastTest
#'
NULL

gammaContrastTest <- function(tab.a, tab.b, ...) UseMethod("gammaContrastTest")



#' @rdname gammaContrastTest
#'
#' @export
#"
gammaContrastTest.divtable <- function(tab.a, tab.b, zero.var.adjust = TRUE,
    n.resample = 10000, method = c("bootstrap", "permute"),
    test.quantiles = c(0.001, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 0.999),
    ...)
{
    stopifnot(inherits(tab.b, 'divtable'))
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
    a.vardist <- list(a = diag(.diversityTest.gower(.diversityTest.distmat(X.a.k))))
    n.a <- sapply(a.vardist, length)
    stopifnot(sum(n.a) == N.a)
    G.a <- length(n.a)
    ans$N.a <- N.a
    ans$G.a <- G.a
    terms.a = .diversityTest.CalcTerms(n.a, a.vardist, zero.var.adjust)
    #cat(sprintf("terms.a$V.p = %f  V.a.tot = %f\n", terms.a$V.p, V.a.tot))

    b.vardist <- list(b = diag(.diversityTest.gower(.diversityTest.distmat(X.b.k))))
    n.b <- unlist(lapply(b.vardist, length))
    stopifnot(sum(n.b) == N.b)
    G.b <- length(n.b)
    ans$N.b <- N.b
    ans$G.b <- G.b
    terms.b = .diversityTest.CalcTerms(n.b, b.vardist, zero.var.adjust)
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
    # Analytic distribution
    ans$df.X2 <- DF
    ans$P.analytic <- pchisq(observed.ln.LR.a.b, df = DF, lower.tail = TRUE)
    # Empirical distribution
    nulldist <- .diversityTest.NullDist(obs = observed.ln.LR.a.b,
                                        n.g = n.a.b, g.vardist = a.b.vardist,
                                        zero.var.adjust, method, n.resample)
    ans$n.resample <- n.resample
    ans$resample.method <- method
    ans$quantiles <- quantile(nulldist, test.quantiles)
    ans$P.empirical <- sum(observed.ln.LR.a.b <= nulldist) / n.resample
    ans$empdist <- nulldist
    ans$a.b.vardist = a.b.vardist
    structure(ans, class = c('diversity_test', 'list'))
}



#' Test for difference in gamma diversity among three datasets
#'
#' @param tab.a First site-by-source table
#' @param tab.b Second site-by-source table
#' @param tab.c Third site-by-source table
#' @param zero.var.adjust Logical, whether to adjust zero-variance groups
#' with minimum value, see \code{\link{BLAHBLAH}}
#' @param n.resample Number of iterations for creation of the null distribution
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
    zero.var.adjust = TRUE, n.resample = 10000,
    method = c("bootstrap", "permute"),
    test.quantiles = c(0.001, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 0.999),
    ...)
{
    stopifnot(inherits(tab.b, 'divtable'))
    stopifnot(inherits(tab.c, 'divtable'))
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
    a.distmat <- .diversityTest.distmat(X.a.k)
    a.vardist <- list(a = diag(.diversityTest.gower(a.distmat)))
    n.a <- sapply(a.vardist, length)
    stopifnot(sum(n.a) == N.a)
    N.a <- sum(n.a)
    G.a <- length(n.a)
    ans$N.a <- N.a
    ans$G.a <- G.a
    terms.a <- .diversityTest.CalcTerms(n.a, a.vardist, zero.var.adjust)
    #cat(sprintf("terms.a$V.p = %f  V.a.tot = %f\n", terms.a$V.p, V.a.tot))

    b.distmat <- .diversityTest.distmat(X.b.k)
    b.vardist <- list(b = diag(.diversityTest.gower(a.distmat)))
    n.b <- sapply(b.vardist, length)
    stopifnot(sum(n.b) == N.b)
    N.b <- sum(n.b)
    G.b <- length(n.b)
    ans$N.b <- N.b
    ans$G.b <- G.b
    terms.b <- .diversityTest.CalcTerms(n.b, b.vardist, zero.var.adjust)
    #cat(sprintf("terms.b$V.p = %f  V.b.tot = %f\n", terms.b$V.p, V.b.tot))

    c.distmat <- .diversityTest.distmat(X.c.k)
    c.vardist <- list(c = diag(.diversityTest.gower(a.distmat)))
    n.c <- sapply(c.vardist, length)
    stopifnot(sum(n.c) == N.c)
    N.c <- sum(n.c)
    G.c <- length(n.c)
    ans$N.c <- N.c
    ans$G.c <- G.c
    terms.c <- .diversityTest.CalcTerms(n.c, c.vardist, zero.var.adjust)
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
    # Analytic distribution
    ans$df.X2 <- DF
    ans$P.analytic <- pchisq(observed.ln.LR.a.b.c, df = DF, lower.tail = TRUE)
    # Empirical distribution
    nulldist <- .diversityTest.NullDist(obs = observed.ln.LR.a.b.c,
                                        n.g = n.a.b.c,
                                        g.vardist = a.b.c.vardist,
                                        zero.var.adjust, method, n.resample)
    ans$n.resample <- n.resample
    ans$resample.method <- method
    ans$quantiles <- quantile(nulldist, test.quantiles)
    ans$P.empirical <- sum(observed.ln.LR.a.b.c <= nulldist) / n.resample
    ans$empdist <- nulldist
    ans$a.b.c.vardist <- a.b.c.vardist
    structure(ans, class = c('diversity_test', 'list'))
}



# ------- Internal functions and constants



# Maximum tolerance for row and column sums for Gower matrix
.diversityTest.epsilon <- 1.0e-12



# Still not quite sure...
.diversityTest.ReverseTerms <- TRUE



# construct a distance matrix from the site x group matrix
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



# Create centroid distances for variances based on
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



# If there is no diversity within a group (ss.g == 0), assign a
# minimum diversity.  The minimum distance prior to dividing by
# (n.g - 1) is 1/(2 * n.g * n.g), so replace 0 diversity with
# half this quantity, divided by (n.g - 1):
#
#        1 / ((4 * n.g * n.g) * (n.g - 1))
#
.diversityTest.ZeroVarAdjust <- function(ss.g, n.g)
{
    nn <- names(ss.g[ss.g == 0]) # names of ss.g==0 elements
    if (length(nn)) # index by names
        ss.g[nn] <- 1 / (4 * n.g[nn] * n.g[nn] * (n.g[nn] - 1))
    ss.g
}



# Calculate terms of the variance, log-likelihood and degrees of freedom
.diversityTest.CalcTerms <- function(n.g, g.vardist, zero.var.adjust = TRUE)
{
    N <- sum(n.g)
    G <- length(n.g)
    V.g <- unlist(lapply(g.vardist, sum)) / (n.g - 1)
    if (zero.var.adjust)
        V.g <- .diversityTest.ZeroVarAdjust(V.g, n.g)
    # ss.pooled
    V.p <- sum((n.g - 1) * V.g) / (N - G)
    term.V.g <- sum((n.g - 1) * log(V.g))
    term.V.p <- (N - G) * log(V.p)
    term.denom <- 1 + ((1 / (3 * (G - 1))) *
                       (sum(1 / (n.g - 1)) - (1 / (N - G))))
    ln.LR <- if (.diversityTest.ReverseTerms)
        (term.V.p - term.V.g) / term.denom
    else
        (term.V.g - term.V.p) / term.denom
    DF <- G - 1
    list(V.g = V.g, V.p = V.p, ln.LR = ln.LR, DF = DF)
}



# Construct null distribution of the variance
.diversityTest.NullDist <- function(obs, n.g, g.vardist,
    zero.var.adjust = TRUE, method = c("bootstrap", "permute"),
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
        terms <- .diversityTest.CalcTerms(n.g, g.vardist, zero.var.adjust)
        nulldist <- c(nulldist, terms$ln.LR)
    }
    sort(nulldist)
}

