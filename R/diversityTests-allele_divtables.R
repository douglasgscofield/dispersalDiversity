#' @include diversity-allele_divtables.R
#' @include diversityTests-divtable.R
NULL



#' Print the result of a allele diversity test (alpha or gamma)
#'
#' Print an object of class \code{allele_diversity_test}, the result of
#' \code{\link{alphaDiversityTest}}, \code{\link{alphaContrastTest}},
#' \code{\link{alphaContrastTest3}}, \code{\link{gammaContrastTest}} and
#' \code{\link{gammaContrastTest3}} when given an argument of class
#' \code{\link{allele_divtables}}.
#'
#' The result of the test printed is the summed effect across all loci;
#' see Sork \emph{et al}. (2015) for details.  If \code{include.loci} is
#' \code{TRUE}, then this is followed by the test result for each locus
#' individually.
#'
#' @param x   Object of class \code{allele_diversity_test}, returned by
#' \code{\link{alphaDiversityTest}}, \code{\link{alphaContrastTest}},
#' \code{\link{alphaContrastTest3}}, \code{\link{gammaContrastTest}} and
#' \code{\link{gammaContrastTest3}} when given an argument of class
#' \code{\link{allele_divtables}}.
#'
#' @param digits  Number of significant digits to use when printing
#' numeric values.  Defaults to the value of the \code{"digits"} option.
#'
#' @param include.loci  If \code{TRUE}, print test results for each
#' individual locus.  Defaults to \code{FALSE}.
#'
#' @return \code{x}, returned invisibly
#'
#' @seealso \code{\link{alphaDiversityTest}}, \code{\link{alphaContrastTest}}, \code{\link{alphaContrastTest3}}, \code{\link{gammaContrastTest}}, \code{\link{gammaContrastTest3}}
#'
#' @references
#'
#' Sork, V. L., Smouse, P. E., Grivet, D. and Scofield, D. G. (In press)
#' Impact of asymmetric male and female gamete dispersal on allelic 
#' diversity and spatial genetic structure in valley oak 
#' (\emph{Quercus lobata} N\'{e}e).  \emph{Evolutionary Ecology}.
#'
#' @examples
#'
#' ## library(readGenalex)  # already loaded as a prerequisite
#' data(Qagr_pericarp_genotypes)  # from readGenalex
#' gt <- createAlleleTables(Qagr_pericarp_genotypes)
#' alpha.test <- alphaDiversityTest(gt)
#' print(alpha.test)
#'
#' @export
#'
print.allele_diversity_test <- function(x, digits = getOption("digits"),
    include.loci = FALSE, ...)
{
    cat(x$method, "\n\n", sep = "")
    cat("data:  ", x$data.name, "\n", sep = "")
    cat("loci:  ", paste(sort(names(x$loci))), "\n", sep = "")
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
    if (include.loci) {
        cat("\n\nTests for individual loci\n\n")
        for (l in sort(names(x$loci)))
            print.diversity_test(x$loci[[l]], ...)
    }
    invisible(x)
}



# Already documented in diversityTests.divtable
#' @rdname alphaDiversityTest
#'
#' @export
#'
alphaDiversityTest.allele_divtables <- function(adt, zero.div.adjust = TRUE, 
    n.resample = 10000, method = c("bootstrap", "permute"),
    test.quantiles = c(0.001, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 0.999),
    ...)
{
    method <- match.arg(method)
    name <- deparse(substitute(adt))
    nm <- names(adt)

    ans <- list(method = "Alpha diversity test, contrast among sites for multiple loci in single data set",
                data.name = name)

    # for the grand test
    ans$observed.ln.LR <- 0
    ans$empdist <- numeric(n.resample)

    # sublists for each locus
    ans$loci <- list()

    for (l in 1:length(adt)) {

        locus <- nm[l]

        this.ans <- list(method = paste("Alpha diversity test, contrast among sites for locus",
                                        locus),
                         data.name = paste(sep = "_", name, locus),
                         locus = locus)

        tab <- adt[[locus]]
        .checkRowSums(tab)
        g.vardist <- .diversityTest.directGowerDiag(tab)
        n.g <- sapply(g.vardist, length)
        N <- sum(n.g)
        G <- length(n.g)
        this.ans$N.samples <- N
        this.ans$N.groups <- G
        terms <- .diversityTest.CalcTerms(n.g, g.vardist, zero.div.adjust)
        this.ans$observed.ln.LR <- terms$ln.LR
        nulldist <- .diversityTest.NullDist(obs = terms$ln.LR, 
                                            n.g = n.g, g.vardist = g.vardist, 
                                            zero.div.adjust = zero.div.adjust, 
                                            method = method, 
                                            n.resample = n.resample)
        this.ans$n.resample <- n.resample
        this.ans$resample.method <- method
        this.ans$quantiles <- quantile(nulldist, test.quantiles)
        this.ans$P <- sum(terms$ln.LR <= nulldist) / n.resample
        this.ans$empdist <- nulldist

        # add the current locus' results to the grand totals
        ans$loci[[ locus ]] <- structure(this.ans, class = c('diversity_test', 'list'))
        ans$observed.ln.LR <- ans$observed.ln.LR + this.ans$observed.ln.LR
        ans$empdist <- ans$empdist + this.ans$empdist
    }

    ans$N.samples <- ans$loci[[ nm[1] ]]$N.samples
    ans$N.groups <- ans$loci[[ nm[1] ]]$N.groups
    ans$n.resample <- n.resample
    ans$resample.method <- method
    ans$quantiles <- quantile(ans$empdist, test.quantiles)
    ans$P <- sum(ans$observed.ln.LR <= ans$empdist) / ans$n.resample
    ans <- structure(ans, class = c('allele_diversity_test', 'list'))
    invisible(ans)
}



# Already documented under alleleContrastTest.divtable
#' @rdname alphaContrastTest
#'
#' @export
#'
alphaContrastTest.allele_divtables <- function(adt.a, adt.b,
    zero.div.adjust = TRUE, n.resample = 10000, 
    method = c("bootstrap", "permute"),
    test.quantiles = c(0.001, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 0.999))
{

    name.a <- deparse(substitute(adt.a))
    name.b <- deparse(substitute(adt.b))
    if (! inherits(adt.b, 'allele_divtables'))
        stop("both ", name.a, " and ", name.b, " must be of class 'allele_divtables'")
    method <- match.arg(method)
    stopifnot(length(adt.a) == length(adt.b) && all(names(adt.a) == names(adt.b)))
    nm <- names(adt.a)

    ans <- list(method = "Alpha diversity test, contrast between two multilocus data sets",
                data.name = paste(sep = ",", name.b, name.b))
    ans$name.a <- name.a
    ans$name.b <- name.b

    # for the grand test
    ans$observed.ln.LR <- 0
    ans$empdist <- numeric(n.resample)

    # sublists for each locus
    ans$loci <- list()

    for (l in 1:length(adt.a)) {

        locus <- nm[l]

        this.ans <- list(method = paste("Alpha diversity test, contrast between 2 datasets for locus", locus),
                         data.name = paste(sep = "_", paste(sep = ",", name.a, name.b), locus),
                         locus = locus)

        tab.a <- adt.a[[locus]]
        .checkRowSums(tab.a)
        a.vardist <- .diversityTest.directGowerDiag(tab.a)
        n.a <- sapply(a.vardist, length)
        N.a <- sum(n.a)
        G.a <- length(n.a)
        this.ans$N.a <- N.a
        this.ans$G.a <- G.a
        terms.a <- .diversityTest.CalcTerms(n.a, a.vardist, zero.div.adjust)

        tab.b <- adt.b[[locus]]
        .checkRowSums(tab.b)
        b.vardist <- .diversityTest.directGowerDiag(tab.b)
        n.b <- sapply(b.vardist, length)
        N.b <- sum(n.b)
        G.b <- length(n.b)
        this.ans$N.b <- N.b
        this.ans$G.b <- G.b
        terms.b <- .diversityTest.CalcTerms(n.b, b.vardist, zero.div.adjust)
        V.a.b.p <- (((N.a - G.a) * terms.a$V.p) + 
                    ((N.b - G.b) * terms.b$V.p)) / (N.a + N.b - G.a - G.b)
        observed.ln.LR.a.b <- ((N.a + N.b - G.a - G.b) * log(V.a.b.p)) - 
                              ((N.a - G.a) * log(terms.a$V.p)) - 
                              ((N.b - G.b) * log(terms.b$V.p))

        # Combine A and B into strata for comparison
        n.a.b <- c(a = N.a, b = N.b)
        a.b.vardist <- list(a = unlist(a.vardist, use.names = FALSE), 
                            b = unlist(b.vardist, use.names = FALSE))
        N <- sum(n.a.b)
        G <- length(n.a.b)
        this.ans$N.samples <- N
        this.ans$N.groups <- G
        this.ans$observed.ln.LR <- observed.ln.LR.a.b
        nulldist <- .diversityTest.NullDist(obs = observed.ln.LR.a.b,
                                            n.g = n.a.b,
                                            g.vardist = a.b.vardist,
                                            zero.div.adjust, method, n.resample)
        this.ans$N.a <- N.a
        this.ans$N.b <- N.b
        this.ans$a.b.vardist <- a.b.vardist
        this.ans$n.resample <- n.resample
        this.ans$resample.method <- method
        this.ans$quantiles <- quantile(nulldist, test.quantiles)
        this.ans$P <- sum(observed.ln.LR.a.b <= nulldist) / n.resample
        this.ans$empdist <- nulldist

        # add the current locus' results to the grand totals
        ans$loci[[ locus ]] <- structure(this.ans, class = c('diversity_test', 'list'))
        ans$observed.ln.LR <- ans$observed.ln.LR + this.ans$observed.ln.LR
        ans$empdist <- ans$empdist + this.ans$empdist
    }

    ans$N.samples <- ans$loci[[ nm[1] ]]$N.samples
    ans$N.groups <- ans$loci[[ nm[1] ]]$N.groups
    ans$N.a <- ans$loci[[ nm[1] ]]$N.a
    ans$N.b <- ans$loci[[ nm[1] ]]$N.b
    ans$n.resample <- n.resample
    ans$resample.method <- method
    ans$quantiles <- quantile(ans$empdist, test.quantiles)
    ans$P <- sum(ans$observed.ln.LR <= ans$empdist) / ans$n.resample
    ans <- structure(ans, class = c('allele_diversity_test', 'list'))
    invisible(ans)
}




# Documentation in gammaContrastTest.divtable
#
#' @rdname gammaContrastTest
#'
#' @export
#'
gammaContrastTest.allele_divtables <- function(adt.a, adt.b,
    zero.div.adjust = TRUE, n.resample = 10000, 
    method = c("bootstrap", "permute"),
    test.quantiles = c(0.001, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 0.999))
{
    name.a <- deparse(substitute(adt.a))
    name.b <- deparse(substitute(adt.b))
    if (! inherits(adt.b, 'allele_divtables'))
        stop("both ", name.a, " and ", name.b, " must be of class 'allele_divtables'")
    method <- match.arg(method)
    stopifnot(length(adt.a) == length(adt.b) && 
              all(names(adt.a) == names(adt.b)))
    nm <- names(adt.a)

    ans <- list(method = "Gamma diversity test, contrast between two multilocus data sets",
                data.name = paste(sep = ",", name.b, name.b))
    ans$name.a <- name.a
    ans$name.b <- name.b

    # for the grand test
    ans$observed.ln.LR <- 0
    ans$empdist <- numeric(n.resample)

    # sublists for each locus
    ans$loci <- list()

    for (l in 1:length(adt.a)) {

        locus <- nm[l]

        this.ans <- list()
        this.ans <- list(method = "Gamma diversity test, contrast between 2 single-locus data sets",
                         data.name = paste(sep = "_", paste(sep = ",", name.a, name.b), locus),
                         locus = locus)

        tab.a <- adt.a[[ nm[l] ]]
        .checkRowSums(tab.a)
        tab.b <- adt.b[[ nm[l] ]]
        .checkRowSums(tab.b)

        X.a.k <- apply(tab.a, 2, sum)
        X.b.k <- apply(tab.b, 2, sum)
        N.a <- sum(X.a.k)
        N.b <- sum(X.b.k)
        R.a.0 <- sum((X.a.k * (X.a.k - 1)) / (N.a * (N.a - 1)))
        R.b.0 <- sum((X.b.k * (X.b.k - 1)) / (N.b * (N.b - 1)))
        V.a.tot <- 1 - R.a.0
        V.b.tot <- 1 - R.b.0
        V.a.b.tot <- ((N.a - 1) * V.a.tot + (N.b - 1) * V.b.tot) / 
                     (N.a + N.b - 2)
        observed.ln.LR.a.b <- ((N.a + N.b - 2) * log(V.a.b.tot)) - 
                              ((N.a - 1) * log(V.a.tot)) - 
                              ((N.b - 1) * log(V.b.tot)) 

        #a.distmat <- .diversityTest.distmat(X.a.k)
        #a.vardist <- list(a = diag(.diversityTest.gower(a.distmat)))
        a.vardist <- .diversityTest.directGowerDiag(X.a.k)
        n.a <- lapply(a.vardist, length)
        stopifnot(sum(n.a) == N.a)
        N.a <- sum(n.a)
        G.a <- length(n.a)
        this.ans$N.a <- N.a
        this.ans$G.a <- G.a
        terms.a <- .diversityTest.CalcTerms(n.a, a.vardist, zero.div.adjust)
        #cat(sprintf("%s  terms.a$V.p = %f  V.a.tot = %f\n", locus, terms.a$V.p, V.a.tot))

        #b.distmat <- .diversityTest.distmat(X.b.k)
        #b.vardist <- list(b = diag(.diversityTest.gower(b.distmat)))
        b.vardist <- .diversityTest.directGowerDiag(X.b.k)
        n.b <- lapply(b.vardist, length)
        stopifnot(sum(n.b) == N.b)
        N.b <- sum(n.b)
        G.b <- length(n.b)
        this.ans$N.b <- N.b
        this.ans$G.b <- G.b
        terms.b <- .diversityTest.CalcTerms(n.b, b.vardist, zero.div.adjust)
        #cat(sprintf("%s  terms.b$V.p = %f  V.b.tot = %f\n", locus, terms.b$V.p, V.b.tot))

        # Combine A and B into stratta for comparison
        n.a.b <- c(a = N.a, b = N.b)
        a.b.vardist <- list(a = unlist(a.vardist, use.names = FALSE), 
                            b = unlist(b.vardist, use.names = FALSE))
        N <- sum(n.a.b)
        G <- length(n.a.b)
        this.ans$N.samples <- N
        this.ans$N.groups <- G
        this.ans$observed.ln.LR <- observed.ln.LR.a.b
        nulldist <- .diversityTest.NullDist(obs = observed.ln.LR.a.b,
                                            n.g = n.a.b,
                                            g.vardist = a.b.vardist,
                                            zero.div.adjust, method, n.resample)
        this.ans$a.b.vardist <- a.b.vardist
        this.ans$n.resample <- n.resample
        this.ans$resample.method <- method
        this.ans$quantiles <- quantile(nulldist, test.quantiles)
        this.ans$P <- sum(observed.ln.LR.a.b <= nulldist) / n.resample
        this.ans$empdist <- nulldist

        # add the current locus' results to the grand totals
        ans$loci[[ locus ]] <- structure(this.ans, class = c('diversity_test', 'list'))
        ans$observed.ln.LR <- ans$observed.ln.LR + this.ans$observed.ln.LR
        ans$empdist <- ans$empdist + this.ans$empdist
    }

    ans$N.samples <- ans$loci[[ nm[1] ]]$N.samples
    ans$N.groups <- ans$loci[[ nm[1] ]]$N.groups
    ans$N.a <- ans$loci[[ nm[1] ]]$N.a
    ans$G.a <- ans$loci[[ nm[1] ]]$G.a
    ans$N.b <- ans$loci[[ nm[1] ]]$N.b
    ans$G.b <- ans$loci[[ nm[1] ]]$G.b
    ans$n.resample <- n.resample
    ans$resample.method <- method
    ans$quantiles <- quantile(ans$empdist, test.quantiles)
    ans$P <- sum(ans$observed.ln.LR <= ans$empdist) / ans$n.resample
    ans <- structure(ans, class = c('allele_diversity_test', 'list'))
    invisible(ans)
}


