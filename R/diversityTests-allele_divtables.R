#' @include diversity-allele_divtables.R
#' @include diversityTests-divtable.R
# For collation, load after the above
NULL



# Already documented in diversityTests.divtable
#
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

    ans <- list()
    ans$name <- name
    ans$names <- nm
    # for the grand test
    ans$observed.ln.LR <- 0
    ans$empdist <- numeric(n.resample)
    # sublists for each locus

    for (l in 1:length(adt)) {

        locus <- nm[l]

        this.ans <- list()
        this.ans$locus <- locus

        tab <- adt[[locus]]

        #g.vardist <- lapply(.diversityTest.distmat(tab), 
        #                    function(x) diag(.diversityTest.gower(x)))  
        g.vardist <- .diversityTest.directGowerDiag(tab)
        n.g <- sapply(g.vardist, length)
        N <- sum(n.g)
        G <- length(n.g)
        this.ans$name <- name
        this.ans$N.samples <- N
        this.ans$N.groups <- G

        terms <- .diversityTest.CalcTerms(n.g, g.vardist, zero.div.adjust)
        PVAL <- pchisq(terms$ln.LR, df = terms$DF, lower.tail = FALSE)
        cat("Bartlett's Test for Heteroscedasticity in Intra-group Variances for", name, "\n")
        cat("Current locus:", locus, "\n")
        cat("For comparison against analytic X^2, see data structure\n")
        #cat("Comparing against X-2 distribution directly:\n")
        #cat(sprintf("N samples = %d  groups = %d  observed.ln.LR = %f  df=(G-1) = %d  P = %g\n", 
        #            N, G, terms$ln.LR, terms$DF, PVAL))
        this.ans$observed.ln.LR <- terms$ln.LR
        this.ans$df.X2 <- terms$DF
        this.ans$P.analytic <- PVAL
        nulldist <- .diversityTest.NullDist(obs = terms$ln.LR, 
                                            n.g = n.g, g.vardist = g.vardist, 
                                            zero.div.adjust = zero.div.adjust, 
                                            method = method, 
                                            n.resample = n.resample)
        PVAL <- sum(terms$ln.LR <= nulldist)/n.resample
        cat("Comparing against bootstrap X-2 distribution:\n")
        q2 <- quantile(nulldist, test.quantiles)
        o <- sprintf("N samples = %d  groups = %d  observed.ln.LR = %f  resamples = %d  boot.X2[%s] = [%s]  P = %g\n", 
            N, G, terms$ln.LR, n.resample, 
            paste(collapse = " ", test.quantiles), 
            paste(collapse = " ", round(q2, 3)),
            PVAL)
        cat(o)
        this.ans$n.resample <- n.resample
        this.ans$resample.method <- method
        this.ans$quants <- q2
        this.ans$P.empirical <- PVAL
        this.ans$empdist <- nulldist

        # add the current locus' results to the grand totals
        ans[[ locus ]] <- this.ans
        ans$observed.ln.LR <- ans$observed.ln.LR + this.ans$observed.ln.LR
        ans$empdist <- ans$empdist + this.ans$empdist
    }

    ans$N.samples <- ans[[ nm[1] ]]$N.samples
    ans$N.groups <- ans[[ nm[1] ]]$N.groups
    ans$n.resample <- n.resample
    ans$resample.method <- method
    cat("Bartlett's Test for Heteroscedasticity in Intra-group Variances for",name,"\n\n")
    cat("Over ALL loci\n")
    cat("Skipping comparison against analytic X^2\n")
    ans$P.empirical <- sum(ans$observed.ln.LR <= ans$empdist)/ans$n.resample
    cat("Contrasting groups, compare summed ln-LR against summed bootstrap X^2 distribution:\n")
    ans$quants <- quantile(ans$empdist, test.quantiles)
    cat(sprintf("N samples = %d  groups = %d  observed.ln.LR = %f  resamples = %d  boot.X2[%s] = [%s]  P = %g\n", 
              ans$N.samples,
              ans$N.groups,
              ans$observed.ln.LR,
              ans$n.resample,
              paste(collapse = " ", test.quantiles),
              paste(collapse = " ", round(ans$quants, 3)),
              ans$P.empirical))
    ####
    class(ans) <- c(class(ans), "allele.diversityTest")
    invisible(ans)
}



# Already documented under alleleContrastTest.divtable
#
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
    .RT <- .diversityTest.ReverseTerms
    .diversityTest.ReverseTerms <- FALSE

    ans <- list()
    ans$name.a <- name.a
    ans$name.b <- name.b
    ans$names <- nm
    # for the grand test
    ans$observed.ln.LR <- 0
    ans$empdist <- numeric(n.resample)
    # sublists for each locus

    for (l in 1:length(adt.a)) {

        locus <- nm[l]

        this.ans <- list()
        this.ans$locus <- locus

        tab.a <- adt.a[[locus]]
        #a.vardist <- lapply(.diversityTest.distmat(tab.a), 
        #                    function(x) diag(.diversityTest.gower(x)))  
        a.vardist <- .diversityTest.directGowerDiag(tab.a)
        n.a <- sapply(a.vardist, length)
        N.a <- sum(n.a)
        G.a <- length(n.a)
        this.ans$name.a <- name.a
        this.ans$name.b <- name.b
        this.ans$N.a <- N.a
        this.ans$G.a <- G.a
        terms.a <- .diversityTest.CalcTerms(n.a, a.vardist, zero.div.adjust)

        tab.b <- adt.b[[locus]]
        #b.vardist <- lapply(.diversityTest.distmat(tab.b), 
        #                    function(x) diag(.diversityTest.gower(x)))  
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
        DF <- G -1
        this.ans$N.samples <- N
        this.ans$N.groups <- G
        PVAL <- pchisq(observed.ln.LR.a.b, df = DF, lower.tail = TRUE)
        cat("Bartlett's Test for Heteroscedasticity in **Inter-group** Alpha Variances between",name.a,"and",name.b,"\n")
        cat("Current locus:", locus, "\n")
        cat("For comparison against analytic X^2, see data structure\n")
        #cat("Comparing against X^2 distribution directly:\n")
        #cat(sprintf("Samples N.a = %d, N.b = %d  groups = %d  observed.ln.LR.a.b = %f  df=(G-1) = %d  P = %g\n", 
        #      N.a, N.b, G, observed.ln.LR.a.b, DF, PVAL))
        this.ans$observed.ln.LR <- observed.ln.LR.a.b
        this.ans$df.X2 <- G - 1
        this.ans$P.analytic <- PVAL
        nulldist <- .diversityTest.NullDist(obs = observed.ln.LR.a.b,
                                            n.g = n.a.b,
                                            g.vardist = a.b.vardist,
                                            zero.div.adjust, method, n.resample)
        PVAL <- sum(observed.ln.LR.a.b <= nulldist)/n.resample
        cat("Contrasting groups, compare against bootstrap X^2 distribution:\n")
        q2 <- quantile(nulldist, test.quantiles)
        o <- sprintf("Samples N.a = %d, N.b = %d  groups = %d  observed.ln.LR.a.b = %f  resamples = %d  boot.X2[%s] = [%s]  P = %g\n\n", 
              N.a, N.b, G, observed.ln.LR.a.b, n.resample, 
              paste(collapse = " ", test.quantiles), 
              paste(collapse = " ", round(q2, 3)), 
              PVAL)
        cat(o)
        this.ans$n.resample <- n.resample
        this.ans$N.a <- N.a
        this.ans$N.b <- N.b
        this.ans$resample.method <- method
        this.ans$quants <- q2
        this.ans$P.empirical <- PVAL
        this.ans$empdist <- nulldist
        this.ans$a.b.vardist <- a.b.vardist

        # add the current locus' results to the grand totals
        ans[[ locus ]] <- this.ans
        ans$observed.ln.LR <- ans$observed.ln.LR + this.ans$observed.ln.LR
        ans$empdist <- ans$empdist + this.ans$empdist

    }
    ans$N.samples <- ans[[ nm[1] ]]$N.samples
    ans$N.groups <- ans[[ nm[1] ]]$N.groups
    ans$N.a <- ans[[ nm[1] ]]$N.a
    ans$N.b <- ans[[ nm[1] ]]$N.b
    ans$n.resample <- n.resample
    ans$resample.method <- method
    cat("Bartlett's Test for Heteroscedasticity in **Inter-group** Alpha Variances between",name.a,"and",name.b,"\n")
    cat("Over ALL loci\n")
    cat("Skipping comparison against analytic X^2\n")
    ans$P.empirical <- sum(ans$observed.ln.LR <= ans$empdist)/ans$n.resample
    cat("Contrasting groups, compare summed ln-LR against summed bootstrap X^2 distribution:\n")
    ans$quants <- quantile(ans$empdist, test.quantiles)
    cat(sprintf("Samples N.a = %d, N.b = %d, groups = %d, observed.ln.LR = %f, resamples = %d, boot.x2[%s] = [%s]   P = %g\n\n\n",
              ans$N.a,
              ans$N.b,
              ans$N.groups,
              ans$observed.ln.LR,
              ans$n.resample,
              paste(collapse = " ", test.quantiles),
              paste(collapse = " ", round(ans$quants, 3)), 
              ans$P.empirical))
    .diversityTest.ReverseTerms <- .RT
    ####
    class(ans) <- c(class(ans), "allele.diversityTest.AlphaContrast")
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

    ans <- list()
    ans$name.a <- name.a
    ans$name.b <- name.b
    ans$names <- nm
    ans$observed.ln.LR <- 0
    ans$empdist <- numeric(n.resample)

    for (l in 1:length(adt.a)) {
        locus <- nm[l]

        this.ans <- list()
        this.ans$locus <- locus

        tab.a <- adt.a[[ nm[l] ]]
        tab.b <- adt.b[[ nm[l] ]]

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
        cat(sprintf("%s  terms.a$V.p = %f  V.a.tot = %f\n", locus, terms.a$V.p, V.a.tot))

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
        cat(sprintf("%s  terms.b$V.p = %f  V.b.tot = %f\n", locus, terms.b$V.p, V.b.tot))

        # Combine A and B into stratta for comparison
        n.a.b <- c(a = N.a, b = N.b)
        a.b.vardist <- list(a = unlist(a.vardist, use.names = FALSE), 
                            b = unlist(b.vardist, use.names = FALSE))
        N <- sum(n.a.b)
        G <- length(n.a.b)
        DF <- G -1
        this.ans$N.samples <- N
        this.ans$N.groups <- G
        PVAL <- pchisq(observed.ln.LR.a.b, df = DF, lower.tail = TRUE)
        cat("Bartlett's Test for Heteroscedasticity in **Inter-group** Gamma Variances between", name.a, "and", name.b, "\n")
        cat("Current locus:", locus, "\n")
        cat("For comparison against analytic X^2, see data structure\n")
        #cat("Comparing against X-2 distribution directly:\n")
        #cat(sprintf("Samples N.a = %d, N.b = %d  groups = %d  observed.ln.LR.a.b = %f  df=(G-1) = %d  P = %g\n", 
        #      N.a, N.b, G, observed.ln.LR.a.b, DF, PVAL))
        this.ans$observed.ln.LR <- observed.ln.LR.a.b
        this.ans$df.X2 <- DF
        this.ans$P.analytic <- PVAL
        nulldist <- .diversityTest.NullDist(obs = observed.ln.LR.a.b,
                                            n.g = n.a.b,
                                            g.vardist = a.b.vardist,
                                            zero.div.adjust, method, n.resample)
        PVAL <- sum(observed.ln.LR.a.b <= nulldist)/n.resample
        cat("Contrasting groups, compare against bootstrap X-2 distribution:\n")
        q2 <- quantile(nulldist, test.quantiles)
        cat(sprintf("Samples N.a = %d, N.b = %d  groups = %d  observed.ln.LR.a.b = %f  resamples = %d  boot.X2[%s] = [%s]  P = %g\n\n", 
              N.a, N.b, G, observed.ln.LR.a.b, n.resample, 
              paste(collapse = " ", test.quantiles), 
              paste(collapse = " ", round(q2, 3)), 
              PVAL))
        this.ans$n.resample <- n.resample
        this.ans$resample.method <- method
        this.ans$quants <- q2
        this.ans$P.empirical <- PVAL
        this.ans$empdist <- nulldist
        this.ans$a.b.vardist <- a.b.vardist

        # add the current locus' results to the grand totals
        ans[[ locus ]] <- this.ans
        ans$observed.ln.LR <- ans$observed.ln.LR + this.ans$observed.ln.LR
        ans$empdist <- ans$empdist + this.ans$empdist

    }
    ans$N.samples <- ans[[ nm[1] ]]$N.samples
    ans$N.groups <- ans[[ nm[1] ]]$N.groups
    ans$N.a <- ans[[ nm[1] ]]$N.a
    ans$G.a <- ans[[ nm[1] ]]$G.a
    ans$N.b <- ans[[ nm[1] ]]$N.b
    ans$G.b <- ans[[ nm[1] ]]$G.b
    ans$n.resample <- n.resample
    ans$resample.method <- method
    cat("Bartlett's Test for Heteroscedasticity in **Inter-group** Gamma Variances between",name.a,"and",name.b,"\n")
    cat("Over ALL loci\n")
    cat("Skipping comparison against analytic X^2\n")
    ans$P.empirical <- sum(ans$observed.ln.LR <= ans$empdist)/ans$n.resample
    cat("Contrasting groups, compare summed ln-LR against summed bootstrap X^2 distribution:\n")
    ans$quants <- quantile(ans$empdist, test.quantiles)
    cat(sprintf("Samples N.a = %d, G.a = %d, N.b = %d, G.b = %d, groups = %d, observed.ln.LR = %f, resamples = %d, boot.x2[%s] = [%s]   P = %g\n\n\n",
              ans$N.a,
              ans$G.a,
              ans$N.b,
              ans$G.b,
              ans$N.groups,
              ans$observed.ln.LR,
              ans$n.resample,
              paste(collapse = " ", test.quantiles),
              paste(collapse = " ", round(ans$quants, 3)),
              ans$P.empirical))
    ####
    class(ans) <- c(class(ans), "allele.diversityTest.GammaContrast")
    invisible(ans)
}


