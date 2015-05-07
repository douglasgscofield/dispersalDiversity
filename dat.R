library(Matrix)

load.dat = function(file = "dat.dat") read.delim(file)

convert.dat = function(d) {
    r = setNames(d$Reads, d$OTU)
    l = split(r, d$Site_Species)
    sstab = do.call("rbind", l)
    names(dimnames(sstab)) = c("Site_Species", "OTU")
    sstab
}


# Maximum tolerance for row and column sums for Gower matrix
alt.diversityTest.epsilon <- 1.0e-12



# Still not quite sure...
alt.diversityTest.ReverseTerms <- TRUE

#' @rdname alphaDiversityTest
#'
#' @export
#'
alt.alphaDiversityTest.divtable <- function(tab, zero.var.adjust = TRUE,
    n.resample = 10000, method = c("bootstrap", "permute"),
    test.quantiles = c(0.001, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 0.999),
    ...)
{
    method <- match.arg(method)
    ans <- list(method = "Alpha diversity test, contrast among sites in single data set")
    ans$data.name <- deparse(substitute(tab))
    g.distmat <- alt.diversityTest.distmat(tab)
    #g.vardist <- lapply(g.distmat, function(x) diag(alt.diversityTest.gower(x)))
    g.vardist <- lapply(g.distmat, alt.diversityTest.gowerDiag)
    n.g <- sapply(g.vardist, length)
    N <- sum(n.g)
    G <- length(n.g)
    ans$N.samples <- N
    ans$N.groups <- G
    terms = alt.diversityTest.CalcTerms(n.g, g.vardist, zero.var.adjust)
    ans$observed.ln.LR <- terms$ln.LR
    # Analytic distribution
    ans$df.X2 <- terms$DF
    ans$P.analytic <- pchisq(terms$ln.LR, df = terms$DF, TRUE)
    # Empirical distribution
    nulldist <- alt.diversityTest.NullDist(obs = terms$ln.LR, n.g = n.g, 
        g.vardist = g.vardist, zero.var.adjust = zero.var.adjust,
        method = method, n.resample=n.resample)
    ans$n.resample <- n.resample
    ans$resample.method <- method
    ans$quantiles <- quantile(nulldist, test.quantiles)
    ans$P.empirical <- sum(terms$ln.LR <= nulldist) / n.resample
    ans$empdist <- nulldist
    structure(ans, class = c('diversity_test', 'list'))
}


# Construct a 0-1 distance matrix from the site x group matrix, each
# entry is 1 of the group is identical and 0 if it is not.  This is
# represented efficiently.
#
alt.diversityTest.distmat <- function(tab, group = dimnames(tab)[[1]],
                                   drop = TRUE)
{
    if (is.null(dim(tab))) {
        dim(tab) <- c(1, length(tab))
        dimnames(tab) <- list(Site = "onedim", Group = names(tab))
    }
    if (dim(tab)[1] > 1 && is.null(group))
        stop("must supply group(s), all groups not supported")
    else if (missing(group) && dim(tab)[1] == 1)
        group <- 1
    G <- dim(tab)[1]
    K <- dim(tab)[2]
    N <- sum(tab)
    N.G <- rowSums(tab)
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



# Return diagonal of matrix of centroid distances for variances based on Gower (1966)
#
# Gower JC. 1966. Some distance properties of latent root and vector
# methods used in multivariate analysis.  Biometrika 53:325-338.
#
alt.diversityTest.gowerDiag <- function(dmat)
{
    if (is.null(dim(dmat))) {
        dim(dmat) <- c(1, length(dmat))
        dimnames(dmat) <- list(Site = "onedim", Group = names(dmat))
    }
    if (! all(dmat == t(dmat)))
        stop("dmat not symmetric")
    d <- -0.5 * dmat
    rd <- rowMeans(d)
    diag(d) + (-rd + -rd) + mean(d)
}



# Create centroid distances for variances based on Gower (1966)
#
# Gower JC. 1966. Some distance properties of latent root and vector
# methods used in multivariate analysis.  Biometrika 53:325-338.
#
alt.diversityTest.gower <- function(dmat)
{
    if (is.null(dim(dmat))) {
        dim(dmat) <- c(1, length(dmat))
        dimnames(dmat) <- list(Site = "onedim", Group = names(dmat))
    }
    if (! all(dmat == t(dmat)))
        stop("dmat not symmetric")
    d <- -0.5 * dmat
    rd <- rowMeans(d)
    gower.mat <- d + outer(-rd, -rd, "+") + mean(d)
    if (! all(abs(rowSums(gower.mat)) <= alt.diversityTest.epsilon))
        stop("abs(rowSums(gower.mat)) > alt.diversityTest.epsilon")
    if (! all(abs(colSums(gower.mat)) <= alt.diversityTest.epsilon))
        stop("abs(colSums(gower.mat)) > alt.diversityTest.epsilon")
    gower.mat
}


# Calculate terms of the variance, log-likelihood and degrees of freedom
alt.diversityTest.CalcTerms <- function(n.g, g.vardist, zero.var.adjust = TRUE)
{
    N <- sum(n.g)
    G <- length(n.g)
    V.g <- sapply(g.vardist, sum) / (n.g - 1)
    if (zero.var.adjust)
        V.g <- alt.diversityTest.ZeroVarAdjust(V.g, n.g)
    # ss.pooled
    V.p <- sum((n.g - 1) * V.g) / (N - G)
    term.V.g <- sum((n.g - 1) * log(V.g))
    term.V.p <- (N - G) * log(V.p)
    term.denom <- 1 + ((1 / (3 * (G - 1))) *
                       (sum(1 / (n.g - 1)) - (1 / (N - G))))
    ln.LR <- if (alt.diversityTest.ReverseTerms)
        (term.V.p - term.V.g) / term.denom
    else
        (term.V.g - term.V.p) / term.denom
    DF <- G - 1
    list(V.g = V.g, V.p = V.p, ln.LR = ln.LR, DF = DF)
}



# Construct null distribution of the variance
alt.diversityTest.NullDist <- function(obs, n.g, g.vardist,
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
        terms <- alt.diversityTest.CalcTerms(n.g, g.vardist, zero.var.adjust)
        nulldist <- c(nulldist, terms$ln.LR)
    }
    sort(nulldist)
}

# If there is no diversity within a group (ss.g == 0), assign a
# minimum diversity.  The minimum distance prior to dividing by
# (n.g - 1) is 1/(2 * n.g * n.g), so replace 0 diversity with
# half this quantity, divided by (n.g - 1):
#
#        1 / ((4 * n.g * n.g) * (n.g - 1))
#
alt.diversityTest.ZeroVarAdjust <- function(ss.g, n.g)
{
    nn <- names(ss.g[ss.g == 0]) # names of ss.g==0 elements
    if (length(nn)) # index by names
        ss.g[nn] <- 1 / (4 * n.g[nn] * n.g[nn] * (n.g[nn] - 1))
    ss.g
}


