#' @include diversity.R
# for collation
NULL

# Using PMI statistics (Grivet et al. 2005, Scofield et al.  2010, 2011) as
# well as alpha, beta, gamma diversity (Scofield et al 2012) calculate allele
# diversity following Sork et al.
#
# @examples
#
# # The workflow to calculate basic allelic diversity statistics:
#
#   dat <- readGenalex("GenAlEx-format-file-of-genotypes.txt")
#   gt <- createAlleleTables(dat)
#   div <- diversityMultilocus(gt)
#
# # For comparing allele diversity between two different samples:
#
#   dat1 <- readGenalex("file-of-genotypes-sample-1.txt")
#   dat2 <- readGenalex("file-of-genotypes-sample-2.txt")
#   gt1 <- createAlleleTables(dat1)
#   gt2 <- createAlleleTables(dat2)
#   alpha.contrast <- allele.alphaContrastTest(gt1, gt2)
#   gamma.contrast <- allele.gammaContrastTest(gt1, gt2)
#
#
# FUNCTIONS
#
# allele.diversity() : 
#    The function calculating diversity for a set of loci.  The single argument
#    is a list produced by createAlleleTables(), and it uses the function
#    allele.diversitySingleLocus().
#
# createAlleleTables() :
#   Take a data.frame of genotypes read by readGenalex(), produce a list of
#   allele count tables used by the other functions.  Each entry of the list is,
#   for each locus, a table of site x allele counts, with row names being the
#   site names, and column names being the names given to the individual
#   alleles.
#
# allele.diversitySingleLocus() : 
#    The single argument is, for a single locus, a table of site x allele
#    counts, with row names being the site names, and column names being the
#    names given to the individual alleles.



#' Convert class \code{'genalex'} object to a list of allele count tables
#'
#' S3 method to convert an object of class \code{'genalex'} to a list of 
#' allele count tables.  This is a generic so that other methods might be
#' written to convert other genetic formats.
#'
#' @examples
#'
#' # The workflow to calculate basic allelic diversity statistics:
#' #
#' # dat <- readGenalex("GenAlEx-format-file-of-genotypes.txt")
#' # gt <- createAlleleTables(dat)
#'
#' @export
#'
#' @name createAlleleTables
#'
NULL

createAlleleTables <- function(x, ...) UseMethod("createAlleleTables")



# TODO: NOTE THAT p AND v are now the name of the dimensions
#
# as.allele_table.genalex?  as.table_list.genalex?  S3 method
#
# remove new.ploidy argument, see orig code below
#
#createAlleleTables.genalex <- function(dat, new.ploidy,
#                                   collapse.alleles = TRUE,
#                                   exclude = c(NA, "0"), quiet = FALSE) {
#    if (! is.genalex(dat))
#        stop("input must be class 'genalex'")
#    if (new.ploidy >= 2 && ! collapse.alleles)
#        stop("Must collapse ploidy when new.ploidy >= 2")
#    if (attr(dat, "ploidy") > new.ploidy)
#        dat <- reduceGenalexPloidy(dat, new.ploidy)

#' @rdname createAlleleTables
#'
#' @export
#'
createAlleleTables.genalex <- function(dat, new.ploidy = 2,
                                       collapse.alleles = TRUE,
                                       exclude = c(NA, "0"), 
                                       quiet = FALSE)
{
    if (new.ploidy >= 2 && ! collapse.alleles)
        stop("Must collapse ploidy when new.ploidy >= 2")
    if (attr(dat, "ploidy") > new.ploidy)
        dat <- reducePloidy(dat, new.ploidy)
    lc <- attr(dat, "locus.columns")
    ln <- attr(dat, "locus.names")
    pop <- attr(dat, "pop.title")
    ans <- list()
    ex <- list()
    for (il in 1:length(lc)) {
        v <- as.vector(unlist(dat[, lc[il]:(lc[il] + new.ploidy - 1)]))
        ex[[ ln[il] ]] <- sum(v %in% exclude)
        p <- rep(dat[[pop]], new.ploidy)
        ans[[ ln[il] ]] <- t(as.matrix(table(v, p, exclude = exclude)))
    }
    if (sum(unlist(ex)) && !quiet)
        cat(sprintf("Excluding %d entries based on 'exclude = c(%s)'\n", 
                    sum(unlist(ex)), paste(collapse = ", ", exclude)))
    class(ans) <- c('allele_tables', 'list')
    return(ans)
}



#' Calculate diversity across a set of allele tables
#'
#' @examples
#'
#' # The workflow to calculate basic allelic diversity statistics:
#' #
#' # dat <- readGenalex("GenAlEx-format-file-of-genotypes.txt")
#' # gt <- createAlleleTables(dat)
#' # div <- diversityMultilocus(gt)
#'
#' @export
#'
#' @name diversityMultilocus
#'
NULL

diversityMultilocus <- function(x, ...) UseMethod("diversityMultilocus")



# TODO: now with method=, the returned values are not correctly named
#
#' @rdname diversityMultilocus
#'
#' @export
#'
diversityMultilocus.allele_tables <- function(lst, method = c("r", "q.nielsen", "q"))
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
    alpha.max <- (N - G) / G
    gamma.max <- (N - 1)
    beta.max <- G
    alpha.scaled <- ((alpha.bar.weighted - 1) / alpha.bar.weighted) * 
                    (alpha.max / (alpha.max - 1))
    gamma.scaled <- ((gamma - 1) / gamma) * (gamma.max / (gamma.max - 1))
    beta.scaled <- ((beta - 1) / beta) * (beta.max / (beta.max - 1))
    # return value
    list(q.gg               = q.gg,
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
         divergence.mean    = divergence.mean
         )
}


#' Calculate allelic diversity for a single locus
#'
#' @param  tab Table of site x allele counts, with row names being the site
#' names and column names being the names given to the individual alleles
#'
#' @return List of diversity estimates for the locus
#'
#' @export diversitySingleLocus
#'
diversitySingleLocus <- function(tab, method = c("r", "q.nielsen", "q"))
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
        alpha.max <- (ploidy * N - G) / G
        ((a - 1) / a) * (alpha.max / (alpha.max - 1))
    }
    body(allele.func.scale.alpha)[[2]][[3]] <- allele.N
    body(allele.func.scale.alpha)[[3]][[3]] <- allele.G
    allele.func.scale.gamma <- function(g, ploidy=2) {
        N <- allele.N
        gamma.max <- N - 1
        ((g - 1) / g) * (N / (N - 1))
    }
    body(allele.func.scale.gamma)[[2]][[3]] <- allele.N
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

    list(allele.table              = tab,
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
         allele.scaled.ploidy2     = allele.scaled.ploidy2
         )
}

