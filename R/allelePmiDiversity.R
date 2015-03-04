#' @include pmiDiversity.R
# for collation
NULL

# Using PMI statistics (Grivet et al. 2005, Scofield et al.  2010, 2011) as
# well as alpha, beta, gamma diversity (Scofield et al 2012) calculate allele
# diversity following Sork et al.
#
# The workflow to calculate basic allelic diversity statistics:
#
#   source("allelePmiDiversity.R")
#   dat = readGenalex("GenAlEx-format-file-of-genotypes.txt")
#   gt = allele.createTableList(dat)
#   div = allele.pmiDiversity(gt)
#
# For comparing allele diversity between two different samples:
#
#   source("allelePmiDiversity.R")
#   source("alleleDiversityTests.R")
#   dat1 = readGenalex("file-of-genotypes-sample-1.txt")
#   dat2 = readGenalex("file-of-genotypes-sample-2.txt")
#   gt1 = allele.createTableList(dat1)
#   gt2 = allele.createTableList(dat2)
#   alpha.contrast = allele.alphaContrastTest(gt1, gt2)
#   gamma.contrast = allele.gammaContrastTest(gt1, gt2)
#
#
# FUNCTIONS
#
# allele.pmiDiversity() : 
#    The function calculating diversity for a set of loci.  The single argument
#    is a list produced by allele.createTableList(), and it uses the function
#    allele.pmiDiversitySingleLocus().
#
# allele.createTableList() :
#   Take a data.frame of genotypes read by readGenalex(), produce a list of
#   allele count tables used by the other functions.  Each entry of the list is,
#   for each locus, a table of site x allele counts, with row names being the
#   site names, and column names being the names given to the individual
#   alleles.
#
# allele.pmiDiversitySingleLocus() : 
#    The single argument is, for a single locus, a table of site x allele
#    counts, with row names being the site names, and column names being the
#    names given to the individual alleles.
#
# CHANGELOG
#
# 0.1   First version, for Sork et al.



#' Convert class \code{'genalex'} object to a list of allele count tables
#'
#' S3 method to convert an object of class \code{'genalex'} to a list of 
#' allele count tables.  This is a generic so that other methods might be
#' written to convert other genetic formats.
#'
#' @export allele.CreateTableList
#'
allele.createTableList <- function(dat, new.ploidy = 2, 
                                   collapse.alleles = TRUE,
                                   exclude = c(NA, "0"), quiet = FALSE) {
    if (! is.genalex(dat))
        stop("input must be class 'genalex'")
    if (new.ploidy >= 2 && ! collapse.alleles)
        stop("Must collapse ploidy when new.ploidy >= 2")
    if (attr(dat, "ploidy") > new.ploidy)
        dat = reduceGenalexPloidy(dat, new.ploidy)
    lc <- attr(dat, "locus.columns")
    ln <- attr(dat, "locus.names")
    pop <- attr(dat, "pop.title")
    ans <- list()
    ex <- list()
    for (il in 1:length(lc)) {
        v <- as.vector(unlist(dat[, lc[il]:(lc[il] + new.ploidy - 1)]))
        ex[[ ln[il] ]] <- sum(v %in% exclude)
        p <- rep(dat[[pop]], new.ploidy)
        ans[[ ln[il] ]] <- t(as.matrix(table(v, p, exclude=exclude)))
    }
    if (sum(unlist(ex)) && !quiet)
        cat(sprintf("Excluding %d entries based on 'exclude=c(%s)'\n", 
                    sum(unlist(ex)), paste(collapse=", ", exclude)))
    return(ans)
}



#' allele.pmiDiversity
#'
#' @export allele.pmiDiversity
#'
allele.pmiDiversity = function(lst) {
  # calculates diversity statistics for a collection of loci, the argument
  # is produced by allele.createTableList()

  # create utility zero vectors and matrices
  pop.locus.0 = matrix(0, nrow=nrow(lst[[1]]), ncol=length(names(lst)), dimnames=list(Pop=rownames(lst[[1]]), Locus=names(lst)))
  pop.pop.0 = matrix(0, nrow=nrow(lst[[1]]), ncol=nrow(lst[[1]]), dimnames=list(Pop=rownames(lst[[1]]), Pop=rownames(lst[[1]])))
  pop.0 = pop.locus.0[, 1]
  locus.0 = pop.locus.0[1, ]
  # accumulator vectors and matrices
  alpha.gk = r.gg = pop.locus.0
  r.gh = pop.pop.0
  overlap.a = divergence.a = R.0.a = r.bar.weighted.a = alpha.bar.weighted.a = locus.0
  r.k = alpha.k = pop.0
  p = list()
  # Go through loci, storing diversity values
  for (a in names(lst)) {
    p[[a]] = allele.pmiDiversitySingleLocus(lst[[a]])
    alpha.gk[, a] = p[[a]]$allele.alpha.g
    r.gg[, a] = p[[a]]$allele.r.gg
    R.0.a[a] = p[[a]]$allele.R.0
    r.bar.weighted.a[a] = p[[a]]$allele.r.bar.weighted
    alpha.bar.weighted.a[a] = p[[a]]$allele.alpha.bar.weighted
    r.gh = r.gh + p[[a]]$allele.r.gh
    overlap.a[a] = p[[a]]$allele.overlap
    divergence.a[a] = p[[a]]$allele.divergence
  }
  # summary variables
  G = length(locus.0)
  N = p[[1]]$allele.N  # number of alleles, so for diploid this will already be 2*N
  K = p[[1]]$allele.G
  # derived variables
  alpha.k[] = apply(alpha.gk, 1, mean)
  r.k[] = apply(r.gg, 1, mean)
  r.bar.weighted = mean(r.bar.weighted.a)
  alpha.bar.weighted = 1 / r.bar.weighted
  R.0 = mean(R.0.a)
  gamma = 1 / R.0
  beta = gamma / alpha.bar.weighted
  # derived overlaps
  r.gh = r.gh / G  # mean of r.gh elements across loci, note both r_gh and r_hg so no '2*' required, verify symmetrical
  overlap = sum(r.gh) / ((K - 1) * sum(r.k))
  divergence = 1 - overlap
  overlap.mean = mean(overlap.a)
  divergence.mean = mean(divergence.a)
  # scaled values, note that N here is already 2*N when diploid
  alpha.max = (N - G) / G
  gamma.max = (N - 1)
  beta.max = G
  alpha.scaled = ((alpha.bar.weighted - 1) / alpha.bar.weighted) * (alpha.max / (alpha.max - 1))
  gamma.scaled = ((gamma - 1) / gamma) * (gamma.max / (gamma.max - 1))
  beta.scaled = ((beta - 1) / beta) * (beta.max / (beta.max - 1))
  # return value
  list(produced.by        = p[[1]]$produced.by,
       r.gg               = r.gg,
       alpha.gk           = alpha.gk,
       r.k                = r.k,
       alpha.k            = alpha.k,
       r.bar.weighted.a   = r.bar.weighted.a,
       R.0.a              = R.0.a,
       r.bar.weighted     = r.bar.weighted,
       alpha.bar.weighted = alpha.bar.weighted,
       gamma              = gamma,
       beta               = beta,
       G                  = G,
       N                  = N,
       K                  = K,
       alpha.scaled       = alpha.scaled,
       gamma.scaled       = gamma.scaled,
       beta.scaled        = beta.scaled,
       r.gh               = r.gh,
       overlap            = overlap,
       divergence         = divergence,
       overlap.mean       = overlap.mean,
       divergence.mean    = divergence.mean)
}


#' Calculate allelic diversity for a single locus
#'
#' @export allele.pmiDiversitySingleLocus
#'
allele.pmiDiversitySingleLocus <- function(tab)
{
  # pmiDiversity() calculates several quantities we use here
  p = pmiDiversity(tab)

  allele.N = p$num.samples
  allele.G = p$num.groups
  allele.r.gg = p$r.gg
  allele.n.g = p$num.samples.group
  allele.n.k = p$num.samples.source
  allele.alpha.g = 1 / allele.r.gg
  allele.r.bar.weighted = sum((allele.n.g - 1) * allele.r.gg) / (allele.N - allele.G)
  allele.alpha.bar.weighted = 1 / allele.r.bar.weighted
  allele.R.0 = sum(allele.n.k * (allele.n.k - 1)) / (allele.N * (allele.N - 1))
  allele.gamma = 1 / allele.R.0
  allele.beta = allele.gamma / allele.alpha.bar.weighted

  # overlap and divergence; recreate C matrix not returned from pmiDiversity()
  allele.r.gh = p$Q.mat; diag(allele.r.gh) = 0 # C matrix here
  allele.overlap = p$r.overlap
  allele.divergence = p$r.divergence

  # functions to scale alpha, beta, gamma, see Sork et al
  # the use of body() replaces N and G in the parse tree with current values
  allele.func.scale.alpha = function(a, ploidy=2) {
    N = allele.N
    G = allele.G
    alpha.max = (ploidy * N - G) / G
    ((a - 1) / a) * (alpha.max / (alpha.max - 1))
  }
  body(allele.func.scale.alpha)[[2]][[3]] = allele.N
  body(allele.func.scale.alpha)[[3]][[3]] = allele.G
  allele.func.scale.gamma = function(g, ploidy=2) {
    N = allele.N
    gamma.max = N - 1
    ((g - 1) / g) * (N / (N - 1))
  }
  body(allele.func.scale.gamma)[[2]][[3]] = allele.N
  allele.func.scale.beta = function(b, ploidy=2) {
    G = allele.G
    beta.max = G
    ((b - 1) / b) * (beta.max / (beta.max - 1))
  }
  body(allele.func.scale.beta)[[2]][[3]] = allele.G
  # scaled to [0, 1] assuming ploidy of 1 and 2
  allele.scaled.ploidy1 = 
      list(allele.alpha.g            = allele.func.scale.alpha(allele.alpha.g, 1),
           allele.alpha.bar.weighted = allele.func.scale.alpha(allele.alpha.bar.weighted, 1),
           allele.gamma              = allele.func.scale.gamma(allele.gamma, 1),
           allele.beta               = allele.func.scale.beta(allele.beta, 1))
  allele.scaled.ploidy2 = 
      list(allele.alpha.g            = allele.func.scale.alpha(allele.alpha.g, 2),
           allele.alpha.bar.weighted = allele.func.scale.alpha(allele.alpha.bar.weighted, 2),
           allele.gamma              = allele.func.scale.gamma(allele.gamma, 2),
           allele.beta               = allele.func.scale.beta(allele.beta, 2))

  list(allele.table              = tab,
       allele.N                  = allele.N,
       allele.G                  = allele.G,
       allele.r.gg               = allele.r.gg,
       allele.n.g                = allele.n.g,
       allele.n.k                = allele.n.k,
       allele.alpha.g            = allele.alpha.g,
       allele.r.bar.weighted     = allele.r.bar.weighted,
       allele.alpha.bar.weighted = allele.alpha.bar.weighted,
       allele.R.0                = allele.R.0,
       allele.r.gh               = allele.r.gh,
       allele.overlap            = allele.overlap,
       allele.divergence         = allele.divergence,
       allele.gamma              = allele.gamma,
       allele.beta               = allele.beta,
       allele.func.scale.alpha   = allele.func.scale.alpha,
       allele.func.scale.gamma   = allele.func.scale.gamma,
       allele.func.scale.beta    = allele.func.scale.beta,
       allele.scaled.ploidy1     = allele.scaled.ploidy1,
       allele.scaled.ploidy2     = allele.scaled.ploidy2)
}

