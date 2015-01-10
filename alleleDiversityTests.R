# alleleDiversityTests.R
#
# Copyright (c) 2014 Douglas G. Scofield, Evolutionary Biology Centre, Uppsala University
#
# douglas.scofield@ebc.uu.se
# douglasgscofield@gmail.com
#
# These statistical tools were developed in collaboration with Peter Smouse
# (Rutgers University) and Victoria Sork (UCLA) and were funded by U.S. National
# Science Foundation awards NSF-DEB-0514956 and NSF-DEB-0516529.
#
# Use as you see fit.  No warranty regarding this code is implied nor should be
# assumed.  Send bug reports etc. to one of the above email addresses.
#
#
#
# Provide functions for performing allele-based diversity tests, following
# Sork et al., extensions of those from Scofield et al. 2012 American
# Naturalist 180(6) 719-732, http://www.jstor.org/stable/10.1086/668202).
#
# Input for tests is a list, one entry per locus, of site x allele counts.
# Each table must have the same format.
#
# (columns) format.
#
# FUNCTIONS PROVIDED 
#
# allele.alphaContrastTest(lst.a, lst.b) : Test whether there is a difference
# in the alpha diversity between two lists of allele diversity datasets
# 
# allele.gammaContrastTest(lst.a, lst.b) : Test whether there is a difference
# in the gamma diversity between two datasets
# 
# CHANGELOG
#
# 0.1   First version, for Sork et al.

.alleleDiversityTestsVersion = 0.1

# This code relies on other files in respository
source("allelePmiDiversity.R")
source("diversityTests.R")

#---------------------------------------------


allele.alphaDiversityTest <- function(lst, 
                               zero.var.adjust=TRUE, 
                               n.resample=10000, 
                               method=c("bootstrap", "permute"))
{
  method <- match.arg(method)
  name = deparse(substitute(tab))
  nm = names(lst)
  q1 <- c(0.5, 0.95, 0.99, 0.999)

  ans <- list()
  ans$name = name
  ans$names = nm
  # for the grand test
  ans$observed.ln.LR = 0
  ans$empdist = numeric(n.resample)
  # sublists for each locus

  for (l in 1:length(lst)) {

    locus = nm[l]

    this.ans = list()
    this.ans$locus = locus

    tab = lst[[locus]]

    g.vardist <- lapply(.diversityTest.distmat(tab), 
                        function(x) diag(.diversityTest.gower(x)))  
    n.g <- unlist(lapply(g.vardist, length))
    N <- sum(n.g)
    G <- length(n.g)
    this.ans$name = name
    this.ans$N.samples <- N
    this.ans$N.groups <- G

    terms = .diversityTest.CalcTerms(n.g, g.vardist, zero.var.adjust)
    PVAL = pchisq(terms$ln.LR, df=terms$DF, lower.tail=FALSE)
    cat("Bartlett's Test for Heteroscedasticity in Intra-group Variances for", name, "\n")
    cat("Current locus:", locus, "\n")
    cat("For comparison against analytic X^2, see data structure\n")
    #cat("Comparing against X-2 distribution directly:\n")
    #cat(sprintf("N samples = %d  groups = %d  observed.ln.LR = %f  df=(G-1) = %d  P = %g\n", 
    #            N, G, terms$ln.LR, terms$DF, PVAL))
    this.ans$observed.ln.LR <- terms$ln.LR
    this.ans$df.X2 <- terms$DF
    this.ans$P.analytic <- PVAL
    nulldist <- .diversityTest.NullDist(obs=terms$ln.LR, 
                                        n.g=n.g, g.vardist=g.vardist, 
                                        zero.var.adjust=zero.var.adjust, 
                                        method=method, n.resample=n.resample)
    PVAL <- sum(terms$ln.LR <= nulldist)/n.resample
    cat("Comparing against bootstrap X-2 distribution:\n")
    q2 <- quantile(nulldist, q1)
    o = sprintf("N samples = %d  groups = %d  observed.ln.LR = %f  resamples = %d  boot.X2[%s] = [%.2f %.2f %.2f %.2f]  P = %g\n", 
        N, G, terms$ln.LR, n.resample, paste(collapse=" ", q1), 
        q2[1], q2[2], q2[3], q2[4], PVAL)
    cat(o)
    this.ans$n.resample <- n.resample
    this.ans$resample.method <- method
    this.ans$quants <- q2
    this.ans$P.empirical <- PVAL
    this.ans$empdist <- nulldist
    this.ans$version <- .diversityTestsVersion

    # add the current locus' results to the grand totals
    ans[[ locus ]] = this.ans
    ans$observed.ln.LR = ans$observed.ln.LR + this.ans$observed.ln.LR
    ans$empdist = ans$empdist + this.ans$empdist
  }

  ans$N.samples = ans[[ nm[1] ]]$N.samples
  ans$N.groups = ans[[ nm[1] ]]$N.groups
  ans$n.resample <- n.resample
  ans$resample.method <- method
  cat("Bartlett's Test for Heteroscedasticity in Intra-group Variances for",name,"\n\n")
  cat("Over ALL loci\n")
  cat("Skipping comparison against analytic X^2\n")
  ans$P.empirical = sum(ans$observed.ln.LR <= ans$empdist)/ans$n.resample
  cat("Contrasting groups, compare summed ln-LR against summed bootstrap X^2 distribution:\n")
  ans$quants = quantile(ans$empdist, q1)
  cat(sprintf("N samples = %d  groups = %d  observed.ln.LR = %f  resamples = %d  boot.X2[%s] = [%.2f %.2f %.2f %.2f]  P = %g\n", 
              ans$N.samples,
              ans$N.groups,
              ans$observed.ln.LR,
              ans$n.resample,
              paste(collapse=" ", q1),
              ans$quants[1], ans$quants[2], ans$quants[3], ans$quants[4],
              ans$P.empirical))
  ans$version = .alleleDiversityTestsVersion
  ####
  class(ans) <- c(class(ans), "allele.diversityTest")
  invisible(ans)
}


#---------------------------------------------


allele.alphaContrastTest = function(lst.a, lst.b,
                             zero.var.adjust=TRUE, 
                             n.resample=10000, 
                             method=c("bootstrap", "permute"))
{
  method <- match.arg(method)
  stopifnot(length(lst.a) == length(lst.b) && all(names(lst.a) == names(lst.b)))
  name.a = deparse(substitute(lst.a))
  name.b = deparse(substitute(lst.b))
  nm = names(lst.a)
  .RT = .diversityTest.ReverseTerms
  .diversityTest.ReverseTerms = FALSE
  q1 <- c(0.5, 0.95, 0.99, 0.999)

  ans <- list()
  ans$name.a = name.a
  ans$name.b = name.b
  ans$names = nm
  # for the grand test
  ans$observed.ln.LR = 0
  ans$empdist = numeric(n.resample)
  # sublists for each locus

  for (l in 1:length(lst.a)) {

    locus = nm[l]

    this.ans = list()
    this.ans$locus = locus

    tab.a = lst.a[[locus]]
    a.vardist <- lapply(.diversityTest.distmat(tab.a), 
                        function(x) diag(.diversityTest.gower(x)))  
    n.a <- unlist(lapply(a.vardist, length))
    N.a <- sum(n.a)
    G.a <- length(n.a)
    this.ans$name.a = name.a
    this.ans$name.b = name.b
    this.ans$N.a <- N.a
    this.ans$G.a <- G.a
    terms.a = .diversityTest.CalcTerms(n.a, a.vardist, zero.var.adjust)

    tab.b = lst.b[[locus]]
    b.vardist <- lapply(.diversityTest.distmat(tab.b), 
                        function(x) diag(.diversityTest.gower(x)))  
    n.b <- unlist(lapply(b.vardist, length))
    N.b <- sum(n.b)
    G.b <- length(n.b)
    this.ans$N.b <- N.b
    this.ans$G.b <- G.b
    terms.b = .diversityTest.CalcTerms(n.b, b.vardist, zero.var.adjust)

    V.a.b.p = (((N.a - G.a) * terms.a$V.p) + ((N.b - G.b) * terms.b$V.p)) / (N.a + N.b - G.a - G.b)
    observed.ln.LR.a.b = ((N.a + N.b - G.a - G.b) * log(V.a.b.p)) - ((N.a - G.a) * log(terms.a$V.p)) - ((N.b - G.b) * log(terms.b$V.p))

    # Combine A and B into strata for comparison
    n.a.b = c(a=N.a, b=N.b)
    a.b.vardist = list(a=unlist(a.vardist, use.names=FALSE), b=unlist(b.vardist, use.names=FALSE))
    N = sum(n.a.b)
    G = length(n.a.b)
    DF = G -1
    this.ans$N.samples <- N
    this.ans$N.groups <- G
    PVAL = pchisq(observed.ln.LR.a.b, df=DF, lower.tail=TRUE)
    cat("Bartlett's Test for Heteroscedasticity in **Inter-group** Alpha Variances between",name.a,"and",name.b,"\n")
    cat("Current locus:", locus, "\n")
    cat("For comparison against analytic X^2, see data structure\n")
    #cat("Comparing against X^2 distribution directly:\n")
    #cat(sprintf("Samples N.a = %d, N.b = %d  groups = %d  observed.ln.LR.a.b = %f  df=(G-1) = %d  P = %g\n", 
    #      N.a, N.b, G, observed.ln.LR.a.b, DF, PVAL))
    this.ans$observed.ln.LR <- observed.ln.LR.a.b
    this.ans$df.X2 <- G - 1
    this.ans$P.analytic <- PVAL
    nulldist = .diversityTest.NullDist(obs=observed.ln.LR.a.b,
                                               n.g=n.a.b,
                                               g.vardist=a.b.vardist,
                                               zero.var.adjust, method, n.resample)
    PVAL <- sum(observed.ln.LR.a.b <= nulldist)/n.resample
    cat("Contrasting groups, compare against bootstrap X^2 distribution:\n")
    q2 <- quantile(nulldist, q1)
    o = sprintf("Samples N.a = %d, N.b = %d  groups = %d  observed.ln.LR.a.b = %f  resamples = %d  boot.X2[%s] = [%.2f %.2f %.2f %.2f]  P = %g\n\n", 
          N.a, N.b, G, observed.ln.LR.a.b, n.resample, paste(collapse=" ", q1), 
          q2[1], q2[2], q2[3], q2[4], PVAL)
    cat(o)
    this.ans$n.resample <- n.resample
    this.ans$N.a = N.a
    this.ans$N.b = N.b
    this.ans$resample.method <- method
    this.ans$quants <- q2
    this.ans$P.empirical <- PVAL
    this.ans$empdist <- nulldist
    this.ans$a.b.vardist = a.b.vardist
    this.ans$version <- .diversityTestsVersion

    # add the current locus' results to the grand totals
    ans[[ locus ]] = this.ans
    ans$observed.ln.LR = ans$observed.ln.LR + this.ans$observed.ln.LR
    ans$empdist = ans$empdist + this.ans$empdist

  }
  ans$N.samples = ans[[ nm[1] ]]$N.samples
  ans$N.groups = ans[[ nm[1] ]]$N.groups
  ans$N.a = ans[[ nm[1] ]]$N.a
  ans$N.b = ans[[ nm[1] ]]$N.b
  ans$n.resample <- n.resample
  ans$resample.method <- method
  cat("Bartlett's Test for Heteroscedasticity in **Inter-group** Alpha Variances between",name.a,"and",name.b,"\n")
  cat("Over ALL loci\n")
  cat("Skipping comparison against analytic X^2\n")
  ans$P.empirical = sum(ans$observed.ln.LR <= ans$empdist)/ans$n.resample
  cat("Contrasting groups, compare summed ln-LR against summed bootstrap X^2 distribution:\n")
  ans$quants = quantile(ans$empdist, q1)
  cat(sprintf("Samples N.a = %d, N.b = %d, groups = %d, observed.ln.LR = %f, resamples = %d, boot.x2[%s] = [%.2f %.2f %.2f %.2f]   P = %g\n\n\n",
              ans$N.a,
              ans$N.b,
              ans$N.groups,
              ans$observed.ln.LR,
              ans$n.resample,
              paste(collapse=" ", q1),
              ans$quants[1], ans$quants[2], ans$quants[3], ans$quants[4],
              ans$P.empirical))
  ans$version <- .alleleDiversityTestsVersion
  .diversityTest.ReverseTerms = .RT
  ####
  class(ans) <- c(class(ans), "allele.diversityTest.AlphaContrast")
  invisible(ans)
}


#---------------------------------------------


allele.gammaContrastTest = function(lst.a, lst.b,
                  zero.var.adjust=TRUE, 
                  n.resample=10000, 
                  method=c("bootstrap", "permute"))
{
  method = match.arg(method)
  stopifnot(length(lst.a) == length(lst.b) && all(names(lst.a) == names(lst.b)))
  name.a = deparse(substitute(lst.a))
  name.b = deparse(substitute(lst.b))
  nm = names(lst.a)
  q1 <- c(0.5, 0.95, 0.99, 0.999)

  ans = list()
  ans$name.a = name.a
  ans$name.b = name.b
  ans$names = nm
  ans$observed.ln.LR = 0
  ans$empdist = numeric(n.resample)

  for (l in 1:length(lst.a)) {

    locus = nm[l]

    this.ans = list()
    this.ans$locus = locus

    tab.a = lst.a[[ nm[l] ]]
    tab.b = lst.b[[ nm[l] ]]

    X.a.k = apply(tab.a, 2, sum)
    X.b.k = apply(tab.b, 2, sum)
    N.a = sum(X.a.k)
    N.b = sum(X.b.k)
    R.a.0 = sum((X.a.k * (X.a.k - 1)) / (N.a * (N.a - 1)))
    R.b.0 = sum((X.b.k * (X.b.k - 1)) / (N.b * (N.b - 1)))
    V.a.tot = 1 - R.a.0
    V.b.tot = 1 - R.b.0
    V.a.b.tot = ((N.a - 1) * V.a.tot + (N.b - 1) * V.b.tot) / (N.a + N.b - 2)
    observed.ln.LR.a.b = ((N.a + N.b - 2) * log(V.a.b.tot)) - ((N.a - 1) * log(V.a.tot)) - ((N.b - 1) * log(V.b.tot)) 

    
    a.vardist <- list(a=diag(.diversityTest.gower(.diversityTest.distmat(X.a.k))))
    n.a <- unlist(lapply(a.vardist, length))
    stopifnot(sum(n.a) == N.a)
    N.a <- sum(n.a)
    G.a <- length(n.a)
    this.ans$N.a <- N.a
    this.ans$G.a <- G.a
    terms.a = .diversityTest.CalcTerms(n.a, a.vardist, zero.var.adjust)
    cat(sprintf("%s  terms.a$V.p = %f  V.a.tot = %f\n", locus, terms.a$V.p, V.a.tot))

    b.vardist <- list(b=diag(.diversityTest.gower(.diversityTest.distmat(X.b.k))))
    n.b <- unlist(lapply(b.vardist, length))
    stopifnot(sum(n.b) == N.b)
    N.b <- sum(n.b)
    G.b <- length(n.b)
    this.ans$N.b <- N.b
    this.ans$G.b <- G.b
    terms.b = .diversityTest.CalcTerms(n.b, b.vardist, zero.var.adjust)
    cat(sprintf("%s  terms.b$V.p = %f  V.b.tot = %f\n", locus, terms.b$V.p, V.b.tot))

    # Combine A and B into stratta for comparison
    n.a.b = c(a=N.a, b=N.b)
    a.b.vardist = list(a=unlist(a.vardist, use.names=FALSE), b=unlist(b.vardist, use.names=FALSE))
    N = sum(n.a.b)
    G = length(n.a.b)
    DF = G -1
    this.ans$N.samples <- N
    this.ans$N.groups <- G
    PVAL = pchisq(observed.ln.LR.a.b, df=DF, lower.tail=TRUE)
    cat("Bartlett's Test for Heteroscedasticity in **Inter-group** Gamma Variances between", name.a, "and", name.b, "\n")
    cat("Current locus:", locus, "\n")
    cat("For comparison against analytic X^2, see data structure\n")
    #cat("Comparing against X-2 distribution directly:\n")
    #cat(sprintf("Samples N.a = %d, N.b = %d  groups = %d  observed.ln.LR.a.b = %f  df=(G-1) = %d  P = %g\n", 
    #      N.a, N.b, G, observed.ln.LR.a.b, DF, PVAL))
    this.ans$observed.ln.LR <- observed.ln.LR.a.b
    this.ans$df.X2 <- DF
    this.ans$P.analytic <- PVAL
    nulldist = .diversityTest.NullDist(obs=observed.ln.LR.a.b,
                                               n.g=n.a.b,
                                               g.vardist=a.b.vardist,
                                               zero.var.adjust, method, n.resample)
    PVAL <- sum(observed.ln.LR.a.b <= nulldist)/n.resample
    cat("Contrasting groups, compare against bootstrap X-2 distribution:\n")
    q2 <- quantile(nulldist, q1)
    cat(sprintf("Samples N.a = %d, N.b = %d  groups = %d  observed.ln.LR.a.b = %f  resamples = %d  boot.X2[%s] = [%.2f %.2f %.2f %.2f]  P = %g\n\n", 
          N.a, N.b, G, observed.ln.LR.a.b, n.resample, paste(collapse=" ", q1), 
          q2[1], q2[2], q2[3], q2[4], PVAL))
    this.ans$n.resample <- n.resample
    this.ans$resample.method <- method
    this.ans$quants <- q2
    this.ans$P.empirical <- PVAL
    this.ans$empdist <- nulldist
    this.ans$a.b.vardist = a.b.vardist
    this.ans$version <- .diversityTestsVersion

    # add the current locus' results to the grand totals
    ans[[ locus ]] = this.ans
    ans$observed.ln.LR = ans$observed.ln.LR + this.ans$observed.ln.LR
    ans$empdist = ans$empdist + this.ans$empdist

  }
  ans$N.samples = ans[[ nm[1] ]]$N.samples
  ans$N.groups = ans[[ nm[1] ]]$N.groups
  ans$N.a = ans[[ nm[1] ]]$N.a
  ans$G.a = ans[[ nm[1] ]]$G.a
  ans$N.b = ans[[ nm[1] ]]$N.b
  ans$G.b = ans[[ nm[1] ]]$G.b
  ans$n.resample = n.resample
  ans$resample.method = method
  cat("Bartlett's Test for Heteroscedasticity in **Inter-group** Gamma Variances between",name.a,"and",name.b,"\n")
  cat("Over ALL loci\n")
  cat("Skipping comparison against analytic X^2\n")
  ans$P.empirical = sum(ans$observed.ln.LR <= ans$empdist)/ans$n.resample
  cat("Contrasting groups, compare summed ln-LR against summed bootstrap X^2 distribution:\n")
  ans$quants = quantile(ans$empdist, q1)
  cat(sprintf("Samples N.a = %d, G.a = %d, N.b = %d, G.b = %d, groups = %d, observed.ln.LR = %f, resamples = %d, boot.x2[%s] = [%.2f %.2f %.2f %.2f]   P = %g\n\n\n",
              ans$N.a,
              ans$G.a,
              ans$N.b,
              ans$G.b,
              ans$N.groups,
              ans$observed.ln.LR,
              ans$n.resample,
              paste(collapse=" ", q1),
              ans$quants[1], ans$quants[2], ans$quants[3], ans$quants[4],
              ans$P.empirical))
  ans$version <- .alleleDiversityTestsVersion
  ####
  class(ans) <- c(class(ans), "allele.diversityTest.GammaContrast")
  invisible(ans)

}


