# diversityTests.R

# Provide functions for performing dispersal diversity tests (Scofield et al
# 2012 American Naturalist 180(6) 719-723,
# http://www.jstor.org/stable/10.1086/668202).  Tests require as input one or
# more tables of counts in sites (rows) X sources (columns) format.

.diversityTestsVersion = "0.3.1"

# Copyright (c) 2012 Douglas G. Scofield, Umeå Plant Science Centre, Umeå, Sweden
#
# douglas.scofield@plantphys.umu.se
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
# FUNCTIONS PROVIDED 
#
# 
# See Scofield et al. in review for methodological details
# 
#
# alphaDiversityTest(tab) : Test for differences in alpha diversity among sites
#                           within a single dataset
# 
# alphaContrastTest(tab.a, tab.b)
#                         : Test whether there is a difference in the alpha
#                           diversity between two datasets
# 
# alphaContrastTest.3(tab.a, tab.b, tab.c)
#                         : Test whether there is a difference in the alpha
#                           diversity among three datasets
# 
# plotAlphaTest(result)   : Plot the list returned from alphaDiversityTest()
#                           or alphaContrastTest() for evaluation
# 
# pairwiseMeanTest(tab)   : Test whether mean pairwise divergence/overlap among 
#                           sites is different from the null espectation
#
# plotPairwiseMeanTest()  : Plot the list returned from the above test for
#                           evaluation
# 
# gammaContrastTest(tab.a, tab.b)
#                         : Test whether there is a difference in the gamma
#                           diversity between two datasets
# 
# gammaContrastTest.3(tab.a, tab.b, tab.c)
#                         : Test whether there is a difference in the gamma
#                           diversity among three datasets
#
#
#
# CHANGELOG
#
# 0.3.1: Modify plotAlphaTest to take n.resample into account
# 0.3: Versioning and collaborators/funding statement.
# 0.2: Add q.nielsen to pairwiseMeanTest.
# 0.1: First release
#
#
# TODO
#
# --- Turn this into an actual R package


source("pmiDiversity.R")

# This code relies on pmiDiversity() and nielsenTransform(), both defined in
# pmiDiversity.R


#---------------------------------------------


alphaDiversityTest <- function(tab, 
                               zero.var.adjust=TRUE, 
                               n.resample=10000, 
                               method=c("bootstrap", "permute"))
{
  method <- match.arg(method)
  ans <- list()
  g.vardist <- lapply(.diversityTest.distmat(tab), 
                      function(x) diag(.diversityTest.gower(x)))  
  n.g <- unlist(lapply(g.vardist, length))
  N <- sum(n.g)
  G <- length(n.g)
  ans$N.samples <- N
  ans$N.groups <- G

  terms = .diversityTest.CalcTerms(n.g, g.vardist, zero.var.adjust)
  PVAL = pchisq(terms$ln.LR, df=terms$DF, 
                lower.tail=ifelse(.diversityTest.ReverseTerms, FALSE, TRUE))
  cat("Bartlett's Test for Heteroscedasticity in Intra-group Variances\n\n")
  cat("Comparing against X-2 distribution directly:\n")
  cat(sprintf("N samples = %d  groups = %d  observed.ln.LR = %f  df=(G-1) = %d  P = %g\n", 
        N, G, terms$ln.LR, terms$DF, PVAL))
  ans$observed.ln.LR <- terms$ln.LR
  ans$df.X2 <- terms$DF
  ans$P.analytic <- PVAL
  nulldist <- .diversityTest.NullDist(obs=terms$ln.LR, 
                                         n.g=n.g, g.vardist=g.vardist, 
                                         zero.var.adjust=zero.var.adjust, 
                                         method=method, n.resample=n.resample)
  PVAL <- sum(terms$ln.LR <= nulldist)/n.resample
  cat("\nComparing against bootstrap X-2 distribution:\n")
  q1 <- if(.diversityTest.ReverseTerms)
    c(0.5, 0.99, 0.999)
  else
    c(0.5, 0.01, 0.001)
  q2 <- quantile(nulldist, q1)
  cat(sprintf("N samples = %d  groups = %d  observed.ln.LR = %f  resamples = %d  boot.X2[%s] = [%.2f %.2f %.2f]  P = %g\n", 
        N, G, terms$ln.LR, n.resample, paste(collapse=" ", q1), 
        q2[1], q2[2], q2[3], PVAL))
  ans$n.resample <- n.resample
  ans$resample.method <- method
  ans$quants <- q2
  ans$P.empirical <- PVAL
  ans$empdist <- nulldist
  ans$version <- .diversityTestsVersion
  ####
  class(ans) <- c(class(ans), "diversityTest")
  invisible(ans)
}


#---------------------------------------------


alphaContrastTest = function(tab.a, tab.b,
                             zero.var.adjust=TRUE, 
                             n.resample=10000, 
                             method=c("bootstrap", "permute"))
{
  method <- match.arg(method)
  .RT = .diversityTest.ReverseTerms
  .diversityTest.ReverseTerms = FALSE
  ans <- list()
  a.vardist <- lapply(.diversityTest.distmat(tab.a), 
                      function(x) diag(.diversityTest.gower(x)))  
  n.a <- unlist(lapply(a.vardist, length))
  N.a <- sum(n.a)
  G.a <- length(n.a)
  ans$N.a <- N.a
  ans$G.a <- G.a
  terms.a = .diversityTest.CalcTerms(n.a, a.vardist, zero.var.adjust)
  # V.a.p = terms.a$V.p

  b.vardist <- lapply(.diversityTest.distmat(tab.b), 
                      function(x) diag(.diversityTest.gower(x)))  
  n.b <- unlist(lapply(b.vardist, length))
  N.b <- sum(n.b)
  G.b <- length(n.b)
  ans$N.b <- N.b
  ans$G.b <- G.b
  terms.b = .diversityTest.CalcTerms(n.b, b.vardist, zero.var.adjust)
  # V.b.p = terms.b$V.p
  V.a.b.p = (((N.a - G.a) * terms.a$V.p) + ((N.b - G.b) * terms.b$V.p)) / (N.a + N.b - G.a - G.b)
  observed.ln.LR.a.b = ((N.a + N.b - G.a - G.b) * log(V.a.b.p)) - ((N.a - G.a) * log(terms.a$V.p)) - ((N.b - G.b) * log(terms.b$V.p))

  # Combine A and B into strata for comparison
  n.a.b = c(a=N.a, b=N.b)
  a.b.vardist = list(a=unlist(a.vardist, use.names=FALSE), b=unlist(b.vardist, use.names=FALSE))
  N = sum(n.a.b)
  G = length(n.a.b)
  DF = G -1
  ans$N.samples <- N
  ans$N.groups <- G
  PVAL = pchisq(observed.ln.LR.a.b, df=DF, lower.tail=TRUE)
  cat("Bartlett's Test for Heteroscedasticity in **Inter-group** Alpha Variances\n\n")
  cat("Comparing against X-2 distribution directly:\n")
  cat(sprintf("Samples N.a = %d, N.b = %d  groups = %d  observed.ln.LR.a.b = %f  df=(G-1) = %d  P = %g\n", 
        N.a, N.b, G, observed.ln.LR.a.b, DF, PVAL))
  ans$observed.ln.LR <- observed.ln.LR.a.b
  ans$df.X2 <- G - 1
  ans$P.analytic <- PVAL
  nulldist = .diversityTest.NullDist(obs=observed.ln.LR.a.b,
                                             n.g=n.a.b,
                                             g.vardist=a.b.vardist,
                                             zero.var.adjust, method, n.resample)
  PVAL <- sum(observed.ln.LR.a.b <= nulldist)/n.resample
  cat("\nContrasting groups, compare against bootstrap X-2 distribution:\n")
  q1 <- c(0.5, 0.01, 0.001)
  q2 <- quantile(nulldist, q1)
  cat(sprintf("Samples N.a = %d, N.b = %d  groups = %d  observed.ln.LR.a.b = %f  resamples = %d  boot.X2[%s] = [%.2f %.2f %.2f]  P = %g\n", 
        N.a, N.b, G, observed.ln.LR.a.b, n.resample, paste(collapse=" ", q1), 
        q2[1], q2[2], q2[3], PVAL))
  ans$n.resample <- n.resample
  ans$resample.method <- method
  ans$quants <- q2
  ans$P.empirical <- PVAL
  ans$empdist <- nulldist
  ans$a.b.vardist = a.b.vardist
  ans$version <- .diversityTestsVersion
  .diversityTest.ReverseTerms = .RT
  ####
  class(ans) <- c(class(ans), "diversityTest.AlphaContrast")
  invisible(ans)

}


#---------------------------------------------


alphaContrastTest.3 = function(tab.a, tab.b, tab.c,
                               zero.var.adjust=TRUE, 
                               n.resample=10000, 
                               method=c("bootstrap", "permute"))
{
  method <- match.arg(method)
  .RT = .diversityTest.ReverseTerms
  .diversityTest.ReverseTerms = FALSE
  ans <- list()

  a.vardist <- lapply(.diversityTest.distmat(tab.a), 
                      function(x) diag(.diversityTest.gower(x)))  
  n.a <- unlist(lapply(a.vardist, length))
  N.a <- sum(n.a)
  G.a <- length(n.a)
  ans$N.a <- N.a
  ans$G.a <- G.a
  terms.a = .diversityTest.CalcTerms(n.a, a.vardist, zero.var.adjust)
  # V.a.p = terms.a$V.p

  b.vardist <- lapply(.diversityTest.distmat(tab.b), 
                      function(x) diag(.diversityTest.gower(x)))  
  n.b <- unlist(lapply(b.vardist, length))
  N.b <- sum(n.b)
  G.b <- length(n.b)
  ans$N.b <- N.b
  ans$G.b <- G.b
  terms.b = .diversityTest.CalcTerms(n.b, b.vardist, zero.var.adjust)
  # V.b.p = terms.b$V.p

  c.vardist <- lapply(.diversityTest.distmat(tab.c), 
                      function(x) diag(.diversityTest.gower(x)))  
  n.c <- unlist(lapply(c.vardist, length))
  N.c <- sum(n.c)
  G.c <- length(n.c)
  ans$N.c <- N.c
  ans$G.c <- G.c
  terms.c = .diversityTest.CalcTerms(n.c, c.vardist, zero.var.adjust)
  # V.c.p = terms.c$V.p

  V.a.b.c.p = (((N.a - G.a) * terms.a$V.p) + ((N.b - G.b) * terms.b$V.p) + ((N.c - G.c) * terms.c$V.p)) / (N.a + N.b + N.c - G.a - G.b - G.c)
  observed.ln.LR.a.b.c = ((N.a + N.b + N.c - G.a - G.b - G.c) * log(V.a.b.c.p)) - ((N.a - G.a) * log(terms.a$V.p)) - ((N.b - G.b) * log(terms.b$V.p)) - ((N.c - G.c) * log(terms.c$V.p))


  # Combine A B C into strata for comparison
  n.a.b.c = c(a=N.a, b=N.b, c=N.c)
  a.b.c.vardist = list(a=unlist(a.vardist, use.names=FALSE), 
                       b=unlist(b.vardist, use.names=FALSE), 
                       c=unlist(c.vardist, use.names=FALSE))
  N = sum(n.a.b.c)
  G = length(n.a.b.c)
  DF = G -1
  ans$N.samples <- N
  ans$N.groups <- G
  PVAL = pchisq(observed.ln.LR.a.b.c, df=DF, lower.tail=TRUE)
  cat("Bartlett's Test for Heteroscedasticity in **Inter-group** Alpha Variances\n\n")
  cat("Comparing against X-2 distribution directly:\n")
  cat(sprintf("Samples N.a = %d, N.b = %d, N.c 0 %d  groups = %d  observed.ln.LR.a.b.c = %f  df=(G-1) = %d  P = %g\n", 
        N.a, N.b, N.c, G, observed.ln.LR.a.b.c, DF, PVAL))
  ans$observed.ln.LR <- observed.ln.LR.a.b.c
  ans$df.X2 <- G - 1
  ans$P.analytic <- PVAL
  nulldist = .diversityTest.NullDist(obs=observed.ln.LR.a.b.c,
                                             n.g=n.a.b.c,
                                             g.vardist=a.b.c.vardist,
                                             zero.var.adjust, method, n.resample)
  PVAL <- sum(observed.ln.LR.a.b.c <= nulldist)/n.resample
  cat("\nContrasting groups, compare against bootstrap X-2 distribution:\n")
  q1 <- c(0.5, 0.01, 0.001)
  q2 <- quantile(nulldist, q1)
  cat(sprintf("Samples N.a = %d, N.b = %d, N.c = %d  groups = %d  observed.ln.LR.a.b.c = %f  resamples = %d  boot.X2[%s] = [%.2f %.2f %.2f]  P = %g\n", 
        N.a, N.b, N.c, G, observed.ln.LR.a.b.c, n.resample, paste(collapse=" ", q1), 
        q2[1], q2[2], q2[3], PVAL))
  ans$n.resample <- n.resample
  ans$resample.method <- method
  ans$quants <- q2
  ans$P.empirical <- PVAL
  ans$empdist <- nulldist
  ans$a.b.c.vardist = a.b.c.vardist
  ans$version <- .diversityTestsVersion
  .diversityTest.ReverseTerms = .RT
  ####
  class(ans) <- c(class(ans), "diversityTest.AlphaContrast.3")
  invisible(ans)

}


#---------------------------------------------


plotAlphaTest <- function(result, 
                          add.analytic=FALSE, 
                          compress.x=TRUE, 
                          do.postscript=FALSE,
                          ...)
{
  if (! "diversityTest" %in% class(result) &&
      ! "diversityTest.AlphaContrast" %in% class(result))
    stop(deparse(substitute(result)), "not a result of alphaDiversityTest() or alphaContrastTest()")
  if (do.postscript)
    eps("pBDT.eps", width=8, height=3)
  par(mar=c(3.5, 3.5, 1.5, 1), mgp=c(2, 0.5, 0))
  plotted.observed.ln.LR <- result$observed.ln.LR
  empdist.range <- diff(range(result$empdist[-result$n.resample]))
  full.range <- diff(range(result$empdist))
  breaks <- 50
  compressed.ratio <- 1.2
  if (compress.x && full.range > empdist.range * compressed.ratio) {
    xlim <- range(result$empdist[-result$n.resample])
    h <- hist(result$empdist[-result$n.resample], 
              breaks=seq(xlim[1], xlim[2], length.out=breaks), 
              plot=FALSE)
    xlim <- xlim * c(1, compressed.ratio)
    plotted.observed.ln.LR <- xlim[2]
  } else {
    h <- hist(result$empdist[-result$n.resample], breaks=50, plot=FALSE)
    xlim <- range(c(h$breaks, result$empdist[result$n.resample]))
  }
  if (add.analytic) {
    analytic.x <- seq(xlim[1], xlim[2], length.out=200)
    analytic.X2 <- dchisq(x=analytic.x, df=result$df.ln.LR)
    ylim <- range(c(0, h$count, analytic.X2))
  } else {
    ylim <- range(c(0, h$count))
  }
  plot(h, xlim=xlim, ylim=ylim, freq=TRUE, 
       xlab=expression("ln(LR) value"), 
       ylab="Frequency", 
       main="", ...)
  lines(rep(plotted.observed.ln.LR, 2), c(ylim[1], ylim[2]*0.5), col="darkgray", lty=1, lwd=2)
  if (add.analytic) {
    lines(analytic.x, analytic.X2, lty=3, lwd=1)
  }
  legend("topright", 
         legend=substitute("Observed ln(LR)" == OX2,
                           list(OX2=sprintf("%.3f", round(result$observed.ln.LR, 3)))),
         bty="n")
  if (do.postscript)
    dev.off()
}


#---------------------------------------------


pairwiseMeanTest <- function(tab, 
                             n.iter=10000, 
                             method=c("r", "q", "q.nielsen"), 
                             statistic=c("divergence", "overlap"))
{
  method <- match.arg(method)
  statistic <- match.arg(statistic)
  pmiresults <- pmiDiversity(tab)
  gammafreq <- apply(tab, 2, sum)
  K <- sum(gammafreq > 0)
  n.g <- apply(tab, 1, sum)
  names.g <- rownames(tab)
  G <- length(n.g)
  K <- length(gammafreq)
  obs <- if (method == "r") {
    if (statistic == "divergence") { pmiresults$r.divergence } else { pmiresults$r.overlap }
  } else if (method == "q") { 
    if (statistic == "divergence") { pmiresults$q.divergence } else { pmiresults$q.overlap }
  } else if (method == "q.nielsen") { 
    if (statistic == "divergence") { pmiresults$q.nielsen.divergence } else { pmiresults$q.nielsen.overlap }
  } else {
    stop("unimplemented method")
  }
  nulldist <- numeric(n.iter)
  for (i in 1:(n.iter - 1)) {
    mat <- c()
    for (g in 1:G) {
      mat <- rbind(mat, t(rmultinom(1, n.g[g], gammafreq)))
    }
    reltab <- sweep(mat, 1, n.g, FUN="/")
    Q.mat <- reltab %*% t(reltab)
    cache.gg <- diag(Q.mat) # q.gg
    diag(Q.mat) <- 0
    if (method == "r") {
      cache.gg <- ((n.g * cache.gg) - 1) / (n.g - 1) # r.gg
    } else if (method == "q.nielsen") {
      cache.gg <- nielsenTransform(cache.gg) # q.nielsen.gg
    }
    if (statistic == "divergence") {
      nulldist[i] <- 1 - sum(Q.mat)/((G - 1)*sum(cache.gg))
    } else {
      nulldist[i] <- sum(Q.mat)/((G - 1)*sum(cache.gg))
    }
  }
  nulldist[n.iter] <- obs
  P.lower <- sum(obs >= nulldist) / n.iter
  P.upper <- sum(obs <= nulldist) / n.iter
  ans <- list(n.cache=G, n.source=K, n.seed=sum(n.g),
              obs=obs, n.iter=n.iter, nulldist=nulldist, 
              P.lower=P.lower, P.upper=P.upper, method=method, statistic=statistic)
  ans$version <- .diversityTestsVersion
  class(ans) <- c(class(ans), "class.PairwiseMeanTest")
  ####
  ans
}


#---------------------------------------------


plotPairwiseMeanTest <- function(result, ...)
{
  if (! "class.PairwiseMeanTest" %in% class(result))
    stop(deparse(substitute(result)), "not a result of pairwiseMeanTest()")
  par(mar=c(2.8, 2.8, 0.5, 0), mgp=c(1.6, 0.4, 0), tcl=-0.25, cex=0.9)
  xlim <- range(c(0, 1, result$nulldist))
  breaks <- seq(xlim[1], xlim[2], length.out=50)
  h <- hist(result$nulldist[-result$n.iter], breaks=breaks, plot=FALSE)
  ylim <- range(c(0, h$counts))
  if (result$method == "q") {
    if (result$statistic == "overlap") {
      plot(h, xlim=xlim, ylim=ylim, freq=TRUE,
          xlab=expression(bar(omega)*minute*" value"), ylab="Frequency", main="", ...)
    } else {
      plot(h, xlim=xlim, ylim=ylim, freq=TRUE,
          xlab=expression(bar(delta)*minute*" value"), ylab="Frequency", main="", ...)
    }
  } else {
    if (result$statistic == "overlap") {
      plot(h, xlim=xlim, ylim=ylim, freq=TRUE,
         xlab=expression(bar(omega)*" value"), ylab="Frequency", main="", ...)
    } else {
      plot(h, xlim=xlim, ylim=ylim, freq=TRUE,
          xlab=expression(bar(delta)*" value"), ylab="Frequency", main="", ...)
    }
  }
  lines(rep(result$obs, 2), c(ylim[1], ylim[2]*0.5), col="black", lty=2, lwd=4)
  if (result$method == "q") {
    if (result$statistic == "overlap") {
      legend("topleft", bty="n", legend=substitute("Observed "*bar(omega)*minute == OBS,
                               list(OBS=round(result$obs, 3))))
    } else {
      legend("topleft", bty="n", legend=substitute("Observed "*bar(delta)*minute == OBS,
                               list(OBS=round(result$obs, 3))))
    }
  } else {
    if (result$statistic == "overlap") {
      legend("topleft", bty="n", legend=substitute("Observed "*bar(omega) == OBS,
                               list(OBS=round(result$obs, 3))))
    } else {
      legend("topleft", bty="n", legend=substitute("Observed "*bar(delta) == OBS,
                               list(OBS=round(result$obs, 3))))
    }
  }
}


#---------------------------------------------


gammaContrastTest = function(tab.a, tab.b,
                  zero.var.adjust=TRUE, 
                  n.resample=10000, 
                  method=c("bootstrap", "permute"))
{
  ans = list()
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
  ans$N.a <- N.a
  ans$G.a <- G.a
  terms.a = .diversityTest.CalcTerms(n.a, a.vardist, zero.var.adjust)
  cat(sprintf("terms.a$V.p = %f  V.a.tot = %f\n", terms.a$V.p, V.a.tot))

  b.vardist <- list(b=diag(.diversityTest.gower(.diversityTest.distmat(X.b.k))))
  n.b <- unlist(lapply(b.vardist, length))
  stopifnot(sum(n.b) == N.b)
  N.b <- sum(n.b)
  G.b <- length(n.b)
  ans$N.b <- N.b
  ans$G.b <- G.b
  terms.b = .diversityTest.CalcTerms(n.b, b.vardist, zero.var.adjust)
  cat(sprintf("terms.b$V.p = %f  V.b.tot = %f\n", terms.b$V.p, V.b.tot))

  # Combine A and B into stratta for comparison
  n.a.b = c(a=N.a, b=N.b)
  a.b.vardist = list(a=unlist(a.vardist, use.names=FALSE), b=unlist(b.vardist, use.names=FALSE))
  N = sum(n.a.b)
  G = length(n.a.b)
  DF = G -1
  ans$N.samples <- N
  ans$N.groups <- G
  PVAL = pchisq(observed.ln.LR.a.b, df=DF, lower.tail=TRUE)
  cat("Bartlett's Test for Heteroscedasticity in **Inter-group** Gamma Variances\n\n")
  cat("Comparing against X-2 distribution directly:\n")
  cat(sprintf("Samples N.a = %d, N.b = %d  groups = %d  observed.ln.LR.a.b = %f  df=(G-1) = %d  P = %g\n", 
        N.a, N.b, G, observed.ln.LR.a.b, DF, PVAL))
  ans$observed.ln.LR <- observed.ln.LR.a.b
  ans$df.X2 <- DF
  ans$P.analytic <- PVAL
  nulldist = .diversityTest.NullDist(obs=observed.ln.LR.a.b,
                                             n.g=n.a.b,
                                             g.vardist=a.b.vardist,
                                             zero.var.adjust, method, n.resample)
  PVAL <- sum(observed.ln.LR.a.b <= nulldist)/n.resample
  cat("\nContrasting groups, compare against bootstrap X-2 distribution:\n")
  q1 <- c(0.5, 0.01, 0.001)
  q2 <- quantile(nulldist, q1)
  cat(sprintf("Samples N.a = %d, N.b = %d  groups = %d  observed.ln.LR.a.b = %f  resamples = %d  boot.X2[%s] = [%.2f %.2f %.2f]  P = %g\n", 
        N.a, N.b, G, observed.ln.LR.a.b, n.resample, paste(collapse=" ", q1), 
        q2[1], q2[2], q2[3], PVAL))
  ans$n.resample <- n.resample
  ans$resample.method <- method
  ans$quants <- q2
  ans$P.empirical <- PVAL
  ans$empdist <- nulldist
  ans$a.b.vardist = a.b.vardist
  ans$version <- .diversityTestsVersion
  ####
  class(ans) <- c(class(ans), "diversityTest.GammaContrast")
  invisible(ans)

}


#---------------------------------------------


gammaContrastTest.3 = function(tab.a, tab.b, tab.c,
                               zero.var.adjust=TRUE, 
                               n.resample=10000, 
                               method=c("bootstrap", "permute"))
{
  ans = list()
  X.a.k = apply(tab.a, 2, sum)
  X.b.k = apply(tab.b, 2, sum)
  X.c.k = apply(tab.c, 2, sum)
  N.a = sum(X.a.k)
  N.b = sum(X.b.k)
  N.c = sum(X.c.k)
  R.a.0 = sum((X.a.k * (X.a.k - 1)) / (N.a * (N.a - 1)))
  R.b.0 = sum((X.b.k * (X.b.k - 1)) / (N.b * (N.b - 1)))
  R.c.0 = sum((X.c.k * (X.c.k - 1)) / (N.c * (N.c - 1)))
  V.a.tot = 1 - R.a.0
  V.b.tot = 1 - R.b.0
  V.c.tot = 1 - R.c.0
  V.a.b.c.tot = ((N.a - 1) * V.a.tot + (N.b - 1) * V.b.tot + (N.c - 1) * V.c.tot) / (N.a + N.b + N.c - 3)
  observed.ln.LR.a.b.c = ((N.a + N.b + N.c - 3) * log(V.a.b.c.tot)) - ((N.a - 1) * log(V.a.tot)) - ((N.b - 1) * log(V.b.tot)) - ((N.c - 1) * log(V.c.tot)) 

  
  a.vardist <- list(a=diag(.diversityTest.gower(.diversityTest.distmat(X.a.k))))
  n.a <- unlist(lapply(a.vardist, length))
  stopifnot(sum(n.a) == N.a)
  N.a <- sum(n.a)
  G.a <- length(n.a)
  ans$N.a <- N.a
  ans$G.a <- G.a
  terms.a = .diversityTest.CalcTerms(n.a, a.vardist, zero.var.adjust)
  #cat(sprintf("terms.a$V.p = %f  V.a.tot = %f\n", terms.a$V.p, V.a.tot))

  b.vardist <- list(b=diag(.diversityTest.gower(.diversityTest.distmat(X.b.k))))
  n.b <- unlist(lapply(b.vardist, length))
  stopifnot(sum(n.b) == N.b)
  N.b <- sum(n.b)
  G.b <- length(n.b)
  ans$N.b <- N.b
  ans$G.b <- G.b
  terms.b = .diversityTest.CalcTerms(n.b, b.vardist, zero.var.adjust)
  #cat(sprintf("terms.b$V.p = %f  V.b.tot = %f\n", terms.b$V.p, V.b.tot))

  c.vardist <- list(c=diag(.diversityTest.gower(.diversityTest.distmat(X.c.k))))
  n.c <- unlist(lapply(c.vardist, length))
  stopifnot(sum(n.c) == N.c)
  N.c <- sum(n.c)
  G.c <- length(n.c)
  ans$N.c <- N.c
  ans$G.c <- G.c
  terms.c = .diversityTest.CalcTerms(n.c, c.vardist, zero.var.adjust)
  #cat(sprintf("terms.c$V.p = %f  V.c.tot = %f\n", terms.c$V.p, V.c.tot))

  # Combine A B C into stratta for comparison
  n.a.b.c = c(a=N.a, b=N.b, c=N.b)
  a.b.c.vardist = list(a=unlist(a.vardist, use.names=FALSE), 
                       b=unlist(b.vardist, use.names=FALSE),
                       c=unlist(c.vardist, use.names=FALSE))
  N = sum(n.a.b.c)
  G = length(n.a.b.c)
  DF = G -1
  ans$N.samples <- N
  ans$N.groups <- G
  PVAL = pchisq(observed.ln.LR.a.b.c, df=DF, lower.tail=TRUE)
  cat("Bartlett's Test for Heteroscedasticity in **Inter-group** Gamma Variances\n\n")
  cat("Comparing against X-2 distribution directly:\n")
  cat(sprintf("Samples N.a = %d, N.b = %d, N.c = %d  groups = %d  observed.ln.LR.a.b.c = %f  df=(G-1) = %d  P = %g\n", 
        N.a, N.b, N.c, G, observed.ln.LR.a.b.c, DF, PVAL))
  ans$observed.ln.LR <- observed.ln.LR.a.b.c
  ans$df.X2 <- DF
  ans$P.analytic <- PVAL
  nulldist = .diversityTest.NullDist(obs=observed.ln.LR.a.b.c,
                                             n.g=n.a.b.c,
                                             g.vardist=a.b.c.vardist,
                                             zero.var.adjust, method, n.resample)
  PVAL <- sum(observed.ln.LR.a.b.c <= nulldist)/n.resample
  cat("\nContrasting groups, compare against bootstrap X-2 distribution:\n")
  q1 <- c(0.5, 0.01, 0.001)
  q2 <- quantile(nulldist, q1)
  cat(sprintf("Samples N.a = %d, N.b = %d, N.c = %d  groups = %d  observed.ln.LR.a.b = %f  resamples = %d  boot.X2[%s] = [%.2f %.2f %.2f]  P = %g\n", 
        N.a, N.b, N.c, G, observed.ln.LR.a.b.c, n.resample, paste(collapse=" ", q1), 
        q2[1], q2[2], q2[3], PVAL))
  ans$n.resample <- n.resample
  ans$resample.method <- method
  ans$quants <- q2
  ans$P.empirical <- PVAL
  ans$empdist <- nulldist
  ans$a.b.c.vardist = a.b.c.vardist
  ans$version <- .diversityTestsVersion
  ####
  class(ans) <- c(class(ans), "diversityTest.GammaContrast.3")
  invisible(ans)

}



#---------------------------------------------
#---------------------------------------------


# Internal functions and constants


.diversityTest.epsilon <- 1.0e-12
.diversityTest.ReverseTerms <- TRUE

.diversityTest.distmat <- function(tab, group=dimnames(tab)[[1]], drop=TRUE)
{
  if (is.null(dim(tab))) {
    dim(tab) = c(1, length(tab))
    dimnames(tab) = list(Site="onedim", Genotype=names(tab))
  } 
  if (dim(tab)[1] > 1 && is.null(group)) {
    stop("must supply group(s), all groups not supported")
  } else if (missing(group) && dim(tab)[1] == 1) group <- 1
  G <- dim(tab)[1]
  K <- dim(tab)[2]
  N <- sum(tab)
  N.G <- apply(tab, 1, sum)
  D <- list()
  for (g in group) {
    Dmat <- matrix(1, N.G[g], N.G[g])
    n.K <- tab[g,][tab[g,] > 0]
    cum.n.K <- cumsum(n.K)
    for (src in 1:length(n.K)) {
      # which rows/cols to 0
      slice <- (cum.n.K[src] - n.K[src] + 1):cum.n.K[src] 
      Dmat[slice,slice] <- 0
    }
    D[[as.character(g)]] <- Dmat
  }
  ####
  return(if (length(D) == 1 && drop) D[[1]] else D)
}


.diversityTest.gower <- function(dmat)
{
  if (is.null(dim(dmat))) {
    dim(dmat) = c(1, length(dmat))
    dimnames(dmat) = list(Site="onedim", Genotype=names(dmat))
  } 
  if (! all(dmat == t(dmat))) stop("dmat not symmetric")
  d <- -0.5*dmat
  rd <- apply(d, 1, mean)
  if (! all(rd == apply(d,2,mean))) stop("dmat not symmetric")
  gower.mat <- d + outer(-rd, -rd, "+") + mean(d)
  if (! all(abs(rowSums(gower.mat)) <= .diversityTest.epsilon))
    stop("abs(rowSums(gower.mat)) > .diversityTest.epsilon")
  if (! all(abs(colSums(gower.mat)) <= .diversityTest.epsilon))
    stop("abs(colSums(gower.mat)) > .diversityTest.epsilon")
  ####
  gower.mat
}


.diversityTest.ZeroVarAdjust = function(ss.g, n.g)
{
  # adjust for no diversity within a group (ss.g == 0), in this case,
  # the minimum distance possible prior to dividing by (n.g-1) is
  # 1/(2*n.g*n.g), so make our replacement for 0 be half this, divided
  # of course by (n.g-1):  1/((4*n.g*n.g)*(n.g-1))
  nn <- names(ss.g[ss.g == 0]) # names of ss.g==0 elements
  if (length(nn) > 0)
    ss.g[nn] <- 1/(4*n.g[nn]*n.g[nn]*(n.g[nn]-1)) # and index by names
  ####
  ss.g
}


.diversityTest.CalcTerms <- function(n.g, g.vardist, zero.var.adjust=TRUE)
{
  N = sum(n.g)
  G = length(n.g)
  V.g = unlist(lapply(g.vardist, sum)) / (n.g - 1)
  if (zero.var.adjust)
    V.g = .diversityTest.ZeroVarAdjust(V.g, n.g)
  V.p = ss.pooled = sum((n.g - 1) * V.g) / (N - G)
  term.V.g = sum((n.g - 1) * log(V.g))
  term.V.p = (N - G) * log(V.p)
  term.denom = 1 + (1 / (3 * (G - 1))) * (sum(1 / (n.g - 1)) - (1 / (N - G)))
  ln.LR = if (.diversityTest.ReverseTerms) 
    (term.V.p - term.V.g) / term.denom
  else 
    (term.V.g - term.V.p) / term.denom
  DF = G - 1
  ####
  list(V.g=V.g, V.p=V.p, ln.LR=ln.LR, DF=DF)
}



.diversityTest.NullDist = function(obs, n.g, g.vardist, 
                                      zero.var.adjust=TRUE, 
                                      method=c("bootstrap", "permute"),
                                      n.resample=10000)
{
  method <- match.arg(method)
  N = sum(n.g)
  G = length(n.g)
  cum.n.g <- cumsum(n.g)
  all.g.vardist <- unlist(g.vardist, use.names=FALSE)
  nulldist <- obs  # observed
  for (i in 2:n.resample) {
    p <- switch(method,
          "permute" = sample(all.g.vardist),
          "bootstrap" = sample(all.g.vardist, replace=TRUE))
    for (n in names(g.vardist)) {
      # peel off a slice of the distance permutation for each group
      slice <- (cum.n.g[n]-n.g[n]+1):cum.n.g[n]
      g.vardist[[n]] <- p[slice]
    }
    terms = .diversityTest.CalcTerms(n.g, g.vardist, zero.var.adjust)
    nulldist <- c(nulldist, terms$ln.LR)
  }
  ####
  sort(nulldist)
}

