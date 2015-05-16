# library(readGenalex) # should be redundant
library(dispersalDiversity)


#########################################
context("Testing diversity.divtable basic attributes")


t <- matrix(c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3), byrow = TRUE, nrow = 3, ncol = 4,
             dimnames = list(G = c("a","b","c"), K = c("t1", "t2", "t3", "t4")))
m <- as.divtable(t)
d <- diversity(m)


test_that("diversity.divtable gets basic attributes right", {
    expect_equal(d$num.groups, 3)
    expect_equal(d$num.sources, 4)
    expect_equal(d$num.samples, 24)
    expect_equal(d$num.samples.group, setNames(c(4,8,12), c("a","b","c")))
    expect_equal(d$num.sources.group, setNames(c(4,4,4), c("a","b","c")))
    expect_equal(d$num.groups.source, setNames(c(3,3,3,3), c("t1","t2","t3","t4")))
    # what is the correct y.gh
    expect_equal(d$y.gh, 8)
    # what is the corect prop.y.0.gh
    expect_equal(d$prop.y.0.gh, matrix(c(0,1,1,1,0,1,1,1,0), byrow = TRUE, nrow = 3, ncol = 3,
                                       dimnames = list(G=c("a","b","c"),G=c("a","b","c"))))
    # ... q.0.gh
    expect_equal(d$q.0.gh, matrix(c(0,1,1,0,0,1,0,0,0), byrow = TRUE, nrow = 3, ncol = 3,
                                       dimnames = list(G=c("a","b","c"),G=c("a","b","c"))))

    expect_equal(d$Q.mat, matrix(0.25, nrow = 3, ncol = 3,
                                 dimnames = list(G=c("a","b","c"),G=c("a","b","c"))))
})

context("Testing diversity.divtable contents of q list")

test_that("diversity.divtable gets q right for scaled diversity", {
    expect_equal(d$q$q.gg, setNames(c(0.25,0.25,0.25), c("a","b","c")))
    expect_equal(d$q$q.bar.0, 0.25)
    expect_equal(d$q$q.unweighted.mean, 0.25)
    expect_equal(d$q$alpha.g, setNames(c(4,4,4), c("a","b","c")))
    expect_equal(d$q$alpha.unweighted.mean, 4)
    expect_equal(d$q$d.alpha, 4)
    expect_equal(d$q$d.gamma, 4)
    expect_equal(d$q$d.beta, 1)
    expect_equal(d$q$diversity.mat, matrix(c(4,1,1,1,4,1,1,1,4), byrow = TRUE, nrow = 3, ncol = 3,
                                           dimnames = list(G=c("a","b","c"),G=c("a","b","c"))))
    expect_equal(d$q$divergence.mat, matrix(0, byrow = TRUE, nrow = 3, ncol = 3,
                                            dimnames = list(G=c("a","b","c"),G=c("a","b","c"))))
    expect_equal(d$q$overlap.mat, matrix(1, byrow = TRUE, nrow = 3, ncol = 3,
                                         dimnames = list(G=c("a","b","c"),G=c("a","b","c"))))
    expect_equal(d$q$overlap, 1)
    expect_equal(d$q$divergence, 0)
})

context("Testing diversity.divtable contents of r list")

test_that("diversity.divtable gets r 'right' for scaled diversity, a = 1 makes r Inf there", {
    expect_equal(d$r$q.gg, setNames(c(0,0.1428571429,0.1818181818), c("a","b","c")))
    expect_equal(d$r$q.bar.0, 0.16)
    expect_equal(d$r$q.unweighted.mean, 0.1082251082)
    expect_equal(d$r$alpha.g, setNames(c(Inf,7,5.5), c("a","b","c")))
    expect_equal(d$r$alpha.unweighted.mean, Inf)
    expect_equal(d$r$d.alpha, 9.24)
    expect_equal(d$r$d.gamma, 4.6)
    expect_equal(d$r$d.beta, 0.4978354978)
    expect_equal(d$r$diversity.mat, matrix(c(Inf,3.5,2.75,3.5,7,1.54,2.75,1.54,5.5), byrow = TRUE, nrow = 3, ncol = 3,
                                           dimnames = list(G=c("a","b","c"),G=c("a","b","c"))))
    expect_equal(d$r$divergence.mat, matrix(c(0,-2.5,-1.75,-2.5,0,-0.54,-1.75,-0.54,0), byrow = TRUE, nrow = 3, ncol = 3,
                                            dimnames = list(G=c("a","b","c"),G=c("a","b","c"))))
    expect_equal(d$r$overlap.mat, matrix(c(1,3.5,2.75,3.5,1,1.54,2.75,1.54,1), byrow = TRUE, nrow = 3, ncol = 3,
                                         dimnames = list(G=c("a","b","c"),G=c("a","b","c"))))
    expect_equal(d$r$overlap, 2.31)
    expect_equal(d$r$divergence, -1.31)
})

context("Testing diversity.divtable contents of q.nielsen list")

test_that("diversity.divtable gets q.nielsen right for scaled diversity, low-sample-size problems here as well", {
    expect_equal(d$q.nielsen$q.gg, setNames(c(0.1666666667,0.1734693878,0.1942148760), c("a","b","c")))
    expect_equal(d$q.nielsen$q.bar.0, 0.1863198644)
    expect_equal(d$q.nielsen$q.unweighted.mean, 0.1781169768)
    expect_equal(d$q.nielsen$alpha.g, setNames(c(6,5.764705882,5.148936170), c("a","b","c")))
    expect_equal(d$q.nielsen$alpha.unweighted.mean, 5.637880684)
    expect_equal(d$q.nielsen$d.alpha, 5.614287969)
    expect_equal(d$q.nielsen$d.gamma, 4.540772532)
    expect_equal(d$q.nielsen$d.beta, 0.8087886759)
    expect_equal(d$q.nielsen$diversity.mat, matrix(c(6,1.47,1.385496183,1.47,5.764705882,1.359862385,1.385496183,1.359862385,5.148936170),
                                           byrow = TRUE, nrow = 3, ncol = 3,
                                           dimnames = list(G=c("a","b","c"),G=c("a","b","c"))))
    expect_equal(d$q.nielsen$divergence.mat, matrix(c(0,-0.47,-0.3854961832,-0.47,0,-0.3598623853,-0.3854961832,-0.3598623853,0),
                                            byrow = TRUE, nrow = 3, ncol = 3,
                                            dimnames = list(G=c("a","b","c"),G=c("a","b","c"))))
    expect_equal(d$q.nielsen$overlap.mat, matrix(c(1,1.47,1.385496183,1.47,1,1.359862385,1.385496183,1.359862385,1),
                                         byrow = TRUE, nrow = 3, ncol = 3,
                                         dimnames = list(G=c("a","b","c"),G=c("a","b","c"))))
    expect_equal(d$q.nielsen$overlap, 1.403571992)
    expect_equal(d$q.nielsen$divergence, -0.4035719922)
})


#----------------------------------------

