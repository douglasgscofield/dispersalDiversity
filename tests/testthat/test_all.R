# library(readGenalex) # should be redundant
library(dispersalDiversity)


#########################################
context("Setting up some test data")

# simple 'genalex' objects

g1 <- data.frame(a = 11:13, a.2 = 14:16, b = 101:103, b.2 = 104:106)
x1 <- genalex(1:3, "snurf", g1)

g2 <- data.frame(a = 21:23, a.2 = 24:26, b = 201:203, b.2 = 204:206)
x2 <- genalex(4:6, "snirf", g2)

# more complex 'genalex' objects

data(Qagr_adult_genotypes)
data(Qagr_pericarp_genotypes)

# create random matrix
set.seed(42)
n.sites <- 5
n.sources <- 10
n.samples <- 200
## data frame of site-source pairs
t <- data.frame(site = sample(n.sites, n.samples, replace = TRUE),
                source = round(runif(n.samples) * n.sources + 0.5))
m1 <- as.divtable(table(t))
d1 <- diversity(m1)

test_that("divtable and diversity data types", {
    expect_is(m1, "divtable")
    expect_is(m1, "table")
    expect_true(is.matrix(m1))
    expect_equal(names(dimnames(m1)), c("site", "source"))
    expect_equal(sum(m1), 200)
    expect_equal(dim(m1), c(5, 10))
    expect_is(d1, "diversity")
    expect_is(d1, "list")
    expect_true(! is.null(d1[["q"]]))
    expect_true(! is.null(d1[["r"]]))
    expect_true(! is.null(d1[["q.nielsen"]]))
})


#----------------------------------------

context("Testing as.divtable methods")

test_that("m1 (was table) vs. xt# (xtabs)", {
    xt1 <- xtabs(data = t)
    mxt1 <- as.divtable(xt1)
    expect_is(mxt1, "divtable")
    expect_is(mxt1, "table")
    expect_true(is.matrix(mxt1))
    expect_true(is.null(attr(mxt1, "call")))
    expect_equal(m1, mxt1)
    expect_equal(m1, as.matrix(mxt1))

    xt2 <- xtabs(~ site + source, data = t)
    mxt2 <- as.divtable(xt2)
    expect_is(mxt2, "divtable")
    expect_is(mxt2, "table")
    expect_true(is.matrix(mxt2))
    expect_true(is.null(attr(mxt2, "call")))
    expect_equal(m1, mxt2)

    xt3 <- with(t, xtabs(~ site + source))
    mxt3 <- as.divtable(xt3)
    expect_is(mxt3, "divtable")
    expect_is(mxt3, "table")
    expect_true(is.matrix(mxt3))
    expect_true(is.null(attr(mxt3, "call")))
    expect_equal(m1, mxt3)
})


test_that("m1 (was table) vs. xm# (matrix)", {
    xm1 <- as.matrix(m1)
    expect_true(is.matrix(xm1))
    mxm1 <- as.divtable(xm1)
    expect_is(mxm1, "divtable")
    expect_is(mxm1, "table")
    expect_true(is.matrix(mxm1))
    expect_equal(m1, mxm1)
})


#----------------------------------------

context("Testing diversity.table()")

test_that("results with various table formats are identical", {
    expect_equal(d1, diversity(table(t)))
    xt <- xtabs(data = t)
    expect_equal(d1, diversity(xt))
})

# Many more tests of specific values


context("Testing allele_divtables-class and createAlleleTables()")

x1.adt <- createAlleleTables(x1)

test_that("each member of allele_divtables is a divtable and table", {
    expect_true(all(sapply(x1.adt, class)[1, ] == "divtable"))
    expect_true(all(sapply(x1.adt, class)[2, ] == "table"))
    expect_true(all(sapply(x1.adt, ncol) == 6))
    expect_true(all(sapply(x1.adt, nrow) == 1))
})

test_that("allele names and values correct in an allele_divtables", {
    expect_output(x1.adt, "$a", fixed = TRUE)
    expect_output(x1.adt, "$b", fixed = TRUE)
    expect_output(x1.adt, " +alleles")
    expect_output(x1.adt, "pop +11 +12 +13 +14 +15 +16")
    expect_output(x1.adt, "snurf +1 +1 +1 +1 +1 +1")
    expect_output(x1.adt, "pop +101 +102 +103 +104 +105 +106")
})

test_that("missing data reported correctly with quiet = FALSE", {
    expect_output(aal <- createAlleleTables(Qagr_adult_genotypes, quiet = FALSE), "Excluding 92 entries based on")
})

context("Testing as.allele_divtables() as createAlleleTables() synonym")

xx1.adt <- as.allele_divtables(x1)

test_that("object format same with as.allele_divtables()", {
    expect_true(all(sapply(xx1.adt, class)[1, ] == "divtable"))
    expect_true(all(sapply(xx1.adt, class)[2, ] == "table"))
    expect_true(all(sapply(xx1.adt, ncol) == 6))
    expect_true(all(sapply(xx1.adt, nrow) == 1))
})

test_that("allele names and values correct in an allele_divtables", {
    expect_output(xx1.adt, "$a", fixed = TRUE)
    expect_output(xx1.adt, "$b", fixed = TRUE)
    expect_output(xx1.adt, " +alleles")
    expect_output(xx1.adt, "pop +11 +12 +13 +14 +15 +16")
    expect_output(xx1.adt, "snurf +1 +1 +1 +1 +1 +1")
    expect_output(xx1.adt, "pop +101 +102 +103 +104 +105 +106")
})


