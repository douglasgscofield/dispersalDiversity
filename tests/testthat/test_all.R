library(dispersalDiversity)


#########################################
context("Setting up some test data")

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
m1 <- table(t)
d1 <- diversity(m1)

test_that("results with various table formats are identical", {
    expect_equal(d1, diversity(table(t)))
    mat <- do.call(rbind, lapply(split(t, t$site), function(x) table(x$source)))
    expect_equal(d1, diversity(mat))
    xt <- xtabs(data = t)
    expect_equal(d1, diversity(xt))
})


# create "identity matrix", one species per site
# create patchy matrix

context("Testing diversity.table()")
context("Testing diversity.allele_tables()")
