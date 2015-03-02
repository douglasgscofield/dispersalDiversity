TODO for dispersalDiversity
---------------------------

* Where should I get data from for examples?  Perhaps the pericarp data from readGenalex?
* Do I need to do a types x sites data structure, or is that a 2d table? Perhaps provide functions to convert tables and xtabs to this?
* Tests

Completed
---------

* Worked out how to deal with `library(RColorBrewer)` use in `membershipPlot`.  If `RColorBrewer` is available then it is used, otherwise it is not.

Create a random distance matrix
------

```R
n.sites <- 5
n.sources <- 10
n.samples <- 1000
# data frame of site-source pairs
t <- data.frame(site = sample(n.sites, n.samples, replace = TRUE),
                source = round(runif(n.samples) * n.sources + 0.5))

# site-by-source matrix
m1 <- do.call(rbind, lapply(split(t, t$site), function(x) table(x$source)))
# this creates a class c("matrix")

# or use xtabs
m2 <- xtabs(data = t)
# this creates a class c("xtabs","table")
```

Both the above methods create a matrix with rownames of site levels, and column names of source levels.  `xtabs` also names the dims after the variables (`site` and `source`).
