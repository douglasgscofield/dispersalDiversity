TODO for dispersalDiversity
---------------------------

* make this an actual package
* come up with some tests
* submit to CRAN
* do I need to do a types x sites data structure, or is that a 2d table?
* perhaps provide functions to convert tables and xtabs to this?

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
m <- do.call(rbind, lapply(split(t, t$site), function(x) table(x$source)))
```
