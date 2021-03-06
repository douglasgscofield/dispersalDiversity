TODO for dispersalDiversity
---------------------------

* Continue on `diversityMultilocus` function returns match `diversity` function returns, also redo `diversitySingleLocus` returns
* Do we need `method = "permute"`? 
* Update figures with existing code
* continue generalising allele diversity tests
* Clean up pairwise diversity tests, return values and print and plot methods
* Make sites/groups/species terminology consistent
* Consider adding other methods to `createAlleleTables` using the `as.genalex` functionality available to make conversions
* copy chunks of `diversity` docs to `diversityMultilocus` and `diversitySingleLocus`
* Note README includes several examples and more documentation on parameters and functions
* Get data permission from VLS
* Tests

Completed
---------

* Added a `plot.divtable` method which produces an annotated `membershipPlot`
* `membershipPlot` annotation arguments `l1` and `l2` are placed by an `annotation` function which can be redefined by the user
* Weighted means and sample variances of `q.gg` for each estimator and of `q.gh`, all as defined in Scofield *et al*. 2011 *Oecologia*, are now included in the results returned by `diversity` for `divtable` objects
* For all diversity tests, the comparison of the log-likelihood value against an analytic &Chi;<sup>2</sup> distribution is removed as this was not an appropriate test
* Diversity tests now produce an error if any group has just one member, and `as.divtable` produces a warning for the same condition
* `membershipPlot` writes the plot into a PDF file on option (the EPS option has been removed)
* The pooled PMI values returned by `diversity` are documented
* The `divtable` and `allele_divtables` classes have separate documentation
* Export and document the `nielsenTransform` function
* Production of pie plots and the `method` argument are removed from `membershipPlot`
* The Gower distance matrix diagonal generation is much, much faster
* Removed `accum.method` and `distance.file` arguments from `gammaAccum` functions.  The `proximity` method was never well thought through; a proximity-based method could be valuable but not in the way it was implemented here.
* For alpha and gamma diversity tests, all return class `'diversity_test'` object which is handled by common `print.diversity_test` and `plot.diversity_test` methods
* New `as.divtable` generic, with methods `as.divtable.table`, `as.divtable.xtabs`, `as.divtable.matrix` and `as.divtable.data.frame`, the latter of which first converts to `matrix`, then to `divtable`
* Settled on new class `divtable`, which is shared with `table`
* For allelic data, new class `allele_divtables`, which is shared with `list` and is a list of `divtable` objects for each locus
* Quantiles checked are now richer and symmetric
* Renamed `allele.createTableList` to S3 generic and method `createAlleleTables` and `createAlleleTables.genalex`.  Added the synonym generic `as.allele_divtables`.
* Stabilised function interfaces using S3 classes.
* An example of random distance matrix creation is added to the README.
* Reworked return value from `diversity` to return separate lists for `q`, `q.nielsen` and `r`
* Changes function name of `pmiDiversity` to `diversity`
* For data, added 2002 and 2004 granary assignments for *Q. lobata*
* Worked out how to deal with `library(RColorBrewer)` use in `membershipPlot`.  If `RColorBrewer` is available then it is used, with the new option `fill.palette = "Dark2"` selecting the palette.  If it is not available, `rainbow` is used.

