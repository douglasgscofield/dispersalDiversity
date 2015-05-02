TODO for dispersalDiversity
---------------------------

* running diversity(granaries_2002_Qlob) produces NaN for L0049-2002, probably because row sum is 1; handle this better
* Clean up pairwise diversity tests, return values and print and plot methods
* Make sites/groups/species terminology consistent
* Consider adding other methods to `createAlleleTables` using the `as.genalex` functionality available to make conversions
* Define `print` and `plot` methods for `allele_divtables`
* Am I happy with `print.divtable` and `plot.divtable` methods?
* Separate doc for `diversity` class
* finish pooled PMI description in `diversity`
* copy chunks of `diversity` docs to `diversityMultilocus` and `diversitySingleLocus`
* Note README includes several examples and more documentation on parameters and functions
* Incorporate weighted means and variances equations from Scofield et al. 2011
* sort out reverseTerms stuff
* reverseTerms especially with empirical PVAL calculation with pchisq
* make allele diversity function returns match `pmiDiversity` function returns
* expand documentation for `nielsenTransform`
* Get data permission from VLS
* Tests

Completed
---------

* For alpha and gamma diversity tests, all return class `'diversity_test'` object which is handled by common `print.diversity_test` and `plot.diversity_test` methods
* New `as.divtable` generic, with methods `as.divtable.table`, `as.divtable.xtabs`, `as.divtable.matrix` and `as.divtable.data.frame`, the latter of which first converts to `matrix`, then to `divtable`
* Settled on new class `divtable`, which is shared with `table`
* For allelic data, new class `allele_divtables`, which is shared with `list` and is a list of `divtable` objects for each locus
* Quantiles checked are now richer and symmetric
* Renamed `allele.createTableList` to S3 generic and method `createAlleleTables` and `createAlleleTables.genalex`.  Added the synonym generic `as.allele_divtables`.
* Started streamlining with S3 classes.
* Incorporated random distance matrix creation in to README.
* Changes function name of `pmiDiversity` to `diversity`
* For data, added 2002 and 2004 granary assignments for *Q. lobata*
* Worked out how to deal with `library(RColorBrewer)` use in `membershipPlot`.  If `RColorBrewer` is available then it is used, with the new option `fill.palette = "Dark2"` selecting the palette.  If it is not available, `rainbow` is used.
* Reworked return value from `pmiDiversity` to return separate lists for `q`, `q.nielsen` and `r`

