TODO for dispersalDiversity
---------------------------

* Define `print` and `plot` methods for `allele_tables`
* Do we want `table` methods, or define a new class that is the same as a table but is specifically a diversity table?
* Class returned by functions that is plottable by plotAlphaTest... or is this plottable by `htest_boot`
* make `createAlleleTables` return matrices?  or perhaps I have settled on a table as the input class?
* finish pooled PMI description in `diversity`
* copy chunks of `diversity` docs to `diversityMultilocus` and `diversitySingleLocus`
* Note README includes several examples and more documentation on parameters and functions
* Incorporate weighted means and variances equations from Scofield et al. 2011
* sort out reverseTerms stuff
* reverseTerms especially with empirical PVAL calculation with pchisq
* make allele diversity function returns match `pmiDiversity` function returns
* expand documentation for `nielsenTransform`
* Get data permission from VLS
* Do I need to do a types x sites data structure, or is that a 2d table? Perhaps provide functions to convert tables and xtabs to this?
* Tests

Completed
---------

* Quantiles checked are now richer and symmetric
* Renamed `allele.createTableList` to S3 generic and method `createAlleleTables` and `createAlleleTables.genalex`
* Started streamlining with S3 classes.
* Incorporated random distance matrix creation in to README.
* Changes function name of `pmiDiversity` to `diversity`
* For data, added 2002 and 2004 granary assignments for *Q. lobata*
* Worked out how to deal with `library(RColorBrewer)` use in `membershipPlot`.  If `RColorBrewer` is available then it is used, with the new option `fill.palette = "Dark2"` selecting the palette.  If it is not available, `rainbow` is used.
* Reworked return value from `pmiDiversity` to return separate lists for `q`, `q.nielsen` and `r`

