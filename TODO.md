TODO for dispersalDiversity
---------------------------

* Continue to straighten out .default cases for diversityTest*.R
* Sort out test differences
* That wierdness when converting to as.divtable for the Thamnolia diversity stuff
* This wierdness:

    mitchell-clustering> m = do.reads.matrix(subset(dat200, Site == "Iceland"), "Species")
    mitchell-clustering> m
               OTU
    Species     OTU 1 OTU 2 OTU 3 OTU 4 OTU 5 OTU 6 OTU 7
      Cetraria      5 26483    23     0     0     0     0
      Thamnolia    19 12637     0     0     0     0     0
    mitchell-clustering> alphaDiversityTest(m)
    Error in UseMethod("alphaDiversityTest") : 
      no applicable method for 'alphaDiversityTest' applied to an object of class "c('matrix', 'integer', 'numeric')"
    mitchell-clustering> as.divtable(m)
               OTU
    Species     OTU 1 OTU 2 OTU 3 OTU 4 OTU 5 OTU 6 OTU 7
      Cetraria      5 26483    23     0     0     0     0
      Thamnolia    19 12637     0     0     0     0     0
    mitchell-clustering> alphaDiversityTest(as.divtable(m))
    Alpha diversity test, contrast among sites in single data set
    ....

* Update figures with existing code
* straighten out eps and pdf file for membershipPlot
* continue generalising allele diversity tests
* running diversity(granaries_2002_Qlob) produces NaN for L0049-2002, probably because row sum is 1; handle this better
* Clean up pairwise diversity tests, return values and print and plot methods
* Make sites/groups/species terminology consistent
* Consider adding other methods to `createAlleleTables` using the `as.genalex` functionality available to make conversions
* Define `print` and `plot` methods for `allele_divtables`
* Am I happy with `print.divtable` and `plot.divtable` methods?
* finish pooled PMI description in `diversity`
* copy chunks of `diversity` docs to `diversityMultilocus` and `diversitySingleLocus`
* Note README includes several examples and more documentation on parameters and functions
* Incorporate weighted means and variances equations from Scofield et al. 2011
* sort out reverseTerms stuff
* reverseTerms especially with empirical PVAL calculation with pchisq
* make allele diversity function returns match `pmiDiversity` function returns
* Get data permission from VLS
* Tests

Completed
---------

* The `divtable` and `allele_divtables` classes have separate documentation
* Export and document the `nielsenTransform` function
* Production of pie plots and the `method` argument are removed from `membershipPlot`
* Optimisation of Gower distance matrix diagonal generation
* Removed `accum.method` and `distance.file` arguments from `gammaAccum` functions.  The `proximity` method was never well thought through; a proximity-based method could be valuable but not in the way it was implemented here.
* For alpha and gamma diversity tests, all return class `'diversity_test'` object which is handled by common `print.diversity_test` and `plot.diversity_test` methods
* New `as.divtable` generic, with methods `as.divtable.table`, `as.divtable.xtabs`, `as.divtable.matrix` and `as.divtable.data.frame`, the latter of which first converts to `matrix`, then to `divtable`
* Settled on new class `divtable`, which is shared with `table`
* For allelic data, new class `allele_divtables`, which is shared with `list` and is a list of `divtable` objects for each locus
* Quantiles checked are now richer and symmetric
* Renamed `allele.createTableList` to S3 generic and method `createAlleleTables` and `createAlleleTables.genalex`.  Added the synonym generic `as.allele_divtables`.
* Stabilised function interfaces using S3 classes.
* Incorporated random distance matrix creation in to README.
* Changes function name of `pmiDiversity` to `diversity`
* For data, added 2002 and 2004 granary assignments for *Q. lobata*
* Worked out how to deal with `library(RColorBrewer)` use in `membershipPlot`.  If `RColorBrewer` is available then it is used, with the new option `fill.palette = "Dark2"` selecting the palette.  If it is not available, `rainbow` is used.
* Reworked return value from `pmiDiversity` to return separate lists for `q`, `q.nielsen` and `r`

