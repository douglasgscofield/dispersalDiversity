TODO for dispersalDiversity
---------------------------

* Note README includes several examples and more documentation on parameters and functions
* Incorporate weighted means and variances equations from Scofield et al. 2011
* sort out reverseTerms stuff
* rename `allele.createTableList` and all the stuff in allelePmiDiversity.R
* make `allele.createTableList` an S3 generic that operates on class `genalex`
* make allele diversity function returns match `pmiDiversity` function returns
* expand documentation for `pmiDiversity` and `nielsenTransform`
* Get data permission from VLS
* Do I need to do a types x sites data structure, or is that a 2d table? Perhaps provide functions to convert tables and xtabs to this?
* Tests

Completed
---------

* Started streamlining with S3 classes.
* Incorporated random distance matrix creation in to README.
* Changes function name of `pmiDiversity` to `diversity`
* For data, added 2002 and 2004 granary assignments for *Q. lobata*
* Worked out how to deal with `library(RColorBrewer)` use in `membershipPlot`.  If `RColorBrewer` is available then it is used, with the new option `fill.palette = "Dark2"` selecting the palette.  If it is not available, `rainbow` is used.
* Reworked return value from `pmiDiversity` to return separate lists for `q`, `q.nielsen` and `r`

