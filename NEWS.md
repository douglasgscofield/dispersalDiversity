0.9.9000 (development version, branch make_R_package)
------
* Diversity tests return class `diversity_test`, which has `print` and `plot` methods
* Gamma accumulation functions return class `gamma_accum`, which has a `plot` method
* `diversity`, `alleleDiversityTest`, `alleleContrastTest`, and `alleleContrastTest3` have methods for `divtable` and `allele_divtables`
* Lists of allele diversity tables have class `allele_divtables`
* Diversity site-by-source tables have class `divtable`, with methods to convert from `table`, `xtabs`, `matrix`, and `data.frame`
* `diversity` returns separate lists for diversity values calculated using the standard diversity calculations (`q`), classical corrections (`r`), and Nielsen et al. (2003 *Molecular Ecology*) corrections (`q.nielsen`).
* `membershipPlot` with `fill.method = "color"` will use the `RColorBrewer` package if it is available; the default palette is `"Dark2"` and this can be changed with the `fill.palette =` option.  Also, `"colour"` is a `fill.method` synonym for `"color"`.
* The datasets `granaries_2002_Qlob` and `granaries_2004_Qlob` are provided, which include assignments of *Quercus lobata* acorns harvested from acorn woodpecker granaries in 2002 and 2004 to seed source trees.
* Numerous functions have been renamed (e.g., `pmiDiversity` to `diversity`) and data structures reconfigured (e.g., the return value from `diversity`)


0.1
------

* Initial version available via Github