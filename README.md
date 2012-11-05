Dispersal Diversity : Statistics and Tests
==========================================


R functions for calculating dispersal diversity statistics and making
comparisons involving dispersal diversity statistics (Scofield _et al._
[_American Naturalist_](http://www.jstor.org/stable/10.1086/668202).  All
functions take as input a simple data structure: a table of site (rows) by
source (columns) counts.  Though we originally developed the diversity tests to
understand seed dispersal in plant populations, the tests themselves should be
useful for any diversity data (such as biodiversity of plant communities, etc.)
that can be expressed with this same data structure.

The `pmiDiversity.R` and `diversityTest.r` source files are required for
performing diversity tests.  If all that is desired are PMI (Grivet _et al._
2005, Scofield _et al._ 2010, Scofield _et al._  2011) and diversity (Scofield _et
al._ <I>Am Nat</I>) statistics (<i>q<sub>gg</sub></i>,
<i>&alpha;<sub>g</sub></i>, etc.) the source file `pmiDiversity.R` contains the
`pmiDiversity()` function that provides these and this can be used separately.

Put all the source files in the same directory, and within your R session
simply

```R
source("diversityTests.R")
```

Additional source files are provided to perform other tasks.
`plotPairwiseMatrix.R` is available for plotting pairwise divergence/overlap
matrices.  More information is available below.  This file requires the
`pmiDiversity.R` source file to be available within the same directory.  Then
simply

```R
source("plotPairwiseMatrix.R")
```

`gammaAccum.R` is available for collecting &gamma; diversity accumulation
information and plotting this.  More information is avaialble below.  This file
requires the `pmiDiversity.R` source file to be available within the same
directory.  Then simply

```R
source("gammaAccum.R")
```


* * *

These statistical tools were developed in collaboration with Peter Smouse
(Rutgers University) and Victoria Sork (UCLA) and were funded by U.S. National
Science Foundation awards NSF-DEB-0514956 and NSF-DEB-0516529.

* * *

pmiDiversity.R
--------------

Defines the R function `pmiDiversity()` which takes a site-by-source table and
produces statistics for Probability of Maternal Identity aka PMI (Grivet _et al._
2005, Scofield _et al._ 2010, Scofield _et al._ 2011) and dispersal
diversity (Scofield _et al._ _Am Nat_).  Three different PMI and diversity
statistics are calculated:

* <i>q<sub>gg</sub></i>-based, known to be biased (Grivet _et al._ 2005)

* <i>r<sub>gg</sub></i>-based, unbiased but poor performers at low sample sizes
  (Grivet _et al._ 2005, Scofield _et al._ _Am Nat_)

* <i>q<sup>*</sup><sub>gg</sub></i>-based, which apply the transformation
  developed by Nielsen _et al._ (2003) to be unbiased and seem to perform well
(Scofield _et al._ 2010, Scofield _et al._ 2011, Scofield _et al._ _Am Nat_).


diversityTests.R
----------------

Defines several R functions which, like `pmiDiversity()`, take a site-by-source
table (one or more) and test diversity statistics within and among them.  See
(Scofield _et al._ _Am Nat_) for methodological details.  The file `pmiDiversity.R`
(see above) is required to be in the same directory, as it provides functions
used here.

`alphaDiversityTest(tab)`
: Test for differences in &alpha; diversity among sites within a single dataset
 
`alphaContrastTest(tab.a, tab.b)`
: Test whether there is a difference in the &alpha; diversity between two datasets

`alphaContrastTest.3(tab.a, tab.b, tab.c)`
: Test whether there is a difference in the &alpha; diversity among three datasets

`plotAlphaTest(result)`
: Plot the list returned from `alphaDiversityTest()` or `alphaContrastTest()` for evaluation

`pairwiseMeanTest(tab)`
: Test whether mean pairwise divergence/overlap among sites is different from the null espectation

`plotPairwiseMeanTest()`
: Plot the list returned from the above test for evaluation

`gammaContrastTest(tab.a, tab.b)`
: Test whether there is a difference in the &gamma; diversity between two datasets

`gammaContrastTest.3(tab.a, tab.b, tab.c)`
: Test whether there is a difference in the &gamma; diversity among three datasets


plotPairwiseMatrix.R
--------------------

Provides a function for plotting pairwise diversity matrices as returned by the
`pmiDiversity()` function, examples of which can be seen in Figure 4A-C of
Scofield _et al._ <I>Am Nat</I>.

`plotPairwiseMatrix()`
: Create a visual plot of pairwise divergence or overlap values as calculated by
`pmiDiversity()`

For example, with `tab` defined as above, plot the divergence matrix based on
<i>r<sub>gg</sub></i> calculations, labelling the axes "Seed Pool", using the
following code: 

```R
pmiD = pmiDiversity(tab)
plotPairwiseMatrix(pairwise.mat=pmiD$r.divergence.mat, 
                   pairwise.mean=pmiD$r.divergence, 
                   statistic="divergence", 
                   axis.label="Seed Pool")

```


gammaAccum.R
------------

Provides functions for calculating &gamma; accumulation across sites, and
plotting the result, examples of which can be seen in Figure 4D-F of Scofield
_et al._ <I>Am Nat</I>.  The file `pmiDiversity.R` (see above) is required to be in
the same directory, as it provides functions used here.

A typical workflow using these functions would be:

```R
rga.result = runGammaAccum(tab)
plotGammaAccum(rga.result)
```

#### runGammaAccum(tab)

Perform a &gamma; diversity accumulation on the site-by-source data in `tab`.
The result is returned in a list, which may be passed to `plotGammaAccum()` to
plot the result.  Several arguments control the method of accumulation and
value of &gamma; calculated.  Only the defaults have been tested; the others were
developed while exploring the data and must be considered experimental.

`tab` 
: Site-by-source table, same format as that passed to `pmiDiversity()`

`gamma.method` 
: Calculate &gamma; using `"r"` (default), `"q.nielsen"` or `"q"`
method (see paper)

`resample.method` 
: `"permute"` (default) or `"bootstrap"`; whether to resample
sites without (`"permute"`) or with (`"bootstrap"`) replacement

`accum.method` 
: `"random"` (default) or `"proximity"`.  If `proximity` is
used, then `distance.file` must be supplied

`distance.file` 
: A file or data.frame containing three columns of data, with
the header/column names being `pool`, `X`, and `Y`, containing the spatial
locations of the seed pools named in the row names of tab; only used with
`accum.method="proximity"`


#### plotGammaAccum(rga.result)

Create a visual plot of &gamma; accumulation results from `runGammaAccum()`.


#### Additional functions

The following functions typically won't be used separately, use `runGammaAccum()`
instead.

`gammaAccum()`
: Workhorse function for &gamma; accumulation

`gammaAccumStats()`
: Extracts stats from the result of `gammaAccum()`

`runGammaAccumSimple()`
: Wrapper that runs and then returns stats from `gammaAccum()`



References
----------

Scofield, D. G., P. E. Smouse, J. Karubian and V. L. Sork.  Accepted.  Use of
&alpha;, &beta;, and &gamma; diversity measures to characterize seed dispersal by animals.
[_American Naturalist_](http://www.jstor.org/stable/10.1086/668202), [data](http://dx.doi.org/10.5061/dryad.40kq7).

Scofield, D. G., V. R. Alfaro, V. L. Sork, D. Grivet, E. Martinez, J. Papp, A.
R. Pluess _et al._ 2011. Foraging patterns of acorn woodpeckers (<I>Melanerpes
formicivorus</I>) on valley oak (<I>Quercus lobata</I> Née) in two California oak
savanna-woodlands. <I>Oecologia</I> 166:187-196.

Scofield, D. G., V. L. Sork, and P. E. Smouse. 2010. Influence of acorn
woodpecker social behaviour on transport of coast live oak (<I>Quercus agrifolia</I>)
acorns in a southern California oak savanna. <I>Journal of Ecology</I> 98:561-572.

Grivet, D., P. E. Smouse, and V. L. Sork. 2005. A novel approach to an old
problem: tracking dispersed seeds. <I>Molecular Ecology</I> 14:3585-3595.

Nielsen, R., D. R. Tarpy, and H. K. Reeve. 2003. Estimating effective paternity
number in social insects and the effective number of alleles in a population.
<I>Molecular Ecology</I> 12:3157-3164.

