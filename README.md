Dispersal Diversity : Statistics and Tests
------------------------------------------

R functions for calculating dispersal diversity statistics and making
comparisons involving dispersal diversity statistics (Scofield et al. in
review).  All functions take as input a simple data structure: a table of site
(rows) by source (columns) counts.  Though we originally developed the
diversity tests to understand seed dispersal in plant populations, the tests
themselves should be useful for any diversity data (such as biodiversity of
plant communities, etc.) that can be expressed with this same data structure.

All source files are required for performing diversity tests.  If all that is
desired are PMI (Grivet et al. 2005) and diversity (Scofield et al. in review)
statistics (<i>q<sub>gg</sub></i>, <i>&alpha;<sub>g</sub></i>, etc.) the source
file `pmiDiversity.R` contains the `pmiDiversity()` function that provides these
and this can be used separately.

Put all the source files in the same directory, and within your R session
simply

    source("diversityTests.R")


### pmiDiversity.R

Defines the R function `pmiDiversity()` which takes a site-by-source table and
produces statistics for Probability of Maternal Identity aka PMI (Grivet et al.
2005 Mol Ecol, Scofield et al. 2010 J Ecol, Scofield et al. 2011) and dispersal
diversity (Scofield et al. in review).  Three different PMI and diversity
statistics are calculated:

* <i>q<sub>gg</sub></i>-based, known to be biased (Grivet et al. 2005)

* <i>r<sub>gg</sub></i>-based, unbiased but poor performers at low sample sizes
  (Grivet et al. 2005, Scofield et al. in review)

* <i>q<sup>*</sup><sub>gg</sub></i>-based, which apply the transformation
  developed by Nielsen et al. (2003) to be unbiased and seem to perform well
(Scofield et al. 2010, Scofield et al. 2011, Scofield et al. in review).


### diversityTests.R

Defines several R functions which, like `pmiDiversity()`, take a site-by-source
table (one or more) and test diversity statistics within and among them.  See
Scofield et al. in review for methodological details.  The file pmiDiversity.R
(see above) is required to be in the same directory, as it provides functions
used here.

`alphaDiversityTest(tab)`
: Test for differences in alpha diversity among sites within a single dataset
 
`alphaContrastTest(tab.a, tab.b)`
: Test whether there is a difference in the alpha diversity between two datasets

`alphaContrastTest.3(tab.a, tab.b, tab.c)`
: Test whether there is a difference in the alpha diversity among three datasets

`pairwiseMeanTest(tab)`
: Test whether mean pairwise divergence/overlap among sites is different from the null espectation

`plotPairwiseMeanTest()`
: Plot the list returned from the above test for evaluation

`plotAlphaTest(result)`
: Plot the list returned from `alphaDiversityTest()` or `alphaContrastTest()` for evaluation

`gammaContrastTest(tab.a, tab.b)`
: Test whether there is a difference in the gamma diversity between two datasets

`gammaContrastTest.3(tab.a, tab.b, tab.c)`
: Test whether there is a difference in the gamma diversity among three datasets

