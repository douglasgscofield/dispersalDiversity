#' List of divtables holding allele diversity data
#'
#' An object of class \code{allele_divtables} is a list of 
#' \code{\link{divtable}} objects, each representing allele count
#' data for a single genetic locus across sites/groups.  Row names are
#' the site/group names, while column names are the individual alleles.
#' This is the basic genetic data object of the
#' \code{\link{dispersalDiversity}} package, and is accepted by the
#' function \code{\link{diversity}}, which generates descriptive statistics,
#' and the functions \code{\link{alphaDiversityTest}},
#' \code{\link{gammaContrastTest}} and many others that test for differences
#' in the structure of genetic diversity.
#'
#' Assembling an \code{allele_divtables} by hand is labourious.  More
#' typically, genotype data will be in a class \code{\link{genalex}} object
#' and then converted to \code{allele_divtables} using
#' \code{\link{createAlleleTables}}.
#'
# @examples
#'
#' @name allele_divtables-class
#' @aliases allele_divtables
#'
NULL


# TODO: what might be a good way to print and plot allele_divtables?



#' Generate an \code{\link{allele_divtables}} object from a class \code{genalex} object
#'
#' S3 method to convert an object of class \code{genalex} to an object of
#' class \code{\link{allele_divtables}}.  a list of \code{\link{divtable}}
#' objects representing site-specific allele counts.  This is an S3 generic
#' so that other methods might be written to convert other genetic
#' formats.
#'
#' Another option for converting genotypes to \code{\link{allele_divtables}}
#' objects is to convert using \code{\link[readGenalex]{as.genalex}} from the
#' \href{http://cran.r-project.org/web/packages/readGenalex/index.html}{readGenalex}
#' package, if there is a method available to perform the conversion.
#'
#' Although missing alleles may be common in genotypic data, there is no
#' provision in \code{\link{diversity}} and other functions in this package
#' for missing data.  Missing alleles are recognised if they match one of the
#' values in \code{exclude}.  The numbers of missing alleles recognised is
#' reported if \code{quiet = FALSE}.
#'
#' @note \code{as.allele_divtables} is a synonym
#'
#' @param x Object of class \code{genalex} holding genotypes to be converted
#' @param exclude Values in \code{x} that indicate missing alleles, these are
#' excluded from the \code{divtable} entries for each locus
#' @param quiet If \code{TRUE}, report the number of missing alleles excluded
#' @return Object of class \code{\link{allele_divtable}}
#'
#' @examples
#'
#' ## Use genotype data from readGenalex package, already loaded
#' data(Qagr_pericarp_genotypes)
#' pal <- createAlleleTables(Qagr_pericarp_genotypes)
#' str(pal)
#' ## The divtable for the first locus
#' pal[1]
#' ## Some missing data
#' data(Qagr_adult_genotypes)
#' aal <- createAlleleTables(Qagr_adult_genotypes, quiet = FALSE)
#'
#' @export createAlleleTables as.allele_divtables as.allele_divtables.default as.allele_divtables.genalex
#'
#' @aliases createAlleleTables as.allele_divtables as.allele_divtables.default as.allele_divtables.genalex
#'
#' @name createAlleleTables
#'
NULL

createAlleleTables <- function(x, ...) UseMethod("createAlleleTables")



#' @rdname createAlleleTables
#'
#' @export
#'
createAlleleTables.default <- function(x, ...)
{
    stop("Cannot create allele divtables, perhaps ", deparse(substitute(x)),
         " be converted to class 'genalex'.  See '?readGenalex::as.genalex'.")
}



#' @rdname createAlleleTables
#'
#' @export
#'
createAlleleTables.genalex <- function(x, exclude = c(NA, "0"),
    quiet = TRUE, ...)
{
    ploidy <- attr(x, "ploidy")
    lc <- attr(x, "locus.columns")
    ln <- attr(x, "locus.names")
    population <- attr(x, "pop.title")
    ans <- list()
    ex <- list()
    for (il in 1:length(lc)) {
        alleles <- as.vector(unlist(x[, lc[il]:(lc[il] + ploidy - 1)]))
        ex[[ ln[il] ]] <- sum(alleles %in% exclude)
        pop <- rep(x[[population]], ploidy)
        ans[[ ln[il] ]] <- as.divtable(table(pop, alleles, exclude = exclude))
    }
    if (sum(unlist(ex)) && ! quiet)
        cat(sprintf("Excluding %d entries based on 'exclude = c(%s)'\n", 
                    sum(unlist(ex)), paste(collapse = ", ", exclude)))
    structure(ans, class = c('allele_divtables', 'list'))
}



as.allele_divtables <- function(x, ...) UseMethod("as.allele_divtables")



as.allele_divtables.default <- function(x, ...)
    createAlleleTables.default(x, ...)



as.allele_divtables.genalex <- function(x, ...)
    createAlleleTables.genalex(x, ...)
