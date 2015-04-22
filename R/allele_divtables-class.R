#' List of divtables holding allele diversity data
#'
#' An object of class \code{allele_divtables} is a list of 
#' \code{\link{divtable}} objects, each representing allele count
#' data for a single genetic locus across sites/groups.  This is the basic
#' genetic data object of the \code{\link{dispersalDiversity}} package,
#' and is accepted by the function \code{\link{diversity}}, which generates
#' descriptive statistics, and the functions \code{\link{alphaDiversityTest}},
#' \code{\link{gammaContrastTest}} and many others that test for differences
#' in the structure of genetic diversity.
#'
#' Assembling an \code{allele_divtables} by hand is labourious.  More
#' typically, genotype data will be present in a \code{\link{genalex}}
#' object and then converted to \code{allele_divtables} using
#' \code{\link{createAlleleTables}}.  Other genotype formats may be handled
#' if a facility for first converting them to class \code{\link{genalex}}
#' is provided by the \code{\link{readGenalex}} package or another tool.
#'
# @examples
#'
#' @name allele_divtables-class
#' @aliases allele_divtables
#'
NULL


# TODO: what might be a good way to print and plot allele_divtables?



#' Calculate an \code{\link{allele_divtables}} object from a class \code{genalex} object
#'
#' S3 method to convert an object of class \code{genalex} to an object of
#' class \code{\link{allele_divtables}}, a list of \code{\link{divtable}}
#' objects representing site-specific allele counts.
#' Each entry of the list is, for each locus, a
#' table of site-by-allele counts, with row names being the site names and
#' column names being the names given to the individual alleles.  This is a
#' generic so that other methods might be written to convert other genetic
#' formats.
#'
#' @examples
#'
#' # The workflow to calculate basic allelic diversity statistics:
#' #
#' # dat <- readGenalex("GenAlEx-format-file-of-genotypes.txt")
#' # gt <- createAlleleTables(dat)
#'
#' @export
#'
#' @name createAlleleTables
#'
NULL

createAlleleTables <- function(x, ...) UseMethod("createAlleleTables")



createAlleleTables.default <- function(x, ...)
{
    xn <- deparse(substitute(x))
    stop("Cannot create allele divtables.  Can ", xn,
         " be converted to class 'genalex'?  See readGenalex::as.genalex.")
}



# as.allele_table.genalex?  as.table_list.genalex?  S3 method
#
# remove new.ploidy argument, see orig code below
#
#createAlleleTables.genalex <- function(dat, new.ploidy,
#                                   collapse.alleles = TRUE,
#                                   exclude = c(NA, "0"), quiet = FALSE) {
#    if (! is.genalex(dat))
#        stop("input must be class 'genalex'")
#    if (new.ploidy >= 2 && ! collapse.alleles)
#        stop("Must collapse ploidy when new.ploidy >= 2")
#    if (attr(dat, "ploidy") > new.ploidy)
#        dat <- reduceGenalexPloidy(dat, new.ploidy)

#' @rdname createAlleleTables
#'
#' @export
#'
createAlleleTables.genalex <- function(x, new.ploidy = 2,
    collapse.alleles = TRUE, exclude = c(NA, "0"), quiet = FALSE)
{
    if (new.ploidy >= 2 && ! collapse.alleles)
        stop("Must collapse ploidy when new.ploidy >= 2")
    if (attr(dat, "ploidy") > new.ploidy)
        dat <- reducePloidy(dat, new.ploidy)
    lc <- attr(dat, "locus.columns")
    ln <- attr(dat, "locus.names")
    population <- attr(dat, "pop.title")
    ans <- list()
    ex <- list()
    for (il in 1:length(lc)) {
        alleles <- as.vector(unlist(dat[, lc[il]:(lc[il] + new.ploidy - 1)]))
        ex[[ ln[il] ]] <- sum(alleles %in% exclude)
        pop <- rep(dat[[population]], new.ploidy)
        ans[[ ln[il] ]] <- t(as.matrix(table(alleles, pop, exclude = exclude)))
    }
    if (sum(unlist(ex)) && !quiet)
        cat(sprintf("Excluding %d entries based on 'exclude = c(%s)'\n", 
                    sum(unlist(ex)), paste(collapse = ", ", exclude)))
    class(ans) <- c('allele_divtables', 'list')
    return(ans)
}



