#' @include divtable-class.R
NULL

#' List of divtables holding allele diversity data
#'
#' An object of class \code{allele_divtables} is a list of
#' \code{\link{divtable}} objects, each representing sites-by-allele counts
#' data for a single genetic locus.  Row and column names are
#' the site names and individual alleles, respectively.
#' This is the basic data object for analysis of genetic data using the
#' \code{\link{dispersalDiversity}} package.  It is accepted by the
#' function \code{\link{diversity}}, which generates descriptive statistics,
#' and the functions \code{\link{alphaDiversityTest}},
#' \code{\link{gammaContrastTest}} and others that test for differences
#' in the structure of genetic diversity.
#'
#' Assembling an \code{allele_divtables} by hand is labourious.  More
#' typically, genotype data will be in a class \code{\link{genalex}} object
#' and then converted to \code{\link{allele_divtables}} using
#' \code{\link{createAlleleTables}}.  The \code{\link{createAlleleTables}}
#' and \code{\link{as.allele_divtables}} functions will attempt to convert
#' non-\code{\link{genalex}} objects using to \code{\link{genalex}} format
#' using \code{\link[readGenalex]{as.genalex}}, but this may not be successful.
#'
#' @examples
#'
#' ## One possible way to plot an \code{allele_divtables} object, plotting
#' ## each locus as a separate divtable in a two-column format
#'
#' data(Qagr_pericarp_genotypes)
#' pal <- createAlleleTables(Qagr_pericarp_genotypes)
#' par(mfrow = c(round(length(pal) / 2), 2))
#' lapply(names(pal), function(n) plot(pal[[n]], main = n, l2 = NULL, las = 2))
#'
#' @name allele_divtables-class
#'
#' @aliases allele_divtables
#'
NULL



#' Generate an allele_divtables object from a class genalex object
#'
#' S3 method to convert an object of class \code{genalex} to an object of
#' class \code{\link{allele_divtables}}.  a list of \code{\link{divtable}}
#' objects representing sites-by-allele counts.  This is an S3 generic
#' so that other methods might be written to convert other genetic
#' formats.
#'
#' If \code{x} is not of class \code{genalex}, an attempt is made to convert
#' it to class \code{genalex} using \code{\link[readGenalex]{as.genalex}}
#' from the
#' \href{http://cran.r-project.org/web/packages/readGenalex/index.html}{readGenalex}
#' package.  An error will be produced if \code{x} is of a class or format
#' that cannot be converted to class \code{genalex}.
#'
#' Another option for converting genotypes to \code{\link{allele_divtables}}
#' objects is to convert to one of the formats recognised by
#' \code{\link[readGenalex]{as.genalex}}.
#'
#' Although missing alleles may be common in genotypic data, there is no
#' provision in \code{\link{diversity}} and other functions in this package
#' for missing data.  Missing alleles are recognised and excluded if they
#' match one of the values in \code{exclude}.  The numbers of missing alleles
#' recognised is reported if \code{quiet = FALSE}.
#'
#' @note \code{as.allele_divtables} is a synonym, unless \code{x} is of class
#' \code{list}.  If so, if the class of each element of \code{x} is of class
#' \code{\link{divtable}}, then the class of \code{x} is changed to
#' \code{\link{allele_divtable}}.  If \code{x} is of class
#' \code{allele_divtables} it is returned unchanged.
#'
#' @param x Object of class \code{genalex} holding genotypes to be converted,
#' or of a class and format that can be converted to \code{genalex} using
#' \code{\link[readGenalex]{as.genalex}}
#'
#' @param exclude Values in \code{x} that indicate missing alleles, these are
#' excluded from the \code{divtable} entries for each locus
#'
#' @param quiet If \code{TRUE}, report the number of missing alleles excluded
#'
#' @return Object of class \code{\link{allele_divtable}}
#'
#' @examples
#'
#' ## Use genotype data from readGenalex package, already loaded
#' data(Qagr_pericarp_genotypes)
#' pal <- createAlleleTables(Qagr_pericarp_genotypes)
#' str(pal)
#'
#' ## The divtable for the first locus
#' pal[[1]]
#'
#' ## allele_divtables removes and can report missing data
#' data(Qagr_adult_genotypes)
#' aal <- createAlleleTables(Qagr_adult_genotypes, quiet = FALSE)
#'
#' @export createAlleleTables as.allele_divtables as.allele_divtables.default as.allele_divtables.genalex as.allele_divtables.list as.allele_divtables.allele_divtables
#'
#' @aliases createAlleleTables as.allele_divtables as.allele_divtables.default as.allele_divtables.genalex as.allele_divtables.list as.allele_divtables.allele_divtables
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
    if (inherits(x, 'data.frame') || inherits(x, 'loci')) {
        x <- as.genalex(x)
        return(createAlleleTables.genalex(x))
    } else {
        stop("Cannot convert to class 'allele_divtables', perhaps ",
             deparse(substitute(x)), " can be converted to class 'genalex'?",
             " See '?readGenalex::as.genalex'.")
    }
}



#' @rdname createAlleleTables
#'
#' @export
#'
createAlleleTables.allele_divtables <- function(x, ...)
    x



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


#----------------------------------
#
# synonyms, documented and exported in createAlleleTables above

as.allele_divtables <- function(x, ...) UseMethod("as.allele_divtables")

as.allele_divtables.default <- function(x, ...)
    createAlleleTables.default(x, ...)

as.allele_divtables.genalex <- function(x, ...)
    createAlleleTables.genalex(x, ...)

# These are not synonyms but are still documented in createAlleleTables above
#
as.allele_divtables.list <- function(x, ...)
{
    if (all(sapply(x, inherits, 'divtable')))
        structure(x, class('allele_divtables', 'list'))
    else stop(deparse(substitute(x)),
              " cannot be converted to class allele_divtables,",
              " all members must be class 'divtable'")
}

as.allele_divtables.allele_divtables <- function(x, ...)
    x

