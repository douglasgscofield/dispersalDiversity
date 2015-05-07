#' Table holding diversity data
#'
#' An object of class \code{divtable}, a diversity table, holds site (rows) by
#' membership (columns) count data expressed using \code{numeric} values.  Such
#' data could represent numbers of individuals observed of specific species,
#' numbers of seeds observed originating from different parental trees, or, as a
#' component of class \code{\link{allele_divtables}}, counts of alleles observed
#' for a specific genetic locus.  This is the basic data object of the
#' \code{\link{dispersalDiversity}} package, and is accepted by the function
#' \code{\link{diversity}}, which generates descriptive statistics, and the
#' functions \code{\link{alphaDiversityTest}}, \code{\link{gammaContrastTest}}
#' and many others that test for differences in the structure of diversity.
#'
#' A \code{divtable} can be assembled by hand, but more typically, diversity
#' data will be read into a \code{matrix} or \code{data.frame} or assembled
#' into a \code{table} or \code{xtabs} and then converted to a
#' \code{divtable} using \code{as.divtable}.  All functions in the
#' \code{\link{dispersalDiversity}} package that accept \code{divtable} also
#' have methods that accept \code{matrix}, \code{data.frame}, \code{table}
#' or \code{xtabs} and then convert these to \code{divtable} prior to
#' analysis using \code{\link{as.divtable}}.
#'
#' @examples
#'
#' ## Generate divtables, starting with random site-species pairs with a
#' ## log-normal distribution of abundances
#' set.seed(42)
#' t <- data.frame(site = sample(n.sites, n.samples, replace = TRUE),
#'                 source = round(rlnorm(n.samples, .1, 1) + 0.5))
#' head(t)
#' ## Using table
#' as.divtable(table(t))
#' ## Using xtabs
#' as.divtable(xtabs(~ site + species, data = t))
#' ## More complicated, create 'matrix' intermediary
#' as.divtable(do.call(rbind, lapply(split(t, t$site),
#'                                   function(x) table(x$source))))
#'
#' @name divtable-class
#' @aliases divtable
#'
NULL



#' Plot \code{divtable} object using \code{membershipPlot}
#'
#' Plot object of class \code{\link{divtable}} using 
#' \code{\link{membershipPlot}}.  This function specifies
#' \code{method = "bar"} and sets \code{l2} to be the site (row) sums of
#' \code{x}.  If you desire other values for these parameters, simply use
#' \code{\link{membershipPlot}} directly.  
#'
#' @param x     Object of class \code{\link{divtable}}
#' @param l2    \code{l2} argument given to \code{\link{membershipPlot}}
#' for column headings
#' @param \dots Additional arguments passed to \code{\link{membershipPlot}}
#'
#' @export
#'
plot.divtable <- function(x, l2 = paste("N =", rowSums(x)), ...) {
    membershipPlot(x, l2 = l2, ...)
}



#' Convert to class \code{'divtable'}
#'
#' @name as.divtable
#' @export
NULL

as.divtable <- function(x, ...) UseMethod("as.divtable")



#' @rdname as.divtable
#'
#' @export
#'
as.divtable.table <- function(x, ...)
{
    x <- as.numeric(x)
    structure(x, class = c('divtable', 'table'))
}



#' @rdname as.divtable
#'
#' @export
#'
as.divtable.matrix <- function(x, ...)
{
    x <- as.numeric(x)
    structure(as.table(x), class = c('divtable', 'table'))
}



#' @rdname as.divtable
#'
#' @export
#'
as.divtable.data.frame <- function(x, ...)
{
    x <- as.numeric(x)
    structure(as.table(as.matrix(x)), class = c('divtable', 'table'))
}



#' @rdname as.divtable
#'
#' @export
#'
as.divtable.xtabs <- function(x, ...)
{
    x <- as.numeric(x)
    structure(as.table(x), class = c('divtable', 'table'))
}

