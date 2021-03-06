#' Table holding diversity data
#'
#' An object of class \code{divtable}, a diversity table, holds site (rows) by
#' membership (columns) count data expressed using \code{numeric} values.  Such
#' data could represent (for example) numbers of individuals observed of
#' specific species, numbers of seeds observed originating from different
#' parental trees, or, as a component of class \code{\link{allele_divtables}},
#' counts of alleles observed for a specific genetic locus.  This is the basic
#' data object of the \code{\link{dispersalDiversity}} package, and is
#' accepted by the function \code{\link{diversity}}, which generates
#' descriptive statistics, and the functions \code{\link{alphaDiversityTest}},
#' \code{\link{gammaContrastTest}} and many others that test for differences
#' in the structure of diversity.
#'
#' A \code{divtable} can be assembled by hand, but more typically
#' data will be read into a \code{matrix} or \code{data.frame} or assembled
#' into a \code{table} or \code{xtabs} and then converted to a
#' \code{divtable} using \code{as.divtable}.  All functions in the
#' \code{\link{dispersalDiversity}} package that accept \code{divtable} also
#' have a default method that will convert \code{matrix}, \code{data.frame},
#' \code{table} and \code{xtabs} to \code{divtable} using
#' \code{\link{as.divtable}}.
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
#'
#' ## Using xtabs
#' as.divtable(xtabs(~ site + species, data = t))
#'
#' ## More complicated, create 'matrix' intermediary
#' as.divtable(do.call(rbind, lapply(split(t, t$site),
#'                                   function(x) table(x$source))))
#'
#' @name divtable-class
#'
#' @aliases divtable
#'
NULL



#' Plot divtable object using membershipPlot
#'
#' Plot object of class \code{\link{divtable}} using
#' \code{\link{membershipPlot}}.  This function sets \code{l2} to be the
#' row sums of \code{x}.  Other arguments will be passed to
#' \code{\link{membershipPlot}}.
#'
#' @param x     Object of class \code{\link{divtable}}
#'
#' @param l2    \code{l2} argument given to \code{\link{membershipPlot}}
#' for column headings
#'
#' @param \dots Additional arguments passed to \code{\link{membershipPlot}}
#'
#' @export
#'
plot.divtable <- function(x, ..., l2 = paste("N =", rowSums(x)), main = deparse(substitute(x))) {
    membershipPlot(x, ..., l2 = l2, main = main)
}



#' Convert to class divtable
#'
#' @name as.divtable
#'
#' @export
NULL

as.divtable <- function(x, ...) UseMethod("as.divtable")



#' @rdname as.divtable
#'
#' @export
#'
as.divtable.table <- function(x, ...)
{
    mode(x) <- "numeric"
    .divtableCheckRowSums(x)
    structure(x, class = c('divtable', 'table'))
}



#' @rdname as.divtable
#'
#' @export
#'
as.divtable.matrix <- function(x, ...)
{
    mode(x) <- "numeric"
    .divtableCheckRowSums(x)
    structure(as.table(x), class = c('divtable', 'table'))
}



#' @rdname as.divtable
#'
#' @export
#'
as.divtable.data.frame <- function(x, ...)
{
    x <- as.matrix(x)
    mode(x) <- "numeric"
    .divtableCheckRowSums(x)
    structure(as.table(x), class = c('divtable', 'table'))
}



#' @rdname as.divtable
#'
#' @export
#'
as.divtable.xtabs <- function(x, ...)
{
    mode(x) <- "numeric"
    .divtableCheckRowSums(x)
    if (! is.null(attr(x, "call")))
        attr(x, "call") <- NULL
    structure(as.table(x), class = c('divtable', 'table'))
}



#' @rdname as.divtable
#'
#' @export
#'
as.divtable.default <- function(x, ...)
{
    stop(deparse(substitute(x)), " cannot be converted to class divtable, ",
         "must be class table, matrix, xtabs or data.frame")
}



.divtableCheckRowSums <- function(x)
{
    if (any(rowSums(x) == 1)) {
        warning("row(s) ", paste(names(which(rowSums(x) == 1))),
                " contain 1 or fewer items and cannot be included in diversity tests",
                call. = FALSE)
    }
}

