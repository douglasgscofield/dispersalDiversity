#' Table holding diversity data
#'
#' An object of class \code{divtable}, a diversity table, holds
#' site (rows) by membership (columns) count data.  Such data could
#' represent numbers of individuals observed of specific species,
#' numbers of seeds observed originating from different parental trees,
#' or, as a component of class \code{\link{allele_divtables}}, counts
#' of alleles observed for a specific genetic locus.  This is the basic
#' data object of the \code{\link{dispersalDiversity}} package,
#' and is accepted by the function \code{\link{diversity}}, which generates
#' descriptive statistics, and the functions \code{\link{alphaDiversityTest}},
#' \code{\link{gammaContrastTest}} and many others that test for differences
#' in the structure of diversity.
#'
#' A \code{divtable} can be assembled by hand, but more typically, diversity
#' data will be read into a \code{matrix} or \code{data.frame} or assembled
#' into a \code{table} or \code{xtabs} and then converted to a
#' \code{divtable} using \code{as.divtable}.  All functions in the
#' \code{\link{dispersalDiversity}} package that accept \code{divtable} also
#' have methods that accept \code{matrix}, \code{data.frame}, \code{table}
#' or \code{xtabs} and then convert these to \code{divtable} prior to
#' analysis.
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

# Now add plot.divtable which uses membershipPlot. Do we need print.divtable?
