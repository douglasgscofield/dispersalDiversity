#' Transform vector of squared frequencies using the method of Nielsen et al. (2003)
#'
#' Inverses of squared frequencies are commonly used to estimate so-called
#' effective numbers of mates, alleles etc., but these estimates are subject
#' to bias especially at low sample sizes.  The common method of handling
#' low sample sizes (the \code{r} estimates provided here) involve a
#' subtraction of 1 from observed counts and produce no values for 1-sized
#' groups.  This function starts with squared frequencies and group sizes
#' and transforms them into approximately unbiased estimates using the
#' method described by Nielsen \emph{et al}. (2003).
#'
#' @param q.gg Vector of squared frequencies
#'
#' @param n.g Vector with group size for each element of \code{q.gg}
#'
#' @return \code{q.gg} vector with the Nielsen et al. transform applied
#'
#' @references
#'
#' Nielsen, R., Tarpy, D. R. and Reeve, H. K. (2003) Estimating effective
#' paternity number in social insects and the effective number of alleles in
#' a population.  \emph{Molecular Ecology} 12:3157-3164.
#'
#' @export
#'
nielsenTransform <- function(q.gg, n.g)
{
    (((q.gg * (n.g+1) * (n.g - 2)) + (3 - n.g)) / ((n.g - 1)^2))
}

