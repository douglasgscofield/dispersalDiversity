#' Transform vector of squared frequencies using the method of Nielsen et al. (2003)
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
#' @export nielsenTransform
#'
nielsenTransform <- function(q.gg, n.g)
{  
    (((q.gg * (n.g+1) * (n.g - 2)) + (3 - n.g)) / ((n.g - 1)^2))
}

