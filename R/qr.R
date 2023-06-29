#' Clark2DT-Function
#'
#' @details
#'
#' @return what is returned by the function?
#'
#' @param par parameters U, P, N estimated within the Clark2DT function
#' @param x represents the distance to the nearest seed source. Must be numeric
#'
#' @details Mixture of Gaussian nuclei that produces tails that are not quite as long. Maximum at seed tree itself and cannot become 0 at x = 0.
Clark2dt <- function(x, par){
  U <- par[1]
  P <- par[2]
  N <- par[3]
  result <- N * ((P)/ (pi*U * ((1+((x^ 2) / U)) ^ (P + 1))))
  return(result)
}

#' lognormal-Function
#' @param par parameters a, σ, N estimated within the lognormal function
#' @param x represents the distance to the nearest seed source. Must be numeric
#'
#' @details The equation means log η_i is a quadratic function of log x_i with a maximum at log a-2σ^2.

lognormal <- function(x, par) {
  log.a <- par[1]
  σ <- exp(par[2])
  N <- par[3]
  result <- N * exp(-(log(x)-log.a)^2 / (2*σ^2)) / (x^2 * 2*pi * sqrt(2*pi*σ^2))
  return(result)
}


#'Function selection
#' param x represents the distance to the nearest seed source. Must be numeric
#' param par are parameters to be estimated

#S.functions <- function(x, par, fun){
#  fun <- match.arg(fun,
#                   c("Clark2dt",
#                     "lognormal"))
#  f <- get(fun)
#  result <- f(x, par)
#  return(result)
#}



