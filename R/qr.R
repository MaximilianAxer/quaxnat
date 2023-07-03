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

#' Clark2dt.log
#'
#' Clark2dt computes the value of the dispersal function from Clark et al. (1999) multiplied by \eqn{N}.
#'
#' @return Numeric vector of function values multiplied by \eqn{N}.
#'
#' @param par Numeric vector with three elements representing the log-transformed parameters \eqn{u} and \eqn{p} and the scaling \eqn{N}.
#' @param x Numeric vector of distances the the nearest seed source.
#'
#' @details Mixture of Gaussian nuclei that produces tails that are not quite as long. Maximum at seed tree itself and cannot become 0 at x = 0.


Clark2dt.log <- function(x, par){
  a <- exp(par[1])
  P <- exp(par[2])
  N <- par[3]
  result <- N * P / (pi*a^2 * (1+(x/a)^2) ^ (P+1))
  return(result)
}

#' exponential.power.log-Function
#' @param par parameters a, b, N estimated within the exponential.power.log function, with a and b representing the scale and the shape parameter of the function, respectively.
#' @param x represents the distance to the nearest seed source. Must be numeric.
#'
#' @details The equation means a fat-tailed distributions with the maximum at zero. As a result of its flexible shape, the exponential power distribution has been applied in a number of theoretical studies that address dispersal.

exponential.power.log <- function(x, par) {
  a <- exp(par[1])
  b <- exp(par[2])
  N <- par[3]
  N * b / (2*pi*a^2*gamma(2/b)) * exp(-(x/a)^b)
}


#' Weibull.log-Function
#' @param par parameters a, b, N estimated within the Weibull.log function with a and b representing the scale and the shape parameter of the function, respectively.
#' @param x represents the distance to the nearest seed source. Must be numeric.
#'
#' @details The distribution is fat-tailed when b ≤ 1 and thin-tailed otherwise. As for the exponential power function, when b = 2, the Weibull degenerates to the normal distribution, but when b = 1, it does not degenerate to the exponential distribution.

Weibull.log <- function(x, par) {
  a <- exp(par[1])
  b <- exp(par[2])
  N <- par[3]
  N * a^-b * b / (2*pi) * (x/a)^(b-2) * exp(-(x/a)^b)
}

#' Geometric.log-Function
#' @param par parameters a, b, N estimated within the geometric.log function with a and b representing the scale and the shape parameter of the function, respectively.
#' @param x represents the distance to the nearest seed source. Must be numeric.
#'
#' @details The function will behave quite differently from the exponential and Weibull distributions. They show a fat tail, whatever the value of the shape parameter (b), and the distributions become increasingly fat-tailed ?as b declines toward ‘1’.

geometric.log <- function(x, par) {
  a <- exp(par[1])
  b <- exp(par[2])
  N <- par[3]
  N * (b-2)*(b-1) / (2*pi*a) * (1+x/a)^-b
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



