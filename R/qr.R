##############################################################################
#' lognormal.log
#' @param par Numeric vector with three elements representing log-transformed 
#' scale and shape parameters, given by the mean \eqn{a} and standard 
#' deviation \eqn{\sigma} of the underlying normal distribution, and the 
#' scaling \eqn{N}.
#' @param x represents the distance to the nearest seed source. Must be 
#' numeric.
#'
#' @details The spatial dispersal density, representing the probability 
#' density function, divided by \eqn{2\pi x}, of the distance of a seed from 
#' its source is here given by
#' \deqn{f(x) = \frac{1}{2\pi x^2 \sqrt{2\pi\sigma^2}} 
#'   e^{-\frac{1}{2\sigma^2}(\log(x/a))^2},}
#' see Greene and Johnson (1989), Stoyan and Wagner (2001). Thus, the 
#' distance is assumed to have the log-normal distribution such that the 
#' log-distance has a normal distribution with mean \eqn{a} and variance 
#' \eqn{\sigma^2}. Note that \eqn{\log f(x)} is a quadratic function of
#' \eqn{\log x} with a maximum at \eqn{\log a-2\sigma^2}.
#'
#' @references
#' Greene, D.F., Johnson, E.A., 1989. A model of wind dispersal of winged or 
#' plumed seeds. *Ecology* **70**(2), 339–347.
#' \doi{10.2307/1937538}
#'
#' Stoyan, D., Wagner, S., 2001. Estimating the fruit dispersion of 
#' anemochorous forest trees. *Ecol. Modell.* **145**, 35–47.
#' \doi{10.1016/S0304-3800(01)00385-4}

lognormal.log <- function(x, par) {
  log.a <- par[1]
  σ <- exp(par[2])
  N <- par[3]
  result <- N * exp(-(log(x)-log.a)^2 / (2*σ^2)) / (x^2 * 2*pi * sqrt(2*pi*σ^2))
  return(result)
}


##############################################################################
#' Clark2dt.log
#'
#' Clark2dt computes the value of the dispersal function from Clark et al. 
#' (1999) multiplied by \eqn{N}.
#'
#' @return Numeric vector of function values multiplied by \eqn{N}.
#'
#' @param par Numeric vector with three elements representing the log-
#' transformed parameters \eqn{u} and \eqn{p} and the scaling \eqn{N}.
#' @param x Numeric vector of distances the the nearest seed source.
#'
#' @details The spatial dispersal density, representing the probability 
#' density function, divided by \eqn{2\pi x}, of the distance of a seed from 
#' its source is here given by
#' \deqn{f(x) = \frac{p}{\pi u (1+x^2/u)^{p+1}},}
#' see Clark et al. (1999), and Austerlitz et al. (2004) with a different 
#' parameterization (\eqn{b=p+1}). This represents a mixture of Gaussian 
#' densities that produces tails that are not quite as long ???longer???. It 
#' has its maximum at the seed tree itself and cannot become \eqn{0} at 
#' \eqn{x=0}.
#'
#' @references
#' Clark, J.S., Silman, M., Kern, R., Macklin, E. and HilleRisLambers, J. 
#' (1999). Seed dispersal near and far: patterns across temperate and tropical 
#' forests. *Ecology* **80**, 1475–1494. 
#' \doi{10.1890/0012-9658(1999)080[1475:SDNAFP]2.0.CO;2}
#'
#' Austerlitz, F., Dick, C.W., Dutech, C., Klein, E.K., Oddou-Muratorio, S., 
#' Smouse, P.E. and Sork, V.L. (2004). Using genetic markers to estimate the 
#' pollen dispersal curve. *Molecular Ecology* **13**, 937–954. 
#' \doi{https://doi.org/10.1111/j.1365-294X.2004.02100.x}

Clark2dt.log <- function(x, par){
  a <- exp(par[1])
  P <- exp(par[2])
  N <- par[3]
  result <- N * P / (pi*a^2 * (1+(x/a)^2) ^ (P+1))
  return(result)
}


##############################################################################
#' exponential.power.log
#' @param par Numeric vector with three elements representing the log-
#' transformed scale and shape parameters \eqn{a} and \eqn{b} of the dispersal 
#' density, and a the scaling \eqn{N}.
#' @param x represents the distance to the nearest seed source. Must be 
#' numeric.
#'
#' @details The spatial dispersal density, representing the probability 
#' density function, divided by \eqn{2\pi x}, of the distance of a seed from 
#' its source is here given by
#' \deqn{f(x) = \frac{b}{2\pi a^2\Gamma(2/b)} e^{-(x/a)^b},}
#' see Austerlitz et al. (2004). This represents a fat-tailed distribution, 
#' and the function has its maximum at zero. As a result of its flexible 
#' shape, the exponential power distribution has been applied in a number of 
#' theoretical studies that address dispersal.
#'
#' @references
#' Austerlitz, F., Dick, C.W., Dutech, C., Klein, E.K., Oddou-Muratorio, S., 
#' Smouse, P.E. and Sork, V.L. (2004). Using genetic markers to estimate the 
#' pollen dispersal curve. *Molecular Ecology* **13**, 937–954. 
#' \doi{https://doi.org/10.1111/j.1365-294X.2004.02100.x}

exponential.power.log <- function(x, par) {
  a <- exp(par[1])
  b <- exp(par[2])
  N <- par[3]
  N * b / (2*pi*a^2*gamma(2/b)) * exp(-(x/a)^b)
}


##############################################################################
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


##############################################################################
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


##############################################################################
##' Clark2DT-Function
##'
##' @details
##'
##' @return what is returned by the function?
##'
##' @param par parameters U, P, N estimated within the Clark2DT function
##' @param x represents the distance to the nearest seed source. Must be numeric
##'
##' @details Mixture of Gaussian nuclei that produces tails that are not quite as long. Maximum at seed tree itself and cannot become 0 at x = 0.
#Clark2dt <- function(x, par){
#  U <- par[1]
#  P <- par[2]
#  N <- par[3]
#  result <- N * ((P)/ (pi*U * ((1+((x^ 2) / U)) ^ (P + 1))))
#  return(result)
#}
#
##'Function selection
##' param x represents the distance to the nearest seed source. Must be numeric
##' param par are parameters to be estimated
#
#S.functions <- function(x, par, fun){
#  fun <- match.arg(fun,
#                   c("Clark2dt",
#                     "lognormal"))
#  f <- get(fun)
#  result <- f(x, par)
#  return(result)
#}


