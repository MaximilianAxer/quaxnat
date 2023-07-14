##############################################################################
#' Dispersal Kernels For Log-Normal Distance Distributions
#'
#' `lognormal` computes the value, multiplied by \eqn{N}, of a dispersal
#' function based on seeds having a distance with a log-normal distribution
#' from the their source.
#'
#' @return Numeric vector of function values multiplied by \eqn{N}.
#'
#' @param par Numeric vector with three elements representing log-transformed
#' scale and shape parameters, given by the mean \eqn{a} and standard
#' deviation \eqn{\sigma} of the underlying normal distribution, and the
#' scaling \eqn{N}.
#' @param x represents the distance to the nearest seed source. Must be
#' numeric.
#'
#' @details The dispersal kernel, i.e. spatial probability density 
#' of the location of a seed relative to its source, is here given by
#' \deqn{k(x)={\Gamma (d/2) \over 
#'   2\pi ^{d/2}\left\|{x}\right\|^{d}\sqrt{2\pi \sigma ^{2}}}
#'   e^{-{1 \over 2\sigma ^{2}}(\log (\left\|{x}\right\|/a))^{2}},}
#' which corresponds to a probability density of the distance given by
#' \deqn{p(r)={1 \over r\sqrt{2\pi \sigma ^{2}}}
#'   e^{-{1 \over 2\sigma ^{2}}(\log (r/a))^{2}},}
#' where \eqn{d} is the spatial dimension and \eqn{\left\|{\,}\right\|} 
#' denotes the Euclidean norm; see Greene and Johnson (1989), Stoyan and 
#' Wagner (2001) for the planar case. Thus, the distance is assumed to have 
#' the log-normal distribution such that the log-distance has a normal 
#' distribution with mean \eqn{a} and variance \eqn{\sigma^2}. Note that 
#' \eqn{\log k(x)} is a quadratic function of \eqn{\left\|{x}\right\|} with 
#' a maximum at \eqn{\log a-d\sigma^2}.
#'
#' @details The spatial dispersal density, representing the probability
#' density function, divided by \eqn{2\pi x}, of the distance of a seed from
#' its source, is here given by
#' \deqn{f(x) = \frac{1}{2\pi x^2 \sqrt{2\pi\sigma^2}}
#'   e^{-\frac{1}{2\sigma^2}(\log(x/a))^2},}
#' see Greene and Johnson (1989), Stoyan and Wagner (2001). Thus, the
#' distance is assumed to have the log-normal distribution such that the
#' log-distance has a normal distribution with mean \eqn{a} and variance
#' \eqn{\sigma^2}. Note that \eqn{\log f(x)} is a quadratic function of
#' \eqn{\log x} with a maximum at \eqn{\log a-2\sigma^2}.
#'
#' Particularly suitable if the maximum regeneration density is not
#' directly at the seed source (e.g. Janzen-Conell effect)(Nathan et al. 
#' 2012)
#'
#' @references
#' Greene, D.F., Johnson, E.A. (1989). A model of wind dispersal of winged or
#' plumed seeds. *Ecology* **70**(2), 339–347.
#' \doi{10.2307/1937538}
#'
#' Stoyan, D., Wagner, S. (2001). Estimating the fruit dispersion of
#' anemochorous forest trees. *Ecol. Modell.* **145**, 35–47.
#' \doi{10.1016/S0304-3800(01)00385-4}
#'
#' Nathan, R., Klein, E., Robledo‐Arnuncio, J.J., Revilla, E. (2012).
#' Dispersal kernels: review, in Clobert, J., Baguette, M., Benton, T.G., 
#' Bullock, J.M. (eds.), *Dispersal ecology and evolution*, 186–210.
#' \doi{10.1093/acprof:oso/9780199608898.003.0015}

lognormal <- function(x, par) {
  log.a <- par[1]
  σ <- exp(par[2])
  N <- par[3]
  result <- N * exp(-(log(x)-log.a)^2/(2*σ^2)) / (2*pi * sqrt(2*pi*σ^2) * x^2)
  return(result)
}


##############################################################################
#' Dispersal Kernels From 2-Dimensional t Distribution
#'
#' `Clark2dt` computes the value of the dispersal function from Clark et al.
#' (1999) multiplied by \eqn{N}.
#'
#' @return Numeric vector of function values multiplied by \eqn{N}.
#'
#' @param par Numeric vector with three elements representing the
#' log-transformed parameters \eqn{a} and \eqn{p} and the scaling \eqn{N}.
#' @param x Numeric vector of distances to the nearest seed source.
#'
#' @details The spatial dispersal density, representing the probability
#' density function, divided by \eqn{2\pi x}, of the distance of a seed from
#' its source, is here given by
#' \deqn{f(x) = \frac{p}{\pi a^2 (1+(x/a)^2)^{p+1}},}
#' see Clark et al. (1999) and Austerlitz et al. (2004) (with
#' parameterizations \eqn{a=\sqrt{u}} and \eqn{p=b-1}, respectively). This
#' represents a mixture of Gaussian densities that exhibits heavier tails. It
#' has its maximum at zero.
#'
#' @references
#' Clark, J.S., Silman, M., Kern, R., Macklin, E., HilleRisLambers, J.
#' (1999). Seed dispersal near and far: patterns across temperate and tropical
#' forests. *Ecology* **80**, 1475–1494.
#' \doi{10.1890/0012-9658(1999)080[1475:SDNAFP]2.0.CO;2}
#'
#' Austerlitz, F., Dick, C.W., Dutech, C., Klein, E.K., Oddou-Muratorio, S.,
#' Smouse, P.E., Sork, V.L. (2004). Using genetic markers to estimate the
#' pollen dispersal curve. *Molecular Ecology* **13**, 937–954.
#' \doi{10.1111/j.1365-294X.2004.02100.x}

Clark2dt <- function(x, par){
  a <- exp(par[1])
  p <- exp(par[2])
  N <- par[3]
  result <- N * p / (pi*a^2 * (1+(x/a)^2) ^ (p+1))
  return(result)
}


##############################################################################
#' Dispersal Kernels From Exponential Power Family
#'
#' `exponential.power` computes the value, multiplied by \eqn{N}, of a 
#' dispersal kernel from the exponential power family, which includes, as
#' special cases, distance distributions based on normal and exponential
#' distributions.
#' `quax_exponential.power` computes the value, multiplied by \eqn{N}, of a 
#' dispersal kernel from an exponential power family, which includes, as
#' special cases, distance distributions based on normal and exponential
#' distributions.
#'
#' @return Numeric vector of function values multiplied by \eqn{N}.
#'
#' @param par Numeric vector with three elements representing the
#' log-transformed scale and shape parameters \eqn{a} and \eqn{b} of the
#' dispersal density, and a the scaling \eqn{N}.
#' @param x represents the distance to the nearest seed source. Must be
#' numeric.
#'
#' @details The dispersal kernel, i.e. spatial probability density 
#' of the location of a seed relative to its source, is here given by
#' \deqn{k(x)={b\Gamma (d/2) \over 2\pi ^{d/2}a^{d}\Gamma (d/b)}
#'    e^{-(\left\|{x}\right\|/a)^{b}},}
#' which corresponds to a probability density of the distance given by
#' \deqn{p(r)={b \over a^{d}\Gamma (d/b)}r^{d-1}e^{-(r/a)^{b}},}
#' where \eqn{d} is the spatial dimension and \eqn{\left\|{\,}\right\|} 
#' denotes the Euclidean norm; see Bateman (1947), Clark et al. (1998), 
#' Austerlitz et al. (2004), Nathan et al. (2012) for the planar case.
#' 
#' @details The spatial dispersal density, representing the probability
#' density function, divided by \eqn{2\pi x}, of the distance of a seed from
#' its source, is here given by
#' \deqn{f(x) = \frac{b}{2\pi a^2\Gamma(2/b)} e^{-(x/a)^b},}
#' see Bateman (1947), Clark et al. (1998), Austerlitz et al. (2004). This
#' function has its maximum at zero and represents a rather flexible family
#' of distributions including the classical bivariate Gaussian kernels, the
#' kernels based on an exponential distribution of distances, and, for
#' \eqn{b<1}, fat-tailed distributions. For \eqn{b>1}, it shows thin-tailed distributions.
#' It has consequently been applied in a number of theoretical studies that
#' address dispersal (Ribbens et al. 1994; Bullock et al. 2017).
#'
#' @references
#' Bateman, A. (1947). Contamination in seed crops: III. relation with
#' isolation distance. *Heredity* **1**, 303–336.
#' \doi{10.1038/hdy.1947.20}
#'
#' Ribbens, E., Silander Jr, J. A., & Pacala, S. W.  (1994). Seedling 
#' recruitment in forests: calibrating models to predict patterns of tree 
#' seedling dispersion. *Ecology* **75**, 1794-1806.
#' \doi{10.2307/1939638}
#'
#' Clark, J.S., Macklin, E., Wood, L. (1998). Stages and spatial scales of
#' recruitment limitation in southern Appalachian forests. *Ecological
#' Monographs* **68**(2), 213–235.
#' \doi{10.2307/2657201}
#'
#' Clark, J.S. (1998). Why trees migrate so fast: confronting theory with
#' dispersal biology and the paleorecord. *The American Naturalist*
#' **152**(2), 204–224.
#' \doi{10.1086/286162}
#'
#' Austerlitz, F., Dick, C.W., Dutech, C., Klein, E.K., Oddou-Muratorio, S.,
#' Smouse, P.E., Sork, V.L. (2004). Using genetic markers to estimate the
#' pollen dispersal curve. *Molecular Ecology* **13**, 937–954.
#' \doi{10.1111/j.1365-294X.2004.02100.x}
#'
#' #'Bullock, J. M., Mallada González, L., Tamme, R., Götzenberger, L., White, S. M., Pärtel, M., Hooftman, D. A.
#' (2017).  A synthesis of empirical plant dispersal kernels.
#' *Journal of Ecology* **105**, 6-19.
#' \doi{10.1111/1365-2745.12666}
#'
#' Nathan, R., Klein, E., Robledo‐Arnuncio, J.J., Revilla, E. (2012).
#' Dispersal kernels: review, in Clobert, J., Baguette, M., Benton, T.G., 
#' Bullock, J.M. (eds.), *Dispersal ecology and evolution*, 186–210.
#' \doi{10.1093/acprof:oso/9780199608898.003.0015}
#'

exponential.power <- function(x, par) {
  a <- exp(par[1])
  b <- exp(par[2])
  N <- par[3]
  N * b / (2*pi*a^2*gamma(2/b)) * exp(-(x/a)^b)
}

quax_exponential.power <- function(
    r = sqrt(.colSums(x^2,nrow(x),d)),
    x, par, N = 1, d = if (missing(x)) 2 else ncol(x)) {
  a <- exp(par[1])
  b <- exp(par[2])
  N * b * gamma(d/2) / (2*pi^(d/2)*a^d*gamma(d/b)) * exp(-(r/a)^b)
}


##############################################################################
#' Dispersal Kernels From Weibull Family
#'
#' `Weibull` computes the value of the dispersal function from Tufto et al.
#' (1997) multiplied by \eqn{N}.
#'
#' @return Numeric vector of function values multiplied by \eqn{N}.
#'
#' @param par Numeric vector with three elements representing the
#' log-transformed parameters \eqn{a} and \eqn{b} and the scaling \eqn{N}.
#' @param x Numeric vector of distances to the nearest seed source.
#'
#' @details The dispersal kernel, i.e. spatial probability density 
#' of the location of a seed relative to its source, is here given by
#' \deqn{k(x)={b\Gamma (d/2) \over 2\pi ^{d/2}a^{b}}\left\|{x}\right\|^{b-d}
#'   e^{-(\left\|{x}\right\|/a)^{b}},}
#' which corresponds to a probability density of the distance given by
#' \deqn{p(r)={b \over a^{b}}r^{b-1}e^{-(r/a)^{b}},}
#' where \eqn{d} is the spatial dimension and \eqn{\left\|{\,}\right\|} 
#' denotes the Euclidean norm; see Tufto et al. (1997) for the planar case.
#' Thus, the distance is assumed to have the Weibull distribution with scale 
#' parameter \eqn{a} and shape parameter \eqn{b}.
#' 
#' @details The spatial dispersal density, representing the probability
#' density function, divided by \eqn{2\pi x}, of the distance of a seed from
#' its source, is here given by
#' \deqn{f(x) = \frac{a^{-b}b}{2\pi} (\frac{x}{a})^{b-2} e^{-(x/a)^b},}
#' see Tufto et al. (1997), Austerlitz et al. (2004).
#' The distribution is fat-tailed when \eqn{b<1} and
#' thin-tailed otherwise (Nathan et al. 2012).
#' For \eqn{b>1}, the mode of the function is at \eqn{x>1}. In this way, the function approaches the normal distribution.
#'
#' @references
#' Tufto, J., Engen, S., Hindar, K. (1997). Stochastic dispersal processes in
#' plant populations, *Theoretical Population Biology* **52**(1), 16–26.
#' \doi{10.1006/tpbi.1997.1306}
#'
#' Austerlitz, F., Dick, C.W., Dutech, C., Klein, E.K., Oddou-Muratorio, S.,
#' Smouse, P.E., Sork, V.L. (2004). Using genetic markers to estimate the
#' pollen dispersal curve. *Molecular Ecology* **13**, 937–954.
#' \doi{10.1111/j.1365-294X.2004.02100.x}
#'
#' #' Nathan, R., Klein, E., Robledo‐Arnuncio, J.J., Revilla, E. (2012).
#' Dispersal kernels: review, in Clobert, J., Baguette, M., Benton, T.G., 
#' Bullock, J.M. (eds.), *Dispersal ecology and evolution*, 186–210.
#' \doi{10.1093/acprof:oso/9780199608898.003.0015}

Weibull <- function(x, par) {
  a <- exp(par[1])
  b <- exp(par[2])
  N <- par[3]
  N * a^-b * b / (2*pi) * (x/a)^(b-2) * exp(-(x/a)^b)
}


##############################################################################
#' Power-Law Dispersal Kernels
#'
#' `power` computes the value of the dispersal function from (WHERE?)
#' multiplied by \eqn{N}.
#'
#' @return Numeric vector of function values multiplied by \eqn{N}.
#'
#' @param par Numeric vector with three elements representing the
#' log-transformed parameters \eqn{a} and \eqn{b} and the scaling \eqn{N}.
#' @param x represents the distance to the nearest seed source. Must be
#' numeric.
#'
#' @details The dispersal kernel, i.e. spatial probability density 
#' of the location of a seed relative to its source, is here given by
#' \deqn{k(x)={\Gamma (d/2) \over 2\pi ^{d/2}a^{d}\Beta(d,b-d)}
#'   (1+{\left\|{x}\right\| \over a})^{-b},}
#' which corresponds to a probability density of the distance given by
#' \deqn{p(x)={1 \over a^{d}\Beta(d,b-d)}r^{d-1}(1+{r \over a})^{-b},}
#' where \eqn{d} is the spatial dimension and \eqn{\left\|{\,}\right\|} 
#' denotes the Euclidean norm; see Nathan et al. (2012) for the planar case.
#'
#' @details The spatial dispersal density, representing the probability
#' density function, divided by \eqn{2\pi x}, of the distance of a seed from
#' its source, is here given by
#' \deqn{f(x) = \frac{(b-2)(b-1)}{2\pi a^2} (1+\frac{x}{a})^{-b},}
#' see Nathan et al. (2012).
#' (CHANGE THE FOLLOWING, SHOULD BE OUR OWN CHARACTERIZATION:) Austerlitz 
#' et al. (2004) characterize it as follows: The geometric and 2dt families 
#' “will behave quite differently from the exponential and Weibull 
#' distributions. They show a fat tail, whatever the value of the shape 
#' parameter (\eqn{b}), and the distributions become increasingly fat-tailed 
#' as \eqn{b} declines toward ‘1’.”
#'
#' @references
#' Nathan, R., Klein, E., Robledo‐Arnuncio, J.J., Revilla, E. (2012).
#' Dispersal kernels: review, in Clobert, J., Baguette, M., Benton, T.G., 
#' Bullock, J.M. (eds.), *Dispersal ecology and evolution*, 186–210.
#' \doi{10.1093/acprof:oso/9780199608898.003.0015}
#'
#' Austerlitz, F., Dick, C.W., Dutech, C., Klein, E.K., Oddou-Muratorio, S.,
#' Smouse, P.E., Sork, V.L. (2004). Using genetic markers to estimate the
#' pollen dispersal curve. *Molecular Ecology* **13**, 937–954.
#' \doi{10.1111/j.1365-294X.2004.02100.x}

power <- function(x, par) {
  a <- exp(par[1])
  b <- exp(par[2])
  N <- par[3]
  N * (b-2)*(b-1) / (2*pi*a^2) * (1+x/a)^-b
}


##############################################################################
##' Clark2DT-Function
##'
##' @details
##'
##' @return what is returned by the function?
##'
##' @param par parameters U, P, N estimated within the Clark2DT function
##' @param x represents the distance to the nearest seed source. Must be
##' numeric
##'
##' @details Mixture of Gaussian nuclei that produces tails that are not
##' quite as long. Maximum at seed tree itself and cannot become 0 at x = 0.
#Clark2dt <- function(x, par){
#  U <- par[1]
#  P <- par[2]
#  N <- par[3]
#  result <- N * ((P)/ (pi*U * ((1+((x^ 2) / U)) ^ (P + 1))))
#  return(result)
#}
#
##'Function selection
##' param x represents the distance to the nearest seed source. Must be
##' numeric
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


