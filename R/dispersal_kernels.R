# Internal functions:
# `rownorms` returns the Euclidean norms of the rows of a matrix x:
rownorms <- function(x,...) if (length(d<-dim(x))>1L && d[2L]>1)
  sqrt(.rowSums(x^2,d[1L],d[2L],...)) else abs(x)
# `surface` returns the volume (area, length) of the surface of a ball with
# radius r in d-dimensional space:
surface <- function(d,r=1) 2*pi^(d/2)/gamma(d/2) * r^(d-1)


##############################################################################
#' Dispersal kernels for log-normal distance distributions
#'
#' `k_lognormal` computes the value, multiplied by \eqn{N}, of a dispersal
#' kernel based on seeds having a distance with a log-normal distribution
#' from the their source.
#'
#' @param x Numeric matrix of positions \eqn{x} relative to the seed source,
#'  or vector of distances \eqn{\left\|{x}\right\|} to the seed source.
#' @param par Numeric vector with two elements representing log-transformed
#'  scale and shape parameters, given by the median distance \eqn{a} and by
#'  the variance \eqn{b} of the underlying normal distribution.
#' @param N The multiplier \eqn{N}.
#' @param d The spatial dimension.
#'
#' @details The dispersal kernel, i.e. spatial probability density
#' of the location of a seed relative to its source, is here given by
#' \deqn{k(x)={\Gamma (d/2) \over
#'   2\pi ^{d/2}\left\|{x}\right\|^{d}\sqrt{2\pi b}}
#'   e^{-{1 \over 2b}(\log (\left\|{x}\right\|/a))^{2}}
#'   ={\Gamma (d/2)e^{d^{2}b/2} \over 2\pi ^{d/2}a^{d}\sqrt{2\pi b}}
#'   e^{-{1 \over 2b}(\log {\left\|{x}\right\| \over a}+db)^{2}},}
#' which corresponds to a probability density of the distance given by
#' \deqn{p(r)={1 \over r\sqrt{2\pi b}}e^{-{1 \over 2b}(\log (r/a))^{2}}
#'   ={e^{b/2} \over a\sqrt{2\pi b}}
#'   e^{-{1 \over 2b}(\log {r \over a}+b)^{2}},}
#' where \eqn{d} is the spatial dimension, \eqn{\left\|{\,}\right\|}
#' denotes the Euclidean norm and the normalizing constant of the kernel
#' involves the \link[base:gamma]{gamma} function; see Greene and Johnson
#' (1989), Stoyan and Wagner (2001) for the planar case. Thus, the distance
#' is assumed to have the \link[stats:Lognormal]{log-normal distribution}
#' such that the log-distance has a normal distribution with mean
#' \eqn{\log a} and variance \eqn{b}. Here \eqn{\log k(x)} is a quadratic
#' function of \eqn{\log \left\|{x}\right\|} with a maximum at
#' \eqn{\log a-db}, while \eqn{\log p(r)} is a quadratic function of
#' \eqn{\log r} with a maximum at \eqn{\log a-b}.
#'
#' This kernel is particularly suitable if the maximum regeneration density
#' is not directly at the seed source (e.g. Janzen–Connell effect), cf.
#' Nathan et al. (2012).
#'
#' @return Numeric vector of function values \eqn{k(x)} multiplied by \eqn{N}.
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
#'
#' @export
#'
#' @examples
#' k_lognormal(2:5, par=c(0,0), d=2)

k_lognormal <- function(x, par, N=1, d=NCOL(x)) {
  r <- rownorms(x)
  log.a <- par[1]
  b <- exp(par[2])
  N / (surface(d)*sqrt(2*pi*b)) *
    exp(d^2*b/2 - d*log.a - (log(r)-log.a+d*b)^2/(2*b))
  # Alternatives for the last line (for nonzero r only):
  # N / (surface(d)*sqrt(2*pi*b)) * exp(-(log(r)-log.a)^2/(2*b)) / r^d
  # N / surface(d,r) * dlnorm(r,log.a,exp(par[2]/2))
}


##############################################################################
#' Dispersal kernels from spatial t distribution
#'
#' `k_t` computes the value, multiplied by \eqn{N}, of the dispersal kernel
#' from Clark et al. (1999) that represents a multivariate t distribution.
#'
#' @param x Numeric matrix of positions \eqn{x} relative to the seed source,
#'  or vector of distances \eqn{\left\|{x}\right\|} to the seed source.
#' @param par Numeric vector with two elements representing the
#'  log-transformed parameters \eqn{a} and \eqn{b}.
#' @param N The multiplier \eqn{N}.
#' @param d The spatial dimension.
#'
#' @details The dispersal kernel, i.e. spatial probability density
#' of the location of a seed relative to its source, is here given by
#' \deqn{k(x)={\Gamma ((b+d)/2) \over \pi ^{d/2}a^{d}\Gamma (b/2)}
#'   (1+{\left\|{x}\right\|^{2} \over a^{2}})^{-(b+d)/2},}
#' which corresponds to a probability density of the distance given by
#' \deqn{p(r)={2 \over a^{d}B(d/2,b/2)}r^{d-1}
#'   (1+{r^{2} \over a^{2}})^{-(b+d)/2},}
#' where \eqn{d} is the spatial dimension, \eqn{\left\|{\,}\right\|}
#' denotes the Euclidean norm and the normalizing constants involve the
#' \link[base:beta]{beta} and \link[base:gamma]{gamma} functions; see Clark
#' et al. (1999) and Austerlitz et al. (2004) for the planar case (with
#' \eqn{a,b} replaced by \eqn{\sqrt{u},2p} and
#' \eqn{a,2b-d}, respectively). This means the position is
#' \eqn{a \over \sqrt{b}} times a random vector having a standard
#' \eqn{d}-variate t distribution with \eqn{b} degrees of freedom (a standard
#' Gaussian vector divided by \eqn{\sqrt{z/b}}, where \eqn{z} is independent
#' and chi-squared distributed with \eqn{b} degrees of freedom), and the
#' squared distance is \eqn{da^{2} \over b} times a random variable having an
#' \link[stats:FDist]{F distribution} with \eqn{d} and \eqn{b} degrees of
#' freedom.
#'
#' This results from the kernel being defined as a mixture of Gaussian
#' kernels with an inverse variance having a
#' \link[stats:GammaDist]{gamma distribution} with shape parameter
#' \eqn{b\over 2} and inverse scale parameter \eqn{a^{2}\over 2}, which for
#' \eqn{a=1} is a \link[stats:Chisquare]{chi-squared distribution} with
#' \eqn{b} degrees of freedom.
#'
#' The dispersal kernel always has its maximum at zero, and the distance has
#' a fat-tailed distribution for all choices of \eqn{b}.
#'
#' @return Numeric vector of function values \eqn{k(x)} multiplied by \eqn{N}.
#'
#' @references
#' Clark, J.S., Silman, M., Kern, R., Macklin, E., HilleRisLambers, J.
#' (1999). Seed dispersal near and far: patterns across temperate and
#' tropical forests. *Ecology* **80**, 1475–1494.
#' \doi{10.1890/0012-9658(1999)080[1475:SDNAFP]2.0.CO;2}
#'
#' Austerlitz, F., Dick, C.W., Dutech, C., Klein, E.K., Oddou-Muratorio, S.,
#' Smouse, P.E., Sork, V.L. (2004). Using genetic markers to estimate the
#' pollen dispersal curve. *Molecular Ecology* **13**, 937–954.
#' \doi{10.1111/j.1365-294X.2004.02100.x}
#'
#' @export
#'
#' @examples
#' k_t(2:5, par=c(0,0), d=2)

k_t <- function(x, par, N=1, d=NCOL(x)) {
  r <- rownorms(x)
  a <- exp(par[1])
  b <- exp(par[2])
  N * 2 / (surface(d) * a^d * beta(d/2, b/2)) / (1+(r/a)^2)^((b+d)/2)
  # Alternatives for the last line (for nonzero r):
  # s<-b/(d*a^2); N * 2 * s / surface(d) * r^(2-d) * df(s*r^2, d, b)
  # s<-b/(d*a^2); N * 2 * s * r * df(s*r^2, d, b) / surface(d,r)
}


##############################################################################
#' Dispersal kernels from exponential power family
#'
#' `k_exponential_power` computes the value, multiplied by \eqn{N}, of a
#' dispersal kernel from the exponential power family that includes, as
#' special cases, Gaussian kernels and kernels that follow an exponential
#' function of the distance.
#'
#' @param x Numeric matrix of positions \eqn{x} relative to the seed source,
#'  or vector of distances \eqn{\left\|{x}\right\|} to the seed source.
#' @param par Numeric vector with two elements representing the
#'  log-transformed scale and shape parameters \eqn{a} and \eqn{b}.
#' @param N The multiplier \eqn{N}.
#' @param d The spatial dimension.
#'
#' @details The dispersal kernel, i.e. spatial probability density
#' of the location of a seed relative to its source, is here given by
#' \deqn{k(x)={b\Gamma (d/2) \over 2\pi ^{d/2}a^{d}\Gamma (d/b)}
#'    e^{-(\left\|{x}\right\|/a)^{b}},}
#' which corresponds to a probability density of the distance given by
#' \deqn{p(r)={b \over a^{d}\Gamma (d/b)}r^{d-1}e^{-(r/a)^{b}},}
#' where \eqn{d} is the spatial dimension, \eqn{\left\|{\,}\right\|}
#' denotes the Euclidean norm and the normalizing constants involve the
#' \link[base:gamma]{gamma} function; see Bateman (1947), Clark et al.
#' (1998), Austerlitz et al. (2004), Nathan et al. (2012) for the planar
#' case. This means the \eqn{b}th power of the distance has a
#' \link[stats:GammaDist]{gamma distribution} with shape parameter
#' \eqn{d/b} and scale parameter \eqn{a^{b}}.
#'
#' The kernel has its maximum at zero and represents a rather flexible family
#' that includes, for \eqn{b=2} the classical Gaussian kernels and for
#' \eqn{b=1}, kernels decreasing exponentially with the distance. For
#' \eqn{b<1} the distance distribution is fat-tailed in the sense of Kot et
#' al. (1996). Such kernels have consequently been applied in a number of
#' theoretical studies that address dispersal (Ribbens et al. 1994, Bullock
#' et al. 2017).
#'
#' @return Numeric vector of function values \eqn{k(x)} multiplied by \eqn{N}.
#'
#' @references
#' Bateman, A. (1947). Contamination in seed crops: III. relation with
#' isolation distance. *Heredity* **1**, 303–336.
#' \doi{10.1038/hdy.1947.20}
#'
#' Kot, M., Lewis, M.A., van den Driessche, P. (1996). Dispersal Data and the
#' Spread of Invading Organisms. *Ecology* **77(7)**, 2027–2042.
#' \doi{10.2307/2265698}
#'
#' Ribbens, E., Silander Jr, J.A., Pacala, S.W. (1994). Seedling recruitment
#' in forests: calibrating models to predict patterns of tree seedling
#' dispersion. *Ecology* **75**, 1794–1806.
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
#' Bullock, J. M., Mallada González, L., Tamme, R., Götzenberger, L., White,
#' S.M., Pärtel, M., Hooftman, D.A. (2017).  A synthesis of empirical plant
#' dispersal kernels. *Journal of Ecology* **105**, 6–19.
#' \doi{10.1111/1365-2745.12666}
#'
#' Nathan, R., Klein, E., Robledo‐Arnuncio, J.J., Revilla, E. (2012).
#' Dispersal kernels: review, in Clobert, J., Baguette, M., Benton, T.G.,
#' Bullock, J.M. (eds.), *Dispersal ecology and evolution*, 186–210.
#' \doi{10.1093/acprof:oso/9780199608898.003.0015}
#'
#' @export
#'
#' @examples
#' k_exponential_power(2:5, par=c(0,0), d=2)

k_exponential_power <- function(x, par, N=1, d=NCOL(x)) {
  r <- rownorms(x)
  a <- exp(par[1])
  b <- exp(par[2])
  N * b / (surface(d) * a^d * gamma(d/b)) * exp(-(r/a)^b)
  # Alternative for the last line (for nonzero r):
  # N / surface(d) * b * r^(b-d) * dgamma(r^b,d/b,,a^b)
}


##############################################################################
#' Dispersal kernels for Weibull distance distributions
#'
#' `k_weibull` computes the value, multiplied by \eqn{N}, of the dispersal
#' kernel from Tufto et al. (1997) based on seeds having a distance with a
#' Weibull distribution from their source.
#'
#' @param x Numeric matrix of positions \eqn{x} relative to the seed source,
#'  or vector of distances \eqn{\left\|{x}\right\|} to the seed source.
#' @param par Numeric vector with two elements representing the
#'  log-transformed scale and shape parameters \eqn{a} and \eqn{b} of the
#'  distance distribution.
#' @param N The multiplier \eqn{N}.
#' @param d The spatial dimension.
#'
#' @details The dispersal kernel, i.e. spatial probability density
#' of the location of a seed relative to its source, is here given by
#' \deqn{k(x)={b\Gamma (d/2) \over 2\pi ^{d/2}a^{b}}\left\|{x}\right\|^{b-d}
#'   e^{-(\left\|{x}\right\|/a)^{b}},}
#' which corresponds to a probability density of the distance given by
#' \deqn{p(r)={b \over a^{b}}r^{b-1}e^{-(r/a)^{b}},}
#' where \eqn{d} is the spatial dimension, \eqn{\left\|{\,}\right\|}
#' denotes the Euclidean norm and the normalizing constants involve the
#' \link[base:gamma]{gamma} function; see Tufto et al. (1997) for the planar
#' case. Thus, the distance is assumed to have the
#' \link[stats:Weibull]{Weibull distribution} with scale parameter \eqn{a}
#' and shape parameter \eqn{b}. Equivalently, the \eqn{b}th power of the
#' distance has an exponential distribution with scale parameter \eqn{a^{b}}.
#'
#' Consequently, if and only if \eqn{b<1}, the distance distribution has
#' a heavier tail than an exponential distribution, although with tail
#' probabilities still decreasing faster than any power law; it is a
#' fat-tailed distribution in the sense of Kot et al. (1996). The kernel
#' coincides with a Gaussian kernel in the special case \eqn{b=d=2}.
#'
#' @return Numeric vector of function values \eqn{k(x)} multiplied by \eqn{N}.
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
#' Kot, M., Lewis, M.A., van den Driessche, P. (1996). Dispersal Data and the
#' Spread of Invading Organisms. *Ecology* **77(7)**, 2027–2042.
#' \doi{10.2307/2265698}
#'
#' Nathan, R., Klein, E., Robledo‐Arnuncio, J.J., Revilla, E. (2012).
#' Dispersal kernels: review, in Clobert, J., Baguette, M., Benton, T.G.,
#' Bullock, J.M. (eds.), *Dispersal ecology and evolution*, 186–210.
#' \doi{10.1093/acprof:oso/9780199608898.003.0015}
#'
#' @export
#'
#' @examples
#' k_weibull(2:5, par=c(0,0), d=2)

k_weibull <- function(x, par, N=1, d=NCOL(x)) {
  r <- rownorms(x)
  a <- exp(par[1])
  b <- exp(par[2])
  if (!N) return(numeric(length(x))) # Compute 0*Inf as 0.
  N * b / (surface(d) * a^b) * r^(b-d) * exp(-(r/a)^b)
  # Alternative for the last line (for nonzero r):
  # N * dweibull(r,b,a) / surface(d,r)
}


##############################################################################
#' Power-law dispersal kernels
#'
#' `k_power` computes the value, multiplied by \eqn{N}, of a dispersal kernel
#' that follows a power law of a constant \eqn{a} plus the distance.
#'
#' @param x Numeric matrix of positions \eqn{x} relative to the seed source,
#' or vector of distances \eqn{\left\|{x}\right\|} to the seed source.
#' @param par Numeric vector with two elements representing the
#' log-transformed parameters \eqn{a} and \eqn{b}.
#' @param N The multiplier \eqn{N}.
#' @param d The spatial dimension.
#'
#' @details The dispersal kernel, i.e. spatial probability density
#' of the location of a seed relative to its source, is here given by
#' \deqn{k(x)={\Gamma (d/2) \over 2\pi ^{d/2}a^{d}B(d,b)}
#'   (1+{\left\|{x}\right\| \over a})^{-(b+d)},}
#' which corresponds to a probability density of the distance given by
#' \deqn{p(r)={1 \over a^{d}B(d,b)}r^{d-1}(1+{r \over a})^{-(b+d)},}
#' where \eqn{d} is the spatial dimension, \eqn{\left\|{\,}\right\|}
#' denotes the Euclidean norm and the normalizing constants involve the
#' \link[base:beta]{beta} and \link[base:gamma]{gamma} functions; see Nathan
#' et al. (2012) for the planar case (with \eqn{b} replaced by \eqn{b-d}).
#' This means the distance is \eqn{da \over b} times a random variable having
#' an \link[stats:FDist]{F distribution} with \eqn{2d} and \eqn{2b} degrees
#' of freedom. This is a fat-tailed distribution for all choices of the
#' parameter \eqn{b}.
#'
#' @return Numeric vector of function values \eqn{k(x)} multiplied by \eqn{N}.
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
#'
#' @export
#'
#' @examples
#' k_power(2:5, par=c(0,0), d=2)

k_power <- function(x, par, N=1, d=NCOL(x)) {
  r <- rownorms(x)
  a <- exp(par[1])
  b <- exp(par[2])
  N / (surface(d) * a^d * beta(d,b)) * (1+r/a)^(-d-b)
  # Alternatives for the last line (for nonzero r):
  # s<-b/(d*a); N * s * df(s*r, 2*d, 2*b) / surface(d,r)
  # q<-r/(a+r); N * a / surface(d) / r^(d+1) * q^2 * dbeta(q,d,b)
  # N * a / (a+r)^2 * dbeta(r/(a+r),d,b) / surface(d,r)
}
