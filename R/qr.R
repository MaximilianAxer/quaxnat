# Internal functions:
# `rownorms` returns the Euclidean norms of the rows of a matrix x:
rownorms <- function(x,...) if (length(d<-dim(x))>1L && d[2L]>1)
  sqrt(.rowSums(x^2,d[1],d[2],...)) else abs(x)
# `surface` returns the volume (area, length) of the surface of a ball with 
# radius r in d-dimensional space:
surface <- function(d,r=1) 2*pi^(d/2)/gamma(d/2) * r^(d-1)

##############################################################################
#' Dispersal Kernels For Log-Normal Distance Distributions
#'
#' `k_lognormal` computes the value, multiplied by \eqn{N}, of a dispersal
#' kernel based on seeds having a distance with a log-normal distribution
#' from the their source.
#'
#' @return Numeric vector of function values \eqn{k(x)} multiplied by \eqn{N}.
#'
#' @param x Numeric matrix of positions \eqn{x} relative to the seed source, 
#' or vector of distances \eqn{\left\|{x}\right\|} to the seed source.
#' @param par Numeric vector with two elements representing log-transformed
#' scale and shape parameters, given by the mean \eqn{a} and standard
#' deviation \eqn{\sigma} of the underlying normal distribution.
#' @param N The multiplier \eqn{N}.
#' @param d The spatial dimension.
#'
#' @details The dispersal kernel, i.e. spatial probability density 
#' of the location of a seed relative to its source, is here given by
#' \deqn{k(x)={\Gamma (d/2) \over 
#'   2\pi ^{d/2}\left\|{x}\right\|^{d}\sqrt{2\pi \sigma ^{2}}}
#'   e^{-{1 \over 2\sigma ^{2}}(\log (\left\|{x}\right\|/a))^{2}},}
#' which corresponds to a probability density of the distance given by
#' \deqn{p(r)={1 \over r\sqrt{2\pi \sigma ^{2}}}
#'   e^{-{1 \over 2\sigma ^{2}}(\log (r/a))^{2}},}
#' where \eqn{d} is the spatial dimension, \eqn{\left\|{\,}\right\|} 
#' denotes the Euclidean norm and the normalizing constant of the kernel 
#' involves the \link[base:beta]{gamma} function; see Greene and Johnson 
#' (1989), Stoyan and Wagner (2001) for the planar case. Thus, the distance 
#' is assumed to have the \link[stats:Lognormal]{log-normal distribution} 
#' such that the log-distance has a normal distribution with mean \eqn{a} and 
#' variance \eqn{\sigma^2}. Here \eqn{\log k(x)} is a quadratic function of 
#' \eqn{\left\|{x}\right\|} with a maximum at \eqn{\log a-d\sigma^2}.
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

k_lognormal <- function(x, par, N=1, d=NCOL(x)) {
  r <- rownorms(x)
  log.a <- par[1]
  σ <- exp(par[2])
  N / (surface(d)*sqrt(2*pi*σ^2)) * exp(-(log(r)-log.a)^2/(2*σ^2)) / r^d
  # Alternative for the last line:
  # N / surface(d,r) * dlnorm(r,log.a,σ)
}


##############################################################################
#' Dispersal Kernels From Spatial t Distribution
#'
#' `k_t` computes the value of the dispersal function from Clark et al.
#' (1999) multiplied by \eqn{N}.
#'
#' @return Numeric vector of function values \eqn{k(x)} multiplied by \eqn{N}.
#'
#' @param x Numeric matrix of positions \eqn{x} relative to the seed source, 
#' or vector of distances \eqn{\left\|{x}\right\|} to the seed source.
#' @param par Numeric vector with two elements representing log-transformed
#' parameters \eqn{a} and \eqn{p}.
#' @param N The multiplier \eqn{N}.
#' @param d The spatial dimension.
#'
#' @details The dispersal kernel, i.e. spatial probability density 
#' of the location of a seed relative to its source, is here given by
#' \deqn{k(x)={\Gamma (d/2) \over \pi ^{d/2}a^{d}\Beta (d/2,p)}
#'   (1+{\left\|{x}\right\|^{2} \over a^{2}})^{-(p+d/2)},}
#' which corresponds to a probability density of the distance given by
#' \deqn{p(r)={2 \over a^{d}\Beta (d/2,p)}r^{d-1}
#'   (1+{r^{2} \over a^{2}})^{-(p+d/2)},}
#' where \eqn{d} is the spatial dimension, \eqn{\left\|{\,}\right\|} 
#' denotes the Euclidean norm and the normalizing constants involve the 
#' \link[base:beta]{beta} and \link[base:beta]{gamma} functions; see Clark 
#' et al. (1999) and Austerlitz et al. (2004) for the planar case (with 
#' parameterizations \eqn{a=\sqrt{u}} and \eqn{p=b-d/2}, respectively). This 
#' means the squared distance is \eqn{da^{2} \over 2p} times a random 
#' variable having an \link[stats:FDist]{F distribution} with \eqn{d} and 
#' \eqn{2p} degrees of freedom.
#'
#' This dispersal kernel represents a mixture of Gaussian densities (with 
#' inverse variance following a gamma distribution with shape parameter 
#' \eqn{p}) that exhibits heavier tails. It has its maximum at zero.
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

k_t <- function(x, par, N=1, d=NCOL(x)) {
  r <- rownorms(x)
  a <- exp(par[1])
  p <- exp(par[2])
  s<-2*p/(d*a^2); N * 2 * s / surface(d) * r^(2-d) * df(s*r^2, d, 2*p)
  # Alternatives for the last line:
  # s<-2*p/(d*a^2); N * 2 * s * r * df(s*r^2, d, 2*p) / surface(d,r)
  # N * 2 / (surface(d) * a^d * beta(d/2, p)) / (1+(r/a)^2)^(p+d/2)
}


##############################################################################
#' Dispersal Kernels From Exponential Power Family
#'
#' `k_exponential.power` computes the value, multiplied by \eqn{N}, of a 
#' dispersal kernel from the exponential power family, which includes, as 
#' special cases, distance distributions based on normal and exponential 
#' distributions.
#'
#' @return Numeric vector of function values \eqn{k(x)} multiplied by \eqn{N}.
#'
#' @param x Numeric matrix of positions \eqn{x} relative to the seed source, 
#' or vector of distances \eqn{\left\|{x}\right\|} to the seed source.
#' @param par Numeric vector with two elements representing the 
#' log-transformed scale and shape parameters \eqn{a} and \eqn{b}.
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
#' \link[base:beta]{gamma} function; see Bateman (1947), Clark et al. (1998), 
#' Austerlitz et al. (2004), Nathan et al. (2012) for the planar case. This 
#' means the \eqn{b}th power of the distance has a gamma distribution with 
#' shape parameter \eqn{d/b} and scale parameter \eqn{a^{b}}.
#' 
#' This function has its maximum at zero and represents a rather flexible 
#' family of distributions including the classical bivariate Gaussian 
#' kernels, the kernels based on an exponential distribution of distances, 
#' and, for \eqn{b<1}, fat-tailed distributions. For \eqn{b>1}, it shows 
#' thin-tailed distributions. It has consequently been applied in a number of 
#' theoretical studies that address dispersal (Ribbens et al. 1994; 
#' Bullock et al. 2017).
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
#' #'Bullock, J. M., Mallada González, L., Tamme, R., Götzenberger, L., White, 
#' #'S. M., Pärtel, M., Hooftman, D. A.
#' (2017).  A synthesis of empirical plant dispersal kernels.
#' *Journal of Ecology* **105**, 6-19.
#' \doi{10.1111/1365-2745.12666}
#'
#' Nathan, R., Klein, E., Robledo‐Arnuncio, J.J., Revilla, E. (2012).
#' Dispersal kernels: review, in Clobert, J., Baguette, M., Benton, T.G., 
#' Bullock, J.M. (eds.), *Dispersal ecology and evolution*, 186–210.
#' \doi{10.1093/acprof:oso/9780199608898.003.0015}
#'

k_exponential.power <- function(x, par, N=1, d=NCOL(x)) {
  r <- rownorms(x)
  a <- exp(par[1])
  b <- exp(par[2])
  N * b / (surface(d) * a^d * gamma(d/b)) * exp(-(r/a)^b)
  # Alternative for the last line:
  # N / surface(d) * b * r^(b-d) * dgamma(r^b,d/b,,a^b)
}


##############################################################################
#' Dispersal Kernels For Weibull Distance Distributions
#'
#' `k_weibull` computes the value, multiplied by \eqn{N}, of the dispersal 
#' function from Tufto et al. (1997) based on seeds having a distance with a 
#' Weibull distribution from the their source.
#'
#' @return Numeric vector of function values \eqn{k(x)} multiplied by \eqn{N}.
#'
#' @param x Numeric matrix of positions \eqn{x} relative to the seed source, 
#' or vector of distances \eqn{\left\|{x}\right\|} to the seed source.
#' @param par Numeric vector with two elements representing the 
#' log-transformed scale and shape parameters \eqn{a} and \eqn{b} of the 
#' distance distribution.
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
#' \link[base:beta]{gamma} function; see Tufto et al. (1997) for the planar 
#' case. Thus, the distance is assumed to have the 
#' \link[stats:Weibull]{Weibull distribution} with scale parameter \eqn{a} 
#' and shape parameter \eqn{b}. Equivalently, the \eqn{b}th power of the 
#' distance has an exponential distribution with scale parameter \eqn{a^{b}}.
#' 
#' The distribution is fat-tailed when \eqn{b<1} and
#' thin-tailed otherwise (Nathan et al. 2012).
#' For \eqn{b>1}, the mode of the function is at \eqn{x>1}. In this way, the 
#' function approaches the normal distribution.
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

k_weibull <- function(x, par, N=1, d=NCOL(x)) {
  r <- rownorms(x)
  a <- exp(par[1])
  b <- exp(par[2])
  N * dweibull(r,b,a) / surface(d,r)
  # Alternative for the last line:
  # N * b / (surface(d) * a^b) * r^(b-d) * exp(-(r/a)^b)
}


##############################################################################
#' Power-Law Dispersal Kernels
#'
#' `k_power` computes the value of the dispersal function from (WHERE?)
#' multiplied by \eqn{N}.
#'
#' @return Numeric vector of function values \eqn{k(x)} multiplied by \eqn{N}.
#'
#' @param x Numeric matrix of positions \eqn{x} relative to the seed source, 
#' or vector of distances \eqn{\left\|{x}\right\|} to the seed source.
#' @param par Numeric vector with two elements representing the 
#' log-transformed parameters \eqn{a} and \eqn{p} of the distance 
#' distribution.
#' @param N The multiplier \eqn{N}.
#' @param d The spatial dimension.
#'
#' @details The dispersal kernel, i.e. spatial probability density 
#' of the location of a seed relative to its source, is here given by
#' \deqn{k(x)={\Gamma (d/2) \over 2\pi ^{d/2}a^{d}\Beta(d,p)}
#'   (1+{\left\|{x}\right\| \over a})^{-(p+d)},}
#' which corresponds to a probability density of the distance given by
#' \deqn{p(x)={1 \over a^{d}\Beta(d,p)}r^{d-1}(1+{r \over a})^{-(p+d)},}
#' where \eqn{d} is the spatial dimension, \eqn{\left\|{\,}\right\|} 
#' denotes the Euclidean norm and the normalizing constants involve the 
#' \link[base:beta]{beta} and \link[base:beta]{gamma} functions; see Nathan 
#' et al. (2012) for the planar case (with the parameterization \eqn{p=b-d}). 
#' This means the distance is \eqn{da \over p} times a random variable 
#' having an \link[stats:FDist]{F distribution} with \eqn{2d} and 
#' \eqn{2p} degrees of freedom.

#' Austerlitz 
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

k_power <- function(x, par, N=1, d=NCOL(x)) {
  r <- rownorms(x)
  a <- exp(par[1])
  p <- exp(par[2])
  N / (surface(d) * a^d * beta(d,p)) * (1+r/a)^(-d-p)
  # Alternatives for the last line:
  # s<-p/(d*a); N * s * df(s*r, 2*d, 2*p) / surface(d,r)
  # q<-r/(a+r); N * a / surface(d) / r^(d+1) * q^2 * dbeta(q,d,p)
  # N * a / (a+r)^2 * dbeta(r/(a+r),d,p) / surface(d,r)
}


##############################################################################
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
