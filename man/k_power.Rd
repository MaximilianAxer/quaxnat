% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dispersal_kernels.R
\name{k_power}
\alias{k_power}
\title{Power-law dispersal kernels}
\usage{
k_power(x, par, N = 1, d = NCOL(x))
}
\arguments{
\item{x}{Numeric matrix of positions \eqn{x} relative to the seed source,
or vector of distances \eqn{\left\|{x}\right\|} to the seed source.}

\item{par}{Numeric vector with two elements representing the
log-transformed parameters \eqn{a} and \eqn{b}.}

\item{N}{The multiplier \eqn{N}.}

\item{d}{The spatial dimension.}
}
\value{
Numeric vector of function values \eqn{k(x)} multiplied by \eqn{N}.
}
\description{
\code{k_power} computes the value, multiplied by \eqn{N}, of a dispersal kernel
that follows a power law of a constant \eqn{a} plus the distance.
}
\details{
The dispersal kernel, i.e. spatial probability density
of the location of a seed relative to its source, is here given by
\deqn{k(x)={\Gamma (d/2) \over 2\pi ^{d/2}a^{d}B(d,b)}
  (1+{\left\|{x}\right\| \over a})^{-(b+d)},}
which corresponds to a probability density of the distance given by
\deqn{p(r)={1 \over a^{d}B(d,b)}r^{d-1}(1+{r \over a})^{-(b+d)},}
where \eqn{d} is the spatial dimension, \eqn{\left\|{\,}\right\|}
denotes the Euclidean norm and the normalizing constants involve the
\link[base:beta]{beta} and \link[base:gamma]{gamma} functions; see Nathan
et al. (2012) for the planar case (with \eqn{b} replaced by \eqn{b-d}).
This means the distance is \eqn{da \over b} times a random variable having
an \link[stats:FDist]{F distribution} with \eqn{2d} and \eqn{2b} degrees
of freedom. This is a fat-tailed distribution for all choices of the
parameter \eqn{b}.
}
\examples{
k_power(2:5, par=c(0,0), d=2)
}
\references{
Nathan, R., Klein, E., Robledo‐Arnuncio, J.J., Revilla, E. (2012).
Dispersal kernels: review, in Clobert, J., Baguette, M., Benton, T.G.,
Bullock, J.M. (eds.), \emph{Dispersal ecology and evolution}, 186–210.
\doi{10.1093/acprof:oso/9780199608898.003.0015}

Austerlitz, F., Dick, C.W., Dutech, C., Klein, E.K., Oddou-Muratorio, S.,
Smouse, P.E., Sork, V.L. (2004). Using genetic markers to estimate the
pollen dispersal curve. \emph{Molecular Ecology} \strong{13}, 937–954.
\doi{10.1111/j.1365-294X.2004.02100.x}
}
