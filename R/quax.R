##############################################################################
#' Estimating potential regeneration densities by quantile regression
#'
#' `quax` estimates parameters of a spatial dispersal kernel that describes
#' the regeneration potential as the \eqn{\tau }th quantile of the
#' regeneration density. Here \eqn{\tau } is between 0 and 1, with typical
#' values close to 1 representing the situation that the full regeneration
#' potential is realized only at a small fraction of all sites.
#'
#' @param ... Vector of positions \eqn{x_{1},...,x_{n}} or distances to the
#'  seed source as required by the specific dispersal kernel. Optionally,
#'  further arguments passed to `optim`, to the default method or to the kernel.
#' @param y Vector of observed values \eqn{y_{1},...,y_{n}} of the regeneration
#'  density of the inventory plot.
#' @param tau Numeric value between 0 and 1. Specifies the quantile \eqn{\tau }
#'  used in the regression.
#' @param fun Function representing the dispersal kernel \eqn{k_{\theta }},
#'  multiplied by \eqn{N}, that is assumed for the regeneration potential.
#'  Values allowed are \code{\link{k_lognormal}}, \code{\link{k_t}},
#'  \code{\link{k_power}}, \code{\link{k_weibull}},
#'  \code{\link{k_exponential_power}} or a custom function with nonnegative
#'  values whose parameters include, in addition to the arguments in `...` not
#'  consumed by `optim` or the default method, the scaling factor `N` and the
#'  spatial dimension `d` (see Examples). The default, `k_lognormal`, is to
#'  fit a model with log-normal distance distributions.
#' @param dim The spatial dimension, by default equal to 2.
#' @param weights Numeric vector of optional non-negative weights \eqn{w_i} of
#'  the observations in the estimation procedure. Default is 1.
#' @param par Numeric vector of initial values for the parameter vector
#'  \eqn{\theta }.
#' @param formula A formula of the form `y ~ x`.
#' @param data,subset,na.action,offset For the formula interface:
#'  Further arguments passed to `model.frame` (along with `weights`).
#'
#' @details The function estimates the parameters \eqn{N} and \eqn{\theta }
#' of the regeneration potential \eqn{Nk_{\theta }} by minimizing
#' \deqn{\displaystyle \sum _{i=1}^{n}w_{i}\rho _{\tau }(y_{i}-Nk_{\theta }
#'   (x_{i})),}
#' where \eqn{\rho _{\tau }(u)=\int _{0}^{u}\tau -\mathbf{1}_{s<0}\;ds=\bigl\{
#'   \begin{smallmatrix}u\tau & \text{if }u\geq 0\\u(\tau -1)&\text{if }u<0
#'   \end{smallmatrix}}
#' (Koenker and Bassett 1978, Chapter 6.6 in Koenker 2005). The preceding
#' line, after subtracting the same expression for \eqn{N=0} and substituting
#' \eqn{s=y_{i}-tk_{\theta }(x_{i})} in the integral, becomes
#' \eqn{\int _{0}^{N}\sum _{i=1}^{n}w_{i}k_{\theta }(x_{i})(
#'   \mathbf{1}_{y_{i}<tk_{\theta }(x_{i})}-\tau )\;dt},
#' and any \eqn{N} such that the last integrand is \eqn{\leq 0} for \eqn{t<N}
#' and \eqn{\geq 0} for \eqn{t>N}, which can always be found as the integrand
#' is increasing in \eqn{t}, minimizes this integral. The integrand being the
#' difference of the sum of \eqn{w_{i}k_{\theta }(x_{i})} over the \eqn{i}
#' with \eqn{y_{i}<tk_{\theta }(x_{i})} and \eqn{\tau } times the sum over
#' all \eqn{i}, with relevant terms for nonzero \eqn{k_{\theta }(x_{i})},
#' this means that the estimate of \eqn{N} for a given vector \eqn{\theta }
#' can be computed as a \eqn{\tau }th quantile. This is implemented as an
#' inner, nested minimization, the result of which is minimized in
#' \eqn{\theta } using \code{\link[stats:optim]{optim}}.
#'
#' This is a rather naive approach to quantile regression that appears to
#' work reasonably well for scaled dispersal kernels \eqn{Nk_{\theta }} as
#' considered here, see Appendix A in Axer et al. (2021). For general
#' quantile regression problems the more sophisticated procedure
#' \code{\link[quantreg:nlrq]{nlrq}} in the package `quantreg`, based on
#' Koenker and Park (1996), is expected to provide better results.
#'
#' In particular, `quax` is subject to the usual numerical issues inherent in
#' optimization: It can get stuck in a local minimum or altogether miss a
#' minimum if the initial values (as specified by the argument `par`) are too
#' far off or if the objective function exhibits bad behavior. Problems can
#' further arise in the dispersal kernels if parameter values passed on a log
#' scale become too large. It is therefore recommended to visually check the
#' results (see Examples). Also, the `optim` arguments `method` and `control`
#' can be added in `...` to select and tune the optimization algorithm, but
#' note that the objective function is usually not differentiable.
#'
#' See Koenker (2005) for a detailed exposition of quantile regression.
#'
#' @return An objcet of class \code{quax} containing the estimated function,
#'  including an attribute `o` containing the results of `optim`.
#'  Generic functions with methods defined for \code{quax} objects invoke these
#'  methods; see \code{\link{summary.quax}} for an example.
#'
#' @references
#' Koenker, R., Bassett, G. (1978). Regression quantiles. *Econometrica*
#' **46**(1), 33–50.
#' \doi{10.2307/1913643}
#'
#' Axer, M., Schlicht, R., Wagner, S. (2021). Modelling potential density of
#' natural regeneration of European oak species (*Quercus robur* L., *Quercus
#' petraea* (Matt.) Liebl.) depending on the distance to the potential seed
#' source: Methodological approach for modelling dispersal from inventory
#' data at forest enterprise level. *Forest Ecology and Management* **482**,
#' 118802.
#' \doi{10.1016/j.foreco.2020.118802}
#'
#' Koenker, R., Park, B.J. (1996). An interior point algorithm for nonlinear
#' quantile regression. *Journal of Econometrics* **71**(1–2), 265–283.
#' \doi{10.1016/0304-4076(96)84507-6}
#'
#' Koenker, R. (2005). Quantile regression. Cambridge University Press.
#' \doi{10.1017/CBO9780511754098}
#'
#' @seealso
#' Function \code{\link[quantreg:nlrq]{nlrq}} in the package
#' \code{\link[quantreg:nlrq]{quantreg}}.
#'
#' @rdname quax
#' @export
#'
#' @examples
#' ## Prepare artificial data:
#' set.seed(0)
#' r <- rgamma(200, shape=2, scale=150)
#' simulated.data <- data.frame(distance = r, density =
#'   rpois(length(r), k_lognormal(r, par=c(6,0), N=1000000, d=2)))
#' plot(density ~ distance, simulated.data)
#'
#' ## Run quax function:
#' f1 <- quax(x = simulated.data$distance, y = simulated.data$density,
#'   tau = 0.9, fun = k_lognormal)
#' summary(f1)
#' curve(f1(x), add=TRUE)
#'
#' ## Do the same using formula interface:
#' f1 <- quax(density ~ distance, simulated.data,
#'   tau = 0.9, fun = k_lognormal)
#' summary(f1)
#' #quantreg::nlrq(density ~ k_lognormal(distance,c(log.a,log.b),N=N,d=2),
#' #  simulated.data, start = c(log.a=6,log.b=0,N=1e6), tau = 0.9)  # similar
#'
#' ## Use another quantile:
#' f2 <- quax(density ~ distance, simulated.data,
#'   tau = 0.99, fun = k_lognormal)
#' summary(f2)
#' curve(f2(x), add=TRUE, lwd=0)
#'
#' ## Show effect of weights:
#' f3 <- quax(density ~ distance, simulated.data,
#'   tau = 0.9, fun = k_lognormal, weights = distance)
#' summary(f3)
#' curve(f3(x), add=TRUE, lty=3)
#'
#' ## Compare various dispersal models:
#' fun <- c("k_lognormal","k_t","k_weibull","k_power","k_exponential_power")
#' for (i in seq_along(fun))
#'   curve(quax(density ~ distance, simulated.data,
#'     tau = 0.9, fun = get(fun[i]), weights = distance)(x),
#'     add=TRUE, col=i, lty=3)
#' legend("topright", fun, col=seq_along(fun), lty=3)
#'
#' ## Use positions in computation:
#' simulated.data$position <- r *
#'   (\(a) cbind(cos(a),sin(a))) (rnorm(length(r)))
#' f3 <- quax(density ~ position, simulated.data,
#'   tau = 0.9, fun = k_lognormal, weights = distance)
#' summary(f3)
#'
#' ## Show problems with bad initial values and try another parameterization:
#' curve(quax(density ~ distance, simulated.data, par = c(log.a=0,log.b=0),
#'   tau = 0.99, fun = k_lognormal)(x), add=TRUE, lty=2)
#' curve(quax(density ~ distance, simulated.data, par = c(a=1,b=1),
#'   tau = 0.99, fun = function(x,par,N,d) if (any(par<=0)) rep(NA,NROW(x))
#'     else k_lognormal(x,log(par),N,d))(x), add=TRUE, lty=2)
#'
#' ## Use custom variant of lognormal model that includes a shift:
#' plot(simulated.data$position)
#' f4 <- quax(density ~ position, simulated.data,
#'   tau = 0.9, par = c(8, 1, 0, 0),
#'   fun = function(x, par, N, d)
#'     k_lognormal(x - rep(par[-(1:2)],each=NROW(x)), par[1:2], N, d)
#' )
#' summary(f4)

quax <- function(...) UseMethod("quax")


#' @rdname quax
#' @importFrom stats optim
#' @export

quax.default <- function(..., y, tau, fun=k_lognormal,
    weights=1, dim=2, par=c(log.a=8, log.b=1)) {

  # Define inner optimization function optN:
  optN <- function(par, y, tau, fun, w, ...,
      gr, method, lower, upper, control, hessian # optional arguments not
  ) {                                            #   to be passed on to fun

    # Compute N that minimizes objective function:
    k <- fun(..., par=par, N=1); stopifnot(length(k)==length(y))
    z <- c(0, y / k)                             # extra 0 (with mass 0) to
    o <- order(z, na.last=NA)                    #   avoid empty o
    s <- cumsum(c(0,w*k)[o])
    N <- z[o[s >= tau*s[length(s)]][1L]]

    # Return value of objective function, attaching N as an attribute:
    res <- y - N*k
    obj <- sum(w * res * (tau - 0.5 + 0.5*sign(res)))
    attr(obj,"N") <- N
    obj
  }

  # Search for parameter vector that minimizes optN:
  o <- optim(par, optN, y=y, tau=tau, fun=fun, w=weights, d=dim, ...)

  # Compute N for that parameter vector:
  obj <- optN(par=o$par, y=y, tau=tau, fun=fun, w=weights, d=dim, ...)

  # Set up and return result:
  formals(fun)$par <- o$par
  formals(fun)$N <- attr(obj,"N")
  formals(fun)$d <- dim
  attr(fun,"o") <- o
  class(fun) <- "quax"
  fun
}


#' @rdname quax
#' @export

quax.formula <- function(formula, data, tau, fun=k_lognormal,
    subset, weights, na.action, offset, ...) {
  cl <- match.call(expand.dots=FALSE)
  cl <- cl[c(1L, match(c("formula","data","subset","weights",
    "na.action","offset"), names(cl), 0L))]
  cl[[1L]] <- quote(stats::model.frame)
  cl$formula[[3L]] <- substitute(0+x, list(x=formula[[3L]]))
  mf <- eval.parent(cl)
  x <- stats::model.matrix(attr(mf,"terms"), mf)
  y <- mf[[1L]]
  o <- stats::model.offset(mf)
  w <- stats::model.weights(mf)
  quax.default(x=x, y=if (is.null(o)) y else y-o,
    tau=tau, fun=fun, weights=if (is.null(w)) 1 else w, ...)
}


##############################################################################
#' Summarizing quantile regression fits of potential regeneration densities
#'
#' This function is the `summary` method for class `quax` objects as returned
#' by \code{\link{quax}}.
#'
#' @param object The function returned by \code{\link{quax}}.
#' @param ... not in use here
#'
#' @details The `value` component of the result can be used to compare the
#'  quality of the fit of different dispersal kernels for the same quantile
#'  to the same data.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{`coefficients`}{The parameters of the estimated dispersal kernel.}
#'   \item{`value`}{The attained value of the objective function that is
#'     minimised in the quantile regression.}
#' }
#'
#' @export
#'
#' @examples
#' ## Prepare artificial data:
#' set.seed(0)
#' r <- rgamma(200, shape=2, scale=150)
#' simulated.data <- data.frame(distance = r, density =
#'   rpois(length(r), k_lognormal(r, par=c(6,0), N=1000000, d=2)))
#' plot(density ~ distance, simulated.data)
#'
#' ## Fit a log-normal and a power-law dispersal kernel to the data:
#' f1 <- quax(density ~ distance, simulated.data,
#'   tau = 0.9, fun = k_lognormal)
#' f2 <- quax(density ~ distance, simulated.data,
#'   tau = 0.9, fun = k_power)
#'
#' ## Compare both fits:
#' summary(f1)
#' summary(f2)

summary.quax <- function(object, ...)
  list(coefficients=formals(object)$par, value=attr(object,"o")$value)
