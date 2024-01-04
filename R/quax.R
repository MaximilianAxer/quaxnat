##############################################################################
#' Estimating Potential Regeneration Densities by Quantile Regression
#'
#' `quax` estimates parameters of a spatial dispersal kernel that describes 
#' the regeneration potential as the \eqn{\tau }th quantile of the 
#' regeneration density. Here \eqn{\tau } is between 0 and 1, with typical 
#' values close to 1 representing the situation that the full regeneration 
#' potential is realized only at a small fraction of all sites.
#'
#' @param ... Vector of positions \eqn{x_{1},...,x_{n}} or distances to the 
#' seed source as required by the specific dispersal kernel. Optionally, 
#' further arguments passed to `optim`, to the default method or to the 
#' kernel.
#' @param y Observed values \eqn{y_{1},...,y_{n}} of the regeneration density 
#' of the inventory plot.
#' @param tau Numeric between 0 and 1. Specifies the quantile \eqn{\tau } 
#' used in the regression.
#' @param fun Dispersal kernel \eqn{k_{\theta }} assumed for the regeneration 
#' potential. Values allowed are \code{\link{k_lognormal}}, 
#' \code{\link{k_t}}, \code{\link{k_power}}, \code{\link{k_weibull}}, 
#' \code{\link{k_exponential_power}} or a custom function (see Examples). The 
#' default, `k_lognormal`, is to fit a model with log-normal distance 
#' distributions.
#' @param dim The spatial dimension, by default equal to 2.
#' @param weights Numeric vector of optional nonnegative weights \eqn{w_i} of 
#' the observations in the estimation procedure. Default is 1.
#' @param par Numeric vector of initial values for the parameter vector 
#' \eqn{\theta }.
#' @param formula A formula of the form `y ~ x`.
#' @param data,subset,na.action,offset For the formula interface: 
#' Further arguments passed to `model.frame` (along with `weights`).
#' @param tol The desired accuracy for the inner optimization (see Details).
#'
#' @details The function estimates the parameters \eqn{N} and \eqn{\theta } 
#' of the regeneration potential \eqn{Nk_{\theta }} by minimizing
#' \deqn{\displaystyle \sum _{i=1}^{n}w_{i}(y_{i}-Nk_{\theta }(x_{i}))\bigl\{
#'   \begin{smallmatrix}\tau \hphantom{-1} &\text{if }y_{i}>Nk_{\theta }
#'   (x_{i})\\ \tau -1&\text{if not}\hphantom{.............}\end{smallmatrix}}
#' (Koenker and Bassett 1978). Due to convexity the minimum in \eqn{N} for a 
#' given vector \eqn{\theta } can always be found by successively shrinking 
#' an interval; this is implemented in an inner, nested minimization using 
#' \code{\link[stats:optimize]{optimize}}, the result of which is minimized 
#' in \eqn{\theta } using \code{\link[stats:optim]{optim}}.
#'
#' This is a rather naive approach to quantile regression that appears to 
#' work reasonably well for scaled dispersal kernels \eqn{Nk_{\theta }} as 
#' considered here, see Appendix A in Axer et al. (2021). For general 
#' quantile regression problems the more sophisticated procedure 
#' \code{\link[quantreg:nlrq]{nlrq}} in the package `quantreg`, based on 
#' Koenker and Park (1996), is expected to provide better results.
#'
##############################################################################
#' In particular, `quax` is subject to the usual numerical issues inherent in 
#' optimization: It can get stuck in a local minimum or altogether miss a 
#' minimum if the initial values (as specified by the argument `par`) are too 
#' far off or if the objective function exhibits bad behavior. Problems can 
#' further arise in the dispersal kernels if the parameter values, which are 
#' passed on a log scale, become too large. It is therefore recommended to 
#' visually check the results (see Examples). Also, the `optim` arguments 
#' `method` and `control` can be added in `...` to select and tune the 
#' optimization algorithm, but note that the objective function is usually not 
#' differentiable.
#'
#' See Koenker (2005) for a detailed exposition of quantile regression. 
#'
#' @return The estimated function, including an attribute `o` containing the 
#' results of `optim`. The returned function has the class attribute 
#' \code{quax}, so generic functions with methods defined for \code{quax} 
#' objects invoke these methods; see \code{\link{summary.quax}} for an 
#' example.
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
#' ## Use custom variant of lognormal model that includes a shift:
#' plot(simulated.data$position)
#' f4 <- quax(density ~ position, simulated.data,
#'   tau = 0.9, par = c(8, 1, 0, 0),
#'   fun = function(x, par, N, d)
#'     k_lognormal(x - rep(par[-(1:2)],each=NROW(x)), par[1:2], N, d)
#' )
#' summary(f4)
#'
#' @export
#' @rdname quax

quax <- function(...) UseMethod("quax")

#' @export
#' @rdname quax
#' @importFrom stats optim optimize

quax.default <- function(..., y, tau, fun=k_lognormal,
    dim=2, weights=1, par=c(log.a=8, log.b=1), tol=1e-50) {

  # Initialize hint for upper bound of N (will be modified by optN):
  Nmax <- 2*sum(y)

  # Define inner optimization function optN:
  optN <- function(par, tol, ...,
      gr, method, lower, upper, control, hessian # optional arguments not 
  ) {                                            #   to be passed on to S

    # Define objective function S:
    S <- function(N, par, y, tau, w, fun, ...) {
      res <- y - fun(..., par=par, N=N)
      sum(w * res * (tau - 0.5 + 0.5*sign(res)))
    }

    # Adjust upper bound of N:
    s <- S(N=0, par=par, ...)
    while (S(N=Nmax, par=par, ...) < s) Nmax <<- 2*Nmax

    # Find N that minimizes S:
    oN <- optimize(S, c(0, Nmax), par=par, ..., tol=tol)

    # Save updated hint for upper bound:
    Nmax <<- (Nmax+oN$minimum)/2

    # Return value of S, attaching N as an attribute:
    obj <- oN$objective
    attr(obj,"N") <- oN$minimum
    obj
  }

  # Search for parameter vector that minimizes optN:
  o <- optim(par, optN, tol=tol, y=y,
    tau=tau, fun=fun, d=dim, w=weights, ...)

  # Compute N for that parameter vector:
  obj <- optN(par=o$par, tol=tol, y=y,
    tau=tau, fun=fun, d=dim, w=weights, ...)

  # Set up and return result:
  formals(fun)$par <- o$par
  formals(fun)$N <- attr(obj,"N")
  formals(fun)$d <- dim
  attr(fun,"o") <- o
  class(fun) <- "quax"
  fun
}

#' @export
#' @rdname quax

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
#' Summarizing Quantile Regression Fits of Potential Regeneration Densities
#'
#' This function is the `summary` method for class `quax` objects as returned
#' by \code{\link{quax}}.
#'
#' @param object The function returned by \code{\link{quax}}.
#' @param ... not in use here
#'
#' @return A list with the following components:
#' \describe{
#'   \item{`coefficients`}{The parameters of the estimated dispersal kernel.}
#'   \item{`value`}{The attained value of the objective function that is
#'     minimised in the quantile regression.}
#' }
#'
#' @details The `value` component of the result can be used to compare the
#' quality of the fit of different dispersal kernels for the same quantile
#' to the same data.
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
#'
#' @export

summary.quax <- function(object, ...)
  list(coefficients=formals(object)$par, value=attr(object,"o")$value)
