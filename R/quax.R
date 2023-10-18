#' @import stats
NULL

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
#' further arguments passed to `optim` (see Details).
#' @param y Observed values \eqn{y_{1},...,y_{n}} of the regeneration density 
#' of the inventory plot.
#' @param tau Numeric between 0 and 1. Specifies the quantile \eqn{\tau } 
#' used for the quantile regression.
#' @param fun Dispersal kernel \eqn{k_{\theta }} assumed for the regeneration 
#' potential. Values allowed are `k_lognormal`, `k_t`, `k_power`, 
#' `k_weibull`, `k_exponential.power` or a custom function (see Examples). 
#' The default, `k_lognormal`, is to fit a model with log-normal distance 
#' distributions.
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
#' \code{\link[stats:optimize]{optimize}}, the result which is minimized in 
#' \eqn{\theta } using \code{\link[stats:optim]{optim}}.
#'
#' @return The estimated function, including an attribute `o` containing the 
#' results of `optim`.
#'
#' @references
#' Koenker, R., Bassett, G. (1978). Regression Quantiles. *Econometrica*
#' **46(1)**, 33–50.
#' \doi{10.2307/1913643}
#'
#' @examples
#' ## Prepare artificial example data:
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
#' fun <- c("k_lognormal","k_t","k_weibull","k_power","k_exponential.power")
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

quax.default <- function(..., y, tau, fun=k_lognormal,
    dim=2, weights=1, par=c(log.a=8, log.b=1), tol=1e-50) {
  Nmax <- 2*sum(y)  # (will be modified by optN)
  optN <- function(par, nframe, tol, ...,
      gr, method, lower, upper, control, hessian # (grab optim arguments)
  ) {
    S <- function(N, par, y, tau, w, fun, ...) {
      res <- y - fun(..., par=par, N=N)
      sum(w * res * (tau - 0.5 + 0.5*sign(res)))
    }
    s <- S(N=0, par=par, ...)
    while (S(N=Nmax, par=par, ...) < s) Nmax <<- 2*Nmax
    oN <- optimize(S, c(0, Nmax), par=par, ..., tol=tol)
    Nmax <<- (Nmax+oN$minimum)/2
    obj <- oN$objective
    attr(obj,"N") <- oN$minimum
    obj
  }
  o <- optim(par, optN, nframe=nframe, tol=tol,
    y=y, tau=tau, fun=fun, d=dim, w=weights, ...)
  obj <- optN(par=o$par, nframe=nframe, tol=tol,
    y=y, tau=tau, fun=fun, d=dim, w=weights, ...)
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
  x <- model.matrix(attr(mf,"terms"), mf)
  y <- mf[[1L]]
  o <- model.offset(mf)
  w <- model.weights(mf)
  quax.default(x=x, y=if (is.null(o)) y else y-o,
    tau=tau, fun=fun, weights=if (is.null(w)) 1 else w, ...)
}


#'summary.quax
#'
#'@description The function for printing the summary of the quantile regression.
#'
#' @return The quality criterion is the value that is minimised in the quantile regression. It represents the weighted sum of absolute residuals.
#' @details The function return a list including the estimated parameters for the quantile regression for the specific distribution function.#'
#'         The estimated function, including an attribute `o` containing the results of `optim`.
#' @examples #create dataframe
#' simulated.data <- data.frame(distance = rlnorm(200, meanlog = 5, sdlog = 1), density = rep(0:10,200))
#'
#'# run quax function
#'f1 <- quax(x = simulated.data$distance, y = simulated.data$density, tau = 0.9, fun = k_lognormal)
#'
#'# run summary.quax
#' summary.quax(f1)
#'
#' @export

summary.quax <- function(f)
  list(coef=formals(f)$par, value= attr(f,"o")$value)
