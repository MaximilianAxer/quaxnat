#' Estimation of Potential Regeneration Densities by Quantile Regression
#'
#' `quax` estimates parameters of a spatial dispersal kernel that describes the potential regeneration density.
#'
#' @param x Numeric vector giving the distance to the nearest seed source for the inventory plot.
#' @param y Observed regeneration density of the inventory plot.
#' @param tau Numeric between 0 and 1. Specifies the quantile used for the quantile regression.
#' @param fun Function assumed for the quantile regression of the regeneration potential. Values allowed are `k_lognormal`, `k_t`, `k_power`, `k_weibull`, `k_exponential.power` or a custom function. The default is to fit a lognormal model.
#' @param weights Numeric vector of weights of the observations in the estimation procedure. Default is 0.
#' @param par Numeric vector of initial values for the parameters to be optimized over, exluding the first parameter `N`.
#' @param ... Further arguments passed to `optim`.
#' @param formula A formula of the form `y ~ x`.
#' @param data,subset,weights,na.action,offset For the formula interface: Further arguments passed to `model.frame` (along with `weights`).
#' @param tol The desired accuracy for the inner optimization, see `optimize`.
#' @details The function return a list including the estimated parameters for the quantile regression for the specific distribution function. ... The global minimum in \eqn{N} for given other parameters can always be found with any desired precision, usually in a small number of steps, by successively shrinking an interval. We realize this as an inner, nested minimization in an internal function `optN`.
#'
#' @return The estimated function, including an attribute `o` containing the results of `optim`.
#' @examples #create dataframe
#' simulated.data <- data.frame(distance = rlnorm(200, meanlog = 5, sdlog = 1), density = rep(0:10,200))
#'
#'# run quax function
#'quax(x = simulated.data$distance, y = simulated.data$density, tau = 0.9, fun = k_lognormal)


#' @rdname quax
quax <- function(...) UseMethod("quax")

#' @rdname quax
quax.default <- function(..., y, tau, fun=k_lognormal,
    dim=2, weights=1, par=c(log.a=8, log.b=1), tol=1e-50) {
  Nmax <- 2*sum(y)  # (will be modified by optN)
  optN <- function(par, nframe, tol, ...,
      gr, method, lower, upper, control, hessian # (grab optim arguments)
  ) {
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

summary.quax <- function(f)
  list(coef=formals(f)$par, value= attr(f,"o")$value)
