#' Estimating Potential Regeneration Densities by Quantile Regression
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

#' @rdname quax
quax <- function(...) UseMethod("quax")

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
