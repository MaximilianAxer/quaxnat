#' optN
#' @description  For given a,σ the global minimum in N can always be found with any desired precision, usually in a small number of steps, by successively shrinking an interval. We realize this as an inner, nested minimization in the function optN.
#'
#'
#' @inheritParams S
#'
#'
optN <- function(x, y, tau, fun, par, w,
                 Nmax = 1000000) {
  s <- S(x = x,
         y = y,
         tau = tau,
         par = c(par, 0),
         w = w,
         fun = fun)
  while (S(x = x,
           y = y,
           tau = tau,
           par = c(par, Nmax),
           w = w,
           fun = fun) < s) { # muss ich jetzt hier den ganzen Bumms eintragen? (x, y, par, tau, fun)
    Nmax <- 100*Nmax
  }
  optimize(function(N) S(x = x,
                         y = y,
                         tau = tau,
                         par = c(par, N),
                         w = w,
                         fun = fun),
           c(0, Nmax),
           tol=1e-50) # muss ich jetzt hier den ganzen Bumms eintragen? (x, y, par, tau, fun)
}
#optN(par = 1:3, x = 1, y = 1, tau = 0.9, fun = "Clark2dt") # dauert sehr lange???


#'SN
#' @inheritParams optN
#'
SN <- function(x, y, tau, fun, w, par, Nmax = 1000000) { # Search for the minimum of the reduced objective
  oN <- optN(x=x, y=y, tau=tau, fun=fun, w=w, par=par, Nmax=Nmax)
  cat("par =", par, "objective =", oN$objective, "N =", oN$minimum, "\n"); flush.console()
  oN$objective
}


#'quax
#'
#'@description The function for estimating the potential regeneration density.
#'
#'
#' @param x Numeric vector giving the distance to the nearest seed source for the inventory plot.
#' @param y Observed regeneration density of the inventory plot.
#' @param tau Numeric between 0 and 1. Specifies the quantile used for the quantile regression. {0;1}
#' @param fun Function assumed for the quantile regression of the regeneration potential. Values allowed are: "exponential", "weibull", "gamma", "Clark2dt", "lognormal", "half-normal". The default is to fit an lognormal model.
#' @param w Numeric vector of weights of the observations in the estimation procedure.
#' @param par Numeric vector of initial values for the parameters to be optimized over, exluding the first parameter `N`.
#' @param ... Further arguments passed to `optim`.
#' @details The function return a list including the estimated parameters for the quantile regression for the specific distribution function.
#'
#' @return The estimated function, including an attribute `o` containing the results of `optim`.

quax <-  function(x, y, tau, fun=lognormal, w=1, par=c(a=8, b=1), ...) {
  o <- optim(par, SN, x = x, y = y, tau = tau, fun = fun, w = w, ...)
  oN <- optN(x=x, y=y, tau=tau, fun=fun, w=w, par=o$par, Nmax = 1000000)
  formals(fun)$par <- c(o$par, N=oN$minimum)
  attr(fun,"o") <- o
  fun
}

#the function assumed for the dispersal distance distribution.
###        Values allowed are: "exponential", "weibull", "gamma", "pearson",
###        "lognormal", "half-normal", "half-cauchy", "half-t", or "user"
###        (a user-defined fonction).
