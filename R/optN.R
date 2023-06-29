#' optN
#' @description  For given a,σ the global minimum in N can always be found with any desired precision, usually in a small number of steps, by successively shrinking an interval. We realize this as an inner, nested minimization in the function optN.
#'
#'
#' @inheritParams S
#'
#'
optN <- function(x, y, tau, fun,par,
                 Nmax = 1000000) {
  s <- S(x = x,
         y = y,
         tau = tau,
         par = c(par, 0),
         fun = fun)
  while (S(x = x,
           y = y,
           tau = tau,
           par = c(par, Nmax),
           fun = fun) < s) { # muss ich jetzt hier den ganzen Bumms eintragen? (x, y, par, tau, fun)
    Nmax <- 100*Nmax
  }
  optimize(function(N) S(x = x,
                         y = y,
                         tau = tau,
                         par = c(par, N),
                         fun = fun),
           c(0, Nmax),
           tol=1e-50) # muss ich jetzt hier den ganzen Bumms eintragen? (x, y, par, tau, fun)
}
#optN(par = 1:3, x = 1, y = 1, tau = 0.9, fun = "Clark2dt") # dauert sehr lange???


#'SN
#' @inheritParams optN
#'
SN <- function(x, y, tau, fun,par, Nmax = 1000000) { # Search for the minimum of the reduced objective
  o <- optN(x, y, tau, fun,par, Nmax = 1000000)
  cat("par =", par, "objective =", o$objective, "N =", o$minimum, "\n"); flush.console()
  o$objective
}


#'quax
#'
#'@description The function for estimating the potential regeneration density.
#'
#'
#' @param x numeric vector giving the distance to the nearest seed source for the inventory plot.
#' @param y Observed regeneration density of the inventory plot.
#' @param tau Represents the quantile, which is used for the quantile regression. {0;1}
#' @param fun function assumed for the quantile regression of the regeneration potential.Values allowed are: "exponential", "weibull", "gamma", "Clark2dt", "lognormal", "half-normal". The default is to fit an lognormal model.
#' @details The function return a list including the estimated parameters for the quantile regression for the specific distribution function.
#'
#' @return  The parameter values are given. In the next step the N optimization is done.

quax <-  function(x, y, tau, fun=lognormal) {
  optim(c(8, 1), SN, x = x, y = y, tau = tau, fun = fun, control=list(reltol=0))
}

#the function assumed for the dispersal distance distribution.
###        Values allowed are: "exponential", "weibull", "gamma", "pearson",
###        "lognormal", "half-normal", "half-cauchy", "half-t", or "user"
###        (a user-defined fonction).
