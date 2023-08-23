#'objective function for minimization
#'
#'
#'@param x represents the distance to the nearest seed source. Must be numeric
#'@param fun dispersal function, which is used
#'@param par are parameters to be estimated
#'@param tau determines the quantile
#'@param y response: regeneration density in corresponding unit (per ha)
#'
#'@description  The parameters are estimated by minimizing a sum of weighted absolute errors, where y_1,...,y_n are the observed regeneration numbers and x_1,...,x_n,η_1,...,η_n are defined as above.
#' While common regression approaches assume that the average error in the model is zero, the errors are asymmetrically weighted for quantile regression (see Fahrmeir et al., 2013, Ch. 10).
#' The objective function is defined and the sum of absolute weighted errors is minimized.
#' \deqn{\sum_{i=1}^n x_i w_i |y_i - \eta_i|}
#'@details ...



S <- function(N, par, y, tau, w, fun, ...) {
  res <- y - if (N) fun(..., par=par, N=N) else 0  # (compute 0*Inf as 0)
  s <- tau - 0.5 + 0.5*sign(res)
  result <- sum(w * s * res)
  return(result)
}

#S(y = 1,  x = 1, par = 1:3, tau = 0.99,  fun = "l")
