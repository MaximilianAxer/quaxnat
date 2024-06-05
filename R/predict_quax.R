##############################################################################
#' Prediction of Potential Regeneration Densites
#'
#' @description Prediction of the potential regeneration density as a function
#' of the distance to the nearest seed tree.
#'
#'
#'
#' @param distmap An object of the Distmap function. A SpatRaster with distances
#'  to the nearest seed tree is used for the prediction of the potential
#'  regeneration densities.
#'
#' @param quax A quax object is used for the prediction. This is a parameterised
#'  dispersal function using quantile regression.
#'
#' @details A SpatRaster is created with the same resolution as the input
#' raster, defined by the study area. The potential regeneration density is
#' calculated and given for each raster cell.
#'
#' @return  The regeneration density is given on the same scale
#' (e.g. numbers per hectare) as in the input data.
#'
#' @examples
#' ## Prepare artificial data:
#' set.seed(0)
#' r <- rgamma(200, shape=2, scale=150)
#' simulated.data <- data.frame(distance = r, density = rpois(length(r),
#' k_lognormal(r, par=c(6,0), N=1000000, d=2)))
#'
#' ## Run quax function:
#' f1 <- quax(x = simulated.data$distance, y = simulated.data$density,
#'           tau = 0.9, fun = k_lognormal)
#'
#' ## Create raster data set
#'  rr <- terra::rast(
#'  matrix(sample(0:10, 20 * 20, replace = TRUE),
#'         nrow = 20, ncol = 20))
#'
#' ## Compute distance for prediction area
#' distance <- Distmap(raster = rr, treespecies = "10")
#'
#' ## Prediction
#' p <- predict_quax(distmap = distance, quax = f1)
#' terra::plot(p)
#' @export


predict_quax <- function(distmap, quax) {
  prediction <- quax(terra::values(distmap))
  terra::values(distmap) <- prediction
  return(distmap)
}

