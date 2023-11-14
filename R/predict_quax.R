##############################################################################
#' Prediction of Potential Regeneration Densites
#'
#' @description Prediction of the potential regeneration density as a function of the distance to the nearest seed tree.
#'
#'
#'
#' @param distmap An object of the Distmap-Function. A SpatRaster with distances to the nearest seed tree is used for the prediction of the potential regeneration densities.
#' @param quax An quax object is used for the prediction. This is a parameterised dispersal function using quantile regression.
#'
#' @details A SpatRaster is created with the same resolution as the input raster, defined by the study area. The potential regeneration density is calculated and given for each raster cell.
#'
#' @return  The regeneration density is given in N/ha.
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
#' #Raster data set
#' rr <- terra::rast(
#'  matrix(sample(0:10, 20 * 20, replace = TRUE),
#'         nrow = 20, ncol = 20))
#'
#' #compute distance for prediction area
#' distance <- Distmap(fe_raster = rr, treespecies = "10")
#'
#' #prediction
#' predict_quax(distmap = distance, quax = f1)
#' 
#' @export


predict_quax <- function(distmap, quax) {
  prediction <- quax(terra::values(distmap))
  terra::values(distmap) <- prediction
  return(distmap)
}

