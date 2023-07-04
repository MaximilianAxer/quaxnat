#' predict_quax
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

predict_quax <- function(distmap, quax) {
  prediction <- quax(distmap)
  return(prediction)
}
