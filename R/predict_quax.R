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
#' @examples
#' #quax-Object
#' f <- quax(Dgl_B0 ~ distance_dgl, VJ_pot, subset=Dgl_B0>0, weights=distance_dgl, tau=0.95, fun=lognormal)
#'
#' #Raster data set
#' r1 <- terra::rast(nrows = 100, ncols = 100, res = 1)
#' rr <- setValues(r1, c(0,1,5,4,7,9,13))
#' plot(rr)
#'
#' #compute distance for prediction area
#' distance <- Distmap(fe_raster = rr, treespecies = "13")
#'
#' #prediction
#'predict_quax(distmap = distance, quax = f)
#'

predict_quax <- function(distmap, quax) {
  prediction <- quax(terra::values(distmap))
  terra::values(distmap) <- prediction
  return(distmap)
}
