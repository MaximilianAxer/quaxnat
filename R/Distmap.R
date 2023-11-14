##############################################################################
#' Distance map to nearest seed source
#'
#' @description Creation of a distance map for the study area. The distance to
#' the nearest seed source is calculated for every raster cell.
#'
#'
#'
#'
#' @param fe_raster Remote sensing raster dataset with tree species classification of specific tree species and tree species groups.
#' @param treespecies Represents the numerical value by which the tree species of interest was encoded in the raster dataset.
#'
#' @details A SpatRaster is created with the same resolution as the input raster, defined by the study area.
#'
#' @return  The distance is given in m.
#'
#'@examples
#' ## Raster data set
#' rr <- terra::rast(
#'  matrix(sample(0:10, 20 * 20, replace = TRUE),
#'         nrow = 20, ncol = 20))
#'
#' ## compute distance for prediction area
#' distance <- Distmap(fe_raster = rr, treespecies = "10")
#' terra::plot(distance)
#'
#' @export

#Distmap <- function(fe_raster, fe_geom, treespecies){
#  window <- terra::segregate(terra::crop(fe_raster, sf::st_buffer(sf::st_union(fe_geom), dist = 1000)))[[treespecies]]
#  distmap <- terra::distance(window, target = 0)
#  return(distmap)
#}
#fe_geom Geodata set representing the study area. This can be a polygon or point dataset. It describes the outer boundary of the study area. A buffer of 1000 m is placed around the Bbox to possibly take into account seed trees outside the study area


Distmap <- function(fe_raster, treespecies){
  window <- terra::segregate(fe_raster)[[treespecies]]
  distmap <- terra::distance(window, target = 0)
  return(distmap)
}
