##############################################################################
#' Distance map to nearest seed source
#'
#' @description Creation of a distance map for the study area. The distance to
#' the nearest seed source is calculated for every raster cell.
#'
#'
#'
#'
#' @param raster Remote sensing raster data set with tree species
#' classification of specific tree species and tree species groups.
#' @param treespecies Represents the numerical value by which the tree species
#'  of interest is encoded in the raster data set.
#'
#' @details A SpatRaster is created with the same resolution as the input raster,
#'  defined by the study area.
#'
#' @return  The distance is given in m.
#'
#'@examples
#' ## Create raster data set
#'  rr <- terra::rast(
#'  matrix(sample(0:10, 20 * 20, replace = TRUE),
#'         nrow = 20, ncol = 20))
#'
#' ## Compute distance for study area
#' distance <- Distmap(raster = rr, treespecies = "10")
#'
#' ## Plot the Distmap
#' terra::plot(distance)
#'
#' @export

#Distmap <- function(raster, geom, treespecies){
#  window <- terra::segregate(terra::crop(raster, sf::st_buffer(sf::st_union(geom), dist = 1000)))[[treespecies]]
#  distmap <- terra::distance(window, target = 0)
#  return(distmap)
#}
#geom Geodata set representing the study area. This can be a polygon or point dataset. It describes the outer boundary of the study area. A buffer of 1000 m is placed around the Bbox to possibly take into account seed trees outside the study area


Distmap <- function(raster, treespecies){
  window <- terra::segregate(raster)[[treespecies]]
  distmap <- terra::distance(window, target = 0)
  return(distmap)
}
