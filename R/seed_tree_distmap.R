##############################################################################
#' Distance map to nearest seed source
#'
#' @description Creation of a distance map for the study area. The distance to
#'  the nearest seed source is calculated for every raster cell.
#'
#' @param raster Raster data set with tree species classification of specific
#'  tree species and tree species groups.
#' @param species Represents the numerical value by which the tree species of
#'  interest is encoded in the raster data set.
#'
#' @return A SpatRaster object containing the distances to seed source.
#'  The object has the same resolution and extent as the input raster.
#'
#' @export
#'
#' @examples
#' ## Create raster data set
#'  rr <- terra::rast(
#'  matrix(sample(0:10, 20 * 20, replace = TRUE),
#'         nrow = 20, ncol = 20))
#'
#' ## Compute distance for study area
#' distance <- seed_tree_distmap(raster = rr, species = "10")
#'
#' ## Plot the seed_tree_distmap
#' terra::plot(distance)
#'

#seed_tree_distmap <- function(raster, geom, species){
#  window <- terra::segregate(terra::crop(raster, sf::st_buffer(sf::st_union(geom), dist = 1000)))[[species]]
#  distmap <- terra::distance(window, target = 0)
#  return(distmap)
#}
#geom Geodata set representing the study area. This can be a polygon or point dataset. It describes the outer boundary of the study area. A buffer of 1000 m is placed around the Bbox to possibly take into account seed trees outside the study area


seed_tree_distmap <- function(raster, species){
  window <- terra::segregate(raster)[[species]]
  distmap <- terra::distance(window, target = 0)
  return(distmap)
}
