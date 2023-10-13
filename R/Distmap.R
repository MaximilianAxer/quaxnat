#' Distmap
#'
#' @description Creation of a distance map for the study area.
#'
#'
#'
#' @param fe_raster Remote sensing raster dataset with tree species classification of specific tree species and tree species groups.
#' @param fe_geom Geodata set representing the study area. This can be a polygon or point dataset. It describes the outer boundary of the study area. A buffer of 1000 m is placed around the Bbox to possibly take into account seed trees outside the study area
#' @param treespecies Represents the numerical value by which the tree species of interest was encoded in the raster dataset.
#'
#' @details A SpatRaster is created with the same resolution as the input raster, defined by the study area.
#'
#' @return  The distance is given in m.
#'
#'#' @examples
#'
#' #Raster data set
#' r1 <- terra::rast(nrows = 100, ncols = 100, res = 1)
#' rr <- setValues(r1, c(0,1,5,4,7,9,13))
#' plot(rr)
#'
#' #compute distance for prediction area
#' distance <- Distmap(fe_raster = rr, treespecies = "13")

Distmap <- function(fe_raster, fe_geom, treespecies){
  window <- terra::segregate(terra::crop(fe_raster, sf::st_buffer(sf::st_union(fe_geom), dist = 1000)))[[treespecies]]
  distmap <- terra::distance(window, target = 0)
  return(distmap)
}


