##############################################################################
#' Extracting distances to nearest seed source for point data
#'
#' @description Extracts a distance for the inventory plots. The distance to the
#'  nearest seed source is used for the analysis of the regeneration potential.
#'
#' @param raster Remote sensing raster dataset with tree species
#' classification of specific tree species and tree species groups.
#'
#' @param geom Geodata set representing the study area. This can be a polygon
#'  or point dataset. This describes the outer boundary of the study area.
#'  A buffer of 1000 m is placed around the Bbox to possibly take into account
#'  seed trees outside the study area
#'
#' @param species Represents the numerical value by which the tree species
#' of interest was encoded in the raster dataset.
#'
#' @details For each inventory plot a distance to the nearest seed source of a
#'  specific tree species is given.
#'
#' @return  The distance is given in m.
#'
#'@examples
#' ## Create raster data set
#' set.seed(2023)
#' rr <- terra::rast(
#' matrix(sample(0:10, 20 * 20, replace = TRUE),
#'         nrow = 20, ncol = 20))
#'
#' ## Create vector data set
#' vec <- terra::vect(rbind(c(5,10), c(5,15)))
#'
#' ## Extract distance for the inventory plot
#' extract_dist(raster=rr, geom=vec, species=10)
#'
#' @export



extract_dist <- function(raster, geom, species){
  extract <- terra::extract(Distmap(raster, species), geom)[,2]
  return(extract)
}

# Old Implementation
#extract_dist <- function(raster, geom, species){
#  extract <- terra::extract(Distmap(raster, geom, species), geom)[,2]
#  return(extract)
#}
#

