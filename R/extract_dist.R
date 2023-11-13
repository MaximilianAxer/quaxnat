##############################################################################
#' Extracting distances to nearest seed source for point data
#'
#' @description Extracts a distance for the inventory plots. The distance to the nearest seed source is used for the analysis of the regeneration potential.
#'
#'
#
#'
#' @param fe_raster Remote sensing raster dataset with tree species classification of specific tree species and tree species groups.
#' @param fe_geom Geodata set representing the study area. This can be a polygon or point dataset. This describes the outer boundary of the study area. A buffer of 1000 m is placed around the Bbox to possibly take into account seed trees outside the study area
#' @param treespecies Represents the numerical value by which the tree species of interest was encoded in the raster dataset.
#'
#' @details For each inventory plot a distance to the nearest seed source is given.
#'
#' @return  The distance is given in m.
#'



extract_dist <- function(fe_raster, fe_geom, treespecies){
  extract <- terra::extract(Distmap(fe_raster, fe_geom, treespecies), fe_geom)[,2]
    return(extract)
}


