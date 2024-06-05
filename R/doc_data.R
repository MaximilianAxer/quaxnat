#' Regeneration densities at inventory plots and potential dispersal distances
#' to nearest seed trees
#'
#' A dataset containing the regeneration densities of beech, oak and Douglas fir
#' of the inventory plots and the distance to the nearest conspecific nearest
#' seed tree.

#' @format A data frame with 484 rows and 7 variables
#' \itemize{
#'   \item id. An identifier for each inventory plot as an integer
#'   \item distance_beech. Distance in m from the plot to the nearest beech (0--3206.57)
#'   \item distance_oak. Distance in m from the plot to the nearest oak (0--1481.2)
#'   \item distance_dgl. Distance in m from the plot to the nearest Douglas fir (0--1807)
#'   \item oak_regen. Regeneration density of oak (N/ha) of the plot (0--30)
#'   \item beech_regen. Regeneration density of beech (N/ha) of the plot (0--30)
#'   \item douglas_regen. Regeneration density of Douglas fir (N/ha) of the plot (0--30)
#' }
#'
#' @docType data
#' @keywords datasets
#' @name regeneration
#' @usage data(regeneration)
"regeneration"
