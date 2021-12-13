#' Check for missing road junctions
#'
#' Check if some road endings are close enough to be considered as being potentially part of a junction.
#' Road endings that are actually touching are not considered. Fixing can be done by running \link{st_snap_lines}
#' with the same \code{tolerance} value and the same road network.
#'
#' @param roads  multiples lines (\code{sf} format)
#' @param tolerance numeric (distance unit). Tolerance value used to check if some road endings should have been snapped together.
#'
#' @return Polygon object (\code{sf} format) corresponding to circles of \code{tolerance} radius around the locations where
#' junctions might be missing.
#' @export
#' @examples
#' library(sf)
#'
#' network <- system.file("extdata", "road_network.gpkg", package="MFFProads")
#'
#' roads <- sf::st_read(network, "invalid_topology")
#' 
#' junctions <- st_check_junctions(roads)
#'
#' par(mfrow = c(1,2))
#' for (k in 1:2)
#' {
#'   junction <- junctions[k,]
#'   bbox <- st_bbox(junction)
#'   plot(st_geometry(junctions), col = "cyan", xlim = c(bbox[c("xmin", "xmax")]), ylim = c(bbox[c("ymin", "ymax")]))
#'   plot(st_geometry(roads), col = "red", lwd = 3, add = TRUE)
#'   title(main = k,
#'         sub = paste0("\nConnected: ", junction[["nb_connected"]],
#'                      "\nUnconnected: ", junction[["nb_unconnected"]]))
#' }
st_check_junctions = function(roads, tolerance = 8)
{
  end   <- lwgeom::st_endpoint(roads)
  start <- lwgeom::st_startpoint(roads)
  ends  <- c(start, end)
  
  group_pts <- sf::st_is_within_distance(ends, ends, tolerance) |>
    lapply(sort) |>
    Filter(f = function(x) { length(x) > 1 }) |>
    unique() |>
    lapply(function(x) { ends[x] })
  
  if (length(group_pts) == 0)
  {
    warning("No junctions found using current tolerance value.", .call = FALSE)
    return(NULL)
  }

  mean_pts <- group_pts |>
    lapply(function(x) { sf::st_point(colMeans(sf::st_coordinates(x))) }) |>
    sf::st_sfc(crs = sf::st_crs(roads))
  
  junctions <- group_pts |>
    lapply(function(x) {
      distances <- as.numeric(sf::st_distance(x, x)[1,])
      actual_junc <- sum(distances == 0)
      potential_junc <- length(distances) - actual_junc
      c(actual_junc, potential_junc)
    }) |>
    do.call(what = rbind)

  has_potential <- junctions[,2] > 0
  
  buffer_junctions <- data.frame(
      nb_connected   = junctions[has_potential, 1],
      nb_unconnected = junctions[has_potential, 2]
    ) |>
    sf::st_sf(mean_pts[has_potential]) |>
    sf::st_buffer(tolerance)
  
  return(buffer_junctions)
}


#' Check for crossing roads
#'
#' Check if some roads are crossing elsewhere than at their ends. These occurences should be fixed as roads crossing
#' each other should be split at their intersection point in a valid topological network.
#'
#' @param roads  multiples lines (\code{sf} format)
#'
#' @return Point object (\code{sf} format) corresponding to intersection points of crossing roads.
#' @export
#' @examples
#' library(sf)
#'
#' network <- system.file("extdata", "road_network.gpkg", package="MFFProads")
#'
#' roads <- sf::st_read(network, "invalid_topology")
#' 
#' crossings <- st_check_crossings(roads)
#'
#' plot(st_geometry(roads), col = "red") # Inspected roads
#' plot(st_geometry(crossings), col = "blue", add = TRUE) # Crossings to fix
st_check_crossings = function(roads)
{
  geom <- sf::st_geometry(roads)
  
  group_cross <- sf::st_crosses(geom, geom) |>
    mapply(1:nrow(roads), FUN = function(x, i) { c(i, x) }) |>
    Filter(f = function(x) { length(x) > 1})
  
  if (length(group_cross) == 0) return(NULL)
  
  # Intersection points are cast to MULTIPOINT to allow
  # extraction of coordinates of a single geometry type
  # (no POINT and MULTIPOINT mixing)
  coords <- group_cross |>
    lapply(function(x) {
      sf::st_intersection(geom[x[1]], geom[x[-1]]) |>
      sf::st_cast("MULTIPOINT") |>
      sf::st_coordinates()
    }) |>
    do.call(what = rbind)
  
  unique <- !duplicated(coords[,1:2])
  
  pts_crossings <- coords[unique, 1:2] |>
    as.data.frame() |>
    sf::st_as_sf(coords = c(1,2), crs = sf::st_crs(roads))
  
  return(pts_crossings)
}
