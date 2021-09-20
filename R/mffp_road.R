#' Align and measure a road from Lidar data
#'
#' From a road (line), extracts the line with a buffer from the point cloud and recomputes
#' the actual positioning of the road and compute some metrics for the road such as its width
#' or drivable width as well as its state (exists, no longer exists)
#'
#' @param road a single line (sf format) used a reference to search and measure road
#' @param roads multiples lines (sf format) used a reference to search and measure roads
#' @param ctg a non-normalized \link[lidR:LAScatalog-class]{LAScatalog} object from lidR package
#' @param dtm RasterLayer storing the DTM with a resolution of at least of 1 m. Can be computed
#' with \link[lidR:grid_terrain]{grid_terrain}.
#' @param param a list of many parameters. See \link{mffproads_default_parameters}.
#' @param water a set of polygons (sf format) of water bodies. This is used to mask the water bodies
#' so they cannot be mistaken as a drivable surfaces. Not mandatory but can help. It also allows to
#' detect bridges.
#' @param confidence numeric. The confidence you have on the location of the reference road. 1 means
#' that the road is 100\% a ground truth. This will skip the relocation step. High values mean that your
#' are confident that reference road is correct and this will help the algorithm. Low values leave more
#' freedom to the algorithm but it becomes also more prone to errors. However this parameter is not
#' very sensitive. Default is 0.6
#'
#' @return a list with a several sf objects including stuff for debugging in \code{[["DEBUG"]]} +
#' stuff for end user in \code{[["OUTPUT"]]}
#' @export
#' @examples
#' library(lidR)
#' library(sf)
#' library(raster)
#'
#' dir  <- system.file("extdata", "", package="MFFProads")
#' road <- system.file("extdata", "road_971487.gpkg", package="MFFProads")
#' dtm  <- system.file("extdata", "dtm_1m.tif", package="MFFProads")
#' ctg  <- readLAScatalog(dir)
#' road <- st_read(road, quiet = TRUE)
#' dtm  <- raster(dtm)
#'
#' plot(dtm, col = gray.colors(25, 0, 1))
#' plot(ctg, add = TRUE)
#' plot(st_geometry(road), add = TRUE, col = "red")
#'
#' res <- measure_road(ctg, road, dtm)
#' res
#'
#' plot(st_geometry(road), col = "red") # Inaccurate road track
#' plot(st_geometry(res), col = "blue", add = TRUE) # Corrected road track
#'
#' mapview::mapview(list(road, res),
#'   layer.name = c("Inaccurate", "Corrected"),
#'   color = c("red", "blue"), map.type = "Esri.WorldImagery")
#'
#' \dontrun{
#' # DEBUG
#'
#' options(MFFProads.debug.finding = TRUE)
#' options(MFFProads.debug.measuring = FALSE)
#' options(MFFProads.debug.metrics = FALSE)
#' res <- measure_road(ctg, road, dtm)
#'
#' options(MFFProads.debug.finding = FALSE)
#' options(MFFProads.debug.measuring = TRUE)
#' options(MFFProads.debug.metrics = FALSE)
#' res <- measure_road(ctg, road, dtm)
#'
#' options(MFFProads.debug.finding = FALSE)
#' options(MFFProads.debug.measuring = FALSE)
#' options(MFFProads.debug.metrics = TRUE)
#' res <- measure_road(ctg, road, dtm)
#' }
#' @useDynLib MFFProads, .registration = TRUE
#' @import data.table
measure_road = function(ctg, road, dtm, water = NULL, confidence = 0.7, param = mffproads_default_parameters)
{
  if (sf::st_geometry_type(road) != "LINESTRING") stop("Expecting LINESTRING geometry for 'road'", call. = FALSE)
  if (nrow(road) > 1) stop("Expecting a single LINESTRING", call. = FALSE)
  if (!methods::is(ctg, "LAScatalog")) stop("Expecting a LAscatalog", call. = FALSE)
  lidR:::assert_is_a_number(confidence)
  lidR:::assert_all_are_in_closed_range(confidence, 0, 1)
  if (!is.null(water)) { if (any(!sf::st_geometry_type(water) %in% c("MULTIPOLYGON", "POLYGON"))) stop("Expecting POLYGON geometry type for 'water'", call. = FALSE) }

  param$search$confidence <- confidence

  # This is the metrics we will estimate on the road. We generate a default output in case we should exit early
  new_road <- road
  new_road$ROADWIDTH     <- NA
  new_road$DRIVABLEWIDTH <- NA
  new_road$RIGHTOFWAY    <- NA
  new_road$PABOVE05      <- NA
  new_road$PABOVE2       <- NA
  new_road$SINUOSITY     <- NA
  new_road$ROADWIDTH     <- NA
  new_road$CONDUCTIVITY  <- NA
  new_road$SCORE         <- NA
  new_road$STATE         <- 0

  # reorder the columns so outputs are consistent even if exiting early
  geom <- attr(new_road, "sf_column")
  names <- names(new_road)
  names <- names[names != geom]
  names <- append(names, geom)
  data.table::setcolorder(new_road, names)

  # Exit early for loops because does not work (yet?)
  p1 <- lwgeom::st_startpoint(road)
  p2 <- lwgeom::st_endpoint(road)
  d  <- as.numeric(sf::st_distance(p1,p2)[1,1])
  if (d < 2) {
    warning("Loop roads are not yet supported. Process aborted.", call. = FALSE)
    return(new_road)
  }

  # Query the roads in a collection of files
  las <- extract_road(ctg, road, param)

  # Exit early. This should never happen
  if (lidR::is.empty(las)) {
    warning("No point found.", call. = FALSE)
    return(new_road)
  }

  # If we can't assume that the road is correctly positioned we need to recompute its
  # location accurately
  if (confidence < 1)
    new_road <- road_relocate(las, road, dtm, water, param)
  else
    new_road$CONDUCTIVITY <- 1

  # We now have an accurate road (hopefully). We can make measurement on it
  param$extraction$road_buffer <- 30
  segment_metrics <- road_measure(las, new_road, param)
  segment_metrics <- sf::st_as_sf(segment_metrics, coords = c("xroad", "yroad"), crs = sf::st_crs(las))

  # We can also improve the coarse measurement given by least cost path
  if (confidence < 1)
  {
    spline <- adjust_spline(segment_metrics)
    spline <- sf::st_simplify(spline, dTolerance = 1)
    spline <- sf::st_set_crs(spline, sf::st_crs(las))
    sf::st_geometry(new_road) <- spline
  }

  # Aggregate metrics for the whole road from each segment
  metrics <- road_metrics(new_road, segment_metrics)
  res <- road_state(segment_metrics, metrics$CONDUCTIVITY, param)
  metrics[["SCORE"]] <- res[[2]]
  metrics[["STATE"]] <- res[[1]]

  # Merge the tables of attributes
  original_geometry <- sf::st_geometry(road)
  new_geometry <- sf::st_geometry(new_road)
  attribute_table <- cbind(sf::st_drop_geometry(road), metrics)
  attribute_table[["geom"]] <- new_geometry

  if (attribute_table[["STATE"]] > 2)
    attribute_table[["geom"]] <- original_geometry

  new_road <- sf::st_as_sf(attribute_table)
  return(new_road)
}

#' @export
#' @rdname measure_road
measure_roads = function(ctg, roads, dtm, water = NULL, confidence = 0.7, param = mffproads_default_parameters)
{
  i <- 1:nrow(roads)
  res <- lapply(i, function(x)
  {
    cat("Road", i, "of", nrow(roads), "\n")
    measure_road(ctg, roads[x,], dtm, water, confidence, param)
  })

  do.call(rbind, res)
}
