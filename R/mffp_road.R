#' Align and measure a road from lidar data
#'
#' From a road (line), extracts the line with a buffer from the point cloud and recomputes
#' the actual positioning of the road and compute some metrics for the road such as its width
#' or drivable width as well as its state (exists, no longer exists)
#'
#' @param road a single line (sf format) used as reference to search and measure road
#' @param roads multiples lines (sf format) used as reference to search and measure roads
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
#' very sensitive.
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
  if (!lidR::is.indexed(ctg)) message("No spatial index for LAS/LAZ files in this collection.")

  param$search$confidence <- confidence

  # This is the metrics we will estimate on the road. We generate a default output in case we should exit early
  new_road <- road
  new_road$ROADWIDTH     <- NA
  new_road$DRIVABLEWIDTH <- NA
  new_road$RIGHTOFWAY    <- NA
  new_road$PABOVE05      <- NA
  new_road$PABOVE2       <- NA
  new_road$SHOULDERS     <- NA
  new_road$SINUOSITY     <- NA
  new_road$CONDUCTIVITY  <- NA
  new_road$SCORE         <- NA
  new_road$STATE         <- 0

  # reorder the columns so outputs are consistent even if exiting early
  ngeom <- attr(new_road, "sf_column")
  names <- names(new_road)
  names <- names[names != ngeom]
  names <- append(names, ngeom)
  data.table::setcolorder(new_road, names)

  # Cut the road is too long or lopp
  len <- as.numeric(sf::st_length(road))
  if (len < 50) {
    warning("Too short road to compute anything.", call. = FALSE)
    return(new_road)
  }

  cut <- floor(len/param[["extraction"]][["road_max_len"]])
  if (cut > 0) { message(sprintf("Long road detected. Splitting the roads in %d chunks of %d m to process.", cut+1, round(sf::st_length(road)/(cut+1)))) }
  if (cut == 0 && st_is_loop(road)) { message(sprintf("Loop detected. Splitting the roads in 2 chunks of %d m to process.", round(sf::st_length(road)/2,0))) ; cut = 1 }

  if (cut > 0)
  {
    cuts  <- seq(0,1, length.out = cut+2)
    from  <- cuts[-length(cuts)]
    to    <- cuts[-1]
    roads <- lapply(seq_along(from), function(i) { lwgeom::st_linesubstring(road, from[i], to[i]) })
    roads <- do.call(rbind, roads)
    res   <- measure_roads(ctg, roads, dtm, water, confidence, param)
    geom  <- st_merge_line(res)
    new_road$ROADWIDTH     <- mean(res$ROADWIDTH)
    new_road$DRIVABLEWIDTH <- mean(res$ROADWIDTH)
    new_road$RIGHTOFWAY    <- mean(res$RIGHTOFWAY)
    new_road$PABOVE05      <- mean(res$PABOVE05)
    new_road$PABOVE2       <- mean(res$PABOVE2)
    new_road$SHOULDERS     <- mean(res$SHOULDERS)
    new_road$SINUOSITY     <- sinuosity(new_road)
    new_road$ROADWIDTH     <- mean(res$ROADWIDTH)
    new_road$CONDUCTIVITY  <- mean(res$CONDUCTIVITY)
    new_road$SCORE         <- road_score(new_road, param)
    new_road$STATE         <- get_state(new_road$SCORE)
    if (new_road$STATE < 3) sf::st_geometry(new_road) <- geom
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
  {
    res <- road_relocate(las, road, dtm, water, param)

    if (res$CONDUCTIVITY == 0)
    {
      warning("Impossible to travel to the end of the road. Road does not exist.", call. = FALSE)
      new_road$ROADWIDTH     <- 0
      new_road$DRIVABLEWIDTH <- 0
      new_road$RIGHTOFWAY    <- NA
      new_road$PABOVE05      <- 100
      new_road$PABOVE2       <- 100
      new_road$SHOULDERS     <- 0
      new_road$SINUOSITY     <- NA
      new_road$CONDUCTIVITY  <- 0
      new_road$SCORE         <- 0
      new_road$STATE         <- 4
      return(new_road)
    }

    new_road <- res
  }
  else
  {
    new_road$CONDUCTIVITY <- 1
  }

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
  metrics[["SCORE"]] <- road_score(metrics, param)
  metrics[["STATE"]] <- get_state(metrics[["SCORE"]])

  # Merge the tables of attributes
  original_geometry <- sf::st_geometry(road)
  new_geometry <- sf::st_geometry(new_road)
  attribute_table <- cbind(sf::st_drop_geometry(road), metrics)
  attribute_table[[ngeom]] <- new_geometry

  if (attribute_table[["STATE"]] > 2)
    attribute_table[[ngeom]] <- original_geometry

  new_road <- sf::st_as_sf(attribute_table)
  return(new_road)
}

#' @export
#' @rdname measure_road
measure_roads = function(ctg, roads, dtm, water = NULL, confidence = 0.7, param = mffproads_default_parameters)
{
  i <- 1:nrow(roads)
  res <- lapply(i, function(j)
  {
    cat("Road", j, "of", nrow(roads), "\n")
    measure_road(ctg, roads[j,], dtm, water, confidence, param)
  })

  do.call(rbind, res)
}
