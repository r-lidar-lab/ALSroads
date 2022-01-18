#' Align and measure a road from lidar data
#'
#' From a reference road (line), extracts the line with a buffer from the point cloud and computes
#' the exact positioning of the road (realignment). Then, using new the accurate  shape, computes
#' road metrics including its width, its drivable width, sinuosity as well as its state in four
#' classes (1. operating, 2. maybe operating, 3. maybe decommissioned, 4. decommissioned). The
#' function \link{st_snap_lines} allows to post-process the output to fix minor inaccuracies and
#' reconnect the roads than may no longer be connected because each road is processed independently.
#'
#' @param road a single line (sf format) used as reference to search and measure road.
#' @param roads multiples lines (sf format) used as reference to search and measure roads
#' @param ctg a non-normalized \link[lidR:LAScatalog-class]{LAScatalog} object from lidR package
#' @param dtm RasterLayer storing the DTM with a resolution of at least of 1 m. Can be computed
#' with \link[lidR:grid_terrain]{grid_terrain}.
#' @param param a list of many parameters. See \link{mffproads_default_parameters}.
#' @param water a set of polygons (sf format) of water bodies. This is used to mask the water bodies
#' so they cannot be mistaken as a drivable surfaces. Not mandatory but can help. It also allows to
#' detect bridges above water.
#' @param ... unused

#' @return An sf object similar to the input with additional attributes and an updated geometry. If
#' the class is 3 or 4 the original geometry is preserved to prevent adding more error. The new attributes
#' are ROADWITH, DRIVABLEWIDTH, PERCABOVEROAD (percentage of points between 0.5 and 5 meter above the road)
#' SHOULDERS (average number of shoulders found), SINUOSITY, CONDUCTIVITY (conductivity per linear meters)
#' SCORE (a road state score) and CLASS (4 classes derived from the SCORE). See references
#'
#' @references Roussel, J-R, Achim A (2022). Correction, update, and enhancement of vectorial forestry
#' road map using ALS data, a pathfinder and seven metrics.
#'
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
#' road <- st_read(road, "original", quiet = TRUE)
#' dtm  <- raster(dtm)
#'
#' # Voluntarily add more error to the road
#' crs <- st_crs(road)
#' st_geometry(road) <- st_geometry(road) + st_sfc(st_point(c(-8, 0)))
#' st_crs(road) <- crs
#'
#' plot(dtm, col = gray.colors(25, 0, 1))
#' plot(ctg, add = TRUE)
#' plot(st_geometry(road), add = TRUE, col = "red")
#'
#' res <- measure_road(ctg, road, dtm)
#' res
#' poly <- sf::st_buffer(res, res$ROADWIDTH/2)
#'
#' plot(st_geometry(road), col = "red") # Inaccurate road track
#' plot(st_geometry(res), col = "blue", add = TRUE) # Corrected road track
#'
#' domain <- "https://servicesmatriciels.mern.gouv.qc.ca:443"
#' path <- "/erdas-iws/ogc/wmts/Inventaire_Ecoforestier/Inventaire_Ecoforestier/default/"
#' tiles <- "GoogleMapsCompatibleExt2:epsg:3857/{z}/{y}/{x}.jpg"
#' url <- paste0(domain, path, tiles)
#' m = mapview::mapview(list(road, poly),
#'   layer.name = c("Inaccurate", "Corrected"),
#'   color = c("red", "blue"), map.type = "Esri.WorldImagery")
#' leaflet::addTiles(m@map, url)
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
measure_road = function(ctg, road, dtm, water = NULL, param = mffproads_default_parameters, ...)
{
  # Plenty of checks before to run anything
  dots <- list(...)
  lidR::opt_progress(ctg) <- getOption("MFFProads.debug.verbose")
  geometry_type_road <- sf::st_geometry_type(road)
  if (geometry_type_road != "LINESTRING") stop(glue::glue("Expecting LINESTRING geometry for 'road' but found {geometry_type_road} geometry instead."), call. = FALSE)
  if (nrow(road) > 1) stop("Expecting a single LINESTRING", call. = FALSE)
  if (!methods::is(ctg, "LAScatalog")) stop("Expecting a LAScatalog", call. = FALSE)
  if (!is.null(water)) { if (any(!sf::st_geometry_type(water) %in% c("MULTIPOLYGON", "POLYGON"))) stop("Expecting POLYGON geometry type for 'water'", call. = FALSE) }
  if (sf::st_is_longlat(road)) stop("Expecting a projected CRS for 'road' but found geographic CRS instead.", call. = FALSE)
  if (param[["extraction"]][["road_max_width"]] > param[["extraction"]][["road_buffer"]]) stop("'road_max_width' parameter must be smaller than 'road_buffer' parameter", call. = FALSE)

  # No warning when called from measure_roads because warning are thrown once only by measure_roads
  if (!isFALSE(dots$Windex))
  {
    # If the collection is not indexed we throw a warning and even an error if the density is high
    # Without spatial indexation the road extraction is terribly long
    if (!lidR::is.indexed(ctg))
    {
      d <- density(ctg)
      if (d < 5)
        message("No spatial index for LAS/LAZ files in this collection.")
      else if (d < 10)
        warning("No spatial index for LAS/LAZ files in this collection.", call. = FALSE)
      else
        stop("No spatial index for LAS/LAZ files in this collection.", call. = FALSE)
    }
  }

  if (getOption("MFFProads.debug.progress")) cat("Progress: ")

  # Check for weird roads. We already found a road with a 90 degrees node at the end of the road
  # and it broke the output
  angles <- st_angles(road)
  if (any(angles > 90))
  {
    if (any(angles[c(1, length(angles))] > 90))
      warning("Sharp turn (< 90 degrees) at one or both ends of the input road. This is weird and may lead to invalid outputs.", call. = FALSE)
    else
      warning("Sharp turn (< 90 degrees) between two consecutive segments of the input road. This is weird and may lead to invalid outputs.", call. = FALSE)
  }

  # This is the metrics we will estimate on the road. We generate a default output in case we should exit early
  new_road <- road
  new_road$ROADWIDTH     <- NA
  new_road$DRIVABLEWIDTH <- NA
  #new_road$RIGHTOFWAY    <- NA
  new_road$PERCABOVEROAD <- NA
  #new_road$PABOVE2       <- NA
  new_road$SHOULDERS     <- NA
  new_road$SINUOSITY     <- NA
  new_road$CONDUCTIVITY  <- NA
  new_road$SCORE         <- NA
  new_road$CLASS         <- 0

  # reorder the columns so outputs are consistent even if exiting early
  ngeom <- attr(new_road, "sf_column")
  names <- names(new_road)
  names <- names[names != ngeom]
  names <- append(names, ngeom)
  data.table::setcolorder(new_road, names)


  # Check the distance between start and end points of the road. If they are too close (e.g. closer than
  # the size of the buffer) we need to reduce the size of the buffer otherwise we are likely to do
  # something wrong.
  p1 <- lwgeom::st_startpoint(road)
  p2 <- lwgeom::st_endpoint(road)
  d  <- as.numeric(sf::st_distance(p1, p2))
  if (d < 2 * param[["extraction"]][["road_buffer"]] & !st_is_loop(road))
    param[["extraction"]][["road_buffer"]] <- (d / (2 * param[["extraction"]][["road_buffer"]] + 5) * 2 * param[["extraction"]][["road_buffer"]])/2

  # Get the units to display informative messages
  dist_unit  <- sf::st_crs(road)$units
  length_min <- 4*param[["extraction"]][["section_length"]]


  # Return the original geometry without computing anything if the road is too short to compute anything
  len <- as.numeric(sf::st_length(road))
  if (len < length_min) {
    warning(glue::glue("Road too short (< {length_min} {dist_unit}) to compute anything. Original road returned."), call. = FALSE)
    verbose("Done\n") ; cat("\n")
    return(new_road)
  }

  # Cut the road in subsections if it is too long or it is a loop
  cut <- floor(len/param[["extraction"]][["road_max_len"]])
  if (cut > 0) { message(sprintf("Long road detected. Splitting the roads in %d chunks of %d %s to process.", cut+1, round(sf::st_length(road)/(cut+1)), dist_unit)) }
  if (cut == 0 && st_is_loop(road)) { message(sprintf("Loop detected. Splitting the roads in 2 chunks of %d %s to process.", round(sf::st_length(road)/2,0), dist_unit)) ; cut = 1 }

  if (cut > 0)
  {
    # If we need to cut the road, the road is spitted and we recursively send each piece to this function
    cuts  <- seq(0,1, length.out = cut+2)
    from  <- cuts[-length(cuts)]
    to    <- cuts[-1]
    roads <- lapply(seq_along(from), function(i) { lwgeom::st_linesubstring(road, from[i], to[i]) })
    roads <- do.call(rbind, roads)
    res   <- measure_roads(ctg, roads, dtm, water, param)
    geom  <- st_merge_line(res)
    new_road$ROADWIDTH     <- mean(res$ROADWIDTH)
    new_road$DRIVABLEWIDTH <- mean(res$ROADWIDTH)
    #new_road$RIGHTOFWAY    <- mean(res$RIGHTOFWAY)
    new_road$PERCABOVEROAD <- mean(res$PERCABOVEROAD)
    #new_road$PABOVE2       <- mean(res$PABOVE2)
    new_road$SHOULDERS     <- mean(res$SHOULDERS)
    new_road$SINUOSITY     <- sinuosity(new_road)
    new_road$ROADWIDTH     <- mean(res$ROADWIDTH)
    new_road$CONDUCTIVITY  <- mean(res$CONDUCTIVITY)
    new_road$SCORE         <- road_score(new_road, param)
    new_road$CLASS         <- get_class(new_road$SCORE)
    if (!is.na(new_road$CLASS) && new_road$CLASS < 3) sf::st_geometry(new_road) <- geom
    verbose("Done\n") ; cat("\n")
    return(new_road)
  }

  # Query the roads in a collection of files
  las <- extract_road(ctg, road, param)

  # Exit early. This should never happen
  if (lidR::is.empty(las)) {
    warning("No point found.", call. = FALSE)
    verbose("Done\n") ; cat("\n")
    return(new_road)
  }

  # If we can't assume that the road is correctly positioned we need to recompute its
  # location accurately
  if (param[["constraint"]][["confidence"]] < 1)
  {
    # This is maybe the most important function of the code. It takes the point cloud, the original
    # road, the dtm, water bodies and param to draw a new accurate line
    res <- least_cost_path(las, road, dtm, water, param)

    # If the conductivity of the result is 0 it means that the path finder was not able to reach the
    # point B. We assume the road does not exist.
    if (res$CONDUCTIVITY == 0)
    {
      warning("Impossible to travel to the end of the road. Road does not exist.", call. = FALSE)
      new_road$ROADWIDTH     <- 0
      new_road$DRIVABLEWIDTH <- 0
      #new_road$RIGHTOFWAY    <- NA
      new_road$PERCABOVEROAD <- 100
      #new_road$PABOVE2       <- 100
      new_road$SHOULDERS     <- 0
      new_road$SINUOSITY     <- NA
      new_road$CONDUCTIVITY  <- 0
      new_road$SCORE         <- 0
      new_road$CLASS         <- 4
      verbose("Done\n") ; cat("\n")
      return(new_road)
    }

    # The new road geometry is very short. It is likely a bug that may arise for curved and super short roads
    # The case is handled to avoid failure but in practice it should not happen for regular road. It is an
    # exception
    if (as.numeric(sf::st_length(res)) < length_min)
    {
      warning(glue::glue("The computed road is too short (< {length_min} {dist_unit}) to compute anything. Original road returned."), call. = FALSE)
      new_road$ROADWIDTH     <- NA
      new_road$DRIVABLEWIDTH <- NA
      #new_road$RIGHTOFWAY    <- NA
      new_road$PERCABOVEROAD <- NA
      #new_road$PABOVE2       <- NA
      new_road$SHOULDERS     <- NA
      new_road$SINUOSITY     <- NA
      new_road$CONDUCTIVITY  <- NA
      new_road$SCORE         <- 0
      new_road$CLASS         <- 4
      verbose("Done\n") ; cat("\n")
      return(new_road)
    }

    new_road <- res
  }
  else
  {
    new_road$CONDUCTIVITY <- 1
  }

  # We now have an accurate road (hopefully). We can make measurement on it. This step extracts the
  # width profiles of the road, the percentage of points and relocate more accurately the centerline.
  slice_metrics <- road_measure(las, new_road, param)

  # I don't remember what it is. I guess it is an hidden debug options
  if (isFALSE(dots$reconstruct_line))
    return(slice_metrics)

  # Smooth the centerline using a spline adjustment
  if (param[["constraint"]][["confidence"]] < 1 && nrow(slice_metrics) > 4L)
  {
    spline <- adjust_spline(slice_metrics)
    spline <- sf::st_simplify(spline, dTolerance = 1)
    spline <- sf::st_set_crs(spline, sf::st_crs(las))
    sf::st_geometry(new_road) <- spline
  }

  # Aggregate metrics for the whole road from each segment. We have one metric per 10 meter slices
  # that are average to return a aggregated metrics
  metrics <- road_metrics(new_road, slice_metrics)
  metrics[["SCORE"]] <- road_score(metrics, param)
  metrics[["CLASS"]] <- get_class(metrics[["SCORE"]])

  # Merge the tables of attributes
  original_geometry <- sf::st_geometry(road)
  new_geometry <- sf::st_geometry(new_road)
  attribute_table <- cbind(sf::st_drop_geometry(road), metrics)
  attribute_table[[ngeom]] <- new_geometry

  # If the class of the road is 3 or 4 it is supposed to do not exist. Thus the centerline found is
  # likely to be irrelevant. We put back the original geometry to avoid doing worst than the original
  if (attribute_table[["CLASS"]] > 2)
    attribute_table[[ngeom]] <- original_geometry

  new_road <- sf::st_as_sf(attribute_table)

  # For a redrawn centerline, check if the ends of the
  # new road are suspiciously close to the edge of the caps
  if (attribute_table[["CLASS"]] %in% c(1,2))
  {
    start_ori <- lwgeom::st_startpoint(road)
    start_new <- lwgeom::st_startpoint(new_road)
    start_diff <- as.numeric(sf::st_distance(start_ori, start_new))

    end_ori <- lwgeom::st_endpoint(road)
    end_new <- lwgeom::st_endpoint(new_road)
    end_diff <- as.numeric(sf::st_distance(end_ori, end_new))

    if (max(start_diff, end_diff) > param[["extraction"]][["road_buffer"]] * 0.8) {
      warning("Road within 20% of the edge of an end cap radius at one or both ends. The computed road may have taken a shortcut through the woods.", call. = FALSE)
    }
  }

  verbose("Done\n") ; cat("\n")

  return(new_road)
}

#' @export
#' @rdname measure_road
measure_roads = function(ctg, roads, dtm, water = NULL, param = mffproads_default_parameters)
{
  if (!lidR::is.indexed(ctg))
  {
    d <- density(ctg)
    if (d < 5)
      message("No spatial index for LAS/LAZ files in this collection.")
    else if (d < 10)
      warning("No spatial index for LAS/LAZ files in this collection.", call. = FALSE)
    else
      stop("No spatial index for LAS/LAZ files in this collection.", call. = FALSE)
  }

  i <- 1:nrow(roads)
  res <- lapply(i, function(j)
  {
    if (getOption("MFFProads.debug.verbose") | getOption("MFFProads.debug.progress")) cat("Road", j, "of", nrow(roads), " ")
    measure_road(ctg, roads[j,], dtm, water, param, Windex = FALSE)
  })

  do.call(rbind, res)
}
