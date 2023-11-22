#' Align and measure a road from lidar data
#'
#' From a reference road (spatial line), extracts the line with a buffer from the point cloud and computes
#' the exact positioning of the road (realignment). Then, using new the accurate shape, computes
#' road metrics including its width, its drivable width, its sinuosity as well as its state in four
#' classes. The function \link{st_snap_lines} allows to post-process the output to fix minor inaccuracies
#'  and reconnect the roads that may no longer be connected because each road is processed independently.
#'
#' @param road a single linestring (sf format) used as reference to search and measure the road.
#' @param roads multiple lines (sf format) used as reference to search and measure the roads
#' @param ctg a non-normalized \link[lidR:LAScatalog-class]{LAScatalog} object from lidR package
#' @param dtm RasterLayer storing the DTM with a resolution of at least of 1 m. Can be computed
#' with \link[lidR:grid_terrain]{grid_terrain}. It can be missing if a conductivity layer is provided
#' @param conductivity RasterLayer storing the pre-computed conductivity. It can be NULL in this case
#' if will be computed on the fly but the layer can be pre-computed with \link{rasterize_conductivity}
#' @param param a list of many parameters. See \link{alsroads_default_parameters}.
#' @param water a set of spatial polygons (sf format) of water bodies. This is used to mask the water
#' bodies so they cannot be mistaken as a drivable surfaces. Not mandatory but can help. It also allows
#' to detect bridges above water.
#' @param ... unused

#' @return An sf object similar to the input with additional attributes and an updated geometry. If
#' the class is 3 or 4 the original geometry is preserved to prevent adding more error. The new attributes
#' are ROADWITH, DRIVABLEWIDTH, PERCABOVEROAD (percentage of points between 0.5 and 5 meter above the road)
#' SHOULDERS (average number of shoulders found), SINUOSITY, CONDUCTIVITY (conductivity per linear meters)
#' SCORE (a road state score) and CLASS (4 classes derived from the SCORE). See references
#'
#' @references Roussel, J.-R., Bourdon, J.-F., Morley, I. D., Coops, N. C., & Achim, A. (2022).
#' Correction , update , and enhancement of vectorial forestry road maps using ALS data a pathfinder
#' and seven metrics. International Journal of Applied Earth Observation and Geoinformation, 114(September),
#' 103020. https://doi.org/10.1016/j.jag.2022.103020
#' @export
#' @examples
#' library(lidR)
#' library(sf)
#' library(raster)
#'
#' dir  <- system.file("extdata", "", package="ALSroads")
#' road <- system.file("extdata", "j5gr_centerline_971487.gpkg", package="ALSroads")
#' dtm  <- system.file("extdata", "j5gr_dtm.tif", package="ALSroads")
#' ctg  <- readLAScatalog(dir)
#' road <- st_read(road, "original", quiet = TRUE)
#' dtm  <- raster(dtm)
#'
#' # Voluntarily add more error to the road
#' crs <- st_crs(road)
#' st_geometry(road) <- st_geometry(road) + st_sfc(st_point(c(-8, 0)))
#' st_crs(road) <- crs
#'
#' plot(dtm, col = gray(1:50/50))
#' plot(ctg, add = TRUE)
#' plot(st_geometry(road), add = TRUE, col = "red")
#'
#' res <- measure_road(ctg, road, dtm = dtm)
#' res
#' poly <- sf::st_buffer(res, res$ROADWIDTH/2)
#'
#' plot(dtm, col = gray(1:50/50))
#' plot(st_geometry(road), col = "red", add = TRUE) # Inaccurate road track
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
#' conductivity <- system.file("extdata", "j5gr_conductivity.tif", package="ALSroads")
#' conductivity <- raster(conductivity)
#' plot(conductivity, col = viridis::viridis(50))
#'
#' res <- measure_road(ctg, road, conductivity = conductivity)
#'
#' plot(st_geometry(road), col = "red") # Inaccurate road track
#' plot(st_geometry(res), col = "blue", add = TRUE) # Corrected road track
#' }
#' @useDynLib ALSroads, .registration = TRUE
#' @import data.table
measure_road = function(ctg, centerline, dtm = NULL, conductivity = NULL, water = NULL, param = alsroads_default_parameters, ...)
{
  # Plenty of checks before to run anything
  dots <- list(...)
  lidR::opt_progress(ctg) <- getOption("ALSroads.debug.verbose")
  geometry_type_road <- sf::st_geometry_type(centerline)
  if (geometry_type_road != "LINESTRING") stop(glue::glue("Expecting LINESTRING geometry for 'centerline' but found {geometry_type_road} geometry instead."), call. = FALSE)
  if (nrow(centerline) > 1) stop("Expecting a single LINESTRING", call. = FALSE)
  if (!methods::is(ctg, "LAScatalog")) stop("Expecting a LAScatalog", call. = FALSE)
  if (!is.null(water)) { if (any(!sf::st_geometry_type(water) %in% c("MULTIPOLYGON", "POLYGON"))) stop("Expecting POLYGON geometry type for 'water'", call. = FALSE) }
  if (sf::st_is_longlat(centerline)) stop("Expecting a projected CRS for 'centerline' but found geographic CRS instead.", call. = FALSE)
  if (param[["extraction"]][["road_max_width"]] > param[["extraction"]][["road_buffer"]]) stop("'road_max_width' parameter must be smaller than 'road_buffer' parameter", call. = FALSE)
  if (is.null(dtm) & is.null(conductivity)) stop("'dtm' and 'conductivity' cannot be both NULL", call. = FALSE)
  if (!is.null(dtm) & !is.null(conductivity)) stop("'dtm' or 'conductivity' must be NULL", call. = FALSE)
  crs <- sf::st_crs(centerline)

  # If the collection is not indexed we throw a warning and even an error if the density is high
  # Without spatial indexation the centerline extraction is terribly long
  # (No warning when called from measure_roads because warning are thrown once only by measure_roads)
  if (!isFALSE(dots$Windex)) alert_no_index(ctg)

  # Display progress
  if (getOption("ALSroads.debug.progress")) cat("Progress: ")

  # Check for weird roads. We already found a road with a 90 degrees node
  # at the end of the road and it broke the output
  warn_weird_road(centerline)

  # We need to handle loop in special way so we define a bool
  is_loop <- st_is_loop(centerline)

  # If the confidence on the location of road is 1 we do not need to relocate
  # We know that the input is already perfect and we only need to measure the road
  relocate <- param[["constraint"]][["confidence"]] < 1

  # Get the units to display informative messages
  dist_unit  <- crs$units

  # We generate a default output in case we should exit early the function.
  # Basically a road with the original geometry and NA metrics
  new_road <- road_class0(centerline)

  # Check the distance between start and end points of the road. If they are too close (e.g. closer than
  # the size of the buffer) we need to reduce the size of the buffer otherwise we are likely to do
  # something wrong.
  if (!is_loop)
  {
    bu <- param[["extraction"]][["road_buffer"]]
    p1 <- lwgeom::st_startpoint(centerline)
    p2 <- lwgeom::st_endpoint(centerline)
    d  <- as.numeric(sf::st_distance(p1, p2))
    if (d < 2 * param[["extraction"]][["road_buffer"]])
      param[["extraction"]][["road_buffer"]] <- (d / (2 * bu + 5) * 2 * bu)/2
  }

  # Return the original geometry without computing anything if the centerline is too short to compute anything
  length_min <- 4*param[["extraction"]][["section_length"]]
  len <- as.numeric(sf::st_length(centerline))
  if (len < length_min)
  {
    warning(glue::glue("Road too short (< {length_min} {dist_unit}) to compute anything. Original road returned."), call. = FALSE)
    verbose("Done\n") ; cat("\n")
    return(new_road)
  }

  # Estimate if the road must be split in subsections. This have to goals
  # - It is an optimization to avoid processing giant conductivity raster from looong road
  # - It allows to handle loop roads cases
  cut <- floor(len/param[["extraction"]][["road_max_len"]])
  if (cut > 0) { message(sprintf("Long road detected. Splitting the roads in %d chunks of %d %s to process.", cut+1, round(len/(cut+1)), dist_unit)) }
  if (cut == 0 & is_loop) { message(sprintf("Loop detected. Splitting the roads in 2 chunks of %d %s to process.", round(len/2,0), dist_unit)) ; cut = 1 }

  # If a split is required do the cut
  if (cut > 0)
  {
    # If we need to cut the road, the road is spitted and we recursively send each
    # piece recursively into this function
    cuts  <- seq(0,1, length.out = cut+2)
    from  <- cuts[-length(cuts)]
    to    <- cuts[-1]
    roads <- lapply(seq_along(from), function(i) { lwgeom::st_linesubstring(centerline, from[i], to[i]) })
    roads <- do.call(rbind, roads)
    res   <- measure_roads(ctg, roads, dtm, conductivity, water, param)
    geom  <- st_merge_line(res)

    # Because we split the line we have multiple independent results
    # We recombine the multiple output in a single one
    new_road$ROADWIDTH     <- mean(res$ROADWIDTH)
    new_road$DRIVABLEWIDTH <- mean(res$ROADWIDTH)
    new_road$PERCABOVEROAD <- mean(res$PERCABOVEROAD)
    new_road$SHOULDERS     <- mean(res$SHOULDERS)
    new_road$SINUOSITY     <- sinuosity(new_road)
    new_road$ROADWIDTH     <- mean(res$ROADWIDTH)
    new_road$CONDUCTIVITY  <- mean(res$CONDUCTIVITY)
    new_road$SCORE         <- road_score(new_road, param)
    new_road$CLASS         <- get_class(new_road$SCORE)
    if (!is.na(new_road$CLASS) && new_road$CLASS < 3) sf::st_geometry(new_road) <- geom
    verbose("Done\n") ; cat("\n")

    # We exit the function. The next code being the regular case we no splitting
    new_road <- rename_sf_column(new_road, centerline)

    return(new_road)
  }

  # Query the roads in a collection of files
  las <- extract_road(ctg, centerline, param)

  # Exit early. This should never happen
  if (lidR::is.empty(las))
  {
    warning("No point found.", call. = FALSE)
    verbose("Done\n") ; cat("\n")
    return(new_road)
  }

  # If we can't assume that the road is correctly positioned we need to recompute its
  # location accurately
  if (relocate)
  {
    # This is maybe the most important function of the code. It takes the point cloud, the original
    # road, the dtm or conductivity layer, water bodies and param to draw a new accurate line
    res <- least_cost_path(las, centerline, dtm, conductivity, water, param)

    if (sf::st_geometry_type(res) == "MULTILINESTRING")
      stop("Internal error. A MULTILINESTRING has been returned. Please report", call. = FALSE)

    # If the conductivity of the result is 0 it means that the path finder was not able to reach the
    # point B. We assume the road does not exist.
    if (res$CONDUCTIVITY == 0)
    {
      warning("Impossible to travel to the end of the road. Road does not exist.", call. = FALSE)
      new_road <- road_class4(centerline)
      verbose("Done\n") ; cat("\n")
      return(new_road)
    }

    # The new centerline geometry is very short? It is likely a bug that may arise for curved
    # and super short roads. The case is handled to avoid failure but in practice it should
    # not happen for regular road. It is an exception that must be handled
    if (as.numeric(sf::st_length(res)) < length_min)
    {
      warning(glue::glue("The computed road is too short (< {length_min} {dist_unit}) to compute anything. Original road returned."), call. = FALSE)
      new_road <- road_class0(centerline)
      new_road$CLASS <- 4
      verbose("Done\n") ; cat("\n")
      return(new_road)
    }

    # Eventually, if the result exist and is not too short the new road is the
    # one returned by least_cost_path
    new_road <- res
  }
  else
  {
    new_road$CONDUCTIVITY <- 1
  }

  # We now have an accurate road (hopefully). We can make measurement on it.
  # This step extracts the width profiles of the road, the percentage of points
  # and relocate more accurately the centerline.
  slice_metrics <- road_measure(las, new_road, param)

  # Smooth the centerline using a spline adjustment
  if (relocate && nrow(slice_metrics) > 4L)
  {
    spline <- adjust_spline(slice_metrics)
    spline <- sf::st_simplify(spline, dTolerance = 1)
    spline <- sf::st_set_crs(spline, crs)
    sf::st_geometry(new_road) <- spline
  }

  # Aggregate metrics for the whole road from each segment. We have one metric each 10 meter slices
  # that are average to return a aggregated metrics
  metrics <- road_metrics(new_road, slice_metrics)
  metrics[["SCORE"]] <- road_score(metrics, param)
  metrics[["CLASS"]] <- get_class(metrics[["SCORE"]])

  # Merge the tables of attributes
  ngeom <- attr(new_road, "sf_column")
  original_geometry <- sf::st_geometry(centerline)
  new_geometry <- sf::st_geometry(new_road)
  attribute_table <- cbind(sf::st_drop_geometry(centerline), metrics)
  attribute_table[[ngeom]] <- new_geometry
  new_road <- sf::st_as_sf(attribute_table)
  sf::st_crs(new_road) <- sf::NA_crs_
  sf::st_crs(new_road) <- crs

  # Hidden option for JF Bourdons
  keep_class = dots$keep_class
  if (is.null(keep_class))
    keep_class = 2L
  else
    stopifnot(is.numeric(keep_class), length(keep_class) == 1)

  if (new_road$CLASS > keep_class)
  {
    sf::st_geometry(new_road) <- sf::st_geometry(centerline)
  }


  verbose("Done\n") ; cat("\n")
  new_road <- rename_sf_column(new_road, centerline)
  return(new_road)
}

#' @export
#' @rdname measure_road
measure_roads = function(ctg, roads, dtm, conductivity = NULL, water = NULL, param = alsroads_default_parameters)
{
  alert_no_index(ctg)

  i <- 1:nrow(roads)
  res <- lapply(i, function(j)
  {
    if (getOption("ALSroads.debug.verbose") | getOption("ALSroads.debug.progress")) cat("Road", j, "of", nrow(roads), " ")
    tryCatch(
    {
      measure_road(ctg, roads[j,], dtm, conductivity, water, param, Windex = FALSE)
    },
    error = function(e)
    {
      warning(paste0("Error in road ", i, ": NULL returned.\nThe error was: ", e))
      return(NULL)
    })
  })

  do.call(rbind, res)
}

alert_no_index <- function(ctg)
{
  is_copc = substr(ctg$filename, nchar(ctg$filename)-8, nchar(ctg$filename))
  if (all(is_copc == ".copc.laz"))
  {
    if (utils::packageVersion("rlas") >= "1.7.0")
      return(invisible())
    else
      message(paste0("copc files are supported using package rlas >= 1.7.0. Currently installed: ", utils::packageVersion("rlas")))
  }

  if (!lidR::is.indexed(ctg))
  {
    d <- lidR::density(ctg)
    if (d < 5)
    {
      message("No spatial index for LAS/LAZ files in this collection.")
      return(invisible())
    }
    else if (d < 10)
    {
      warning("No spatial index for LAS/LAZ files in this collection.", call. = FALSE)
      return(invisible())
    }
    else
    {
      stop("No spatial index for LAS/LAZ files in this collection.")
    }
  }
}

road_class0 <- function(centerline)
{
  new_road <- centerline
  new_road$ROADWIDTH     <- NA
  new_road$DRIVABLEWIDTH <- NA
  new_road$PERCABOVEROAD <- NA
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

  new_road <- rename_sf_column(new_road, centerline)

  return(new_road)
}

road_class4 <- function(centerline)
{
  new_road <- centerline
  new_road$ROADWIDTH     <- 0
  new_road$DRIVABLEWIDTH <- 0
  new_road$PERCABOVEROAD <- 100
  new_road$SHOULDERS     <- 0
  new_road$SINUOSITY     <- NA
  new_road$CONDUCTIVITY  <- 0
  new_road$SCORE         <- 0
  new_road$CLASS         <- 4

  # reorder the columns so outputs are consistent even if exiting early
  ngeom <- attr(new_road, "sf_column")
  names <- names(new_road)
  names <- names[names != ngeom]
  names <- append(names, ngeom)
  data.table::setcolorder(new_road, names)

  new_road <- rename_sf_column(new_road, centerline)

  return(new_road)
}

warn_weird_road <- function(centerline)
{
  angles <- st_angles(centerline)
  if (any(angles > 90))
  {
    if (any(angles[c(1, length(angles))] > 90))
      warning("Sharp turn (< 90 degrees) at one or both ends of the input road. This is weird and may lead to invalid outputs.", call. = FALSE)
    else
      warning("Sharp turn (< 90 degrees) between two consecutive segments of the input road. This is weird and may lead to invalid outputs.", call. = FALSE)
  }
}

rename_sf_column <- function(x,as)
{
  # Ensure the sf_colum is the same than the input
  current <- attr(x, "sf_column")
  name    <- attr(as, "sf_column")
  names(x)[names(x) == current] = name
  sf::st_geometry(x) = name
  x
}

