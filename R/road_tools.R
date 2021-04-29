#' Align and measure a road from Lidar data
#'
#' From a road (line), extracts the line with a buffer from the point cloud and recomputes
#' the actual positioning of the road and compute some metrics for the road such as its width
#' or drivable width as well as its state (exists, no longer exists)
#'
#' @param road a single line (sf format)
#' @param ctg a non-normalized LAScatalog object from lidR package
#' @param param a list of many parameters
#' @param water a set of polygons (sf format) of water bodies. This can help the method for road
#' very close to river or lakes
#' @param relocate bool. If the road is mis-positioned it will
#' try to retrieve its actual location before to compute some metrics
#'
#' @return a list with a several sf objects including stuff for debugging in \code{[["DEBUG"]]} +
#' stuff for end user in \code{[["OUTPUT"]]}
#' @export
#' @examples
#' library(lidR)
#' library(sf)
#'
#' dir  <- system.file("extdata", "", package="MFFProads")
#' file <- system.file("extdata", "road_971487.gpkg", package="MFFProads")
#' ctg  <- readLAScatalog(dir)
#' road <- st_read(file, quiet = TRUE)
#'
#' plot(ctg)
#' plot(st_geometry(road), add = TRUE, col = "red")
#'
#' res <- measure_road(road, ctg, mffproads_default_parameters, relocate = TRUE)
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
#' res <- measure_road(road, ctg, mffproads_default_parameters, relocate = TRUE)
#'
#' options(MFFProads.debug.finding = FALSE)
#' options(MFFProads.debug.measuring = TRUE)
#' options(MFFProads.debug.metrics = FALSE)
#' res <- measure_road(road, ctg, mffproads_default_parameters, relocate = TRUE)
#'
#' options(MFFProads.debug.finding = FALSE)
#' options(MFFProads.debug.measuring = FALSE)
#' options(MFFProads.debug.metrics = TRUE)
#' res <- measure_road(road, ctg, mffproads_default_parameters, relocate = TRUE)
#' }
#' @useDynLib MFFProads, .registration = TRUE
#' @import data.table
measure_road = function(road, ctg, param, water = NULL, relocate = FALSE)
{
  if (sf::st_geometry_type(road) != "LINESTRING") stop("Expecting LINESTRING geometry for 'road'", call. = FALSE)
  if (nrow(road) > 1) stop("Expecting a single LINESTRING", call. = FALSE)
  if (!methods::is(ctg, "LAScatalog")) stop("Expecting a LAscatalog", call. = FALSE)
  lidR:::assert_is_a_bool(relocate)
  if (!is.null(water)) { if (any(!sf::st_geometry_type(water) %in% c("MULTIPOLYGON", "POLYGON"))) stop("Expecting POLYGON geometry type for 'water'", call. = FALSE) }

  # This is the metrics we will estimate on the road. We generate a default output in case we should exit early
  SCORE <- NA
  new_road <- road
  new_road$ROADWIDTH     <- NA
  new_road$DRIVABLEWIDTH <- NA
  new_road$RIGHTOFWAY    <- NA
  new_road$SINUOSITY     <- NA
  new_road$ROADWIDTH     <- NA
  new_road$STATE         <- NA
  new_road$SCORE         <- NA

  # The roads that are too short are unlikely to be well measured. Instead of returning
  # poor results I prefer to return NA with a state 0.
  if (sf::st_length(road) < units::as_units(75, "m"))
  {
    warning("Road too short to be reliably measured", call. = FALSE)
    return(new_road)
  }

  # Query the roads in a collection of files
  las <- extract_road(ctg, road, param)

  if (lidR::is.empty(las))
  {
    warning("No point found.", call. = FALSE)
    return(new_road)
  }

  # Classify water
  if (!is.null(water))
  {
    if (nrow(water) > 0L)
    {
        Classification <- NULL
        inwater <- Z<- NULL
        water <- sf::st_buffer(water, 10)
        water <- suppressWarnings(methods::as(water, "Spatial"))
        las <- lidR::merge_spatial(las, water, "inwater")
        Zwater <- las$Z[las$inwater]
        breaks <- seq(min(Zwater), max(Zwater), 0.5)
        d <- findInterval(Zwater, breaks)
        d <- table(d)
        Zwater <- breaks[which.max(d)]
        las@data[inwater == TRUE & Z < Zwater + 1, Classification := lidR::LASWATER]
    }
    else
    {
      water = NULL
    }
  }

  # If we can't assume that the road is correctly positioned we need to recompute its
  # location accurately
  find_scores <- NULL
  if (relocate)
  {
    segment_locations_ <- road_relocate(las, road, param)
    find_scores <- segment_locations_$find_score
    SCORE <- stats::median(find_scores, na.rm = TRUE)
    new_road <- adjust_spline(segment_locations_)
  }

  # From the correctly located road, compute the metrics for each location
  if (relocate) param$extraction$road_buffer = param$extraction$road_buffer/2
  segment_metrics_ <- road_measure(las, new_road, param)

  # Aggregate metrics for the whole road from each segment
  metrics <- road_metrics(new_road, segment_metrics_, find_scores)
  metrics$SCORE <- SCORE
  metrics$STATE <- road_state(segment_metrics_, param, find_scores)

  # Merge the tables of attributes
  original_geometry <-  sf::st_geometry(road)
  new_geometry <- sf::st_geometry(new_road)
  attribute_table <- cbind(sf::st_drop_geometry(road), metrics)
  attribute_table$geom <- if (metrics$STATE > 2)  original_geometry else new_geometry
  new_road <- sf::st_as_sf(attribute_table)
  sf::st_crs(new_road) <- sf::st_crs(road)

  if (getOption("MFFProads.debug") == FALSE)
    return(new_road)
  else
    return(list(new_road = new_road,
                raw_location = segment_locations_))
}

#' Generic function to loop through all sections of a road and compute either some metrics
#' or the actual location of the road if the reference one is mispositioned
#'
#' @param road a single spatial line in sf format
#' @param las a non-normalized LAS point cloud coming with a spatial index
#' @param dtm a DTM corresponding to the LAS point cloud
#' @param param a list of many parameters
#' @param mode can be 'find' or 'measure' to work either in search mode or measurement mode
#'
#' @return a data.table with one row per segment and on column per metrics
#' @noRd
road_generic_loop = function(road, las, param, mode = 'measure')
{
  if (is.null(las@index[["quadtree"]]))  stop("The point cloud is not spatially indexed", call. = TRUE)
  mode <- match.arg(mode, choices = c("measure", "find"))

  # xc is the location of the road. In search mode it is equal to "find" because
  # we are searching it. Otherwise it is NULL because we know that the road is at 0.
  # Indeed we extracted the point cloud based on the centre line.
  xc <- if (mode == 'find') "find" else NULL

  # Retrieve start and end points of the road
  coords_road <- sf::st_coordinates(road)
  start <- coords_road[1,1:2]
  end <- coords_road[nrow(coords_road),1:2]

  # Split the path in n sections of ~ same length
  path_lenght <- as.numeric(sf::st_length(road))
  points <- sf::st_sample(road, round(path_lenght/param[["extraction"]][["section_length"]]), type =  "regular")
  points <- sf::st_cast(points, "POINT")
  points <- sf::st_as_sf(points)

  # Retrieve the "temporal" position along the line so we can order the points later
  CC  <- sf::st_coordinates(points)
  CCX <- CC[,1]
  CCY <- CC[,2]
  dx  <- c(CCX, end[1]) - c(start[1], CCX)
  dy  <- c(CCY, end[2]) - c(start[2], CCY)
  dd  <- sqrt(dx*dx+dy*dy)
  dist <- cumsum(dd[-length(dd)])
  pdist <- dist/dist[length(dist)]
  points$DISTANCE <- path_lenght*pdist

  # previous_xc records the last two positions of the road. In 'measure' mode it is NULL
  # because we now that the road it at 0 but in search mode we used the last two finding
  # to help getting consistent and robust positioning. It is initialized to Inf because
  # we don't know were is the road.
  prev_xc = if (mode == "find") c(Inf, Inf) else NULL

  # We can now loop through each segment between two consecutive points
  n <- nrow(CC)
  ids <- 1:(n-1)
  .segment_metrics <- vector("list", n)
  for (i in ids)
  {
    # update progress estimation
    if (mode == "find")
      cat("Searching the road true position ", round(i/length(ids)*100,0), "%\r", sep = "")
    else
      cat("Computing the road metrics ", round(i/length(ids)*100,0), "%\r", sep = "")

    utils::flush.console()

    # Extraction of the segment
    p1 <- CC[i,1:2]
    p2 <- CC[i+1,1:2]
    dx <- p1[1] - p2[1]
    dy <- p1[2] - p2[2]
    angle  <- atan(dy/dx)
    if (angle < 0) angle <- angle + pi
    center <- (p1+p2)/2
    height <- sqrt(sum((p1-p2)^2))
    width <- param[["extraction"]][["road_buffer"]]
    xmin <- center[1] - width/2
    xmax <- center[1] + width/2
    ymin <- center[2] - height/2
    ymax <- center[2] + height/2
    las_segment <- clip_orectangle_with_index(las, xmin, ymin, xmax, ymax, angle-pi/2)

    # Normalize the segment. In practice both normalized and raw point cloud will be used
    # but we can unnormalized on-the-fly and thus later we will pass only the normalized
    # /!\ Normalization can fail because we never know what point cloud we get for a specific
    # segment. No ground point, no point, not enough ground points and so on. Few error
    # along the road are allowed with the try-catch block
    nlas_segment <- tryCatch(
    {
       lidR::normalize_height(las_segment, lidR::tin(), Wdegenerated = FALSE, na.rm = TRUE)
    },
    error = function(e)
    {
      f <- tempfile(fileext = ".las")
      lidR::writeLAS(las_segment, f)
      message(glue::glue("Normalization impossible in segment {i}. Segment skiped with error : {e} "))
      message(glue::glue("The LAS objects that caused the failure has been saved in {f}."))
      return(NULL)
    })

    if (is.null(nlas_segment)) next

    # We need to put put the point cloud along a single axis maintaining a valid LAS
    # object. Code adapted from lidR::clip_transect
    rot <- matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), ncol = 2)
    coords <- as.matrix(lidR:::coordinates(nlas_segment))
    coords[,1] <- coords[,1] - center[1]
    coords[,2] <- coords[,2] - center[2]
    coords <- coords %*% rot
    X <- coords[,1]
    Y <- coords[,2]
    lidR:::fast_quantization(X, las@header@PHB[["X scale factor"]], 0)
    lidR:::fast_quantization(Y, las@header@PHB[["Y scale factor"]], 0)
    nlas_segment@data[["X"]] <- X
    nlas_segment@data[["Y"]] <- Y
    nlas_segment@header@PHB[["X offset"]] <- 0
    nlas_segment@header@PHB[["Y offset"]] <- 0
    nlas_segment <- lidR::las_update(nlas_segment)
    data.table::setattr(nlas_segment, "rotation", rot)
    data.table::setattr(nlas_segment, "offset", center)

    # Computation of the metrics. The output is either the location of the road segment +
    # a location score if we are in search mode or the metrics (state, width, ...) of the
    # road segment.
    # /!\ Here again there are many rare cases were it could fail such as water only + bridge,
    # missing points. Hard to find each specific case. Few error along the road are allowed
    # with the try-catch block
    m <- tryCatch(
    {
      if (mode == "find")
        m <- segment_location_metrics(nlas_segment, param, prev_xc)
      else
        m <- segment_road_metrics(nlas_segment, param)

      m
    },
    error = function(e)
    {
      f <- tempfile(fileext = ".las")
      nlas_segment <- lidR::add_lasattribute(nlas_segment, name = "Zref", desc = "Absolute Elvation")
      lidR::writeLAS(nlas_segment, f)
      message(glue::glue("Computation impossible in segment {i}. segment_*_metrics() failed with error : {e} "))
      message(glue::glue("The LAS objects that caused the failure has been saved in {f}. segment_metrics() called with prev_xc = {cat(prev_xc)}"))
      return(NULL)
    })

    if (is.null(m)) next

    # We add some additional metrics to know the "temporal" location of this segment
    # and record the original position of the road + we updated privious locations
    m$distance_to_start <- points[["DISTANCE"]][i]
    if (mode == "find")
    {
      m$xref <- center[1]
      m$yref <- center[2]
      prev_xc[2] <- prev_xc[1]
      prev_xc[1] <- m$xc
    }

    .segment_metrics[[i]] <- m
  }

  cat("\n")

  return(data.table::rbindlist(.segment_metrics))
}

#' Compute the "average" metrics for the road from all the metrics of each segment
#'
#' @param road a single spatial line in sf format
#' @param segment_metrics the metric of each segment as outputted by road_generic_loop(..., mode = "measure")
#'
#' @return An updated road with the metrics including the average state, width, drivable width, ...
#' @noRd
road_metrics =  function(road, segment_metrics, find_scores = NULL)
{
  # Average percentage of point above x
  avg_percentage_above_05 = round(mean(segment_metrics[["pzabove05"]]), 1)
  avg_percentage_above_2 = round(mean(segment_metrics[["pzabove2"]]), 1)

  # Average width
  road_width <- round(mean(segment_metrics[["road_width"]]),1)

  # Average drivable width.
  drivable_width <- round(mean(segment_metrics[["drivable_width"]]), 1)

  # Measure of the right of way not available yet
  right_of_way <- NA_real_

  # Sinuosity on XY and on Z
  z <- segment_metrics[["zroad"]]
  x <- segment_metrics[["distance_to_start"]]
  S <- round(sinuosity(road), 2)

  road_metrics = data.frame(
    ROADWIDTH = road_width,
    DRIVABLEWIDTH = drivable_width,
    RIGHTOFWAY = right_of_way,
    PABOVE05 = avg_percentage_above_05,
    PABOVE2 = avg_percentage_above_2,
    SINUOSITY = S)

  if (getOption("MFFProads.debug.metrics")) plot_road_metrics(road_metrics, segment_metrics, find_scores)

  return(road_metrics)
}

road_state = function(segment_metrics, param, find_scores = NULL)
{
  n = 2
  pzabove05 = segment_metrics$pzabove05
  drivable_width = segment_metrics$drivable_width

  # The overall state is the median state
  # Estimate the existence of a road in the segment based on percentage of veg points and width
  p <- param[["state"]][["percentage_veg_thresholds"]]
  road_exist <- ifelse(
    pzabove05 > p[3],
    0,
    ifelse (
      pzabove05 >= p[1] & pzabove05 <= p[3],
      100-(pzabove05-p[1])/((p[3]-p[1])/100),
      100)
    )

  p <- param[["state"]][["drivable_width_thresholds"]]
  drivable_exist <- ifelse(
    drivable_width < p[1],
    0,
    ifelse(
      drivable_width >= p[1] & drivable_width <= p[2],
      (drivable_width-1)/(p[2]-p[1])*100,
      100)
  )

  score_exist = 0
  if (!is.null(find_scores))
  {
    n <- 3
    p <- param[["state"]][["score_thresholds"]]
    score_exist <- ifelse(
      find_scores < p[1],
      0,
      ifelse(
        find_scores >= p[1] & find_scores <= p[2],
        (find_scores-1)/(p[2]-p[1])*100,
        100)
    )
  }

  pexist <- (mean(road_exist, na.rm = TRUE) + mean(drivable_exist, na.rm = TRUE) + mean(score_exist, na.rm = TRUE))/n
  state <- 5 - as.integer(cut(pexist, breaks = c(-1,20,40,70,101)))
  state <- round(stats::median(state, 0))

  return(state)
}

road_relocate = function(las, road, param)
{
  distance_to_start <- NULL
  segment_locations <- road_generic_loop(road, las, param, mode = 'find')
  data.table::setorder(segment_locations, distance_to_start)
  segment_locations <- sf::st_as_sf(segment_locations, coords = c("xroad", "yroad"))
}

road_measure = function(las, road, param)
{
  road_generic_loop(road, las, param, mode = 'measure')
}
# get_zline = function(las, roads)
# {
#   if (is(las, "LAS"))
#   {
#     x <- sf::st_cast(roads, "POINT")
#     dtm = grid_metrics(lasfilterground(las), ~list(Z = mean(Z)), 5)
#     sp <- as(x, "Spatial")
#     return(list(raster::extract(dtm, sp)+1))
#   }
#
#   f = function(i)
#   {
#     x <- sf::st_cast(roads[i,], "POINT")
#     dtm = grid_metrics(lasfilterground(las[[i]]), ~list(Z = mean(Z)), 5)
#     sp <- as(x, "Spatial")
#     raster::extract(dtm, sp)+1
#   }
#
#   return(lapply(1:length(las), f))
# }


# ~~~~ THIS IS NEEDED FOR SPLITTING THE LINES (DISABLED FOR THE MOMENT)
# Computing stable moving mode for the road state
#segment_metrics_[['stable_state']] <- stable_moving_mode(as.integer(ma(segment_metrics_[["state"]])))

# ~~~~ THIS IS NEEDED FOR SPLITTING THE LINES (DISABLED FOR THE MOMENT)
# if (relocate)
# {
#   # At this stage we have a road shape but some section might not have been found because the
#   # road does not actually exists as a whole or no longer exist locally. If the road has been
#   # relocated but some section are not found then the road shape is likely a non sense.
#   # The solution is to restore the original road so at least the output is not worst
#   # than the original
#
#   # If some section of the road cannot be found (states 3+) restore the original
#   # road for these sections. It will be better than nothing...
#   v <- rle(segment_metrics_[["stable_state"]])
#   n <- length(v$lengths)
#   p <- segment_metrics_[["distance_to_start"]][cumsum(v$length)]
#   p <- c(0, p)
#
#   subroads = vector("list", n)
#   for (i in 1:n)
#   {
#     start = p[i]
#     end = p[i+1]
#     sub = spline_as_points[spline_as_points$distance_to_start >= start & spline_as_points$distance_to_start <= end,]
#
#     if (v$values[i] > 2)
#     {
#       subroads[[i]] = sf::st_as_sf(sf::st_sfc(sf::st_linestring(cbind(sub$xref, sub$yref))))
#     }
#     else
#     {
#       subroads[[i]] = sf::st_as_sf(sf::st_sfc(sf::st_linestring(sf::st_coordinates(sub))))
#     }
#   }
#
#   tmp = do.call(rbind, subroads)
#   road_found_as_line = sf::st_as_sf(sf::st_sfc(sf::st_linestring(sf::st_coordinates(tmp)[,-3])))
# }

# Slit the path where the state change
#road_found_splitted <- cut_line(road_found_as_line, at = segment_metrics_[["stable_state"]], metrics = segment_metrics_)

# ~~~~ THIS IS NEEDED FOR SPLITTING THE LINES (DISABLED FOR THE MOMENT)
# For further processing we need a 'temporal' dimension for this road
#CC  <- sf::st_coordinates(spline_as_points)
#CCX <- as.numeric(CC[,1])
#CCY <- as.numeric(CC[,2])
#dx  <- CCX[-1] - CCX[-length(CCX)]
#dy  <- CCY[-1] - CCY[-length(CCY)]
#dd  <- sqrt(dx*dx+dy*dy)
#spline_as_points[["distance_to_start"]] = c(0,cumsum(dd))
