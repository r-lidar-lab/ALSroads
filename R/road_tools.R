#' Extracts sequentially slices perpendicularly to the road
#'
#' Extracts sequentially slices perpendicularly to the road and makes operation on each slice
#' to return profiles along the road
#'
#' @noRd
road_measure = function(las, road, param)
{
  if (is.null(las@index[["quadtree"]]))  stop("The point cloud is not spatially indexed", call. = TRUE)

  # Retrieve start and end points of the road
  coords_road <- sf::st_coordinates(road)
  start <- coords_road[1,1:2]
  end <- coords_road[nrow(coords_road),1:2]

  # Split the path in n sections of ~ same length
  path_lenght <- as.numeric(sf::st_length(road))
  each <- param[["extraction"]][["section_length"]]/path_lenght
  at <- round(seq(0,1, each),4)
  at[length(at)] = 1
  points <- sf::st_line_sample(road, sample = at)
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

  # We can now loop through each segment between two consecutive points
  n <- nrow(CC)
  ids <- 1:(n-1)
  .segment_metrics <- vector("list", n)
  print_at = round(seq(1, n, length.out = 20)) # prints every ~5%
  for (i in ids)
  {
    if (i %in% print_at)
    {
      verbose("Computing road metrics... ", round(i/length(ids)*100,0), "%\r", sep = "")
      utils::flush.console()
    }

    # Extraction of the segment
    p1 <- CC[i,1:2]
    p2 <- CC[i+1,1:2]
    dx <- p1[1] - p2[1]
    dy <- p1[2] - p2[2]
    angle  <- atan(dy/dx)
    if (angle < 0) angle <- angle + pi
    center <- (p1+p2)/2
    height <- sqrt(sum((p1-p2)^2))
    width <- param[["extraction"]][["road_max_width"]]
    xmin <- center[1] - width/2
    xmax <- center[1] + width/2
    ymin <- center[2] - height/2
    ymax <- center[2] + height/2
    las_slice <- clip_orectangle_with_index(las, xmin, ymin, xmax, ymax, angle-pi/2)

    # Normalize the segment. In practice both normalized and raw point cloud will be used
    # but we can unnormalized on-the-fly and thus later we will pass only the normalized version
    # /!\ Normalization can fail because we never know what point cloud we get for a specific
    # segment. No ground point, no point, not enough ground points and so on. Few error
    # along the road are allowed with the try-catch block
    nlas_slice <- tryCatch(
    {
      lidR::normalize_height(las_slice, lidR::tin(), Wdegenerated = FALSE, na.rm = TRUE)
    },
    error = function(e)
    {
      if (isTRUE(getOption("MFFProads.debug")))
      {
        f <- tempfile(fileext = ".las")
        lidR::writeLAS(las_slice, f)
        message(glue::glue("Normalization impossible in segment {i}. Segment skiped with error : {e} "))
        message(glue::glue("The LAS objects that caused the failure has been saved in {f}."))
      }
      return(NULL)
    })

    if (is.null(nlas_slice)) next

    # We need to rotate the point cloud along a single axis maintaining a valid LAS object.
    # This allow to work in 2D (see figure 7)
    nlas_slice <- las_rotate(nlas_slice, angle, center)

    # We now have a 2D-like points cloud of a slice perpendicular to the road. We can perform measurements
    # such as road width, percentage of points and so on...
    # /!\ Here again there are many rare cases were it could fail such as water only + bridge,
    # missing points. Hard to find each specific case. Few error along the road are allowed
    # with the try-catch block
    m <- tryCatch(
    {
      slice_metrics(nlas_slice, param)
    },
    error = function(e)
    {
      if (isTRUE(getOption("MFFProads.debug")))
      {
        f <- tempfile(fileext = ".las")
        nlas_slice <- lidR::add_lasattribute(nlas_slice, name = "Zref", desc = "Absolute Elvation")
        lidR::writeLAS(nlas_slice, f)
        cat("\n")
        cat(glue::glue("Computation impossible in segment {i}. segment_road_metrics() failed with error : {e}\n"))
        cat(glue::glue("The LAS objects that caused the failure has been saved in {f}.\n"))
        cat("\n")
      }
      return(NULL)
    },
    warning = function(e)
    {
      if (isTRUE(getOption("MFFProads.debug")))
      {
        f <- tempfile(fileext = ".las")
        nlas_slice <- lidR::add_lasattribute(nlas_slice, name = "Zref", desc = "Absolute Elvation")
        lidR::writeLAS(nlas_slice, f)
        cat("\n")
        cat(glue::glue("Computation warning in segment {i}. segment_road_metrics() failed with error : {e}\n"))
        cat(glue::glue("The LAS objects that caused the failure has been saved in {f}.\n"))
        cat("\n")
      }
    })

    if (is.null(m)) next

    # We add some additional metrics to know the "temporal" location of this segment
    # and record the original position of the road + we updated previous locations
    m$distance_to_start <- points[["DISTANCE"]][i]
    .segment_metrics[[i]] <- m
  }

  verbose("Computing road metrics... 100%\n")

  segment_metrics <- data.table::rbindlist(.segment_metrics)
  segment_metrics <- sf::st_as_sf(segment_metrics, coords = c("xroad", "yroad"), crs = sf::st_crs(las))

  # Insert back end points because road_measure discards the ending of the road
  # This piece of code should go in road_measure I guess
  end <- segment_metrics[nrow(segment_metrics),]
  end$distance_to_start <- sf::st_length(road)
  end <- sf::st_set_geometry(end, lwgeom::st_endpoint(road))
  end <- sf::st_set_crs(end, sf::NA_crs_)
  end <- sf::st_set_crs(end, sf::st_crs(las))
  segment_metrics <- rbind(segment_metrics, end)
  return(segment_metrics)
}

#' Compute the "average" metrics for the road from all the metrics of each segment
#'
#' @param road a single spatial line in sf format
#' @param segment_metrics the metric of each segment as outputted by road_generic_loop(..., mode = "measure")
#'
#' @return An updated road with the metrics including the average state, width, drivable width, ...
#' @noRd
road_metrics =  function(road, segment_metrics)
{
  # Average percentage of point above x
  avg_percentage_above_05 <- round(mean(segment_metrics[["pzabove05"]], na.rm = TRUE), 1)
  avg_percentage_above_2 <- round(mean(segment_metrics[["pzabove2"]], na.rm = TRUE), 1)
  avg_shoulders <- round(mean(segment_metrics[["number_accotements"]], na.rm = TRUE)/2*100, 1)

  # Average width
  road_width <- round(mean(segment_metrics[["road_width"]], na.rm = TRUE), 1)

  # Average drivable width.
  drivable_width <- round(mean(segment_metrics[["drivable_width"]], na.rm = TRUE), 1)

  # Measure of the right of way not available yet
  right_of_way <- NA_real_

  # Sinuosity on XY and on Z
  z <- segment_metrics[["zroad"]]
  x <- segment_metrics[["distance_to_start"]]
  S <- round(sinuosity(road), 2)

  road_metrics = data.frame(
    ROADWIDTH = road_width,
    DRIVABLEWIDTH = drivable_width,
    #RIGHTOFWAY = right_of_way,
    PERCABOVEROAD = avg_percentage_above_05,
    SHOULDERS = avg_shoulders,
    SINUOSITY = S,
    CONDUCTIVITY = road$CONDUCTIVITY)

  if (getOption("MFFProads.debug.metrics")) plot_road_metrics(road, road_metrics, segment_metrics)

  return(road_metrics)
}

road_score = function(metrics, param)
{
  P     <- metrics$PERCABOVEROAD
  W     <- metrics$DRIVABLEWIDTH
  S     <- metrics$SHOULDERS
  Sigma <- metrics$CONDUCTIVITY

  # The overall state is the median state
  # Estimate the existence of a road in the segment based on percentage of veg points and width
  p <- param[["state"]][["percentage_veg_thresholds"]]
  road_exist <- activation(P, p, "piecewise-linear",  asc = FALSE)

  p <- param[["state"]][["drivable_width_thresholds"]]
  drivable_exist <- activation(W, p, "piecewise-linear",  asc = TRUE)

  p <- param[["state"]][["conductivity_thresholds"]]
  score_exist <- activation(Sigma, p, "piecewise-linear",  asc = TRUE)

  p <- param[["state"]][["shoulder_thresholds"]]
  shoulder_exist <- activation(S, p, "piecewise-linear",  asc = TRUE)

  pexist <- (road_exist + drivable_exist + score_exist + shoulder_exist) / 4*100

  verbose("Estimating the state of the road...\n")
  verbose("   - Estimated probability based on vegetation (P):", round(road_exist,1), "\n")
  verbose("   - Estimated probability based on road size (W):", round(drivable_exist,1), "\n")
  verbose("   - Estimated probability based on conductivity (sigma):", round(score_exist,1), "\n")
  verbose("   - Estimated probability based on road shoulders (S):", round(shoulder_exist,1), "\n")
  verbose("   - Estimated probability:", round(pexist,1), "\n")

  return(round(pexist,1))
}

get_class =  function(score)
{
  5 - as.integer(cut(score, breaks = c(-1,25,50,75,101), right = FALSE))
}

las_rotate <- function(las, angle, center)
{
  xfactor <- las@header@PHB[["X scale factor"]]
  yfactor <- las@header@PHB[["Y scale factor"]]
  rot <- matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), ncol = 2)
  coords <- as.matrix(lidR:::coordinates(las))
  coords[,1] <- coords[,1] - center[1]
  coords[,2] <- coords[,2] - center[2]
  coords <- coords %*% rot
  X <- coords[,1]
  Y <- coords[,2]
  lidR:::fast_quantization(X, xfactor, 0)
  lidR:::fast_quantization(Y, yfactor, 0)
  las@data[["X"]] <- X
  las@data[["Y"]] <- Y
  las@header@PHB[["X offset"]] <- 0
  las@header@PHB[["Y offset"]] <- 0
  las <- lidR::las_update(las)
  data.table::setattr(las, "rotation", rot)
  data.table::setattr(las, "offset", center)
  return(las)
}


