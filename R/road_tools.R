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
  avg_percentage_above_05 = round(mean(segment_metrics[["pzabove05"]], na.rm = TRUE), 1)
  avg_percentage_above_2 = round(mean(segment_metrics[["pzabove2"]], na.rm = TRUE), 1)
  avg_shoulders = round(mean(segment_metrics[["number_accotements"]], na.rm = TRUE)/2*100, 1)

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
    RIGHTOFWAY = right_of_way,
    PABOVE05 = avg_percentage_above_05,
    PABOVE2 = avg_percentage_above_2,
    SHOULDERS = avg_shoulders,
    SINUOSITY = S,
    CONDUCTIVITY = road$CONDUCTIVITY)

  if (getOption("MFFProads.debug.metrics")) plot_road_metrics(road, road_metrics, segment_metrics)

  return(road_metrics)
}

road_score = function(road, param)
{
  pzabove05 = road$PABOVE05
  drivable_width = road$DRIVABLEWIDTH
  embankement = road$SHOULDERS
  score = road$CONDUCTIVITY

  # The overall state is the median state
  # Estimate the existence of a road in the segment based on percentage of veg points and width
  p <- param[["state"]][["percentage_veg_thresholds"]]
  road_exist <- ifelse(
    pzabove05 > p[2],
    0,
    ifelse (
      pzabove05 >= p[1] & pzabove05 <= p[2],
      100-(pzabove05-p[1])/((p[2]-p[1])/100),
      100)
    )

  p <- param[["state"]][["drivable_width_thresholds"]]
  drivable_exist <- ifelse(
    drivable_width < p[1],
    0,
    ifelse(
      drivable_width >= p[1] & drivable_width <= p[2],
      (drivable_width-p[1])/(p[2]-p[1])*100,
      100)
  )

  p <- param[["state"]][["conductivity_thresholds"]]
  score_exist <- ifelse(
    score < p[1],
    0,
    ifelse(
      score >= p[1] & score <= p[2],
      (score-p[1])/(p[2]-p[1])*100,
      100)
  )

  p <- param[["state"]][["shoulder_thresholds"]]
  embamkement_exist <- ifelse(
    embankement < p[1],
    0,
    ifelse(
      embankement >= p[1] & embankement <= p[2],
      (embankement-p[1])/(p[2]-p[1])*100,
      100)
  )

  pexist <- (road_exist + drivable_exist + score_exist + embamkement_exist) / 4

  verbose("Estimating the state of the road...\n")
  verbose("   - Estimated probability based on vegetation:", round(road_exist,1), "\n")
  verbose("   - Estimated probability based on road size:", round(drivable_exist,1), "\n")
  verbose("   - Estimated probability based on conductivity:", round(score_exist,1), "\n")
  verbose("   - Estimated probability based on road shoulders:", round(embamkement_exist,1), "\n")
  verbose("   - Estimated probability:", round(pexist,1), "\n")

  return(round(pexist,1))
}

get_state =  function(score)
{
  5 - as.integer(cut(score, breaks = c(-1,25,50,75,101)))
}

road_relocate = function(las, road, dtm, water, param)
{
  poly2 <- sf::st_buffer(road, param$extraction$road_buffer/2)

  dtm <- raster::crop(dtm, poly2)
  dtm <- raster::mask(dtm, poly2)

  res <- round(raster::res(dtm)[1], 2)
  if (res < 1)
    dtm <- raster::aggregate(dtm, fact = 1/res, fun = mean)
  else if (res > 1)
    stop("The DTM must have a resolution of 1 m or less.")

  # Compute high resolution conductivity map
  conductivities <- grid_conductivity(las, road, dtm, water)
  conductivity <- conductivities$conductivity
  #plot(conductivity,  col = viridis::inferno(20), main = "High res conductivity")
  #plot(chemin, add = T, col = "yellow")

  # Compute low resolution conductivity with mask
  conductivity <- mask_conductivity(conductivity, road, param)
  #plot(conductivity,  col = viridis::inferno(20), main = "Low res conductivity")
  #plot(chemin, add = T, col = "red")

  # Compute start and end points
  AB <- start_end_points(road, param)
  A  <- AB$A
  B  <- AB$B
  #points(A, col = "red")
  #points(B, col = "red")

  # Compute the transition
  trans <- transition(conductivity)

  # Find the path
  verbose("Computing least cost path...\n")
  path <- find_path(trans, road, A, B, param)

  return(path)
}

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
  print_at = round(seq(1, 526, length.out = 10)) # prints every ~10%
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
      if (isTRUE(getOption("MFFProads.debug")))
      {
        f <- tempfile(fileext = ".las")
        lidR::writeLAS(las_segment, f)
        message(glue::glue("Normalization impossible in segment {i}. Segment skiped with error : {e} "))
        message(glue::glue("The LAS objects that caused the failure has been saved in {f}."))
      }
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
    m <- tryCatch({ segment_road_metrics(nlas_segment, param)},
      error = function(e)
      {
        if (isTRUE(getOption("MFFProads.debug")))
        {
          f <- tempfile(fileext = ".las")
          nlas_segment <- lidR::add_lasattribute(nlas_segment, name = "Zref", desc = "Absolute Elvation")
          lidR::writeLAS(nlas_segment, f)
          message(glue::glue("Computation impossible in segment {i}. segment_*_metrics() failed with error : {e} "))
          message(glue::glue("The LAS objects that caused the failure has been saved in {f}. segment_metrics()"))
        }
        return(NULL)
      })

    if (is.null(m)) next

    # We add some additional metrics to know the "temporal" location of this segment
    # and record the original position of the road + we updated previous locations
    m$distance_to_start <- points[["DISTANCE"]][i]
    .segment_metrics[[i]] <- m
  }

  verbose("Computing road metrics... 100%\n")

  return(data.table::rbindlist(.segment_metrics))
}
