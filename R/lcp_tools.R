#' Computation of the least cost path
#'
#' The function computes the conductivity layers, compute the global conductivity, mask the waterbodies,
#' add ending caps, generates location of point A and B, computes the least cost path and returns
#' the new centerline geometry.
#'
#' @param las LAS. the point cloud extracted with a buffer from the centerline
#' @param centerline sf, the road
#' @param DTM,water DTM and water body mask
#' @param list. parameters
#'
#' @noRd
least_cost_path = function(las, centerline, dtm, water, param)
{
  # Compute the limit polygon and mask the DTM
  hold <- sf::st_buffer(centerline, param[["extraction"]][["road_buffer"]])
  dtm  <- raster::crop(dtm, hold)
  dtm  <- raster::mask(dtm, hold)

  # If the DTM has a resolution more than 1 m, aggregate to 1 m. No need of a 50 cm DTM. If the resolution
  # is less than 1 meter, abort.
  res <- round(raster::res(dtm)[1], 2)
  if (res < 1)
    dtm <- raster::aggregate(dtm, fact = 1/res, fun = mean)
  else if (res > 1)
    stop("The DTM must have a resolution of 1 m or less.")

  # Compute high resolution conductivity map (1 m). The function is capable to return all intermediate
  # conductivity layer but it is for debugging and illustration purpose. We only need the final conductivity
  sigma <- grid_conductivity(las, centerline, dtm, water, param, return_all = FALSE)$conductivity

  # Compute low resolution conductivity with mask
  sigma <- mask_conductivity(sigma, centerline, param)

  # Compute start and end points
  AB <- start_end_points(centerline, param)
  A  <- AB$A
  B  <- AB$B

  # Compute the transition
  trans <- transition(sigma)

  # Find the path
  verbose("Computing least cost path...\n")
  path <- find_path(trans, centerline, A, B, param)

  return(path)
}

#' Compute the conductivity raster sigma
#'
#' Compute the conductivity rasters sigma_s, sigma_e, sigma_r and so on and the final conductivity layer
#' sigma
#'
#' @noRd
grid_conductivity <- function(las, centerline, dtm, water = NULL, param, return_all = FALSE)
{
  use_intensity <- "Intensity" %in% names(las@data)
  display <- getOption("MFFProads.debug.finding")

  verbose("Computing conductivity maps...\n")

  nlas <- suppressMessages(lidR::normalize_height(las, dtm, na.rm = TRUE))

  # Handle bridge case:
  # If a road crosses a water body it builds a bridge, i.e. a polygon in which we will
  # force a conductivity of 1 later
  bridge = NULL
  if (!is.null(water) && length(water) > 0)
  {
    id <- NULL
    water <- sf::st_geometry(water)
    bbox <- suppressWarnings(sf::st_bbox(las))
    bbox <- sf::st_set_crs(bbox, sf::st_crs(water))
    water <- sf::st_crop(water, bbox)
    bridge <- sf::st_intersection(sf::st_geometry(centerline), water)
    if (length(bridge) > 0) bridge <- sf::st_buffer(bridge, 5)
  }

  # Terrain metrics using the raster package (slope, roughness)
  terrain <- raster::terrain(dtm, opt = c("slope","roughness"), unit = "degrees")
  slope <- terrain$slope
  roughness <- terrain$roughness

  # Slope-based conductivity
  s <- param$conductivity$s
  sigma_s   <- dtm
  sigma_s[] <- activation(slope[], s, "piecewise-linear", asc = FALSE)

  if (display) raster::plot(sigma_s, col = viridis::viridis(25), main = "Conductivity slope")
  verbose("   - Slope conductivity map \n")

  # Roughness-based conductivity
  r <- param$conductivity$r
  sigma_r   <- dtm
  sigma_r[] <- activation(roughness[], r, "thresholds", asc = FALSE)

  if (display) raster::plot(sigma_r, col = viridis::viridis(3), main = "Conductivity roughness")
  verbose("   - Roughness conductivity map\n")

  # Edge-based conductivity
  e    <- param$conductivity$e
  slop <- raster::as.matrix(slope)
  sobl <- sobel(slop)
  sigma_e <- dtm
  sigma_e[] <- sobl
  sigma_e[] <- activation(sigma_e[], e, "thresholds", asc = FALSE)

  if (display) raster::plot(sigma_e, col = viridis::viridis(3), main = "Conductivity Sobel edges")
  verbose("   - Sobel conductivity map\n")

  # Intensity-based conductivity
  if (use_intensity)
  {
    q <- param$conductivity$q

    # Detect outliers of intensity and change their value. This is not perfect but avoid many troubles
    Intensity <- NULL
    outliers <- as.integer(stats::quantile(las$Intensity, probs = 0.98))
    nlas@data[Intensity > 2L*outliers, Intensity := 2L*outliers]

    # Switch Z and Intensity trick to use fast lidR internal function
    Z <- nlas$Z
    nlas@data[["Z"]] <-  nlas@data[["Intensity"]]
    imax <- lidR:::rasterize_fast(nlas, dtm, 0, "max")
    imin <- lidR:::rasterize_fast(nlas, dtm, 0, "min")
    irange <- imax - imin
    nlas@data[["Z"]] <- Z

    if (display) raster::plot(irange, col = viridis::inferno(10), main = "Intensity range")

    th <- stats::quantile(irange[], probs = q, na.rm = TRUE)
    sigma_i <- dtm
    sigma_i[] <- activation(irange[], th, "thresholds", asc = FALSE)

    if (display) raster::plot(sigma_i, col = viridis::viridis(3), main = "Conductivity intensity")
    verbose("   - Intensity conductivity map\n")
  }
  else
  {
    sigma_i   <- dtm
    sigma_i[] <- 0
  }

  # CHM-based conductivity
  h <- param$conductivity$h
  chm <- lidR::grid_canopy(nlas, dtm, lidR::p2r())
  sigma_chm <- dtm
  sigma_chm[] <- activation(chm[], h, "thresholds", asc = FALSE)

  if (display) raster::plot(sigma_chm, col = viridis::inferno(3), main = "Conductivity CHM")
  verbose("   - CHM conductivity map\n")

  # Lowpoints-based conductivity
  # Check the presence/absence of lowpoints
  z1 <- 1
  z2 <- 3
  th <- 1

  tmp <- lidR::filter_poi(nlas, Z > z1, Z < z2)
  lp  <- lidR::grid_density(tmp, dtm)*(raster::res(dtm)[1]^2)
  sigma_lp <- dtm
  sigma_lp[] <- activation(lp[], th, "thresholds", asc = FALSE)

  if (display) raster::plot(sigma_lp, col = viridis::inferno(3), main = "Bottom layer")
  verbose("   - Bottom layer conductivity map\n")

  # Density-based conductivity
  # Notice that the paper makes no mention of smoothing
  q <- param$conductivity$d

  gnd <- lidR::filter_ground(nlas)
  d   <- lidR::grid_density(gnd, dtm)
  M   <- matrix(1,3,3)
  M[2,2] <- 2
  d <- raster::focal(d, M, mean, padValue = NA, na.rm = T, pad = T)

  val <- d[]
  val <- val[val > 0]
  th  <- stats::quantile(val, probs = q)
  sigma_d <- dtm
  sigma_d[] <- activation(d[], th, "thresholds")

  if (display)  raster::plot(sigma_d, col = viridis::inferno(3), main = "Conductivity density")
  verbose("   - Density conductivity map\n")


  # Final conductivity sigma
  max_coductivity <- 1 * 1 * 1 * (2 * 1 + 1 + 1 + as.numeric(use_intensity))
  sigma <- sigma_s *sigma_lp * sigma_e * (2 * sigma_d + sigma_chm + sigma_r + sigma_i)
  sigma <- sigma/max_coductivity

  verbose("   - Global conductivity map\n")

  if (display) raster::plot(sigma, col = viridis::inferno(15), main = "Conductivity 1m")

  # The output is sigma but we can also return everything to illustrate the paper
  if (!return_all)
  {
    out <- raster::stack(sigma)
    names(out) <- "conductivity"
  }
  else
  {
    out <- raster::stack(slope, roughness, sigma_s, sigma_r, sigma_e, chm, sigma_chm, d, sigma_d, irange, sigma_i, sigma_lp, sigma)
    names(out) <- c("slope", "roughness", "conductivity_slope", "conductivity_roughness", "conductivity_edge", "chm", "conductivity_chm", "density", "conductivity_density", "intensity", "conductivity_intensity", "conductivity_bottom", "conductivity")
  }

  if (!is.null(water) && length(water) > 0)
    out <- raster::mask(out, sf::as_Spatial(water), inverse = TRUE)

  if (length(bridge) > 0)
  {
    cells <- raster::cellFromPolygon(out$conductivity, sf::as_Spatial(bridge))
    cells <- unlist(cells)
    tmp   <- out$conductivity
    tmp[cells] <- 1
    out$conductivity <- tmp

    if (display) raster::plot(out$conductivity, col = viridis::inferno(15), main = "Conductivity 1m with bridge")
  }

  return(out)
}

#' Aggregate the conductivity and modify some pixels
#'
#' Aggregate the conductivity  to a resolution of two meters and update some pixels by masking
#' the map with the bounding polygon, multiplying by a distance to road cost factor (not described
#' in the paper) and add terminal caps conductive pixels to allow driving from the point A and B
#' further apart the road
#' @noRd
mask_conductivity <- function(conductivity, centerline, param)
{
  verbose("Computing conductivity masks...\n")

  # Aggregation and boundary masking
  hull <- sf::st_buffer(centerline, param$extraction$road_buffer)
  conductivity <- raster::aggregate(conductivity, fact = 2, fun = mean, na.rm = TRUE)
  conductivity <- raster::mask(conductivity, hull)

  if (getOption("MFFProads.debug.finding")) raster::plot(conductivity, col = viridis::inferno(15), main = "Conductivity 2m")
  verbose("   - Aggregation and masking\n")

  # Penalty factor based on distance-to-road.
  # Actually small penalty since this part was not very successful. Not described in the paper.
  # Could maybe be removed. Yet it might be useful if putting more constraints
  p <- sf::st_buffer(centerline, 1)
  f <- fasterize::fasterize(p, conductivity)
  f <- raster::distance(f)
  f <- raster::mask(f, hull)
  fmin <- min(f[], na.rm = T)
  fmax <- max(f[], na.rm = T)
  target_min <- 1-param[["constraint"]][["confidence"]]
  f <- (1-(((f - fmin) * (1 - target_min)) / (fmax - fmin)))
  conductivity <- f*conductivity

  if (getOption("MFFProads.debug.finding")) raster::plot(f, col = viridis::viridis(25), main = "Distance factor")
  verbose("   - Road rasterization and distance factor map\n")

  # Set a conductivity of 1 in the caps and 0 on the outer half ring link in figure 6
  # We could use raster::cellFromPolygon but it is slow. This workaround using lidR is complex butfast.
  # Maybe using terra we could simplify the code.
  caps <- make_caps(centerline, param)
  xy <- raster::xyFromCell(conductivity, 1: raster::ncell(conductivity))
  xy <- as.data.frame(xy)
  xy$z <- 0
  names(xy) <- c("X", "Y", "Z")
  xy <- lidR::LAS(xy, lidR::LASheader(xy))
  lidR::projection(xy) <- sf::st_crs(caps$caps)

  res <- !is.na(lidR:::point_in_polygons(xy, caps$caps))
  conductivity[res] <- 1
  res <- !is.na(lidR:::point_in_polygons(xy, caps$shields))
  conductivity[res] <- 0

  if (getOption("MFFProads.debug.finding")) raster::plot(conductivity, col = viridis::inferno(15), main = "Conductivity with end caps")
  verbose("   - Add full conductivity end blocks\n")

  return(conductivity)
}

start_end_points = function(centerline, param)
{
  poly1 <- sf::st_buffer(centerline, param$extraction$road_buffer, endCapStyle = "FLAT")
  start <- lwgeom::st_startpoint(centerline)
  end   <- lwgeom::st_endpoint(centerline)
  start <- sf::st_buffer(start, param$extraction$road_buffer)
  end   <- sf::st_buffer(end, param$extraction$road_buffer)
  start <- sf::st_difference(start, poly1)
  end   <- sf::st_difference(end, poly1)
  A     <- sf::st_centroid(start)
  B     <- sf::st_centroid(end)
  A     <- sf::st_coordinates(A)
  B     <- sf::st_coordinates(B)
  return(list(A = A, B = B))
}

transition <- function(conductivity)
{
  verbose("Computing graph map...\n")

  trans <- gdistance::transition(conductivity, transitionFunction = mean, directions = 8)
  verbose("   - Transition graph\n")

  trans <- gdistance::geoCorrection(trans)
  verbose("   - Geocorrection graph\n")

  return(trans)
}

find_path = function(trans, centerline, A, B, param)
{
  caps <- make_caps(centerline, param)$caps
  trans@crs <- methods::as(sf::NA_crs_, "CRS") # workaround to get rid of rgdal warning

  cost <- gdistance::costDistance(trans, A, B)

  if (is.infinite(cost))
  {
    verbose("    - Impossible to reach the end of the road\n")
    path <- sf::st_geometry(centerline)
    path <- sf::st_as_sf(path)
    path$CONDUCTIVITY <- 0
    return(path)
  }

  path <- gdistance::shortestPath(trans, A, B, output = "SpatialLines") |> suppressWarnings()
  path <- sf::st_as_sf(path)
  len  <- sf::st_length(path)
  path <- sf::st_simplify(path, dTolerance = 3)
  path <- sf::st_set_crs(path, sf::NA_crs_)
  path <- sf::st_set_crs(path, sf::st_crs(centerline))
  path <- sf::st_difference(path, caps)
  path$CONDUCTIVITY <- round(as.numeric(len/cost),2)

  if (getOption("MFFProads.debug.finding")) plot(sf::st_geometry(path), col = "red", add = T, lwd = 2)

  return(path)
}

sobel <- function(img)
{
  # define horizontal and vertical Sobel kernel
  Shoriz <- matrix(c(1, 2, 1, 0, 0, 0, -1, -2, -1), nrow = 3)
  Svert <- t(Shoriz)
  nas <- is.na(img)
  img[nas] <- 0

  # get horizontal and vertical edges
  imgH <- EBImage::filter2(img, Shoriz)
  imgV <- EBImage::filter2(img, Svert)

  # combine edge pixel data to get overall edge data
  hdata <- EBImage::imageData(imgH)
  vdata <- EBImage::imageData(imgV)
  edata <- sqrt(hdata^2 + vdata^2)
  edata[nas] <- NA
  edata
}

make_caps <- function(centerline, param)
{
  XY <- sf::st_coordinates(centerline)[,1:2]
  n <- nrow(XY)
  buf <- param[["extraction"]][["road_buffer"]]

  angles <- st_angles(centerline)
  start_angle <- angles[1]
  end_angle <- angles[length(angles)]

  ii <- 2
  if (start_angle > 90) ii <- 3

  jj <- n-1
  if (end_angle > 90) jj <- n-2

  start <- sf::st_sfc(sf::st_linestring(XY[c(1,ii),]), crs = sf::st_crs(centerline))
  end <- sf::st_sfc(sf::st_linestring(XY[c(jj,n),]), crs = sf::st_crs(centerline))

  poly1 <- sf::st_geometry(sf::st_buffer(start, buf, endCapStyle = "FLAT"))
  poly2 <- sf::st_geometry(sf::st_buffer(end,   buf, endCapStyle = "FLAT"))
  poly3 <- sf::st_geometry(sf::st_buffer(centerline, buf))

  A <- lwgeom::st_startpoint(centerline)
  B <- lwgeom::st_endpoint(centerline)

  cap_A <- sf::st_buffer(A, buf)
  cap_B <- sf::st_buffer(B, buf)

  shield_A <- sf::st_buffer(A, buf - 4)
  shield_B <- sf::st_buffer(B, buf - 4)
  shield <- sf::st_union(shield_A, shield_B)

  caps_A <- sf::st_difference(cap_A, poly1)
  caps_A <- sf::st_cast(caps_A, "POLYGON")
  if (length(caps_A) > 1)
    caps_A <- caps_A[which.max(sf::st_area(caps_A))]

  caps_B <- sf::st_difference(cap_B, poly2)
  caps_B <- sf::st_cast(caps_B, "POLYGON")
  if (length(caps_B) > 1)
    caps_B <- caps_B[which.max(sf::st_area(caps_B))]

  caps <- c(caps_A, caps_B)
  caps <- sf::st_union(caps)
  shield <- sf::st_difference(caps, shield)

  sf::st_crs(caps) <- sf::st_crs(centerline)
  sf::st_crs(shield) <- sf::st_crs(centerline)

  return(list(caps = caps, shields = shield))
}
