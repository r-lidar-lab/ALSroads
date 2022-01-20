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
