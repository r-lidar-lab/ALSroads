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
least_cost_path = function(las, centerline, dtm, conductivity, water, param)
{
  display <- getOption("ALSroads.debug.finding")

  if (is.null(dtm) & is.null(conductivity)) stop("'dtm' and 'conductivity' cannot be both NULL", call. = FALSE)
  if (!is.null(dtm) & !is.null(conductivity)) stop("'dtm' or 'conductivity' must be NULL", call. = FALSE)

  # Compute the limit polygon where we actually need to work
  # This allows to crop the raster (DTM or conductivty) and extract
  # only the part we are working with
  hold <- sf::st_buffer(centerline, param[["extraction"]][["road_buffer"]] + 2)

  # Clip and mask the conductivity (if provided)
  # If the conductivity has a resolution more than 2 m, aggregate to 2 m.
  # If the resolution is less than 2 meter, abort.
  if (!is.null(conductivity))
  {
    conductivity <- raster::crop(conductivity, hold)
    conductivity <- raster::mask(conductivity, hold)

    res <- round(raster::res(conductivity)[1], 2)

    if (res > 2)
      stop("The conductivity must have a resolution of 1 m or less.")

    if (res < 2)
      conductivity <- raster::aggregate(conductivity, fact = 2/res, fun = mean, na.rm = TRUE)
  }

  # If no conductivity is provided it means that the DTM is. We compute the conductivity
  if (is.null(conductivity))
  {
    verbose("Computing conductivity maps...\n")
    conductivity <- rasterize_conductivity(las, dtm = dtm, param = param, return_all = FALSE, return_stack = FALSE)
  }

  # Handle bridge case:
  # If a road crosses a water body it builds a bridge, i.e. a polygon in which we will
  # force a conductivity of 1 later
  bridge = NULL
  if (!is.null(water) && length(sf::st_geometry(water)) > 0)
  {
    id <- NULL
    water <- sf::st_geometry(water)
    bbox <- suppressWarnings(sf::st_bbox(las))
    bbox <- sf::st_set_crs(bbox, sf::st_crs(water))
    water <- sf::st_crop(water, bbox)
    bridge <- sf::st_intersection(sf::st_geometry(centerline), water)
    if (length(bridge) > 0) bridge <- sf::st_buffer(bridge, 5)
    if (length(water)  > 0) conductivity <- raster::mask(conductivity, sf::as_Spatial(water), inverse = TRUE)
  }

  if (length(bridge) > 0)
  {
    cells <- raster::cellFromPolygon(conductivity, sf::as_Spatial(bridge))
    cells <- unlist(cells)
    tmp   <- conductivity
    tmp[cells] <- 1
    conductivity <- tmp

    if (display) raster::plot(conductivity, col = viridis::inferno(15), main = "Conductivity 1m with bridge")
  }

  # Compute low resolution conductivity with mask
  conductivity <- mask_conductivity(conductivity, centerline, param)

  # Compute start and end points
  AB <- start_end_points(centerline, param)
  A  <- AB$A
  B  <- AB$B

  # Find the path
  verbose("Computing least cost path...\n")
  path <- find_path(conductivity, centerline, A, B, param)

  return(path)
}

#' Mask the conductivity and modify some pixels
#'
#' Update some pixels by masking the map with the bounding polygon, multiplying by a distance to
#' road cost factor (not described in the paper) and add terminal caps conductive pixels
#'  to allow driving from the point A and B further apart the road
#' @noRd
mask_conductivity <- function(conductivity, centerline, param)
{
  verbose("Computing conductivity masks...\n")

  # Boundary masking
  hull <- sf::st_buffer(centerline, param$extraction$road_buffer)

  # Fix an issue for road at the very edge of a catalog
  bb_hull <- sf::st_bbox(hull)
  bb_cond <- sf::st_bbox(conductivity)
  if (bb_hull[1] < bb_cond[1] | bb_hull[2] < bb_cond[2] | bb_hull[3] > bb_cond[3] | bb_hull[4] > bb_cond[4])
    conductivity = raster::extend(conductivity, bb_hull)

  conductivity <- raster::mask(conductivity, hull)

  if (getOption("ALSroads.debug.finding")) raster::plot(conductivity, col = viridis::inferno(15), main = "Conductivity 2m")
  verbose("   - Masking\n")

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

  if (getOption("ALSroads.debug.finding")) raster::plot(f, col = viridis::viridis(25), main = "Distance factor")
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

  if (getOption("ALSroads.debug.finding")) raster::plot(conductivity, col = viridis::inferno(15), main = "Conductivity with end caps")
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

find_path = function(conductivity, centerline, A, B, param)
{
  caps <- make_caps(centerline, param)$caps

  trans <- transition(conductivity)
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

  # Trim vertex inside caps
  coords <- st_coordinates(path)[,1:2]
  val <- raster::extract(conductivity, coords)
  path <- sf::st_linestring(coords[val != 1,])

  path <- sf::st_sfc(path)
  path <- sf::st_sf(path, crs = sf::st_crs(centerline))
  path <- sf::st_simplify(path, dTolerance = 3)
  path$CONDUCTIVITY <- round(as.numeric(len/cost),2)

  if (getOption("ALSroads.debug.finding")) plot(sf::st_geometry(path), col = "red", add = T, lwd = 2)

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

  shield_A <- sf::st_buffer(A, buf - 5)
  shield_B <- sf::st_buffer(B, buf - 5)
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
