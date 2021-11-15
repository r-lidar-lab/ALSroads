grid_conductivity <- function(las, road, dtm, water = NULL)
{
  use_intensity <- "Intensity" %in% names(las@data)

  verbose("Computing conductivity maps...\n")

  nlas <- suppressMessages(lidR::normalize_height(las, dtm, na.rm = TRUE))

  # Handle bridge case
  bridge = NULL
  if (!is.null(water) && length(water) > 0)
  {
    id <- NULL
    water <- sf::st_geometry(water)
    bbox <- suppressWarnings(sf::st_bbox(las))
    bbox <- sf::st_set_crs(bbox, sf::st_crs(water))
    water <- sf::st_crop(water, bbox)
    bridge <- sf::st_intersection(sf::st_geometry(road), water)
    if (length(bridge) > 0) bridge <- sf::st_buffer(bridge, 5)
  }

  terrain <- raster::terrain(dtm, opt = c("slope","roughness"), unit = "degrees")

  dt <- system.time({
  smin <- 5
  smax <- 20
  a <- -1/(smax-smin)
  b <- -a*smax
  slope <- terrain$slope
  conductivity_slope <- slope
  conductivity_slope[] <- a*slope[] + b
  conductivity_slope[ slope < smin] <- 1
  conductivity_slope[ slope > smax] <- 0

  if (getOption("MFFProads.debug.finding"))
    raster::plot(conductivity_slope, col = viridis::viridis(25), main = "Conductivity slope")
  })
  verbose("   - Slope conductivity map in ", round(dt[3],2), "s \n", sep = "")


  dt <- system.time({
  rough <- terrain$roughness
  conductivity_rough <- rough
  conductivity_rough[ rough < 0.2 ] <- 1
  conductivity_rough[ rough >= 0.2 ] <- 1/2
  conductivity_rough[ rough > 0.3] <- 0

  if (getOption("MFFProads.debug.finding"))
    raster::plot(conductivity_rough, col = viridis::viridis(3), main = "Conductivity roughness")

  })
  verbose("   - Roughness conductivity map in ", round(dt[3],2), "s \n", sep = "")

  dt <- system.time({
  edge <- slope
  tmp <- raster::as.matrix(edge)
  tmp <- sobel(tmp)
  edge[]<- tmp
  #raster::plot(edge, col = viridis::viridis(30), main = "Sobel")

  conductivity_edge <- edge
  conductivity_edge[ edge < 15 ] <- 1
  conductivity_edge[ edge >= 15 ] <- 1/2
  conductivity_edge[ edge > 40 ] <- 0

  if (getOption("MFFProads.debug.finding"))
    raster::plot(conductivity_edge, col = viridis::viridis(3), main = "Conductivity Sobel edges")
  })
  verbose("   - Sobel conductivity map in ", round(dt[3],2), "s \n", sep = "")

  if (use_intensity)
  {
  dt <- system.time({
  Intensity <- NULL
  q <- as.integer(stats::quantile(las$Intensity, probs = 0.98))
  nlas@data[Intensity > 2L*q, Intensity := 2L*q]

  Z = nlas$Z
  nlas@data[["Z"]] <-  nlas@data[["Intensity"]] # trick to use fast C_rasterize


  if (utils::packageVersion("lidR") < "4.0.0")
  {
    lay  <- lidR:::rOverlay(nlas, dtm, buffer = 0)
    imax <- lidR:::C_rasterize(nlas, lay, FALSE, 1L)
    imin <- lidR:::C_rasterize(nlas, lay, FALSE, 2L)
  }
  else
  {
    imax <- lidR:::rasterize_fast(nlas, dtm, 0, "max")
    imin <- lidR:::rasterize_fast(nlas, dtm, 0, "min")
  }

  irange <- imax - imin
  nlas@data[["Z"]] <- Z

  if (getOption("MFFProads.debug.finding"))
    raster::plot(irange, col = viridis::inferno(10), main = "Intensity range")

  val <- irange[]
  th <- stats::quantile(val, probs = c(0.1, 0.25, 0.25), na.rm = TRUE)
  conductivity_intensity <- irange
  conductivity_intensity[ irange >= 0 ] <- 1
  conductivity_intensity[ irange > th[1] ] <- 1/2
  conductivity_intensity[ irange > th[2] ] <- 1/4
  conductivity_intensity[ irange > th[3]] <- 0

  if (getOption("MFFProads.debug.finding"))
    raster::plot(conductivity_intensity, col = viridis::viridis(3), main = "Conductivity intensity")
  })
  verbose("   - Intensity conductivity map in ", round(dt[3],2), "s \n", sep = "")
  }
  else
  {
    conductivity_intensity <- rOverlay(nlas, dtm, buffer = 0)
    conductivity_intensity[] <- 0
  }

  dt <- system.time({
  chm <- lidR::grid_canopy(nlas, dtm, lidR::p2r())
  #raster::plot(chm, col = height.colors(50))

  conductivity_chm <- chm
  conductivity_chm[ chm < 0.5 ] <- 1
  conductivity_chm[ chm >= 0.5 ] <- 1/2
  conductivity_chm[ chm > 1 ] <- 0

  if (getOption("MFFProads.debug.finding"))
    raster::plot(conductivity_chm, col = viridis::inferno(3), main = "Conductivity CHM")
  })
  verbose("   - CHM conductivity map in ", round(dt[3],2), "s \n", sep = "")


  dt <- system.time({
  tmp <- lidR::filter_poi(nlas, Z > 1, Z < 3)
  d12 <- lidR::grid_density(tmp, dtm)*(raster::res(dtm)[1]^2)
  d12[d12 <= 1] <- 0
  d12[d12 >= 1] <- 1
  d12 <- 1-d12

  if (getOption("MFFProads.debug.finding"))
    raster::plot(d12, col = viridis::inferno(3), main = "Bottom layer")
  })
  verbose("   - Bottom layer conductivity map in ", round(dt[3],2), "s \n", sep = "")


  dt <- system.time({
  gnd <- lidR::filter_ground(nlas)
  d <- lidR::grid_density(gnd, dtm)
  M <- matrix(1,3,3)
  M[2,2] <- 2
  d <- raster::focal(d,M, mean, padValue = NA, na.rm = T, pad = T)

  val <- d[]
  val <- val[val > 0]
  th <- stats::quantile(val, probs = c(0.25, 0.75, 0.95))
  conductivity_density <- d
  conductivity_density[ d == 0 ] <- 0
  conductivity_density[ d > th[1] ] <- 1/4
  conductivity_density[ d > th[2] ] <- 1/2
  conductivity_density[ d > th[3]] <- 1

  if (getOption("MFFProads.debug.finding"))
    raster::plot(conductivity_density, col = viridis::inferno(3), main = "Conductivity density")
  })
  verbose("   - Density conductivity map in ", round(dt[3],2), "s \n", sep = "")


  dt <- system.time({
  max_coductivity <- 1 * 1 * 1 * (2 * 1 + 1 + 1 + as.numeric(use_intensity))
  conductivity_all <- conductivity_slope * d12 * conductivity_edge * (2 * conductivity_density + conductivity_chm + conductivity_rough + conductivity_intensity)
  conductivity_all <- conductivity_all/max_coductivity
  })
  verbose("   - Global conductivity map in ", round(dt[3],2), "s \n", sep = "")

  if (getOption("MFFProads.debug.finding"))
    raster::plot(conductivity_all, col = viridis::inferno(15), main = "Conductivity 1m")

  s <- raster::stack(slope,
                     rough,
                     conductivity_slope,
                     conductivity_rough,
                     conductivity_edge,
                     chm,
                     conductivity_chm,
                     d,
                     conductivity_density,
                     irange,
                     conductivity_intensity,
                     d12,
                     conductivity_all)
  names(s) <- c("slope",
                "roughness",
                "conductivity_slope",
                "conductivity_roughness",
                "conductivity_edge",
                "chm",
                "conductivity_chm",
                "density",
                "conductivity_density",
                "intensity",
                "conductivity_intensity",
                "conductivity_bottom",
                "conductivity")

  if (!is.null(water) && length(water) > 0)
    s <- raster::mask(s, sf::as_Spatial(water), inverse = TRUE)

  if (length(bridge) > 0)
  {
    cells =  raster::cellFromPolygon(s$conductivity, sf::as_Spatial(bridge))
    cells = unlist(cells)
    tmp = s$conductivity
    tmp[cells] = 1
    s$conductivity = tmp

    if (getOption("MFFProads.debug.finding"))
      raster::plot(s$conductivity, col = viridis::inferno(15), main = "Conductivity 1m with bridge")
  }

  return(s)
}

mask_conductivity <- function(conductivity, road, param)
{
  verbose("Computing conductivity masks...\n")

  poly1 <- sf::st_buffer(road, param$extraction$road_buffer/2, endCapStyle = "FLAT")
  hull  <- sf::st_buffer(road, param$extraction$road_buffer/2)
  caps  <- make_caps(road, param)
  poly4 <- sf::st_buffer(road, 1)

  dt <- system.time({
  conductivity <- raster::aggregate(conductivity, fact = 2, fun = mean, na.rm = TRUE)
  conductivity <- raster::mask(conductivity, hull)

  if (getOption("MFFProads.debug.finding"))
    raster::plot(conductivity, col = viridis::inferno(15), main = "Conductivity 2m")
  })
  verbose("   - Aggregation and masking in ", round(dt[3],2), "s \n", sep = "")

  dt <- system.time({
  # Confidence factor
  f <- fasterize::fasterize(poly4, conductivity)
  f <- raster::distance(f)

  f <- raster::mask(f, hull)
  fmin <- min(f[], na.rm = T)
  fmax <- max(f[], na.rm = T)
  target_min <- 1-param[["constraint"]][["confidence"]]
  f <- (1-(((f - fmin) * (1 - target_min)) / (fmax - fmin)))

  if (getOption("MFFProads.debug.finding"))
    raster::plot(f, col = viridis::viridis(25), main = "Distance factor")
  })
  verbose("   - Road rasterization and distance factor map in ", round(dt[3],2), "s \n", sep = "")


  dt <- system.time({
    # could use raster::cellFromPolygon but is slow
    conductivity <- f*conductivity

    xy <- raster::xyFromCell(conductivity, 1: raster::ncell(conductivity))
    xy <- as.data.frame(xy)
    xy$z <- 0
    names(xy) <- c("X", "Y", "Z")
    xy <- lidR::LAS(xy, lidR::LASheader(xy), crs = sf::st_crs(caps$caps))

    if (utils::packageVersion("lidR") < "4.0.0")
    {

      res <- lidR:::C_in_polygon(xy, sf::st_as_text(sf::st_geometry(caps$caps)), 1)
      conductivity[res] <- 1
      res <- lidR:::C_in_polygon(xy, sf::st_as_text(sf::st_geometry(caps$shields)), 1)
      conductivity[res] <- 0
    }
    else
    {
      res <- !is.na(lidR:::point_in_polygons(xy, caps$caps))
      conductivity[res] <- 1
      res <- !is.na(lidR:::point_in_polygons(xy, caps$shields))
      conductivity[res] <- 0
    }

    if (getOption("MFFProads.debug.finding"))
      raster::plot(conductivity, col = viridis::inferno(15), main = "Conductivity with end caps")
  })
  verbose("   - Add full conductivity end blocks in ", round(dt[3],2), "s \n", sep = "")


  return(conductivity)
}

start_end_points = function(road, param)
{
  poly1 <- sf::st_buffer(road, param$extraction$road_buffer/2, endCapStyle = "FLAT")
  start <- lwgeom::st_startpoint(road)
  end   <- lwgeom::st_endpoint(road)
  start <- sf::st_buffer(start, param$extraction$road_buffer/2)
  end <- sf::st_buffer(end, param$extraction$road_buffer/2)
  start <- sf::st_difference(start, poly1)
  end <- sf::st_difference(end, poly1)
  A <- sf::st_centroid(start)
  B <- sf::st_centroid(end)
  A <-sf::st_coordinates(A)
  B <- sf::st_coordinates(B)
  return(list(A = A, B = B))
}

transition <- function(conductivity)
{
  verbose("Computing graph map...\n")

  dt <- system.time({
  trans <- gdistance::transition(conductivity, transitionFunction = mean, directions = 8)
  })
  verbose("   - Transition graph in ", round(dt[3],2), "s \n", sep = "")

  dt <- system.time({
  trans <- gdistance::geoCorrection(trans)
  })
  verbose("   - Geocorrection graph in ", round(dt[3],2), "s \n", sep = "")

  return(trans)
}

find_path = function(trans, road, A, B, param)
{
  caps  <- make_caps(road, param)$caps
  trans@crs <- methods::as(sf::NA_crs_, "CRS") # workaround to get rid of rgdal warning

  cost <- gdistance::costDistance(trans, A, B)

  if (is.infinite(cost))
  {
    verbose("    - Impossible to reach the end of the road\n")
    path <- sf::st_geometry(road)
    path <- sf::st_as_sf(path)
    path$CONDUCTIVITY <- 0
    return(path)
  }

  path <- gdistance::shortestPath(trans, A, B, output = "SpatialLines") |> suppressWarnings()
  path <- sf::st_as_sf(path)
  len  <- sf::st_length(path)
  path <- sf::st_simplify(path, dTolerance = 3)
  path <- sf::st_set_crs(path, sf::NA_crs_) |> sf::st_set_crs(sf::st_crs(road))
  path <- sf::st_difference(path, caps)
  #path$cost <- cost
  #path$cost_per_unit <- cost/len
  path$CONDUCTIVITY <- round(as.numeric(len/cost),2)

  if (getOption("MFFProads.debug.finding"))
    plot(sf::st_geometry(path), col = "red", add = T, lwd = 2)

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

make_caps <- function(road, param)
{
  XY <- sf::st_coordinates(road)[,1:2]
  n <- nrow(XY)
  buf <- param[["extraction"]][["road_buffer"]]/2

  angles <- st_angles(road)
  start_angle <- angles[1]
  end_angle <- angles[length(angles)]

  ii <- 2
  if (start_angle > 90) ii <- 3

  jj <- n-1
  if (end_angle > 90) jj <- n-2

  start <- sf::st_sfc(sf::st_linestring(XY[c(1,ii),]), crs = sf::st_crs(road))
  end <- sf::st_sfc(sf::st_linestring(XY[c(jj,n),]), crs = sf::st_crs(road))

  poly1 <- sf::st_geometry(sf::st_buffer(start, buf, endCapStyle = "FLAT"))
  poly2 <- sf::st_geometry(sf::st_buffer(end,   buf, endCapStyle = "FLAT"))
  poly3 <- sf::st_geometry(sf::st_buffer(road, buf))

  A <- lwgeom::st_startpoint(road)
  B <- lwgeom::st_endpoint(road)

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

  sf::st_crs(caps) <- sf::st_crs(road)
  sf::st_crs(shield) <- sf::st_crs(road)

  return(list(caps = caps, shields = shield))
}

