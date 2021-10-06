grid_conductivity <- function(las, road, dtm, water = NULL)
{
  cat("Computing conductivity maps...\n")

  nlas <- suppressMessages(lidR::normalize_height(las, dtm, na.rm = TRUE))

  # Handle bridge case
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
  cat("   - Slope conductivity map in ", round(dt[3],2), "s \n", sep = "")


  dt <- system.time({
  rough <- terrain$roughness
  conductivity_rough <- rough
  conductivity_rough[ rough < 0.2 ] <- 1
  conductivity_rough[ rough >= 0.2 ] <- 1/2
  conductivity_rough[ rough > 0.3] <- 0

  if (getOption("MFFProads.debug.finding"))
    raster::plot(conductivity_rough, col = viridis::viridis(3), main = "Conductivity roughness")

  })
  cat("   - Roughness conductivity map in ", round(dt[3],2), "s \n", sep = "")

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
  cat("   - Sobel conductivity map in ", round(dt[3],2), "s \n", sep = "")

  dt <- system.time({
  q <- stats::quantile(las$Intensity, probs = 0.98)
  nlas@data[Intensity > 2*q, Intensity := 2*q]

  Z = nlas$Z
  nlas@data[["Z"]] <-  nlas@data[["Intensity"]] # trick to use fast C_rasterize
  irange <- lidR:::rOverlay(nlas, dtm, buffer = 0)
  imax   <- lidR:::C_rasterize(nlas, irange, FALSE, 1L)
  imin   <- lidR:::C_rasterize(nlas, irange, FALSE, 2L)
  irange[] <- imax - imin
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
  cat("   - Intensity conductivity map in ", round(dt[3],2), "s \n", sep = "")

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
  cat("   - CHM conductivity map in ", round(dt[3],2), "s \n", sep = "")


  dt <- system.time({
  tmp <- lidR::filter_poi(nlas, Z > 1, Z < 3)
  d12 <- lidR::grid_density(tmp, dtm)*(raster::res(dtm)[1]^2)
  d12[d12 <= 1] <- 0
  d12[d12 >= 1] <- 1
  d12 <- 1-d12

  if (getOption("MFFProads.debug.finding"))
    raster::plot(d12, col = viridis::inferno(3), main = "Bottom layer")
  })
  cat("   - Bottom layer conductivity map in ", round(dt[3],2), "s \n", sep = "")


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
  cat("   - Density conductivity map in ", round(dt[3],2), "s \n", sep = "")


  dt <- system.time({
  max_coductivity <- 1 * 1 * 1 * (2 * 1 + 1 + 1 + 1)
  conductivity_all <- conductivity_slope * d12 * conductivity_edge * (2 * conductivity_density + conductivity_chm + conductivity_rough + conductivity_intensity)
  conductivity_all <- conductivity_all/max_coductivity
  })
  cat("   - Global conductivity map in ", round(dt[3],2), "s \n", sep = "")

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
  cat("Computing conductivity masks...\n")

  dt <- system.time({
  poly1 <- sf::st_buffer(road, param$extraction$road_buffer/2, endCapStyle = "FLAT")
  poly2 <- sf::st_buffer(road, param$extraction$road_buffer/2)
  sf::st_agr(poly1) <- "constant"
  sf::st_agr(poly2) <- "constant"
  poly3 <- sf::st_difference(poly2, poly1)
  poly4 <- sf::st_buffer(road, 1)

  conductivity <- raster::aggregate(conductivity, fact = 2, fun = mean, na.rm = TRUE)
  conductivity <- raster::mask(conductivity, poly2)

  if (getOption("MFFProads.debug.finding"))
    raster::plot(conductivity, col = viridis::inferno(15), main = "Conductivity 2m")
  })
  cat("   - Aggregation and masking in ", round(dt[3],2), "s \n", sep = "")

  dt <- system.time({
  # Confidence factor
  f <- fasterize::fasterize(poly4, conductivity)
  f <- raster::distance(f)

  f <- raster::mask(f, poly2)
  fmin <- min(f[], na.rm = T)
  fmax <- max(f[], na.rm = T)
  target_min <- 1-param$search$confidence
  f <- (1-(((f - fmin) * (1 - target_min)) / (fmax - fmin)))

  if (getOption("MFFProads.debug.finding"))
    raster::plot(f, col = viridis::viridis(25), main = "Distance factor")
  })
  cat("   - Road rasterization and distance factor map in ", round(dt[3],2), "s \n", sep = "")


  dt <- system.time({
    # could use raster::cellFromPolygon but is slow
    conductivity <- f*conductivity
    xy <- raster::xyFromCell(conductivity, 1: raster::ncell(conductivity))
    xy <- as.data.frame(xy)
    xy$z <- 0
    names(xy) <- c("X", "Y", "Z")
    xy <- lidR::LAS(xy, lidR::LASheader(xy))
    res <- lidR:::C_in_polygon(xy, sf::st_as_text(sf::st_geometry(poly3)), 1)
    conductivity[res] <- 1

    if (getOption("MFFProads.debug.finding"))
      raster::plot(conductivity, col = viridis::inferno(15), main = "Conductivity with end caps")
  })
  cat("   - Add full conductivity end blocks in ", round(dt[3],2), "s \n", sep = "")


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
  cat("Computing graph map...\n")

  dt <- system.time({
  trans <- gdistance::transition(conductivity, transitionFunction = mean, directions = 8)
  })
  cat("   - Transition graph in ", round(dt[3],2), "s \n", sep = "")

  dt <- system.time({
  trans <- gdistance::geoCorrection(trans)
  })
  cat("   - Geocorrection graph in ", round(dt[3],2), "s \n", sep = "")

  return(trans)
}

find_path = function(trans, road, A, B, param)
{
  poly1 <- sf::st_buffer(road, param$extraction$road_buffer/2, endCapStyle = "FLAT")
  poly2 <- sf::st_buffer(road, param$extraction$road_buffer/2)
  sf::st_agr(poly1) <- "constant"
  sf::st_agr(poly2) <- "constant"
  poly3 <- sf::st_geometry(sf::st_difference(poly2, poly1))
  trans@crs <- methods::as(sf::NA_crs_, "CRS") # workaround to get rid of rgdal warning

  cost <- gdistance::costDistance(trans, A, B)

  if (is.infinite(cost))
  {
    cat("    - Impossible to reach the end of the road\n")
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
  path <- sf::st_difference(path, poly3)
  #path$cost <- cost
  #path$cost_per_unit <- cost/len
  path$CONDUCTIVITY <- round(as.numeric(len/cost),2)

  if (getOption("MFFProads.debug.finding"))
    plot(sf::st_geometry(path), col = "red", add = T)

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
