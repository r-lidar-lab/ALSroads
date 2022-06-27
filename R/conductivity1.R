#' Compute the conductivity raster sigma
#'
#' Compute the conductivity raster sigma from lidar data according to Roussel et al. 2020
#' (see references)
#'
#' @param las an object of class LAS or LAScatalog from lidR
#' @param dtm RasterLayer. If NULL is provided a DTM is computed on the fly. But if a DTM is already
#' available it can be given to the function.
#' @param water sf POLYGON to mask water bodies
#' @param param a list of many parameters. See \link{alsroads_default_parameters}.
#' @param ... ignored
#'
#' @examples
#' library(lidR)
#' library(raster)
#' dir  <- system.file("extdata", "", package="ALSroads")
#' dtm  <- system.file("extdata", "j5gr_dtm.tif", package="ALSroads")
#' ctg  <- readLAScatalog(dir)
#' dtm  <- raster(dtm)
#' las  <- readLAS(ctg$filename[1])
#'
#' sigma <- rasterize_conductivity(las, dtm = dtm)
#' plot(sigma, col = viridis::inferno(30))
#'
#' \donttest{
#' sigma <- rasterize_conductivity(ctg, dtm = dtm)
#' plot(sigma, col = viridis::viridis(30))
#' }
#' @return a RasterLayer or SpatRaster
#' @export
rasterize_conductivity1 <- function(las, dtm = NULL, water = NULL, param = alsroads_default_parameters, ...)
{
  UseMethod("rasterize_conductivity1", las)
}

#' @export
rasterize_conductivity1.LAS <- function(las, dtm = NULL, water = NULL, param = alsroads_default_parameters, ...)
{
  use_intensity <- "Intensity" %in% names(las)
  display <- getOption("ALSroads.debug.finding")
  dots <- list(...)
  return_all <- isTRUE(dots$return_all)
  return_stack <- isTRUE(dots$return_stack)
  no_aggregate <- isTRUE(dots$no_aggregate)
  pkg <- if (is.null(dtm)) getOption("lidR.raster.default") else lidR:::raster_pkg(dtm)

  # # Check if some flightline are not aligned wich each other. This is really harmful
  # # for the conductivity in the CHM. It is is the case we must realign and recompute a DTM
  #offsets <- flightlines_z_misalignment_matrix(las)
  # invalid = offsets$count < 100
  # offsets$offsets[invalid] = NA
  # offsets$count[invalid] = 0
  # dz <-  na.omit(as.numeric(offsets$offsets))
  # if (any(abs(dz) > 0.07))
  # {
  #   warning("Detection of misaligned flightlines. Automatic realignment and recomputation of the DTM will be performed", call. = FALSE)
  #   las <- flightlines_z_realignment(las, offsets)
  #   dtm <- lidR::rasterize_terrain(las, 1, lidR::tin(), pkg = "raster")
  # }

  # Extract the DTM
  if (is.null(dtm))
  {
    dtm <- lidR::rasterize_terrain(las, 1, lidR::tin(), pkg = "raster")
  }
  else if (lidR:::is_raster(dtm))
  {
    res <- round(raster::res(dtm)[1], 2)
    if (res > 1) stop("The DTM must have a resolution of 1 m or less.")

    dtm <- raster::crop(dtm, lidR::extent(las))

    if (res < 1)
      dtm <- raster::aggregate(dtm, fact = 1/res, fun = mean)
  }
  else
  {
    stop("dtm must be a RasterLayer or must be NULL")
  }

  mask = NULL
  if (!is.null(water) && length(sf::st_geometry(water)) > 0)
  {
    id <- NULL
    water <- sf::st_geometry(water)
    bbox <- suppressWarnings(sf::st_bbox(las))
    bbox <- sf::st_set_crs(bbox, sf::st_crs(water))
    mask <- sf::st_crop(water, bbox)
    if (length(mask) > 0)
    {
      las = lidR::classify_poi(las, lidR::LASWATER, roi = water)
      mask <- sf::as_Spatial(mask)
    }
    else
    {
      mask = NULL
    }
  }

  # plot(dtm, col = gray(1:30/30))

  # Force to use raster
  if (lidR:::raster_pkg(dtm) == "terra")
    dtm <- raster::raster(dtm)

  template2m <- raster::aggregate(dtm, fact = 2)
  template2m[] <- 0

  nlas <- lidR::normalize_height(las, dtm) |> suppressMessages() |> suppressWarnings()
  nlas@data[Classification == LASWATER & Z > 2, Classification := LASBRIGDE]
  bridge = lidR::filter_poi(nlas, Classification == LASBRIGDE)
  bridge = sf::st_coordinates(bridge, z = FALSE)


  # Terrain metrics using the raster package (slope, roughness)
  slope <- terra::terrain(dtm, opt = c("slope"), unit = "degrees")
  if (!is.null(mask)) slope <- raster::mask(slope, mask, inverse = T)

  #plot(slope, col = gray(1:30/30))

  # Slope-based conductivity
  s <- param$conductivity$s
  #sigma_s   <- dtm
  #sigma_s[] <- activation(slope[], s, "piecewise-linear", asc = FALSE)
  sigma_s <- activation(slope, s, "piecewise-linear", asc = FALSE)
  #sigma_s = activation2(slope)
  sigma_s = raster::aggregate(sigma_s, fact = 2, fun = mean)

  if (display) raster::plot(sigma_s, col = viridis::viridis(25), main = "Conductivity slope")
  verbose("   - Slope conductivity map \n")

  # # Roughness-based conductivity
  # r <- param$conductivity$r
  # f = function(z)
  # {
  #   x = matrix(z,3,3)
  #
  #   d1 = abs(x[1,2] - (x[1,1] + x[1,3])/2)
  #   d2 = abs(x[3,2] - (x[3,1] + x[3,3])/2)
  #   d3 = abs(x[2,1] - (x[1,1] + x[1,3])/2)
  #   d4 = abs(x[2,3] - (x[3,1] + x[3,3])/2)
  #   d5 = abs(x[2,2] - (x[1,1] + x[3,3])/2)
  #   d6 = abs(x[2,2] - (x[1,3] + x[3,1])/2)
  #
  #   sd(c(d1,d2,d3,d4,d5,d6))
  # }
  #
  # rough = dtm |> terra::rast() |> terra::focal(matrix(1,3,3), f) |> raster::raster()
  # if (!is.null(mask)) rough <- raster::mask(rough, mask, inverse = T)
  # if(display) raster::plot(rough, col = gray(1:30/30))
  # sigma_r = activation(rough, c(0.02, 0.05), "piecewise-linear", asc = FALSE)
  # sigma_r = raster::aggregate(sigma_r, fact = 2)


  #if (display) raster::plot(sigma_r, col = viridis::viridis(25), main = "Conductivity roughness")
  #verbose("   - Roughness conductivity map\n")

  # Edge-based conductivity
  e    <- param$conductivity$e
  sigma_e <- sobel.RasterLayer(slope)
  #plot(sigma_e, col = gray(1:30/30))
  #sigma_e <- activation(sigma_e, e, "thresholds", asc = FALSE)
  sigma_e = activation2(sigma_e)
  sigma_e = raster::aggregate(sigma_e, fact = 2, fun = mean)

  if (display) raster::plot(sigma_e, col = viridis::viridis(20), main = "Conductivity Sobel edges")
  verbose("   - Sobel conductivity map\n")

  # Intensity-based conductivity
  sigma_i   <- template2m
  sigma_i[] <- 0
  if (use_intensity)
  {
    irange = intensity_range_by_flightline(las, template2m)
    irange = terra::focal(irange, matrix(1,3,3), mean, na.rm = T)
    if (!is.null(mask)) irange <- raster::mask(irange, mask, inverse = T)

    if (display) raster::plot(irange, col = heat.colors(20), main = "Intensity range")

    #th <- stats::quantile(irange[], probs = q, na.rm = TRUE)
    th <- c(0.25,0.35)
    #sigma_i <- template2m
    sigma_i <- activation(irange, th, "piecewise-linear", asc = FALSE)
    #sigma_i = activation2(irange)

    if (display) raster::plot(sigma_i, col = viridis::viridis(20), main = "Conductivity intensity")
    verbose("   - Intensity conductivity map\n")
  }

  # CHM-based conductivity
  #h <- c(1,1.5)
  #chm = grid_canopy(nlas, template2m, p2r())
  #chm <- chm_by_flightline(nlas, template2m)
  #if (display) plot(chm, col = height.colors(25))
  #sigma_h <- chm

  #if (!is.null(mask)) chm <- raster::mask(chm, mask, inverse = T)
  #sigma_h <- activation(chm, h, "piecewise-linear", asc = FALSE)
  #sigma_h <- activation2(chm)
  #h <- param$conductivity$h
  #chm = rasterize_canopy(nlas, template2m, p2r())
  #sigma_h <- activation(chm, h, "piecewise-linear", asc = FALSE)
  #plot(chm, col = height.colors(25))

  #if (display) raster::plot(sigma_h, col = viridis::inferno(25), main = "Conductivity CHM")
  #verbose("   - CHM conductivity map\n")

  # Lowpoints-based conductivity
  # Check the presence/absence of lowpoints
  #th <- 2
  #sigma_lp = low_points_by_flightline(nlas, template2m, th)
  #if (!is.null(mask)) sigma_lp <- raster::mask(sigma_lp, mask, inverse = T)

  #if (display) raster::plot(sigma_lp, col = viridis::inferno(20), main = "Number of low point")
  #verbose("   - Bottom layer conductivity map\n")

  # Density-based conductivity
  # Notice that the paper makes no mention of smoothing
  q <- param$conductivity$d

  d <- density_gnd_by_flightline(las, template2m, drop_angles = 0)
  d = terra::focal(d, matrix(1,3,3), mean, na.rm = T)
  th <- mean(d[], na.rm = T)
  d[is.na(d)] = 0

  if (!is.null(mask)) d <- raster::mask(d, mask, inverse = T, updatevalue = 0)
  if (display) raster::plot(d, col = viridis::inferno(15), main = "Density of ground points")


  th <- c(th*1.1, th*2)
  #val <- d[]
  #val <- val[val > 0]
  #th  <- stats::quantile(val, probs = q)
  #sigma_d <- template2m
  sigma_d <- activation(d, th, "piecewise-linear")
  #sigma_d = 1-activation2(d)

  if (display)  raster::plot(sigma_d, col = viridis::inferno(25), main = "Conductivity density")
  verbose("   - Density conductivity map\n")


  # Final conductivity sigma
  #alpha = param$conductivity$alpha
  #alpha$i = alpha$i * as.numeric(use_intensity)
  max_coductivity <- 2 + as.numeric(use_intensity)
  sigma <- (sigma_e > 0.2)  * (sigma_s > 0.6) * (sigma_d + sigma_i + sigma_e)
  sigma <- sigma/max_coductivity
  sigma[is.na(sigma)] <- 0.1 # lakes
  cells = raster::cellFromXY(sigma, bridge)
  sigma[cells] = 0.75

  if (display) raster::plot(sigma, col = viridis::inferno(25), main = "Conductivity 2 m", maxpixels = 1e6)
  verbose("   - Global conductivity map\n")

  # The output is sigma but we can also return everything to illustrate the paper
  if (!return_all & return_stack) {
    out <- raster::stack(sigma)
    names(out) <- "conductivity"
  } else if (!return_all & !return_stack) {
    out <- sigma
  } else {
    out <- raster::stack(slope, roughness, sigma_s, sigma_e, chm, sigma_h, d, sigma_d, irange, sigma_i, sigma_lp, sigma)
    names(out) <- c("slope", "roughness", "conductivity_slope", "conductivity_roughness", "conductivity_edge", "chm", "conductivity_chm", "density", "conductivity_density", "intensity", "conductivity_intensity", "conductivity_bottom", "conductivity")
  }

  if (display) raster::plot(out, col = viridis::inferno(15), main = "Conductivity 2m")

  if (pkg == "terra")
    out <- terra::rast(out)

  return(out)
}

#' @export
rasterize_conductivity1.LAScluster = function(las, dtm = NULL, water = NULL, param = alsroads_default_parameters, ...)
{
  x <- lidR::readLAS(las)
  if (lidR::is.empty(x)) return(NULL)

  sigma <- rasterize_conductivity1(x, dtm, water, param, ...)
  sigma <- lidR:::raster_crop(sigma, lidR::st_bbox(las))
  return(sigma)
}

#' @export
rasterize_conductivity1.LAScatalog = function(las, dtm = NULL, water = NULL, param = alsroads_default_parameters, ...)
{
  # Enforce some options
  if (lidR::opt_select(las) == "*") lidR::opt_select(las) <- "xyzciap"

  # Compute the alignment options including the case when res is a raster/stars/terra
  alignment <- lidR:::raster_alignment(2)

  if (lidR::opt_chunk_size(las) > 0 && lidR::opt_chunk_size(las) < 2*alignment$res)
    stop("The chunk size is too small. Process aborted.", call. = FALSE)

  lidR::opt_chunk_buffer(las) <- 50

  # Processing
  options <- list(need_buffer = TRUE, drop_null = TRUE, raster_alignment = alignment, automerge = TRUE)
  output  <- lidR::catalog_apply(las, rasterize_conductivity1, dtm = dtm, water = water, param = param, ..., .options = options)
  return(output)
}

intensity_range_by_flightline <- function(las, res)
{
  .N <- PointSourceID <- NULL

  res <- terra::rast(res)
  ids <- unique(las$PointSourceID)

  if (length(ids) == 1)
  {
    ans = rasterize_intensityrange(las, res)
    return(raster::raster(ans))
  }

  ans <- vector("list", 0)
  for (i in ids) {
    psi <- lidR::filter_poi(las, PointSourceID == i)
    q = quantile(psi$Intensity, probs = 0.95)
    psi@data[Intensity>q, Intensity := q]
    ans[[as.character(i)]] <- rasterize_intensityrange(psi, res)
  }

  ans =terra::rast(ans)
  #plot(ans, col = heat.colors(50))
  ans = terra::stretch(ans, maxv = 1)
  #plot(ans, col = heat.colors(50))
  ans <- terra::tapp(ans, 1, fun = max, na.rm = TRUE)
  #plot(ans, col = heat.colors(50))
  return(raster::raster(ans))
}

density_gnd_by_flightline <- function(las, res, drop_angles = 3)
{
  .N <- PointSourceID <- NULL

  res <- terra::rast(res)
  ids <- unique(las$PointSourceID)
  gnd = filter_ground(las)
  angle_max = max(abs(las$ScanAngleRank))
  gnd = filter_poi(gnd, abs(ScanAngleRank) < angle_max-drop_angles)

  if (length(ids) == 1)
  {
    ans = rasterize_density(gnd, res)
    return(raster::raster(ans))
  }

  ans <- vector("list", 0)
  for (i in ids) {
    psi <- lidR::filter_poi(gnd, PointSourceID == i)
    d <- rasterize_density(psi, res)
    d[d == 0] = NA
    q = quantile(d[], probs = 0.95, na.rm = TRUE)
    d[d>q] = q
    ans[[as.character(i)]] = d
  }

  ans = terra::rast(ans)
  #plot(ans, col = heat.colors(50))
  ans = terra::stretch(ans, maxv = 1)
  #plot(ans, col = heat.colors(50))
  ans <- terra::tapp(ans, 1, fun = max, na.rm = TRUE)
  #plot(ans, col = heat.colors(50))
  return(raster::raster(ans))
}

density_lp_by_flightline <- function(las, res)
{
  .N <- PointSourceID <- NULL

  res <- terra::rast(res)
  ids <- unique(las$PointSourceID)
  z1  <- 0.5
  z2  <- 3
  tmp <- lidR::filter_poi(las, Z > z1, Z < z2)

  if (length(ids) == 1)
  {
    ans = rasterize_density(gnd, res)
    return(raster::raster(ans))
  }

  ans <- vector("list", 0)
  for (i in ids)
  {
    psi <- lidR::filter_poi(tmp, PointSourceID == i)
    d <- rasterize_density(psi, res)
    d[d == 0] = NA
    q = quantile(d[], probs = 0.95, na.rm = TRUE)
    d[d>q] = q
    ans[[as.character(i)]] = d
  }

  ans = terra::rast(ans)
  #plot(ans, col = heat.colors(50))
  ans = terra::stretch(ans, maxv = 1)
  #plot(ans, col = heat.colors(50))
  ans <- terra::tapp(ans, 1, fun = max, na.rm = TRUE)
  #plot(ans, col = heat.colors(50))
  return(raster::raster(ans))
}

low_points_by_flightline <- function(nlas, res, th = 2)
{
  .N <- PointSourceID <- NULL

  z1 <- 0.5
  z2 <- 3
  res <- terra::rast(res)
  ids <- unique(nlas$PointSourceID)

  if (length(ids) == 1)
  {
    tmp <- lidR::filter_poi(nlas, Z > z1, Z < z2)
    ans = lidR:::rasterize_fast(tmp, res, 0, "count", pkg = "terra")
    return(raster::raster(ans > th))
  }

  ans <- vector("list", 0)
  for (i in ids) {
    psi <- lidR::filter_poi(nlas, PointSourceID == i, Z > z1, Z < z2)
    d <- lidR:::rasterize_fast(psi, res, 0, "count", pkg = "terra")
    d[is.na(d)] = 0
    ans[[as.character(i)]] = d <= 3
  }

  ans = terra::rast(ans)
  #plot(ans, col = heat.colors(50))
  ans <- terra::tapp(ans, 1, fun = min, na.rm = TRUE)
  #plot(ans, col = heat.colors(50))
  return(raster::raster(ans))
}

chm_by_flightline <- function(nlas, res)
{
  .N <- PointSourceID <- NULL

  res <- terra::rast(res)
  ids <- unique(nlas$PointSourceID)

  if (length(ids) == 1)
  {
    ans = rasterize_canopy(nlas, res, p2r())
    return(raster::raster(ans))
  }

  ans <- vector("list", 0)
  for (i in ids) {
    psi <- lidR::filter_poi(nlas, PointSourceID == i)
    d <- rasterize_canopy(psi, res, lidR::p2r())
    ans[[as.character(i)]] = d
  }

  ans = terra::rast(ans)
  ans = terra::mask(ans, terra::vect(water), inverse = T)
  #plot(ans, col = height.colors(50))
  ans = terra::stretch(ans, minv = 0)
  #plot(ans, col = height.colors(50))
  ans <- terra::tapp(ans, 1, fun = min, na.rm = TRUE)
  #plot(ans, col = height.colors(50))
  ans =  ans - min(ans[], na.rm = TRUE)
  #plot(ans, col = height.colors(50))

  return(raster::raster(ans))
}


rasterize_intensityrange <- function(las, res)
{
  # Detect outliers of intensity and change their value. This is not perfect but avoid many troubles
  Intensity <- NULL
  outliers <- as.integer(stats::quantile(las$Intensity, probs = 0.98))
  las@data[Intensity > outliers, Intensity := outliers]

  # Switch Z and Intensity trick to use fast lidR internal function
  Z <- las[["Z"]]
  las@data[["Z"]] <-  las@data[["Intensity"]]
  imax <- lidR:::rasterize_fast(las, res, 0, "max")
  imin <- lidR:::rasterize_fast(las, res, 0, "min")
  irange <- imax - imin
  return(irange)
}

