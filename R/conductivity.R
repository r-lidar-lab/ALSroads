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
rasterize_conductivity <- function(las, dtm = NULL, water = NULL, param = alsroads_default_parameters, ...)
{
  UseMethod("rasterize_conductivity", las)
}

#' @export
rasterize_conductivity.LAS <- function(las, dtm = NULL, water = NULL, param = alsroads_default_parameters, ...)
{
  use_intensity <- "Intensity" %in% names(las)
  display <- getOption("ALSroads.debug.finding")
  dots <- list(...)
  return_all <- isTRUE(dots$return_all)
  return_stack <- isTRUE(dots$return_stack)
  no_aggregate <- isTRUE(dots$no_aggregate)
  pkg <- if (is.null(dtm)) getOption("lidR.raster.default") else lidR:::raster_pkg(dtm)

  # Check if some flightline are not aligned wich each other. This is really harmful
  # for the conductivity in the CHM. It is is the case we must realign and recompute a DTM
  offsets <- flightlines_z_misalignment_matrix(las)
  dz <-  na.omit(as.numeric(offsets$offsets))
  if (any(abs(dz) > 0.07))
  {
    warning("Detection of misaligned flightlines. Automatic realignment and recomputation of the DTM will be performed", call. = FALSE)
    las <- flightlines_z_realignment(las, offsets)
    dtm <- lidR::rasterize_terrain(las, 1, lidR::tin(), pkg = "raster")
  }

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
    if (length(mask) > 0) mask <- sf::as_Spatial(mask) else mask = NULL
  }

  # plot(dtm, col = gray(1:30/30))

  # Force to use raster
  if (lidR:::raster_pkg(dtm) == "terra")
    dtm <- raster::raster(dtm)

  template2m <- raster::aggregate(dtm, fact = 2)
  template2m[] <- 0

  nlas <- lidR::normalize_height(las, dtm) |> suppressMessages() |> suppressWarnings()

  # Terrain metrics using the raster package (slope, roughness)
  terrain   <- terra::terrain(dtm, opt = c("slope","roughness"), unit = "degrees")
  slope     <- terrain$slope
  roughness <- terrain$roughness
  #plot(slope, col = gray(1:30/30))
  #plot(roughness, col = gray(1:30/30))

  # Slope-based conductivity
  s <- param$conductivity$s
  #sigma_s   <- dtm
  #sigma_s[] <- activation(slope[], s, "piecewise-linear", asc = FALSE)
  if (!is.null(mask)) slope <- raster::mask(slope, mask, inverse = T)
  sigma_s <- activation(slope, s, "piecewise-linear", asc = FALSE)
  #sigma_s = activation2(slope)
  sigma_s = raster::aggregate(sigma_s, fact = 2, fun = mean)

  if (display) raster::plot(sigma_s, col = viridis::viridis(25), main = "Conductivity slope")
  verbose("   - Slope conductivity map \n")

  # Roughness-based conductivity
  r <- param$conductivity$r
  #sigma_r   <- dtm
  #sigma_r[] <- activation(roughness[], r, "piecewise-linear", asc = FALSE)
  #sigma_r = activation2(roughness)
  #sigma_r = raster::aggregate(sigma_r, fact = 2)

  #if (display) raster::plot(sigma_r, col = viridis::viridis(25), main = "Conductivity roughness")
  #verbose("   - Roughness conductivity map\n")

  # Edge-based conductivity
  e    <- param$conductivity$e
  slop <- raster::as.matrix(slope)
  sobl <- sobel(slop)
  sigma_e <- dtm
  sigma_e[] <- sobl
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

    if (display) raster::plot(irange, col = heat.colors(20), main = "Intensity range")

    #th <- stats::quantile(irange[], probs = q, na.rm = TRUE)
    th <- c(0.1,0.3)
    #sigma_i <- template2m
    sigma_i <- activation(irange, th, "piecewise-linear", asc = FALSE)
    #sigma_i = activation2(irange)

    if (display) raster::plot(sigma_i, col = viridis::viridis(20), main = "Conductivity intensity")
    verbose("   - Intensity conductivity map\n")
  }

  # CHM-based conductivity
  h <- param$conductivity$h
  chm <- lidR::grid_canopy(nlas, template2m, lidR::p2r())
  #plot(chm, col = height.colors(25))
  #sigma_h <- chm

  if (!is.null(mask)) chm <- raster::mask(chm, mask, inverse = T)
  sigma_h <- activation(chm, h, "piecewise-linear", asc = FALSE)
  #sigma_h <- activation2(chm)

  if (display) raster::plot(sigma_h, col = viridis::inferno(25), main = "Conductivity CHM")
  verbose("   - CHM conductivity map\n")

  # Lowpoints-based conductivity
  # Check the presence/absence of lowpoints
  z1 <- 0.5
  z2 <- 3
  th <- 8

  tmp <- lidR::filter_poi(nlas, Z > z1, Z < z2)
  lp  <- lidR:::rasterize_fast(tmp, template2m, 0, "count", pkg = "raster")
  lp[is.na(lp)] <- 0
  if (!is.null(mask)) lp <- raster::mask(lp, mask, inverse = T)
  #sigma_lp <- template2m
  sigma_lp <- activation(lp, th, "thresholds", asc = FALSE)
  #sigma_lp <- activation2(lp)

  if (display) raster::plot(lp, col = viridis::inferno(20), main = "Number of low point")

  if (display) raster::plot(sigma_lp, col = viridis::inferno(20), main = "Bottom layer")
  verbose("   - Bottom layer conductivity map\n")

  # Density-based conductivity
  # Notice that the paper makes no mention of smoothing
  q <- param$conductivity$d

  gnd <- lidR::filter_ground(nlas)
  d   <- lidR::grid_density(gnd, template2m)
  if (!is.null(mask)) d <- raster::mask(d, mask, inverse = T, updatevalue = 0)
  if (display) raster::plot(d, col = viridis::inferno(10), main = "Density of ground points")

  overlap = pixel_metrics(gnd, ~length(unique(PointSourceID)), res = d, pkg = "raster")
  overlap[is.na(overlap)] = 1
  d = d/overlap
  if (display) raster::plot(d, col = viridis::inferno(10), main = "Density of ground points")

  th <- mean(d[], na.rm = T)
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
  sigma <- (sigma_lp > 0.5) * (sigma_e > 0.3)  * (sigma_s > 0.6) * (sigma_d + sigma_h + sigma_i)
  sigma <- sigma/max_coductivity
  sigma[sigma < 0.1] <- 0.1
  sigma[is.na(sigma)] <- 0.1

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
rasterize_conductivity.LAScluster = function(las, dtm = NULL, water = NULL, param = alsroads_default_parameters, ...)
{
  x <- lidR::readLAS(las)
  if (lidR::is.empty(x)) return(NULL)

  sigma <- rasterize_conductivity(x, dtm, water, param, ...)
  sigma <- lidR:::raster_crop(sigma, lidR::st_bbox(las))
  return(sigma)
}

#' @export
rasterize_conductivity.LAScatalog = function(las, dtm = NULL, water = NULL, param = alsroads_default_parameters, ...)
{
  # Enforce some options
  if (lidR::opt_select(las) == "*") lidR::opt_select(las) <- "xyzci"

  # Compute the alignment options including the case when res is a raster/stars/terra
  alignment <- lidR:::raster_alignment(2)

  if (lidR::opt_chunk_size(las) > 0 && lidR::opt_chunk_size(las) < 2*alignment$res)
    stop("The chunk size is too small. Process aborted.", call. = FALSE)

  lidR::opt_chunk_buffer(las) <- 50

  # Processing
  options <- list(need_buffer = TRUE, drop_null = TRUE, raster_alignment = alignment, automerge = TRUE)
  output  <- lidR::catalog_apply(las, rasterize_conductivity, dtm = dtm, water = water, param = param, ..., .options = options)
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
    raster::raster(ans)
  }

  ans <- vector("list", 0)
  for (i in ids) {
    psi <- lidR::filter_poi(las, PointSourceID == i)
    q = quantile(psi$Intensity, probs = 0.98)
    psi@data[Intensity>q, Intensity := q]
    ans[[as.character(i)]] <- rasterize_intensityrange(psi, res)
  }

  ans =terra::rast(ans)
  #plot(ans, col = heat.colors(50))
  ans = terra::stretch(ans, maxv = 1)
  #plot(ans, col = heat.colors(50))
  ans <- terra::tapp(ans, 1, fun = mean, na.rm = TRUE)
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

