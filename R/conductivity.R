#' Compute the conductivity raster sigma
#'
#' Compute the conductivity raster sigma from lidar data according to Roussel et al. 2020
#' (see references)
#'
#' @param las an object of class LAS or LAScatalog from lidR
#' @param dtm RasterLayer. If NULL is provided a DTM is computed on the fly. But if a DTM is already
#' available it can be given to the function.
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
rasterize_conductivity <- function(las, dtm = NULL, param = alsroads_default_parameters, ...)
{
  UseMethod("rasterize_conductivity", las)
}

#' @export
rasterize_conductivity.LAS <- function(las, dtm = NULL, param = alsroads_default_parameters, ...)
{
  use_intensity <- "Intensity" %in% names(las)
  display <- getOption("ALSroads.debug.finding")
  dots <- list(...)
  return_all <- isTRUE(dots$return_all)
  return_stack <- isTRUE(dots$return_stack)
  no_aggregate <- isTRUE(dots$no_aggregate)
  pkg <- if (is.null(dtm)) getOption("lidR.raster.default") else lidR:::raster_pkg(dtm)

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

  # plot(dtm, col = gray(1:30/30))

  # Force to use raster
  if (lidR:::raster_pkg(dtm) == "terra")
    dtm <- raster::raster(dtm)

  nlas <- lidR::normalize_height(las, dtm) |> suppressMessages() |> suppressWarnings()

  # Terrain metrics using the raster package (slope, roughness)
  terrain   <- terra::terrain(dtm, opt = c("slope","roughness"), unit = "degrees")
  slope     <- terrain$slope
  roughness <- terrain$roughness
  #plot(slope, col = gray(1:30/30))
  #plot(roughness, col = gray(1:30/30))

  # Slope-based conductivity
  s <- param$conductivity$s
  sigma_s   <- dtm
  sigma_s[] <- activation(slope[], s, "piecewise-linear", asc = FALSE)

  if (display) raster::plot(sigma_s, col = viridis::viridis(25), main = "Conductivity slope")
  verbose("   - Slope conductivity map \n")

  # Roughness-based conductivity
  r <- param$conductivity$r
  sigma_r   <- dtm
  sigma_r[] <- activation(roughness[], r, "piecewise-linear", asc = FALSE)

  if (display) raster::plot(sigma_r, col = viridis::viridis(25), main = "Conductivity roughness")
  verbose("   - Roughness conductivity map\n")

  # Edge-based conductivity
  e    <- param$conductivity$e
  slop <- raster::as.matrix(slope)
  sobl <- sobel(slop)
  sigma_e <- dtm
  sigma_e[] <- sobl
  #plot(sigma_e, col = gray(1:30/30))
  sigma_e[] <- activation(sigma_e[], e, "thresholds", asc = FALSE)

  if (display) raster::plot(sigma_e, col = viridis::viridis(3), main = "Conductivity Sobel edges")
  verbose("   - Sobel conductivity map\n")

  # Intensity-based conductivity
  if (use_intensity)
  {
    q <- param$conductivity$q

    # Normalize intensity o different flightlins
    if (data.table::uniqueN(nlas$PointSourceID) > 1)
    {
      u <- nlas@data[, list(N = .N, Imin = min(Intensity), Imean = as.numeric(median(Intensity)), Isd = fsd(Intensity)), by = PointSourceID]
      data = nlas@data[, c("PointSourceID", "Intensity")]
      data[, nIntensity := (Intensity - mean(Intensity))/sd(Intensity), by = PointSourceID]
      data[is.na(nIntensity), nIntensity := 0]

      i <- which.max(u$N)
      Iref = u$Imean[i]
      u$f = u$Imean/Iref
      u[f == 0, f := 1]
      idx <- match(nlas$PointSourceID, u$PointSourceID)
      Icorr = nlas$Intensity/u$f[idx]
      nlas$Intensity <- as.integer(Icorr)
    }

    # Detect outliers of intensity and change their value. This is not perfect but avoid many troubles
    Intensity <- NULL
    outliers <- as.integer(stats::quantile(nlas$Intensity, probs = 0.98))
    nlas@data[Intensity > 2L*outliers, Intensity := 2L*outliers]

    # Switch Z and Intensity trick to use fast lidR internal function
    Z <- nlas$Z
    nlas@data[["Z"]] <-  nlas@data[["Intensity"]]
    imax <- lidR:::rasterize_fast(nlas, dtm, 0, "max", pkg = "raster")
    imin <- lidR:::rasterize_fast(nlas, dtm, 0, "min", pkg = "raster")
    irange <- imax - imin
    nlas@data[["Z"]] <- Z

    if (display) raster::plot(irange, col = heat.colors(20), main = "Intensity range")

    th <- stats::quantile(irange[], probs = q, na.rm = TRUE)
    sigma_i <- dtm
    sigma_i[] <- activation(irange[], th, "piecewise-linear", asc = FALSE)

    if (display) raster::plot(sigma_i, col = viridis::viridis(20), main = "Conductivity intensity")
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
  sigma_h <- dtm
  sigma_h[] <- activation(chm[], h, "piecewise-linear", asc = FALSE)

  if (display) raster::plot(sigma_h, col = viridis::inferno(25), main = "Conductivity CHM")
  verbose("   - CHM conductivity map\n")

  # Lowpoints-based conductivity
  # Check the presence/absence of lowpoints
  z1 <- 0.5
  z2 <- 3
  th <- 2

  tmp <- lidR::filter_poi(nlas, Z > z1, Z < z2)
  lp  <- lidR:::rasterize_fast(tmp, dtm, 0, "count", pkg = "raster")
  lp[is.na(lp)] <- 0
  sigma_lp <- dtm
  sigma_lp[] <- activation(lp[], th, "thresholds", asc = FALSE)

  if (display) raster::plot(lp, col = viridis::inferno(20), main = "Number of low point")
  if (display) raster::plot(sigma_lp, col = viridis::inferno(2), main = "Bottom layer")
  verbose("   - Bottom layer conductivity map\n")

  # Density-based conductivity
  # Notice that the paper makes no mention of smoothing
  q <- param$conductivity$d

  gnd <- lidR::filter_ground(nlas)
  d   <- lidR::grid_density(gnd, dtm)
  M   <- matrix(1,3,3)
  M[2,2] <- 2
  d <- raster::focal(d, M, mean, padValue = NA, na.rm = T, pad = T)

  if (display) raster::plot(d, col = viridis::inferno(10), main = "Density of ground points")

  val <- d[]
  val <- val[val > 0]
  th  <- stats::quantile(val, probs = q)
  sigma_d <- dtm
  sigma_d[] <- activation(d[], th, "piecewise-linear")

  if (display)  raster::plot(sigma_d, col = viridis::inferno(25), main = "Conductivity density")
  verbose("   - Density conductivity map\n")

  # Final conductivity sigma
  alpha = param$conductivity$alpha
  alpha$i = alpha$i * as.numeric(use_intensity)
  max_coductivity <- sum(unlist(alpha))
  sigma <- sigma_s * sigma_lp * sigma_e * (alpha$d * sigma_d + alpha$h * sigma_h + alpha$r * sigma_r + alpha$i * sigma_i)
  sigma <- sigma/max_coductivity

  if (display) raster::plot(sigma, col = viridis::inferno(25), main = "Conductivity 1m")
  verbose("   - Global conductivity map\n")

  # The output is sigma but we can also return everything to illustrate the paper
  if (!return_all & return_stack) {
    out <- raster::stack(sigma)
    names(out) <- "conductivity"
  } else if (!return_all & !return_stack) {
    out <- sigma
  } else {
    out <- raster::stack(slope, roughness, sigma_s, sigma_r, sigma_e, chm, sigma_h, d, sigma_d, irange, sigma_i, sigma_lp, sigma)
    names(out) <- c("slope", "roughness", "conductivity_slope", "conductivity_roughness", "conductivity_edge", "chm", "conductivity_chm", "density", "conductivity_density", "intensity", "conductivity_intensity", "conductivity_bottom", "conductivity")
  }

  if (!no_aggregate)
    out <- raster::aggregate(out, fact = 2, fun = mean, na.rm = TRUE)

  if (display) raster::plot(out, col = viridis::inferno(15), main = "Conductivity 2m")

  if (pkg == "terra")
    out <- terra::rast(out)

  return(out)
}

#' @export
rasterize_conductivity.LAScluster = function(las, dtm = NULL, param = alsroads_default_parameters, ...)
{
  x <- lidR::readLAS(las)
  if (lidR::is.empty(x)) return(NULL)

  sigma <- rasterize_conductivity(x, dtm, param, ...)
  sigma <- lidR:::raster_crop(sigma, lidR::st_bbox(las))
  return(sigma)
}

#' @export
rasterize_conductivity.LAScatalog = function(las, dtm = NULL, param = alsroads_default_parameters, ...)
{
  # Enforce some options
  if (lidR::opt_select(las) == "*") lidR::opt_select(las) <- "xyzci"

  # Compute the alignment options including the case when res is a raster/stars/terra
  alignment <- lidR:::raster_alignment(1)

  if (lidR::opt_chunk_size(las) > 0 && lidR::opt_chunk_size(las) < 2*alignment$res)
    stop("The chunk size is too small. Process aborted.", call. = FALSE)

  # Processing
  options <- list(need_buffer = TRUE, drop_null = TRUE, raster_alignment = alignment, automerge = TRUE)
  output  <- lidR::catalog_apply(las, rasterize_conductivity, dtm = dtm, param = param, ..., .options = options)
  return(output)
}
