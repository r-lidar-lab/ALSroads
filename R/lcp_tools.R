grid_conductivity <- function(las, road, dtm, water = NULL)
{
  cat("Computing conductivity maps... 0%\r", sep = "")
  utils::flush.console()

  nlas <- suppressMessages(lidR::normalize_height(las, dtm, na.rm = TRUE))

  # Handle bridge case
  if (!is.null(water))
  {
    id <- NULL
    water <- sf::st_geometry(water)
    if (length(water) > 0)
    {
      inter <- sf::st_intersects(water, road)
      cnt <- sum(sapply(inter, length))
      if (cnt >= 1)
      {
        Classification <- NULL
        inwater <- Z<- NULL
        las <- lidR::merge_spatial(las, sf::as_Spatial(water))
        if (any(!is.na(las$id)))
        {
          Zwater <- las$Z[!is.na(las$id)]
          breaks <- seq(min(Zwater), max(Zwater), 0.5)
          d <- findInterval(Zwater, breaks)
          d <- table(d)
          Zwater <- breaks[which.max(d)]
          las@data[!is.na(id) & Z > Zwater + 5 & Classification != lidR::LASGROUND, Classification := lidR::LASBRIGDE]

          if (sum(las$Classification == lidR::LASBRIGDE) > 0)
          {
            bridge <- lidR::grid_metrics(las, ~max(Z), dtm, filter = ~Classification == LASBRIGDE)
            Zbridge <- mean(bridge[], na.rm = TRUE)
            bridge <- raster::buffer(bridge, 10)
            nab = is.na(bridge)
            dtm[!nab] = Zbridge

            contour <- raster::rasterToPolygons(bridge, dissolve = T)
            contour <- sf::st_geometry(sf::st_as_sf(contour))
            water <- sf::st_difference(water, contour)

            message("There is a brige")
          }
        }
      }
    }
  }

  terrain <- terrain(dtm, opt = c("slope","roughness"), unit = "degrees")

  smin <- 5
  smax <- 20
  a <- -1/(smax-smin)
  b <- -a*smax
  slope <- terrain$slope
  conductivity_slope <- slope
  conductivity_slope[] <- a*slope[] + b
  conductivity_slope[ slope < smin] <- 1
  conductivity_slope[ slope > smax] <- 0

  cat("Computing conductivity maps... 15%\r", sep = "")
  utils::flush.console()

  if (getOption("MFFProads.debug.finding"))
    raster::plot(conductivity_slope, col = viridis::viridis(25), main = "Conductivity slope")

  rough <- terrain$roughness
  conductivity_rough <- rough
  conductivity_rough[ rough < 0.2 ] <- 1
  conductivity_rough[ rough >= 0.2 ] <- 1/2
  conductivity_rough[ rough > 0.3] <- 0

  if (getOption("MFFProads.debug.finding"))
    raster::plot(conductivity_rough, col = viridis::viridis(3), main = "Conductivity roughness")

  cat("Computing conductivity maps... 30%\r", sep = "")
  utils::flush.console()

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

  cat("Computing conductivity maps... 45%\r", sep = "")
  utils::flush.console()

  q <- stats::quantile(las$Intensity, probs = 0.98)
  irange <- lidR::grid_metrics(las, ~diff(range(Intensity)), dtm, filter = ~Intensity < 2*q)

  if (getOption("MFFProads.debug.finding"))
    raster::plot(irange, col = viridis::inferno(10), main = "Intensity range")

  val <- irange[]
  th <- stats::quantile(val, probs = c(0.1, 0.25, 0.25), na.rm = TRUE)
  conductivity_intensity <- irange
  conductivity_intensity[ irange >= 0 ] <- 1
  conductivity_intensity[ irange > th[1] ] <- 1/2
  conductivity_intensity[ irange > th[2] ] <- 1/4
  conductivity_intensity[ irange > th[3]] <- 0

  cat("Computing conductivity maps... 60%\r", sep = "")
  utils::flush.console()

  if (getOption("MFFProads.debug.finding"))
    raster::plot(conductivity_intensity, col = viridis::viridis(3), main = "Conductivity intensity")


  chm <- lidR::grid_canopy(nlas, dtm, lidR::p2r())
  #raster::plot(chm, col = height.colors(50))

  conductivity_chm <- chm
  conductivity_chm[ chm < 0.5 ] <- 1
  conductivity_chm[ chm >= 0.5 ] <- 1/2
  conductivity_chm[ chm > 1 ] <- 0

  cat("Computing conductivity maps... 75%\r", sep = "")
  utils::flush.console()

  if (getOption("MFFProads.debug.finding"))
    raster::plot(conductivity_chm, col = viridis::inferno(3), main = "Conductivity CHM")

  tmp <- lidR::filter_poi(nlas, Z > 1, Z < 3)
  d12 <- lidR::grid_density(tmp, dtm)*(raster::res(dtm)^2)
  d12[d12 <= 1] <- 0
  d12[d12 >= 1] <- 1
  d12 <- 1-d12
  #raster::plot(d12)

  cat("Computing conductivity maps... 90%\r", sep = "")
  utils::flush.console()

  gnd <- lidR::filter_ground(nlas)
  d <- lidR::grid_density(gnd, dtm)
  M <- matrix(1,3,3)
  M[2,2] <- 2
  d <- raster::focal(d,M, mean, padValue = NA, na.rm = T, pad = T)
  #raster::plot(d, col = viridis::inferno(25))

  val <- d[]
  val <- val[val > 0]
  th <- stats::quantile(val, probs = c(0.25, 0.75, 0.95))
  conductivity_density <- d
  conductivity_density[ d == 0 ] <- 0
  conductivity_density[ d > th[1] ] <- 1/4
  conductivity_density[ d > th[2] ] <- 1/2
  conductivity_density[ d > th[3]] <- 1

  cat("Computing conductivity maps... 100%\r", sep = "")
  utils::flush.console()

  if (getOption("MFFProads.debug.finding"))
    raster::plot(conductivity_density, col = viridis::inferno(3), main = "Conductivity density")

  max_coductivity = 1 * 1 * 1 * (2 * 1 + 1 + 1 + 1)
  conductivity_all = conductivity_slope * d12 * conductivity_edge * (2 * conductivity_density + conductivity_chm + conductivity_rough + conductivity_intensity)
  conductivity_all = conductivity_all/max_coductivity

  if (getOption("MFFProads.debug.finding"))
    raster::plot(conductivity_all, col = viridis::inferno(15), "Conductivity")

  s <- raster::stack(slope,
                     rough,
                     conductivity_slope,
                     conductivity_rough,
                     conductivity_edge,
                     conductivity_chm,
                     conductivity_density,
                     conductivity_intensity,
                     conductivity_all)
  names(s) <- c("slope",
                "roughness",
                "conductivity_slope",
                "conductivity_roughness",
                "conductivity_edge",
                "conductivity_chm",
                "conductivity_density",
                "conductivity_intensity",
                "conductivity")

  if (!is.null(water) && length(water) > 0)
    s <- raster::mask(s, sf::as_Spatial(water), inverse = TRUE)

  cat("\n")

  return(s)
}

mask_conductivity <- function(conductivity, road, param)
{
  cat("Computing conductivity masks... 0%\r", sep = "")
  utils::flush.console()

  poly1 <- sf::st_buffer(road, param$extraction$road_buffer/2, endCapStyle = "FLAT")
  poly2 <- sf::st_buffer(road, param$extraction$road_buffer/2)
  sf::st_agr(poly1) <- "constant"
  sf::st_agr(poly2) <- "constant"
  poly3 <- sf::st_difference(poly2, poly1)

  conductivity <- raster::aggregate(conductivity, fact = 2, fun = mean, na.rm = TRUE)
  conductivity <- raster::mask(conductivity, poly2)

  cat("Computing conductivity masks... 25%\r", sep = "")
  utils::flush.console()

  # Confidence factor
  f <- raster::rasterize(road, conductivity)
  f <- raster::distance(f)

  cat("Computing conductivity masks... 50%\r", sep = "")
  utils::flush.console()

  f <- raster::mask(f, poly2)
  fmin <- min(f[], na.rm = T)
  fmax <- max(f[], na.rm = T)
  target_min <- 1-param$search$confidence
  f <- (1-(((f - fmin) * (1 - target_min)) / (fmax - fmin)))

  cat("Computing conductivity masks... 75%\r", sep = "")
  utils::flush.console()

  conductivity <- f*conductivity
  cells <- raster::cellFromPolygon(conductivity, sf::as_Spatial(poly3))
  conductivity[cells[[1]]] <- 1

  cat("Computing conductivity masks... 100%\r", sep = "")
  utils::flush.console()

  cat("\n")

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

find_path = function(trans, road, A, B, param)
{
  poly1 <- sf::st_buffer(road, param$extraction$road_buffer/2, endCapStyle = "FLAT")
  poly2 <- sf::st_buffer(road, param$extraction$road_buffer/2)
  sf::st_agr(poly1) <- "constant"
  sf::st_agr(poly2) <- "constant"
  poly3 <- sf::st_geometry(sf::st_difference(poly2, poly1))
  trans@crs <- methods::as(sf::NA_crs_, "CRS") # workaround to get rid of rgdal warning
  path <- gdistance::shortestPath(trans, A, B, output = "SpatialLines")
  cost <- gdistance::costDistance(trans, A, B)
  path <- sf::st_as_sf(path)
  len  <- sf::st_length(path)
  path <- sf::st_simplify(path, dTolerance = 3)
  path <- sf::st_set_crs(path, sf::NA_crs_) |> sf::st_set_crs(sf::st_crs(poly3))
  path <- sf::st_difference(path, poly3)
  #path$cost <- cost
  #path$cost_per_unit <- cost/len
  path$SCORE <- round(as.numeric(len/cost),2)
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
