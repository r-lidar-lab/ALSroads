#' Drive along an unknown road
#'
#' Drive along an unknown road in a conductivity raster starting only with a small piece
#' of road segment pointing in the right direction.
#'
#' @param starting_road  line (\code{sf} format)
#' @param conductivity  raster (\code{raster} format)
#' @param fov  numeric. Field of view (degrees) ahead of the search vector.
#' @param sightline  numeric (distance unit). Search distance used to find the next most probable
#' point on the road.
#'
#' @return list, \code{line} being the road found (as \code{sfc}) and \code{cost} a numeric vector
#' of cost at each invidual vertices of the line.
#' @export
#' @examples
#' library(raster)
#' library(sf)
#'
#' conductivity <- system.file("extdata", "j53e_conductivity.tif", package = "MFFProads")
#' starting <- system.file("extdata", "j53e_drived_road.gpkg", package = "MFFProads")
#'
#' conductivity <- raster(conductivity)
#' starting <- st_read(starting, "starting_road", quiet = TRUE)
#'
#' res <- drive_road(starting, conductivity)
#'
#' plot(conductivity, col = viridis::viridis(50))
#' plot(starting, add = TRUE, col = "green", lwd = 3)
#' plot(res, add = TRUE, col = "red", lwd = 2)
drive_road <- function(starting_road, conductivity, fov = 45, sightline = 10)
{
  # The idea is to only extract values that are mostly ahead of the line

  t0 <- Sys.time()

  # Threshold of cost to stop the search.
  cost_max <- sightline * 6

  # Constant definition
  buf_dist <- c("ahead" = 900, "behind" = 100, "side" = 400)  # Could be set as a parameter

  # We need an orthonormal raster
  resolution <- raster::res(conductivity)
  if (resolution[1] != resolution[2]) stop("'conductivity' raster must have the same resolution in both X and Y axis.", call. = FALSE)
  resolution <- resolution[1]

  # Initialization of list of coordinates with starting line
  start_pts <- sf::st_geometry(sf::st_cast(starting_road, "POINT"))
  n <- length(start_pts)
  start <- start_pts[n-1][[1]]
  end   <- start_pts[n][[1]]
  list_coords <- list(start, end)
  heading <- get_heading(start, end)

  sub_aoi_conductivity <- query_conductivity_aoi(conductivity, starting_road)
  resolution <- raster::res(sub_aoi_conductivity)[1]
  #plot(conductivity, col = viridis::viridis(50))
  #plot(raster::extent(sub_aoi_conductivity), col = "red", add = T)
  #plot(starting_road, add = T, lwd = 3, col = "red")

  # Init view angles as a function of the resolution of the raster and the sightline
  angular_resolution <- floor((resolution / sightline) * (180 / pi))
  angles <- seq(angular_resolution, fov, angular_resolution)
  angles <- c(rev(-angles), 0, angles)
  angles_rad <- angles * pi / 180

  # Set penalty coefficient for each view angle
  penalty_at_45_degrees <- 1.5  # Could be set as a parameter
  penalty <- abs(angles) / fov
  penalty <- penalty * (penalty_at_45_degrees - 1) + 1

  # Loop initialisation
  current_cost <- 0
  k <- 2
  cat("Driving the conductivity raster\n")
  while (current_cost <= cost_max)
  {
    if (k %% 2 == 0) cat("", (k-1)*10, " m\r", sep = "")
    flush.console()

    # Compute heading from the two previous points
    p1 <- list_coords[[k-1]]
    p2 <- list_coords[[k]]
    heading <- get_heading(p1, p2)

    # Create all possible ends ahead of the previous point
    X <- p2[1] + sightline * cos(heading + angles_rad)
    Y <- p2[2] + sightline * sin(heading + angles_rad)
    M <- data.frame(X,Y)

    # Generate a very small area of interest to analyse
    start <- sf::st_sfc(p2)
    ends <- sf::st_as_sf(M, coords = c("X", "Y"))
    search_zone <- sf::st_bbox(c(start, sf::st_geometry(ends))) + c(-1,-1,1,1)*resolution
    aoi <- raster::crop(sub_aoi_conductivity, search_zone)
    #plot(aoi, col = viridis::viridis(50))
    #plot(start, add = T, pch = 19, col = "red")
    #plot(ends, pch = 19, add = T, col = "red")

    # Compute the transition matrix for this AOI
    trans <- transition(aoi)

    # Check if some ends fall outside of the AOI conductivity raster
    # If any reload another raster part further ahead.
    val <- raster::extract(aoi, M)
    if (any(is.na(val)))
    {
      #cat("\nReaching bounds of AOI. Loading next AOI\n")

      # Make a newly sub aoi conductivity raster further ahead
      line <- sf::st_cast(c(p1, p2), "LINESTRING")
      new_aoi_poly <- line |>
        st_extend_line(c(buf_dist["ahead"], buf_dist["behind"])) |>
        sf::st_buffer(buf_dist["side"], endCapStyle = "FLAT") |>
        sf::st_sfc() |>
        sf::st_set_crs(sf::st_crs(starting_road))

      sub_aoi_conductivity <- query_conductivity_aoi(conductivity, new_aoi_poly)
      #plot(conductivity, col = viridis::viridis(50))
      #plot(raster::extent(sub_aoi_conductivity), col = "red", add = T)
      #plot(starting_road, add = T, lwd = 3, col = "red")

      # Check again if some ends fall outside of the newly cropped conductivity raster
      val <- raster::extract(sub_aoi_conductivity, M)
      if (any(is.na(val))) {
        warning("Drive stopped early. Edge of conductivity raster has been reached.", call. = FALSE)
        cost_max <- -Inf
      }

      # Re-extract the AOI and recompute the transition
      aoi <- raster::crop(sub_aoi_conductivity, search_zone)
      trans <- transition(aoi)
    }

    # Find cost of each edge point
    if (!is.infinite(cost_max))
    {
      cost <- gdistance::costDistance(trans, sf::as_Spatial(start), sf::as_Spatial(ends)) |> suppressWarnings()
      cost <- as.numeric(cost)
      cost[is.infinite(cost)] <- 2 * max(cost[!is.infinite(cost)])
      cost2 <- ma(cost * penalty)
      ends$cost <- cost2
      #plot(raster::crop(conductivity, raster::extent(ends) + 100), col = viridis::viridis(50))
      #plot(ends, pal = viridis::viridis, pch = 19, add = T)
      #plot(p1, col = "red", pch = 19, add = T)
      #plot(p2, col = "red", pch = 19, add = T)

      # Select point coordinates with the lowest
      # penalty adjusted cost
      idx_min <- which.min(cost2)
      current_cost <- cost[idx_min]
      W <- M[idx_min,] |> as.matrix() |> sf::st_point()

      k <- k + 1
      list_coords[[k]] <- W
    }
  }

  newline <- sf::st_sfc(list_coords) |>
    sf::st_coordinates() |>
    sf::st_linestring() |>
    sf::st_sfc() |>
    sf::st_set_crs(sf::st_crs(starting_road)) |>
    sf::st_simplify(dTolerance = 2)

  tf <- Sys.time()
  dt <- tf-t0
  dt <- round(units::as_units(dt),1)
  dist <- sf::st_length(newline)
  hdt = dt
  units(hdt) <- units::make_units(h)
  units(dist) <- units::make_units(km)
  dist = round(dist, 2)
  speed = round(dist/hdt)
  cat("Processed ended in", dt, units::deparse_unit(dt),
      ": road of", dist, units::deparse_unit(dist),
      "driven at", speed, units::deparse_unit(speed),
      "\n")

  return(newline)
}

#' Find potential secondary roads from main road
#'
#' Tries to find secondary roads (branches) joining the main road using
#' a conductivity raster and the line geometry of the main road.
#'
#' @param road  line (\code{sf} format). Known road.
#' @param conductivity  raster (\code{raster} format). Conductivity raster covering the road.
#'
#' @return lines (\code{sf} format) representing potential starting points of branching roads.
#' @export
#' @examples
#' library(sf)
#' library(terra)
#'
#' road <- system.file("extdata", "j53e_drived_road.gpkg", package="MFFProads")
#' conductivity <- system.file("extdata", "j53e_conductivity.tif", package="MFFProads")
#'
#' road <- st_read(road, "drived_road", quiet = TRUE)
#' conductivity <- rast(conductivity)
#'
#' road_branches <- find_road_branches(road, conductivity)
#'
#' plot(conductivity, col = viridis::viridis(50))
#' plot(road, col = "green", lwd = 3, add = TRUE)
#' plot(road_branches, col = "red", add = TRUE, lwd = 3)
find_road_branches <- function(road, conductivity)
{
  if (is(conductivity, "RasterLayer"))
    conductivity <- terra::rast(conductivity)

  resolution <- terra::res(conductivity)
  if (resolution[1] != resolution[2]) stop("'conductivity' raster must have the same resolution in both X and Y axis.", call. = FALSE)
  resolution <- resolution[1]

  resolution_min <- 2  # Could be set as a parameter
  if (resolution < resolution_min)
  {
    agg_factor <- ceiling(resolution_min / resolution)
    resolution <- resolution * agg_factor
    cat("Step 0: downsample conductivity to", resolution, "m\n")
    conductivity_crop <- terra::aggregate(conductivity, fact = agg_factor, fun = mean, na.rm = TRUE)
  }

  # 3 levels of buffer around the road
  d1 <-3
  d2 <- 20
  d3 <- 50
  buffer0 <- sf::st_buffer(road, d1)
  buffer1 <- sf::st_buffer(road, d2)
  buffer2 <- sf::st_buffer(road, d3)
  search_zone <- sf::st_difference(buffer2, buffer1)

  # Compute the average conductivity of the road
  m_conductivity <- terra::extract(conductivity, terra::vect(buffer0))[[2]]
  m_conductivity <- quantile(m_conductivity, na.rm = T, probs = 0.33)

  # Ensure continuity of the road
  conductivity <- terra::crop(conductivity, terra::vect(buffer2))
  conductivity <- terra::mask(conductivity, terra::vect(buffer0), updatevalue = 1, inverse = TRUE)

  # Inside the search zone raster, only keep pixels with high conductivity
  conductivity_search_zone <- terra::mask(conductivity, terra::vect(buffer2))
  conductivity_search_zone[conductivity_search_zone < m_conductivity] <- NA
  # plot(conductivity, col = viridis::viridis(50))
  # plot(conductivity_search_zone,  col = viridis::viridis(50))

  # Find the biggest clump of pixels. This is the road and its branches (hopefully).
  # (Can be reworked easily to ensure it is really the road)
  conductivity_search_zone_clump <- conductivity_search_zone |> terra::patches(directions = 8)
  tmp <- table(conductivity_search_zone_clump[])
  best_patch <- as.numeric(names(which.max(tmp)))
  conductivity_search_zone_clump[conductivity_search_zone_clump != best_patch] <- NA
  # plot(conductivity_search_zone_clump)

  # Extract the branch from the road my masking the road
  conductivity_search_zone_clump <- terra::mask(conductivity_search_zone_clump, terra::vect(search_zone))
  #plot(conductivity_search_zone_clump)

  # Recompute again patches to remove small clumps
  conductivity_search_zone_clump <- terra::patches(conductivity_search_zone_clump, directions = 8)
  tmp <- table(conductivity_search_zone_clump[])
  valid_patches <- as.numeric(names(tmp[tmp > 20 & tmp < 500]))
  conductivity_search_zone_clump[!conductivity_search_zone_clump[] %in% valid_patches] <- NA
  # plot(conductivity_search_zone_clump, col = lidR::random.colors(25))

  # Convert clumps into point
  pts_clump <- terra::as.points(conductivity_search_zone_clump)
  pts_clump_conductivity <- terra::extract(conductivity, pts_clump)
  pts_clump$conductivity <- pts_clump_conductivity[[2]]
  pts_clump <- sf::st_as_sf(pts_clump)
  pts_clump <- dplyr::group_by(pts_clump, patches) |> dplyr::summarise(avg_conductivity = mean(conductivity))
  sf::st_agr(pts_clump) <- "constant"
  # plot(conductivity, col = viridis::viridis(50))
  # plot(pts_clump["avg_conductivity"], cex = 0.1, pch = 19, add = T)

  # linear model to get the orientation of the patches and remove non linear patches
  r <- numeric(nrow(pts_clump))
  for (i in seq_along(pts_clump$geometry))
  {
    #plot(pts_clump$geometry[i])
    xy <- sf::st_coordinates(pts_clump$geometry[i])[,-3]
    mx <- cov(xy)
    ei <- eigen(mx)
    r[i] <- ei$values[1]/ei$values[2]
  }
  pts_clump$r = abs(r)
  sf::st_agr(pts_clump) <- "constant"
  #plot(conductivity, col = viridis::viridis(50))
  #plot(pts_clump["r"], cex = 0.2, add = T, pal = viridis::inferno(20))
  #plot(pts_clump["r"], cex = 0.2)
  pts_clump <- pts_clump[pts_clump$r > 4,]
  pts_clump <- pts_clump[pts_clump$avg_conductivity > 0.5,]
  #plot(conductivity, col = viridis::viridis(50))
  #plot(pts_clump["avg_conductivity"], cex = 0.25, pch = 19, add = T)

  # Compute the centroid of each clump. This is the end for a LCP
  ends <- sf::st_centroid(pts_clump) |> sf::st_geometry()
  #plot(conductivity,  col = viridis::viridis(50))
  #plot(st_geometry(ends), cex = 1, pch = 19, col = "red", add = T)

  starts <- sf::st_nearest_points(ends, road) |> lwgeom::st_endpoint()
  #plot(sf::st_geometry(start), cex = 1, pch = 19, col = "pink", add = T)

  segments <- vector("list", length(stars))
  for(i in seq_along(starts))
  {
    bbox <- sf::st_bbox(c(starts[i], ends[i]))
    bbox <- sf::st_as_sfc(bbox)
    bbox <- sf::st_buffer(bbox, 10)
    sigma <- raster::crop(conductivity, bbox)
    #plot(sigma,  col = viridis::viridis(50))
    #plot(st_geometry(starts[i]), cex = 1, pch = 19, col = "red", add = T)
    #plot(st_geometry(ends[i]), cex = 1, pch = 19, col = "red", add = T)
    trans <- MFFProads:::transition(raster::raster(sigma))
    path  <- gdistance::shortestPath(trans, sf::as_Spatial(starts[i]), sf::as_Spatial(ends[i]), output = "SpatialLines") |> suppressWarnings()
    cost  <- gdistance::costDistance(trans, sf::as_Spatial(starts[i]), sf::as_Spatial(ends[i]))
    path  <- sf::st_as_sf(path)
    path  <- sf::st_set_crs(path, sf::NA_crs_)
    path  <- sf::st_set_crs(path, sf::st_crs(road))
    path  <- sf::st_simplify(path, dTolerance = 2)
    cost  <- as.numeric(cost/sf::st_length(path))
    path$cost <- cost
    #plot(path, add = T, col = "red")
    segments[[i]] <- path
  }
  segments <- do.call(rbind, segments)
  segments <- segments[segments$cost < 1.5,]
  #plot(conductivity,  col = viridis::viridis(50))
  #plot(segments["cost"], add = T,lwd = 3)

  segments <- sf::st_difference(segments, sf::st_buffer(road, 5))
  #plot(conductivity,  col = viridis::viridis(50))
  #plot(segments, add = T, col = "red",lwd = 3)

  return(sf::st_geometry(segments))
}


#' Conductivity raster from boolean raster
#'
#' Conductivity raster from boolean raster. Can be useful if one has a road raster generated
#' by classification algorithm and whant to extract centerlines from it. The moving window
#' will enhance the center of an object. In cases of roads (narrow elongated shape), the centerline
#' will becone more apparent.
#'
#' @param x  raster (\code{raster} format). Boolean raster (1/0) where all 1 represent a pixel that might be categorized as a road.
#' @param w_max  numeric. Maximum window size to use (start from 3x3). Must be an odd number over 3.
#'
#' @return  raster (\code{raster} format) of conductivity with values ranging from 0 to 1.
#' @export
conductivity_from_bool <- function(x, w_max = 15)
{
  if ((w_max %% 2 == 0) | (w_max <= 3)) stop("'w_max' must be an odd number over 3.", call. = FALSE)

  conductivity <- seq(3, w_max, 2) |>
    lapply(function(w) raster::focal(x, matrix(1, w, w))/w^2) |>
    raster::brick() |>
    raster::calc(fun = mean)

  return(conductivity)
}

# Very similar to st_ends_heading()
get_heading <- function(from, to)
{
  from <- sf::st_coordinates(from)
  to   <- sf::st_coordinates(to)
  M <- rbind(from, to)
  Ax <- M[2,1]
  Ay <- M[2,2]
  Bx <- M[1,1]
  By <- M[1,2]
  heading <- atan2(Ay-By, Ax-Bx)
}

query_conductivity_aoi <- function(conductivity, starting_road)
{
  resolution <- raster::res(conductivity)[1]

  buf_dist <- c("ahead" = 900, "behind" = 100, "side" = 400)  # Could be set as a parameter

  aoi <- starting_road |>
    st_extend_line(c(buf_dist["ahead"], buf_dist["behind"])) |>
    sf::st_buffer(buf_dist["side"], endCapStyle = "FLAT")
  #plot(aoi, add = T, col= "red")
  #plot(starting_road, add = T, lwd = 3)

  conductivity_crop <- raster::crop(conductivity, sf::as_Spatial(aoi))
  #plot(conductivity_crop, col = viridis::viridis(50))
  #plot(starting_road, add = T, lwd = 3, col = "red")

  resolution_min <- 2  # Could be set as a parameter
  if (resolution < resolution_min)
  {
    agg_factor <- ceiling(resolution_min / resolution)
    resolution <- resolution * agg_factor
    cat("Downsampling conductivity to", resolution, "m\n")
    conductivity_crop <- raster::aggregate(conductivity_crop, fact = agg_factor, fun = mean, na.rm = TRUE)
  }

  conductivity_crop[is.na(conductivity_crop)] <- 0

  return(conductivity_crop)
}

#' Full workflow example
#'
#' An example of full workflow
#'
#' @examples
#' library(terra)
#' library(raster)
#' library(sf)
#'
#' conductivity <- system.file("extdata", "j53e_conductivity.tif", package = "MFFProads")
#' starting <- system.file("extdata", "j53e_drived_road.gpkg", package = "MFFProads")#'
#' conductivity <- raster(conductivity)
#' starting <- st_read(starting, "starting_road", quiet = TRUE)
#'
#' road <- drive_road(starting, conductivity)
#'
#' plot(conductivity, col = viridis::viridis(50), legend = FALSE)
#' plot(starting, add = TRUE, col = "green", lwd = 3)
#' plot(road, add = TRUE, col = "red", lwd = 2)
#'
#' road_branches <- find_road_branches(road, conductivity)
#'
#' plot(road_branches, col = "lightblue", add = TRUE, lwd = 3)
#'
#' new_roads <- vector("list", length(road_branches))
#' for (i in seq_along(road_branches))
#'   new_roads[[i]] <- drive_road(road_branches[i], conductivity)
#'
#' new_roads <- do.call(c, new_roads)
#' plot(new_roads, add = TRUE, col = "red", lwd = 2)
NULL
