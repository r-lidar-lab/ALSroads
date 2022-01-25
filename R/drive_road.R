# Very similar to st_ends_heading()
get_heading <- function(from, to)
{
  M <- rbind(from, to)
  Ax <- M[2,1]
  Ay <- M[2,2]
  Bx <- M[1,1]
  By <- M[1,2]
  heading <- atan2(Ay-By, Ax-Bx)
}

#' Drive along an unknown road
#'
#' Attempt to drive along an unknown road in a conductivity raster starting only with a small piece
#' of road segment pointing in the right direction.
#'
#' @param starting_road  line (\code{sf} format)
#' @param conductivity  raster (\code{raster} format)
#' @param fov  numeric. Field of view (degrees) ahead of the search vector.
#' @param radius  numeric (distance unit). Search radius used to find the next most probable point on the road.
#' @param cost_max  numeric. Maximal cost allowed in the conductivity raster for point candidate to continue the search.
#' If no value is provided \code{radius * 6} will be used.
#'
#' @return list, \code{line} being the road found (as \code{sfc}) and \code{cost} a numeric vector
#' of cost at each invidual vertices of the line.
#' @export
#' @examples
#' library(raster)
#' library(sf)
#'
#' conductivity <- system.file("extdata", "drived_conductivity.tif", package = "MFFProads")
#' drived_road <- system.file("extdata", "drived_road.gpkg", package = "MFFProads")
#'
#' conductivity <- raster(conductivity)
#' starting_line <- st_read(drived_road, "starting_road", quiet = TRUE)
#'
#' res <- drive_road(starting_line, conductivity)
#'
#' plot(conductivity, col = viridis::viridis(50))
#' plot(res$line, add = TRUE, col = "red", lwd = 2)
#' plot(starting_line, add = TRUE, col = "green", lwd = 3)
#' plot(res$cost, type = "l")
drive_road <- function(starting_road, conductivity, fov = 45, radius = 10, cost_max = NULL)
{
  if (is.null(cost_max)) cost_max <- radius * 6

  resolution <- raster::res(conductivity)
  if (resolution[1] != resolution[2]) stop("'conductivity' raster must have the same resolution in both X and Y axis.", call. = FALSE)
  resolution <- resolution[1]

  # Initialisation of list of coordinates with starting line
  start_coords <- sf::st_coordinates(starting_road)
  start_coords[,3] <- 0
  list_coords <- list(start_coords[1,], start_coords[2,])

  # Crop conductivity raster
  # The idea is to only extract values that are mostly ahead of the line
  p1 <- list_coords[[1]]
  p2 <- list_coords[[2]]
  heading <- get_heading(p1, p2)

  buf_dist <- c("ahead" = 900, "behind" = 100, "side" = 400)  # Could be set as a parameter

  aoi <- starting_road |>
    st_extend_line(c(buf_dist["ahead"], buf_dist["behind"])) |>
    sf::st_buffer(buf_dist["side"], endCapStyle = "FLAT") |>
    sf::as_Spatial()

  conductivity_crop <- raster::crop(conductivity, aoi)

  resolution_min <- 2  # Could be set as a parameter
  if (resolution < resolution_min)
  {
    agg_factor <- ceiling(resolution_min / resolution)
    resolution <- resolution * agg_factor
    cat("Step 0: downsample conductivity to", resolution, "m\n")
    conductivity_crop <- raster::aggregate(conductivity_crop, fact = agg_factor, fun = mean, na.rm = TRUE)
  }

  if (radius < resolution) stop("Radius must be equal or larger than resolution of 'conductivity' raster.", call. = FALSE)

  # Set possible search angles
  angular_resolution <- floor((resolution / radius) * (180 / pi))
  angles <- seq(angular_resolution, fov, angular_resolution)
  angles <- c(rev(-angles), 0, angles)
  angles_rad <- angles * pi / 180

  # Set penalty coefficient for each search angle
  penalty_at_45_degrees <- 1.5  # Could be set as a parameter
  penalty <- abs(angles) / fov
  penalty <- penalty * (penalty_at_45_degrees - 1) + 1

  # Transition matrix
  cat("Step 1: computation of the transition matrix\n")
  trans <- conductivity_crop |>
    transition()

  # Loop initialisation
  cost_selected <- 0
  k <- 2
  cat("Step 2: driving the conductivity raster\n")
  while (cost_selected <= cost_max)
  {
    if (k %% 5 == 0) cat(" | ", (k-1)*10, " m", sep = "")

    # Compute heading from the two previous points
    p1 <- list_coords[[k-1]]
    p2 <- list_coords[[k]]
    heading <- get_heading(p1, p2)

    # Create all possible ends ahead of the previous point
    X <- p2[1] + radius * cos(heading + angles_rad)
    Y <- p2[2] + radius * sin(heading + angles_rad)

    M <- data.frame(X,Y)
    ends <- sf::st_as_sf(M, coords = c("X", "Y")) |>
      sf::as_Spatial()

    # Check if some ends fall outside of the cropped conductivity raster
    val <- raster::extract(conductivity_crop, M)
    if (any(is.na(val)))
    {
      # Make a newly cropped conductivity raster further ahead
      line <- c(p1[1:2], p2[1:2]) |>
        matrix(ncol = 2, byrow = TRUE) |>
        sf::st_linestring() |>
        sf::st_sfc(crs = sf::st_crs(starting_road))

      aoi <- line |>
        st_extend_line(c(buf_dist["ahead"], buf_dist["behind"])) |>
        sf::st_buffer(buf_dist["side"], endCapStyle = "FLAT") |>
        sf::as_Spatial()

      conductivity_crop <- raster::crop(conductivity, aoi)

      resolution <- raster::res(conductivity)[1]
      resolution_min <- 2  # Could be set as a parameter
      if (resolution < resolution_min)
      {
        agg_factor <- ceiling(resolution_min / resolution)
        resolution <- resolution * agg_factor
        conductivity_crop <- raster::aggregate(conductivity_crop, fact = agg_factor, fun = mean, na.rm = TRUE)
      }

      # Check again if some ends fall outside of the newly cropped conductivity raster
      val <- raster::extract(conductivity_crop, M)
      if (any(is.na(val)))
      {
        warning("Drive stopped early. Edge of conductivity raster has been reached.", call. = FALSE)
        cost_max <- -Inf
      }

      # Generate a new transition matrix
      trans <- conductivity_crop |>
        transition()
    }

    # Find cost of each edge point
    start <- sf::st_point(p2) |>
      sf::st_sfc() |>
      sf::as_Spatial()
    cost <- as.numeric(gdistance::costDistance(trans, start, ends)) |> suppressWarnings()
    cost[is.infinite(cost)] <- 2 * max(cost[!is.infinite(cost)])

    # Adjust cost based on angle penalties
    # and smooth value of adjacent points
    cost2 <- ma(cost * penalty)

    # Select point coordinates with the lowest
    # penalty adjusted cost
    idx_min <- which.min(cost2)
    cost_selected <- cost[idx_min]
    W <- M[idx_min,] |>
      as.matrix() |>
      cbind(cost_selected)

    # Save the least-cost pixel coordinates
    k <- k + 1
    list_coords[[k]] <- W
  }

  m_coords <- do.call(rbind, list_coords)
  newline <- m_coords[-nrow(m_coords), 1:2] |>
    sf::st_linestring() |>
    sf::st_sfc(crs = sf::st_crs(starting_road)) |>
    sf::st_simplify(dTolerance = 2)

  return(list(line = newline, cost = m_coords[,3]))
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
#' library(raster)
#'
#' road <- system.file("extdata", "drived_road.gpkg", package="MFFProads")
#' conductivity <- system.file("extdata", "drived_conductivity.tif", package="MFFProads")
#'
#' road <- st_read(road, "drived_road", quiet = TRUE)
#' conductivity <- raster(conductivity)
#' conductivity <- raster::aggregate(conductivity, fact = 2, fun = mean, na.rm = TRUE)
#'
#' road_branches <- find_road_branches_2(road, conductivity)
#'
#' plot(conductivity, col = viridis::viridis(50))
#' plot(road, col = "green", lwd = 3, add = TRUE)
#' plot(road_branches, col = "red", add = TRUE, lwd = 3)
find_road_branches <- function(road, conductivity)
{
  if (is(conductivity, "RasterLayer"))
    conductivity <- terra::rast(conductivity)

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
  valid_patches <- as.numeric(names(tmp[tmp > 20]))
  conductivity_search_zone_clump[!conductivity_search_zone_clump[] %in% valid_patches] <- NA
  # plot(conductivity_search_zone_clump, col = lidR::random.colors(25))

  # Convert clumps into point
  pts_clump <- terra::as.points(conductivity_search_zone_clump)
  pts_clump_conductivity <- terra::extract(conductivity, pts_clump)
  pts_clump$conductivity <- pts_clump_conductivity[[2]]
  pts_clump <- sf::st_as_sf(pts_clump)
  pts_clump <- dplyr::group_by(pts_clump, patches) |> dplyr::summarise(avg_conductivity = mean(conductivity))
  sf::st_agr(pts_clump) <- "constant"
  # plot(pts_clump["avg_conductivity"], cex = 0.25, pch = 19)

  # linear model to get the orientation of the patches and remove non linear patches
  slopes <- numeric(nrow(pts_clump))
  intercepts <- numeric(nrow(pts_clump))
  rsquare <- numeric(nrow(pts_clump))
  for (i in seq_along(pts_clump$geometry))
  {
    xy <- sf::st_coordinates(pts_clump$geometry[i])[,-3]
    xy <- as.data.frame(xy)
    lm <- stats::lm(Y~X, data = xy)
    co <- coefficients(lm)
    rsquare[i] <- summary(lm)$r.squared
    slopes[i] <- co[2]
    intercepts[i] <- co[1]
  }
  keep <- which(rsquare > 0.5)
  pts_clump <- pts_clump[keep,]
  slopes <- slopes[keep]
  intercepts <- intercepts[keep]
  #plot(pts_clump["avg_conductivity"], cex = 0.25, pch = 19)

  # Compute the centroid of each clump
  pts_clump_centroid <- sf::st_centroid(pts_clump)
  #plot(conductivity,  col = viridis::viridis(50))
  #plot(pts_clump_centroid["avg_conductivity"], cex = 1, pch = 19, pal = viridis::viridis, add = T)

  # compute the start of the road which is at the inersection between the segment of known
  # slope/intercept and the road
  xy <- sf::st_coordinates(pts_clump_centroid)
  segments <- vector("list", nrow(xy))
  for (i in 1:nrow(xy))
  {
    centroid <- xy[i,]
    xstart <- centroid[1] - 80
    xend   <- centroid[1] + 80
    ystart <- slopes[i] * xstart + intercepts[i]
    yend <- slopes[i] * xend + intercepts[i]
    m <- matrix(c(xstart, xend, ystart, yend), 2, 2)
    segments[[i]] <- sf::st_sfc(sf::st_linestring(m))
  }
  segments <- do.call(c, segments)
  segments <- sf::st_set_crs(segments, sf::st_crs(road))
  #plot(segments, add = T, col = "red",lwd = 3)

  ex_road <- st_extend_line(road, 100)   # extent the line to get intersection at the end
  starting_roads_pot <- sf::st_intersection(ex_road, segments)
  starting_roads_pot <- sf::st_cast(starting_roads_pot, "POINT")
  # plot(conductivity,  col = viridis::viridis(50))
  # plot(starting_roads_pot, col = "red", pch = 19, cex = 1, add = TRUE)

  init_segments <- vector("list", length(starting_roads_pot))
  for (i in 1:length(starting_roads_pot))
  {
    a <- starting_roads_pot[i]
    b <- sf::st_geometry(pts_clump_centroid)[i]
    init_segments[[i]] <- sf::st_cast(sf::st_union(a,b),"LINESTRING")
  }
  init_segments <- do.call(c, init_segments)
  # plot(conductivity,  col = viridis::viridis(50))
  # plot(init_segments, col = "red", add = TRUE, lwd = 3)

  return(init_segments)
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
#' @examples
#' library(raster)
#'
#' segmented_road <- system.file("extdata", "segmented_road.tif", package="MFFProads")
#' segmented_road <- raster(segmented_road)
#'
#' conductivity <- conductivity_from_bool(segmented_road)
#' plot(conductivity)
conductivity_from_bool <- function(x, w_max = 15)
{
  if ((w_max %% 2 == 0) | (w_max <= 3)) stop("'w_max' must be an odd number over 3.", call. = FALSE)

  conductivity <- seq(3, w_max, 2) |>
    lapply(function(w) raster::focal(x, matrix(1, w, w))/w^2) |>
    raster::brick() |>
    raster::calc(fun = mean)

  return(conductivity)
}

st_ends_heading <- function(line)
{
  M <- sf::st_coordinates(line)
  i <- c(2, nrow(M) - 1)
  j <- c(1, -1)

  headings <- mapply(i, j, FUN = function(i, j) {
    Ax <- M[i-j,1]
    Ay <- M[i-j,2]
    Bx <- M[i,1]
    By <- M[i,2]
    unname(atan2(Ay-By, Ax-Bx))
  })

  return(headings)
}

st_extend_line <- function(line, distance, end = "BOTH")
{
  if (!(end %in% c("BOTH", "HEAD", "TAIL")) | length(end) != 1) stop("'end' must be 'BOTH', 'HEAD' or 'TAIL'")

  M <- sf::st_coordinates(line)[,1:2]
  keep <- !(end == c("TAIL", "HEAD"))

  ends <- c(1, nrow(M))[keep]
  headings <- st_ends_heading(line)[keep]
  distances <- if (length(distance) == 1) rep(distance, 2) else rev(distance[1:2])

  M[ends,] <- M[ends,] + distances[keep] * c(cos(headings), sin(headings))
  newline <- sf::st_linestring(M)

  # If input is sfc_LINESTRING and not sfg_LINESTRING
  if (is.list(line)) newline <- sf::st_sfc(newline, crs = sf::st_crs(line))

  return(newline)
}
