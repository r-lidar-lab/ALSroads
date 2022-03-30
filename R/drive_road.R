#' Drive along an unknown road
#'
#' Drive along an unknown road in a conductivity raster starting only with a small piece
#' of road segment pointing in the right direction.
#'
#' @param seed  line (\code{sf} format)
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
#' conductivity <- system.file("extdata", "j53e_conductivity.tif", package = "ALSroads")
#' starting <- system.file("extdata", "j53e_drived_road.gpkg", package = "ALSroads")
#'
#' conductivity <- raster(conductivity)
#' starting <- st_read(starting, "seed", quiet = TRUE)
#'
#' res <- drive_road(starting, conductivity)
#'
#' plot(conductivity, col = viridis::viridis(50))
#' plot(starting, add = TRUE, col = "green", lwd = 3)
#' plot(res, add = TRUE, col = "red", lwd = 2)
drive_road <- function(seed, conductivity, existing_network = NULL, fov = 85, sightline = 50, min_conductivity = 0.4, ..., disp = FALSE)
{
  t0 <- Sys.time()

  seed = sf::st_simplify(seed, dTolerance = 5)

  # Threshold of cost to stop the search.
  cost_max <- sightline * 1/min_conductivity

  # Constant definition
  # The idea is to only extract values that are mostly ahead of the line from a potentially very large VRT
  buf_dist <- c("ahead" = 900, "behind" = 100, "side" = 400)  # Could be set as a parameter

  # We need an orthonormal raster
  resolution <- raster::res(conductivity)
  if (resolution[1] != resolution[2]) stop("'conductivity' raster must have the same resolution in both X and Y axis.", call. = FALSE)
  resolution <- resolution[1]

  # Initialization of list of coordinates with starting line
  start_pts <- sf::st_cast(sf::st_geometry(seed), "POINT")
  n <- length(start_pts)
  start <- start_pts[n-1][[1]]
  end   <- start_pts[n][[1]]
  list_coords <- list(start, end)
  list_lines <- list(sf::st_geometry(seed)[[1]])
  list_intersections <- list()
  heading <- get_heading(start, end)

  sub_aoi_conductivity <- query_conductivity_aoi(conductivity, seed)
  if (!raster::inMemory(sub_aoi_conductivity))
    stop("Internal error: sub_aoi_conductivity is not loaded in memory")

  resolution <- raster::res(sub_aoi_conductivity)[1]
  #plot(conductivity, col = viridis::inferno(50))
  #plot(raster::extent(sub_aoi_conductivity), col = "red", add = T)
  #plot(seed, add = T, lwd = 3, col = "red")
  #plot(sub_aoi_conductivity, col = viridis::inferno(50))

  # Mask the existing road to avoid driving a known road
  if (!is.null(existing_network) && length(existing_network) > 0)
  {
    bb <- sf::st_bbox(sub_aoi_conductivity)
    bb <- sf::st_set_crs(bb, sf::st_crs(existing_network))
    mask <- sf::st_crop(existing_network, bb)
    if (length(mask) > 0)
    {
      mask <- sf::st_buffer(mask, dist = 10)
      #plot(mask, add = T, col = "red")
      sub_aoi_conductivity <- raster::mask(sub_aoi_conductivity, sf::as_Spatial(mask), inverse = TRUE, updatevalue = 0)
      #plot(sub_aoi_conductivity, col = viridis::inferno(50))
    }
  }

  # Init view angles as a function of the resolution of the raster and the sightline
  angles_rad <- generate_angles(resolution, sightline, fov)

  # Loop initialisation
  current_cost <- 0
  k <- 2
  overcost = 0
  disp = F
  cat("Driving the conductivity raster\n")
  while (overcost <= 1 & !is.infinite(cost_max))
  {
    if (k %% 2 == 0) cat("", (k-1)*sightline, " m\r", sep = "")
    flush.console()

    # Compute heading from the two previous points
    p1 <- list_coords[[k-1]]
    p2 <- list_coords[[k]]
    heading <- get_heading(p1, p2)

    # Create all possible ends ahead of the previous points
    start <- sf::st_sfc(p2)
    ends <- generate_ends(p2, angles_rad, sightline, heading)
    ends$angle = angles_rad

    # Generate a very small area of interest to analyse
    search_zone <- sf::st_bbox(c(start, sf::st_geometry(ends))) + c(-1,-1,1,1)*resolution
    aoi <- raster::crop(sub_aoi_conductivity, search_zone)
    #aoi <- terra::crop(sub_aoi_conductivity, terra::ext(c(search_zone[1], search_zone[3],search_zone[2],search_zone[4]))) # -10 ms

    pol <- c(start, ends$geometry, start)
    pol <- sf::st_coordinates(pol)
    pol <- sf::st_polygon(list(pol))
    pol <- sf::st_sfc(pol)
    pol <- sf::st_buffer(pol, resolution)
    #aoi <- raster::mask(aoi, sf::as_Spatial(pol))
    #plot(sf::st_as_sfc(search_zone))
    #plot(aoi, col = viridis::viridis(50), add = T)
    #plot(start, add = T, pch = 19, col = "red")
    #plot(ends, pch = 19, add = T, col = "red")
    #plot(pol, add = T, border = "green")

    # Compute the transition matrix for this AOI
    trans <- transition(aoi, geocorrection = TRUE)
    #trans <- transition(raster::raster(aoi), geocorrection = TRUE) # + 20 ms

    # Check if some ends fall outside of the AOI conductivity raster
    # If any reload another raster part further ahead.
    val <- raster::extract(aoi, ends)
    #terra:::extract(aoi, st_coordinets(ends)) # - 2 ms

    if (any(is.na(val)))
    {
      #cat("\nReaching bounds of AOI. Loading next AOI\n")

      # Make a newly sub aoi conductivity raster further ahead
      line <- sf::st_cast(c(p1, p2), "LINESTRING")
      new_aoi_poly <- line |>
        st_extend_line(c(buf_dist["ahead"], buf_dist["behind"])) |>
        sf::st_buffer(buf_dist["side"], endCapStyle = "FLAT") |>
        sf::st_sfc() |>
        sf::st_set_crs(sf::st_crs(seed))

      sub_aoi_conductivity <- query_conductivity_aoi(conductivity, new_aoi_poly)
      #plot(conductivity, col = viridis::inferno(50))
      #plot(raster::extent(sub_aoi_conductivity), col = "red", add = T)
      #plot(seed, add = T, lwd = 3, col = "red")

      if (!is.null(existing_network) && length(existing_network) > 0)
      {
        bb <- sf::st_bbox(sub_aoi_conductivity)
        bb <- sf::st_set_crs(bb, sf::st_crs(existing_network))
        mask <- sf::st_crop(existing_network, bb)

        if (length(mask) > 0)
        {
          mask <- sf::st_buffer(mask, dist = 10)
          #plot(mask, add = T, col = "red")
          sub_aoi_conductivity <- raster::mask(sub_aoi_conductivity, sf::as_Spatial(mask), inverse = TRUE, updatevalue = 0)
          #plot(sub_aoi_conductivity, col = viridis::inferno(50))
        }
      }

      # Check again if some ends fall outside of the newly cropped conductivity raster
      # If yes we are close to the edge of the raster. Try again with half the sightline
      val <- raster::extract(sub_aoi_conductivity, ends)
      if (any(is.na(val)))
      {
        tmp_angles_rad <- generate_angles(resolution, sightline/2, fov)
        ends <- generate_ends(p2, tmp_angles_rad, sightline/2, heading)
        ends$angle <- tmp_angles_rad
      }

      # If we still have NA we reached the border of the conductivity map
      val <- raster::extract(sub_aoi_conductivity, ends)
      if (any(is.na(val)))
      {
        warning("Drive stopped early. Edge of conductivity raster has been reached.", call. = FALSE)
        cost_max <- -Inf
      }

      # Re-extract the AOI and recompute the transition
      aoi <- raster::crop(sub_aoi_conductivity, search_zone)
      trans <- transition(aoi, geocorrection = TRUE)
    }

    # Find cost of each edge point
    if (!is.infinite(cost_max))
    {
      ans <- find_reachable(start, ends, trans, cost_max)
      if (is.null(ans))
      {
        cost_max = -Inf
        break
      }

      idx_main <- ans$idx_main
      idx_other <- ans$idx_other
      cost <- ans$cost
      current_cost <- cost[idx_main]

      ends_other <- ends[idx_other,]
      cost_other <- cost[idx_other]
      ends_other <- ends_other[cost_other < current_cost*1.3,]


      end <- ends[idx_main,]

      W <- sf::st_geometry(end)[[1]]
      L <- gdistance::shortestPath(trans, sf::st_coordinates(start), sf::st_coordinates(end), output = "SpatialLines") |> suppressWarnings()
      L <- sf::st_geometry(sf::st_as_sf(L))
      list_coords[[k+1]] <- W
      list_lines[[k]] <- L[[1]]

      # Protect against infinite loops
      mask <- sf::st_buffer(L, dist = 5, endCapStyle = "FLAT")
      sub_aoi_conductivity <- raster::mask(sub_aoi_conductivity, sf::as_Spatial(mask), inverse = TRUE, updatevalue = 0)
      #sub_aoi_conductivity <- terra::mask(sub_aoi_conductivity, terra::vect(mask), inverse = TRUE, updatevalue = 0) # -20 ms

      if (current_cost > cost_max)
        overcost = overcost + 1
      else
        overcost = 0

      I <- NULL
      if (length(sf::st_geometry(ends_other)) > 0)
      {
        I <- gdistance::shortestPath(trans, sf::st_coordinates(start), sf::st_coordinates(ends_other), output = "SpatialLines") |> suppressWarnings()
        I <- sf::st_geometry(sf::st_as_sf(I))
        I <- sf::st_difference(I,L)
        for(II in I) list_intersections <- c(list_intersections, list(II))
      }

      # >>>>>>>>>>>>>>>>>>>
      rot = function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)
      cntrd = rbind(sf::st_coordinates(start), sf::st_coordinates(end))
      cntrd = c(mean(cntrd[,1]), mean(cntrd[,2]))
      #cntrd = sf::st_point(cntrd) |> sf::st_sfc()

      M = sf::st_coordinates(start)
      M = rotate_pi(M, cntrd)
      start2 = sf::st_point(M) |> sf::st_sfc()

      M = sf::st_coordinates(ends)
      M = rotate_pi(M, cntrd)
      ends2 = M %>% as.data.frame %>% sf::st_as_sf(coords = c(1,2))
      ends2$angle = ends$angle
      ends2 <- ends2[abs(ends2$angle) > 15*pi/180,]
      tmp <- raster::extract(aoi, sf::st_coordinates(ends2))
      ends2 <- ends2[!is.na(tmp),]

      #plot(aoi, col = viridis::inferno(50))
      #plot(start, add = TRUE, pch = 19)
      #plot(ends, col = 'red', add = TRUE, pch = 19)
      #plot(start2, add = TRUE, pch = 19, col = "black")
      #plot(ends2, col = 'red', add = TRUE, pch = 19)

      ans2 = find_reachable(start2, ends2, trans, cost_max)
      I2 <- NULL
      if (!is.null(ans2))
      {
        idx_main2 <- ans2$idx_main
        idx_other2 <- ans2$idx_other
        idx2 <- c(idx_main2, idx_other2)
        cost2 <- ans2$cost
        ends2 <- ends2[idx2,]
        cost2 <- cost2[idx2]
        ends2 <- ends2[cost2 < current_cost*1.3,]
        #plot(cost)
        #plot(ends2, add =T, pch = 19, col = "green")

        if (length(sf::st_geometry(ends2)) > 0)
        {
          I2 <- gdistance::shortestPath(trans, sf::st_coordinates(start2), sf::st_coordinates(ends2), output = "SpatialLines") |> suppressWarnings()
          I2 <- sf::st_geometry(sf::st_as_sf(I2))
          I2 <- sf::st_difference(I2,L)
          I2 <- I2[as.numeric(sf::st_length(I2)) > 0.5 * sightline]
          for(II in I2) list_intersections <- c(list_intersections, list(II))
        }
      }
      k <- k + 1

      if (disp)
      {
        if (length(angles_rad) == length(cost))
        {
          plot(angles_rad, cost, type = "b")
          lines(angles_rad, ma(cost, n = 5), type = "b", col = "red")
          abline(v = angles_rad[c(idx_main, idx_other)])
          abline(h = cost_max)
        }
        raster::plot(raster::crop(sub_aoi_conductivity, raster::extent(ends) + 100), col = viridis::inferno(50))
        ends$cost = cost
        tryCatch({
          plot(ends["cost"], pal = viridis::viridis, pch = 19, add = T, breaks = "quantile")
        }, error = function(x) {
          plot(ends$geometry, col = "red", add = T, pch = 19, cex = 0.5)
        })
        path = gdistance::shortestPath(trans, sf::as_Spatial(start), sf::as_Spatial(ends), output = "SpatialLines") |> suppressWarnings()

        plot(p1, col = "red", pch = 19, add = T)
        plot(p2, col = "red", pch = 19, add = T)
        sp::plot(path, add = T, col = "red")
        plot(L, add = T, col = "red", lwd = 3)
        if(!is.null(I)) plot(I, add = T, col = "blue", lwd = 3)
        if(!is.null(I2)) plot(I2, add = T, col = "cornflowerblue", lwd = 3)
        plot(sf::st_sfc(W), add = T, col = "green", pch = 19, cex = 2)
        print(current_cost)
      }
    }
  }

  list_lines = list_lines[1:(length(list_lines)-overcost)]

  newline <- sf::st_sfc(list_lines) |>
    sf::st_coordinates() |>
    sf::st_linestring() |>
    sf::st_sfc() |>
    sf::st_set_crs(sf::st_crs(seed)) |>
    sf::st_simplify(dTolerance = 2)

  if (length(list_intersections) > 1)
  {
    intersections <- sf::st_sfc(list_intersections) |>
    sf::st_set_crs(sf::st_crs(seed)) |>
    sf::st_simplify(dTolerance = 2)
  }
  else
  {
    intersections <- NULL
  }

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

  return(list(road = newline, seeds = intersections))
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

query_conductivity_aoi <- function(conductivity, seed)
{
  resolution <- raster::res(conductivity)[1]

  buf_dist <- c("ahead" = 900, "behind" = 100, "side" = 400)  # Could be set as a parameter

  aoi <- seed |>
    st_extend_line(c(buf_dist["ahead"], buf_dist["behind"])) |>
    sf::st_buffer(buf_dist["side"], endCapStyle = "FLAT")
  #plot(aoi, add = T, col= "red")
  #plot(seed, add = T, lwd = 3)

  conductivity_crop <- raster::crop(conductivity, sf::as_Spatial(aoi))
  #plot(conductivity_crop, col = viridis::viridis(50))
  #plot(seed, add = T, lwd = 3, col = "red")

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
#' conductivity <- system.file("extdata", "j53e_conductivity.tif", package = "ALSroads")
#' starting <- system.file("extdata", "j53e_drived_road.gpkg", package = "ALSroads")#'
#' conductivity <- raster(conductivity)
#' starting <- st_read(starting, "seed", quiet = TRUE)
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


generate_angles <- function(resolution, sightline, fov)
{
  angular_resolution <- floor(1.5 * (resolution / sightline) * (180 / pi))
  angles <- seq(angular_resolution, fov, angular_resolution)
  angles <- c(rev(-angles), 0, angles)
  angles_rad <- angles * pi / 180
  angles_rad
}

generate_ends <- function(center, angles, sightline, direction)
{
  center <- sf::st_coordinates(center)
  X <- center[1] + sightline * cos(direction + angles)
  Y <- center[2] + sightline * sin(direction + angles)
  M <- data.frame(X,Y)
  sf::st_as_sf(M, coords = c("X", "Y"))
}


local_maxima <- function(x)
{
  # Use -Inf instead if x is numeric (non-integer)
  x = na.omit(x)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}

local_maxima_with_height <- function(x)
{
  x = na.omit(x)
  i = local_maxima(x)
  j = local_maxima(-x)

  if (length(i) == 0 | length(j) == 0)
    return(data.frame(idx = 0, depth = 0))

  #plot(x)
  #abline(v = i, col = "blue")
  #abline(v = j, col = "red")

  max = data.frame(idx = i, max = T)
  min = data.frame(idx = j, max = F)
  df = rbind(min, max)
  df = df[order(df$idx),]

  depths = numeric(length(i))
  k = 1
  for (u in seq_along(df$idx))
  {
    if (df$max[u] == TRUE)
    {
      ldepth = x[df$idx[u]] - x[df$idx[u-1]]
      rdepth = x[df$idx[u]] - x[df$idx[u+1]]
      if (length(ldepth) == 0 || is.na(ldepth)) ldepth = 0
      if (length(rdepth) == 0 || is.na(rdepth)) rdepth = 0
      depths[k] = min(ldepth, rdepth)
      k = k+1
    }
  }

  return(data.frame(idx = i, depth = depths))
}

find_reachable <- function(start, ends, trans, cost_max)
{
  if (length(sf::st_geometry(ends)) < 6) return(NULL)

  cost <- gdistance::costDistance(trans, sf::st_coordinates(start), sf::st_coordinates(ends))
  cost <- as.numeric(cost)
  if (all(is.infinite(cost))) return(NULL)
  cost[is.infinite(cost)] <- 2 * max(cost[!is.infinite(cost)])

  # Select local minima of cost
  smooth <- 5
  minima <- local_maxima_with_height(-ma(cost, n = smooth))
  if (sum(minima$depth > 4) > 0) {
    minima <- minima[minima$depth > 4,]
    minima <- minima$idx + as.integer(smooth/2)
    minima <- sapply(minima, function(i){
      j = (i-3):(i+3)
      return(j[which.min(cost[j])])
    })
  } else {
    minima <- which.min(cost)
  }

  keep = cost[minima] <= cost_max
  if (sum(keep) > 0) minima = minima[keep]

  # Main road
  idx_main <- which.min(cost)


  # Potential intersection
  idx_other <- minima[minima != idx_main]

  return(list(idx_main = idx_main, idx_other = idx_other, cost = cost))
}

rotate_pi <- function(X, C)
{
  M = matrix(c(-1,0,0,1), 2,2)
  X[,1] = X[,1] - C[1]
  X[,2] = X[,2] - C[2]
  X = X %*% M
  X[,1] = X[,1] + C[1]
  X[,2] = X[,2] + C[2]
  X
}
