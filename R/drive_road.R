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

  if (is(conductivity, "RasterLayer")) stop("conductivity must be a SpatRaster")

  trace = sf::st_geometry(seed)

  # Because we want an "overall" direction and we want to drop smooth little imperfections
  seed = sf::st_simplify(seed, dTolerance = 5)

  # Threshold of cost to stop the search.
  cost_max <- sightline * 1/min_conductivity

  # We need an orthogonal raster
  resolution <- terra::res(conductivity)
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

  sub_aoi_conductivity <- query_conductivity_aoi(conductivity, seed, smooth = F)

  # Mask the existing road to avoid driving a known road
  sub_aoi_conductivity <- mask_existing_network(sub_aoi_conductivity, existing_network)

  # Init view angles as a function of the resolution of the raster and the sightline
  angles_rad <- generate_angles(resolution, sightline, fov)

  # Loop initialisation
  current_cost <- 0
  k <- 2
  overcost = 0
  novercost = 0
  dovercost = 0

  cat("Driving the conductivity raster\n")
  while (dovercost <= 250 & !is.infinite(cost_max) )
  {
    cat("", (k-1)*sightline*0.8, " m\r", sep = "")
    flush.console()

    # Compute heading from the two previous points
    p1 <- list_coords[[k-1]]
    p2 <- list_coords[[k]]
    heading <- get_heading(p1, p2)

    # Create all possible ends ahead of the previous points
    start <- sf::st_sfc(p2)
    ends <- generate_ends(p2, angles_rad, sightline, heading)
    ends$angle = angles_rad

    # Generate a very small area of interest corresponding to what can be seen in the line of sight
    search_zone <- sf::st_bbox(c(start, sf::st_geometry(ends)))
    search_zone <- terra::ext(c(search_zone[1], search_zone[3],search_zone[2],search_zone[4])) + 10
    aoi <- terra::crop(sub_aoi_conductivity, search_zone)
    #aoi <- terra::stretch(aoi, maxv = 1)
    p <- sf::st_polygon(list(rbind(sf::st_coordinates(start), sf::st_coordinates(ends), sf::st_coordinates(start))))
    p <- sf::st_sfc(p)
    p <- sf::st_buffer(p, 2)
    p <- terra::vect(p)
    terra::crs(p) <- terra::crs(aoi)
    aoi <- terra::mask(aoi, p)
    #plot(aoi, col = viridis::inferno(50), range = c(0,1))
    #plot(start,add =T, col = "red", pch = 19)
    #plot(ends$geometry,add =T, col = "red", pch = 19)

    # Compute the transition matrix for this AOI
    trans <- transition(aoi, geocorrection = TRUE)

    # Check if some ends fall outside of the AOI
    # If any we  reached the loaded part of the conductivity.
    # Reload another raster part further ahead.
    val <- terra::extract(aoi, sf::st_coordinates(ends))[[1]]
    if (any(is.na(val)))
    {
      # Make a newly sub aoi conductivity raster further ahead
      line <- sf::st_cast(c(p1, p2), "LINESTRING") |> sf::st_sfc()
      line <- sf::st_set_crs(line, sf::st_crs(seed))

      sub_aoi_conductivity <- query_conductivity_aoi(conductivity, line, smooth = F)
      #plot(conductivity, col = viridis::inferno(50))
      #plot(raster::extent(sub_aoi_conductivity), col = "red", add = T)
      #plot(seed, add = T, lwd = 3, col = "red")

      sub_aoi_conductivity <- mask_existing_network(sub_aoi_conductivity, existing_network)

      # Check again if some ends fall outside of the newly cropped conductivity raster
      # If yes we are close to the edge of the raster. Try again with half the sightline
      val <- terra::extract(sub_aoi_conductivity, sf::st_coordinates(ends))[[1]]
      if (any(is.na(val)))
      {
        tmp_angles_rad <- generate_angles(resolution, sightline/2, fov)
        ends <- generate_ends(p2, tmp_angles_rad, sightline/2, heading)
        ends$angle <- tmp_angles_rad
      }

      # If we still have NA we reached the border of the conductivity map
      val <- terra::extract(sub_aoi_conductivity, sf::st_coordinates(ends))[[1]]
      if (any(is.na(val)))
      {
        warning("Drive stopped early. Edge of conductivity raster has been reached.", call. = FALSE)
        cost_max <- -Inf
        break
      }

      # Re-extract the AOI and recompute the transition
      aoi <- terra::crop(sub_aoi_conductivity, search_zone)
      #aoi <- terra::stretch(aoi, maxv = 1)
      trans <- transition(aoi, geocorrection = TRUE)
    }

    # Estimates the cost to reach each end points and returns the main direction and
    # other possible directions (i.e. intersections)
    ans <- find_reachable(start, ends, trans, cost_max)

    # The only way to get infinity values is to reach an existing road
    # We need to reduce the sightling to reduu=ce down the approach speed
    if (any(ans$cost >= 9999))
    {
      tmp_angles_rad <- generate_angles(resolution, sightline/2, fov)
      ends <- generate_ends(p2, tmp_angles_rad, sightline/2, heading)
      ends$angle <- tmp_angles_rad
      ans <- find_reachable(start, ends, trans, cost_max/2)
    }

    if (is.null(ans))
    {
      warning("Driving stopped because not reachable point have been found", call. = FALSE)
      cost_max = -Inf
      break
    }

    idx_main <- ans$minima$idx[1]
    depth_main = ans$minima$depth[1]
    idx_other <- ans$minima$idx[-1]
    cost <- ans$cost
    current_cost <- cost[idx_main]
    if (current_cost < 1.5*cost_max & depth_main > 500) current_cost = min(current_cost,cost_max)

    end <- ends[idx_main,]
    ends_other <- ends[idx_other,]
    cost_other <- cost[idx_other]

    L <- gdistance::shortestPath(trans, sf::st_coordinates(start), sf::st_coordinates(end), output = "SpatialLines")
    L <- sf::st_geometry(sf::st_as_sf(L))
    C <- sf::st_coordinates(L)
    C <- C[1:(nrow(C)*0.8),-3]
    L <- sf::st_sfc(sf::st_linestring(C))
    end <- lwgeom::st_endpoint(L)

    list_coords[[k+1]] <- sf::st_geometry(end)[[1]]
    list_lines[[k]] <- L[[1]][-nrow(L[[1]]),]
    trace <- st_join_linestring(trace, L)

    # Protect against infinite loops
    sub_aoi_conductivity <- mask_passage(sub_aoi_conductivity, trace, 0, 0.95,sf::st_crs(seed))
    #plot(sub_aoi_conductivity, col = viridis::inferno(50))

    if (current_cost > cost_max)
    {
      novercost = novercost + 1
      dovercost = dovercost + sightline
      overcost = overcost + (current_cost - cost_max)
    }
    else
    {
      novercost = 0
      dovercost = 0
      overcost = 0
    }

    I <- NULL
    if (length(sf::st_geometry(ends_other)) > 0)
    {
      I <- gdistance::shortestPath(trans, sf::st_coordinates(start), sf::st_coordinates(ends_other), output = "SpatialLines") |> suppressWarnings()
      I <- sf::st_geometry(sf::st_as_sf(I))
      I <- sf::st_difference(I,L)
      for(II in I)
      {
        if (sf::st_geometry_type(II) == "LINESTRING")
        {
          list_intersections <- c(list_intersections, list(II))
          sub_aoi_conductivity <- mask_passage(sub_aoi_conductivity, II, 0.05, 1, sf::st_crs(seed))
        }
      }
    }

    k <- k + 1

    if (disp)
    {
      if (length(angles_rad) == length(cost))
      {
        plot(angles_rad, cost, type = "l", ylim = c(sightline, max(cost) *1.1), col = "green", log = "y")
        #lines(angles_rad, spike_preserving_smooth(cost, n = 9), col = "darkgreen")
        lines(angles_rad, ma(cost, n = 7), col = "darkgreen")
        abline(v = angles_rad[c(idx_main, idx_other)], col = "darkgreen",lwd = 3)

        abline(h = cost_max, lty = 3)
        abline(h = 1.5*cost_max, lty = 3)
        abline(v = angles_rad[idx_main] + c(-0.2612, 0.2612), lty = 3)
        # if (!is.null(ans2))
        # {
        #   lines(angles_rad2, ans2$cost, col = "lightblue")
        #   lines(angles_rad2, spike_preserving_smooth(ans2$cost, n = 9), col = "blue")
        #   lines(angles_rad2, ma(ans2$cost, n = 7), col = "blue")
        #   abline(v = angles_rad2[c(idx2)], col = "blue")
        #
        # }
      }
      terra::plot(terra::crop(sub_aoi_conductivity, terra::ext(aoi) + 80), col = viridis::inferno(50), main = k, range = c(0,1))
      terra::plot(aoi, add = T, col = viridis::inferno(50))
      plot(terra::ext(aoi), add = T)
      ends$cost = cost
      tryCatch({
        plot(ends["cost"], pal = viridis::viridis, pch = 19, add = T, cex = 0.25, breaks = "quantile")
      }, error = function(x) {
        plot(ends$geometry, col = "red", add = T, pch = 19, cex = 0.5)
      })

      tryCatch({
        paths = gdistance::shortestPath(trans, sf::st_coordinates(start), sf::st_coordinates(ends), output = "SpatialLines") |> suppressWarnings()
        paths = sf::st_as_sf(paths)
        paths$cost = ends$cost
        plot(paths, add = T,  pal = viridis::viridis)
        plot(p1, col = "red", pch = 19, add = T)
        plot(p2, col = "red", pch = 19, add = T)
      }, error = function(x) print(x))

      plot(L, add = T, col = "red", lwd = 4)
      if(!is.null(I)) plot(I, add = T, col = "blue", lwd = 3)
      #if(!is.null(I2)) plot(I2, add = T, col = "cornflowerblue", lwd = 3)
      plot(end, add = T, col = "green", pch = 19, cex = 2)
      cat("cost =", current_cost, "\n")
      cat("overcost =", overcost, "\n")
    }

  }

  list_lines = list_lines[1:(length(list_lines)-novercost)]

  newline <- do.call(rbind, list_lines) |>
    sf::st_linestring() |>
    sf::st_sfc() |>
    sf::st_set_crs(sf::st_crs(seed))


  if (!is.null(existing_network))
  {
    tail_point = lwgeom::st_endpoint(newline)
    distance_to_network = sf::st_distance(existing_network, tail_point)
    if (length(existing_network) > 0 && as.numeric(distance_to_network) < 75)
    {
      u = sf::st_simplify(newline, dTolerance = 2)

      #plot(u, xlim = st_bbox(u)[c(1, 3)] + c(-20,20),  ylim = st_bbox(u)[c(2, 4)] + c(-20,20))
      #plot(existing_network, add = T)
      u = st_extend_line(u, 75, end = "TAIL" )
      #plot(u, xlim = st_bbox(u)[c(1, 3)] + c(-20,20),  ylim = st_bbox(u)[c(2, 4)] + c(-20,20))
      #plot(existing_network, add = T)
      p = sf::st_intersection(u, existing_network)
      if (length(p) == 1)
      {
        M = rbind(sf::st_coordinates(newline)[,1:2], sf::st_coordinates(p)[,1:2])
        newline = sf::st_linestring(M) |> sf::st_sfc() |> sf::st_set_crs(sf::st_crs(seed))
      }

      #plot(newline)
      #plot(existing_network, add = T)
    }
  }

  if (length(list_intersections) >= 1)
  {
    intersections <- sf::st_sfc(list_intersections) |>
      sf::st_set_crs(sf::st_crs(seed))

    p = lwgeom::st_startpoint(intersections)
    d = as.numeric(sf::st_distance(p,newline))
    intersections <- intersections[d < 1]
  }
  else
  {
    intersections <- NULL
  }


  nintersection <- length(intersections)
  len = as.numeric(sf::st_length(newline))
  dintersection = nintersection/(len/1000)

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

  return(list(road = newline, seeds = intersections, dintersection = dintersection))
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

query_conductivity_aoi <- function(conductivity, seed, smooth = TRUE)
{
  resolution <- terra::res(conductivity)[1]

  buf_dist <- c("ahead" = 900, "behind" = 100, "side" = 400)  # Could be set as a parameter

  aoi <- seed |>
    st_extend_line(c(buf_dist["ahead"], buf_dist["behind"])) |>
    sf::st_buffer(buf_dist["side"], endCapStyle = "FLAT")
  #plot(aoi, add = T, col= "red")
  #plot(seed, add = T, lwd = 3)

  conductivity_crop <- terra::crop(conductivity, terra::vect(aoi))
  #plot(conductivity_crop, col = viridis::inferno(50))
  #plot(seed, add = T, lwd = 3, col = "red")

  if (smooth)
  {
    zero = conductivity_crop == 0
    conductivity_crop[zero] = NA
    kernel = matrix(1, 3, 3)
    conductivity_crop <- terra::focal(conductivity_crop, kernel, "mean", na.rm = TRUE)
    conductivity_crop[zero] = 0
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

  if (j[1] > i[1])
    j = c(1,j)

  if (j[length(j)] < i[length(i)])
    j = c(j, length(x))

  #plot(-x, type = "b")
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

  ans = data.frame(idx = i, depth = depths)
  return(ans)
}

find_reachable <- function(start, ends, trans, cost_max)
{
  if (length(sf::st_geometry(ends)) < 6) return(NULL)

  cost <- gdistance::costDistance(trans, sf::st_coordinates(start), sf::st_coordinates(ends))
  cost <- as.numeric(cost)

  if (all(is.infinite(cost))) return(NULL)

  cost[is.infinite(cost)] <- max(9999, max(cost[!is.infinite(cost)]))

  # Select local minima of cost
  smooth <- 9
  scost = -ma(cost, n = smooth) #spike_preserving_smooth(cost, n = smooth)
  scost = ma(scost, n = 5)
  offset = floor(smooth/2) + floor(5/2)
  scost = as.numeric(na.omit(scost))
  minima = pracma::findpeaks(scost, zero = "+")

  if (is.null(minima)) return(NULL)
  if (nrow(minima) == 0) return(NULL)

  depth = minima[,1]
  center = minima[,2] + offset
  left = minima[,3] + offset
  right = minima[,4] + offset
  cost2 = depth

  #plot(seq_along(cost), cost, type = "l", col = "red", log = "y")
  #lines(-scost, type = "l", col = "blue")
  #abline(v = center)
  #abline(v = left, col = "red")
  #abline(v = right, col = "blue")

  for (i in seq_along(center))
  {
    lhs = cost[left[i]]
    rhs = cost[right[i]]
    center[i] = which.min(cost[left[i]:right[i]]) + left[i] - 1
  }

  for (i in seq_along(center))
  {
    lhs = cost[left[i]]
    rhs = cost[right[i]]
    cost2[i] = cost[center[i]]
    depth[i] = min(lhs, rhs) - cost2[i]
  }

  minima = data.frame(idx = center, depth = depth, cost = cost2)

  depth_thresolds = c(20,500)
  factor_thresold = c(1,1.5)
  slope = diff(factor_thresold)/diff(depth_thresolds)
  intercept = factor_thresold[1] - slope * depth_thresolds[1]
  cost_max_multiplier_for_depth = minima$depth*slope + intercept
  cost_max_multiplier_for_depth[cost_max_multiplier_for_depth > 1.5] <- 1.5

  keep = (minima$cost <= cost_max_multiplier_for_depth * cost_max)
  if (sum(keep) > 0)
    minima = minima[keep,]
  else
    minima = minima[which.min(minima$cost),]

  data.table::setorder(minima, cost)

  return(list(minima = minima, cost = cost))
}

rotate_pi <- function(X, C)
{
  X[,1] = -(X[,1] - C[1]) + C[1]
  X[,2] = -(X[,2] - C[2]) + C[2]
  X
}

mask_existing_network <- function(x, existing_network)
{
  if (!is.null(existing_network) && length(existing_network) > 0)
  {
    bb <- sf::st_bbox(x)
    bb <- sf::st_set_crs(bb, sf::st_crs(existing_network))
    mask <- sf::st_crop(existing_network, bb)
    if (length(mask) > 0)
    {
      mask <- sf::st_buffer(mask, dist = 10)
      x <- terra::mask(x, terra::vect(mask), inverse = TRUE, updatevalue = 0)
    }
  }

  return(x)
}

spike_preserving_smooth <- function(x, n = 5)
{
  y = ma(x, n)
  dp <- x-y
  th = -10
  dp[is.na(dp)] = 0
  y[dp < th] = x[dp < th]
  #plot(x, type = "l", ylim = c(0, max(x) *1.1))
  #lines(y, col = "red")
  y
}

mask_passage <- function(raster, lines, from, to, crs)
{
  if (is(lines, "sfg")) lines <- sf::st_sfc(lines)
  mask <- lwgeom::st_linesubstring(lines, from, to, 0.95)
  mask <- sf::st_buffer(mask, dist = 5, endCapStyle = "FLAT")
  mask <- sf::st_set_crs(mask, crs)
  raster <- terra::mask(raster, terra::vect(mask), inverse = TRUE, updatevalue = 0)
  return(raster)
}

st_join_linestring = function(x,y)
{
  X = sf::st_coordinates(x)
  Y = sf::st_coordinates(y)
  M = rbind(X,Y)[,1:2]
  L = sf::st_sfc(sf::st_linestring(M))
  return(L)
}
