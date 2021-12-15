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


# Smooth adjacent costs with a mean moving windows
ma <- function(x, n = 3)
{
  if (anyNA(x)) x <- zoo::na.approx(x, na.rm = FALSE)
  as.numeric(stats::filter(x, rep(1 / n, n), sides = 2))
}


#' Drive along an unknowned road
#'
#' Attempt to drive along an unknowned road in a conductivity raster starting
#' only with a small piece of road segment pointing in the right direction.
#'
#' @param start_line  line (\code{sf} format)
#' @param conductivity  raster (\code{terra} format)
#' @param radius  numeric (distance unit). Search radius used to find the next most probable point on the road.
#' @param fov  numeric. Field of view (degrees) ahead of the search vector.
#' @param cost_max  numeric. Maximal cost allowed in the conductivity raster for point candidate to continue the search.
#'
#' @return list, \code{line} being the road found (as \code{sfc}) and \code{cost} a numeric vector of cost at each invidual vertices of the line.
#' @export
#' @examples
#' library(terra)
#' library(sf)
#'
#' conductivity <- rast("conductivity.tif")
#' start_line <- st_read("start_line.gpkg")
#' 
#' res <- drive_road(start_line, conductivity, cost_max = 60)
#' 
#' raster::plot(conductivity, col = viridis::viridis(50))
#' plot(start_line, add = T, col = "red")
#' raster::plot(res$line, add = TRUE, col = "red", lwd = 2)
#' plot(res$cost, type = "l")
drive_road <- function(start_line, conductivity, radius = 10, fov = 45, cost_max = 50)
{
  # TODO Defaut value of "cost_max" might need to be something else as it is function
  #      of the "radius" parameter and the final resolution of "conductivity"
  resolution <- terra::res(conductivity)[1]

  if (resolution < 2)
  {
    cat("Step 0: downsample conductivity\n")
    agg_factor <- ceiling(2 / resolution)
    conductivity <- terra::aggregate(conductivity, fact = agg_factor, fun = mean, na.rm = TRUE)
    resolution <- terra::res(conductivity)[1]
  }

  angular_resolution = floor((180*resolution)/(pi*radius))
  angles = seq(-fov, fov, angular_resolution) * pi / 180

  # Transition matrix
  cat("Step 1: computation of the transition matrix\n")
  trans <- conductivity |>
    raster::raster() |>
    gdistance::transition(transitionFunction = mean, directions = 8) |>
    gdistance::geoCorrection()

  # Initialisation of list of coordinates with starting line
  start_coords <- sf::st_coordinates(start_line)
  start_coords[,3] <- 0
  list_coords <- list(start_coords[1,], start_coords[2,])

  # Loop initialisation
  still_good <- TRUE
  k <- 2
  cat("Step 2: driving the conductivity raster\n")
  while(still_good)
  {
    if (k %% 5 == 0) cat(" | ", (k-1)*10, " m", sep = "")

    # Compute heading from the two previous points
    p1 <- list_coords[[k-1]]
    p2 <- list_coords[[k]]
    heading <- get_heading(p1, p2)

    # Create all possible ends ahead of the previous point
    X <- p2[1] + radius * cos(heading + angles)
    Y <- p2[2] + radius * sin(heading + angles)

    M <- data.frame(X,Y)
    ends <- sf::st_as_sf(M, coords = c("X", "Y")) |>
      sf::as_Spatial()
    
    # Find cost of each edge point
    start <- sf::st_point(p2) |>
      sf::st_sfc() |>
      sf::as_Spatial()
    cost <- as.numeric(gdistance::costDistance(trans, start, ends))
    cost[is.infinite(cost)] <- 2 * max(cost[!is.infinite(cost)])

    # Add cost penalties based on heading
    # Up to 50% on heading of 45 degrees
    penalty_at_45_degrees <- 1.5
    penalty <- (abs(angles) * 180/pi) / fov
    penalty <- penalty * (penalty_at_45_degrees - 1) + 1
    cost2 <- ma(cost * penalty)

    # Select point coordinates with the lowest
    # penalty adjusted cost
    idx_min <- which.min(cost2)
    cost_selected <- cost[idx_min]
    W <- M[idx_min,] |>
      as.matrix() |>
      cbind(cost_selected)

    # Check if maximal cost criterea is met
    # before saving the least-cost pixel coordinates
    # or stopping the drive
    if(cost_selected > cost_max)
    {
      still_good <- FALSE
    }
    else
    {
      k <- k + 1
      list_coords[[k]] <- W
    }
  }

  m_coords <- do.call(rbind, list_coords)
  newline <- m_coords[,1:2] |>
    sf::st_linestring() |>
    sf::st_sfc(crs = sf::st_crs(start_line))

  return(list(line = newline, cost = m_coords[,3]))
}