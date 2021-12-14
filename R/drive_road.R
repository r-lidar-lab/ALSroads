# Very similar to make_caps()
# "fov" like in "field-of-view
make_fov <- function(pt, heading, radius)
{
  M1 <- sf::st_coordinates(pt)
  M2 <- M1 + radius * c(cos(heading), sin(heading))
  
  line <- rbind(M1, M2) |>
    sf::st_linestring() |>
    sf::st_sfc(crs = sf::st_crs(pt))
  
  rectangle <- sf::st_buffer(line, radius, endCapStyle = "FLAT")
  circle <- sf::st_buffer(pt, radius)
  
  fov <- sf::st_intersection(circle, rectangle)
  
  return(fov)
}


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
#' @param cost_max  numeric. Maximal cost allowed in the conductivity raster for point candidate to continue the search.
#' @param fov  numeric. Field-of-view (degrees) of the search vector.
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
#' res <- drive_road(start_line, conductivity, radius = 10, cost_max = 50)
#' 
#' raster::plot(conductivity, col = viridis::viridis(50))
#' plot(start_line, add = T, col = "red")
#' raster::plot(res$line, add = TRUE, col = "red", lwd = 2)
#' plot(res$cost, type = "l")
drive_road <- function(start_line, conductivity, radius = 10, cost_max = 50, fov = 45)
{
  resolution <- terra::res(conductivity)[1]

  if (resolution < 2)
  {
    cat("Steps 0: downsample conductivity\n")
    conductivity <- terra::aggregate(conductivity, fact = 2, fun = mean, na.rm = TRUE)
  }

  resolution <- terra::res(conductivity)[1]
  angular_resolution = floor((180*resolution)/(pi*radius))
  angles = seq(-fov, fov, angular_resolution) * pi / 180

  # Transiton matrix
  cat("Steps 1: computation of the transition matrix\n")

  trans <- conductivity |>
    as("Raster") |>
    gdistance::transition(transitionFunction = mean, directions = 8) |>
    gdistance::geoCorrection()


  # Initialisation of list of coordinates with starting line
  start_coords <- sf::st_coordinates(start_line)
  start_coords[,3] <- 0
  list_coords <- list(start_coords[1,], start_coords[2,])

  # Loop initialisation
  resolution <- stars::st_dimensions(conductivity)$x$delta
  crs <- sf::st_crs(start_line)
  still_good <- TRUE
  k <- 2
  cat("Steps 2: driving the conductivity raster\n")
  while(still_good)
  {
    if (k %% 5 == 0) cat(" | ", (k-1)*10, " m", sep = "")

    # Compute heading from the two previous points
    p1 <- list_coords[[k-1]]
    p2 <- list_coords[[k]]
    heading <- get_heading(p1, p2)

    # Create search zone which is only ahead of the previous point

    X <- p2[1] + radius * cos(heading+angles)
    Y <- p2[2] + radius * sin(heading+angles)

    # Find cost and heading of each edge point
    start <- sf::st_sfc(sf::st_point(p2))
    start <- sf::as_Spatial(start)

    M <- data.frame(X,Y)
    ends <- sf::st_as_sf(M, coords = c("X", "Y"))
    ends <- sf::as_Spatial(ends)
    cost <- as.numeric(gdistance::costDistance(trans, start, ends))
    cost[is.infinite(cost)] <- 2*max(cost[!is.infinite(cost)])


    # Penalisation de 50% du cout a 45 degres
    penalty_at_45_degrees = 1.5
    penalty =  (abs(angles)*180/pi)/fov
    penalty = penalty * (penalty_at_45_degrees - 1) + 1
    cost2 = cost * penalty_at_45_degrees     ### ICI UNE ERREUR, DEVRAIT ÊTRE   cost2 = cost * penalty
    cost2 = ma(cost2)
    #plot(angles*180/pi, cost2, type = "l")
    #lines(angles*180/pi, cost, col = "red")

    idx_min <- which.min(cost2)
    cost_selected <- cost[idx_min]
    W <- as.matrix(M[idx_min,])
    W <- cbind(W, cost_selected)

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
  newline <- sf::st_linestring(m_coords[,1:2])
  newline <- sf::st_sfc(newline, crs = crs)

  return(list(line = newline, cost = m_coords[,3]))
}