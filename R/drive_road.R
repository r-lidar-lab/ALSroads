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


#' Drive along an unknowned road
#'
#' Attempt to drive along an unknowned road in a conductivity raster starting
#' only with a small piece of road segment pointing in the right direction.
#'
#' @param start_line  line (\code{sf} format)
#' @param conductivity  line (\code{stars} format)
#' @param radius  numeric (distance unit). Search radius used to find the next most probable point on the road.
#' @param cost_max  numeric. Maximal cost allowed in the conductivity raster for point candidate to continue the search.
#'
#' @return list, \code{line} being the road found (as \code{sfc}) and \code{cost} a numeric vector of cost at each invidual vertices of the line.
#' @export
#' @examples
#' library(stars)
#' library(sf)
#'
#' conductivity <- read_stars("conductivity.tif")
#' start_line <- st_read("start_line.gpkg")
#' 
#' res <- drive_road(start_line, conductivity, radius = 10, cost_max = 50)
#' 
#' raster::plot(as(conductivity, "Raster"))
#' raster::plot(res$line, add = TRUE)
#' boxplot(res$cost)
drive_road <- function(start_line, conductivity, radius = 10, cost_max = 50)
{
  # Transiton matrix
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
  cat("Steps: 1")
  while(still_good)
  {
    cat("..", k, sep = "")

    # Compute heading from the two previous points
    p1 <- list_coords[[k-1]][1:2]
    p2 <- list_coords[[k]][1:2]
    heading <- get_heading(p1, p2)


    # Create search zone which is only ahead of the previous point
    pt <- sf::st_sfc(sf::st_point(p2), crs = crs)
    fov_cap <- make_fov(pt, heading, radius)

    ext <- sf::st_bbox(fov_cap)
    ext <- c(floor(ext[1:2]), ceiling(ext[3:4])) # only works because raster resolution is 1 m !!
    pts_ext <- sf::st_make_grid(ext,
                                cellsize = resolution,
                                offset = ext[1:2] - resolution/2,
                                crs = crs,
                                what = "centers")

    pts_cap <- pts_ext |>
      sf::st_intersection(fov_cap) |>
      stars::st_extract(x = conductivity)


    # Find all points at edge of search radius
    is_edge_cells <- as.numeric(sf::st_distance(pt, pts_cap)) > (radius - resolution)
    pts_edge <- pts_cap[is_edge_cells,]


    # Find cost and heading of each edge point
    pt_sp <- sf::as_Spatial(pt)
    pts_edge_sp <- sf::as_Spatial(sf::st_geometry(pts_edge))
    cost <- gdistance::costDistance(trans, pt_sp, pts_edge_sp)
    headings <- sapply(sf::st_geometry(pts_edge), function(x) get_heading(p2, sf::st_coordinates(x)))


    # Find weighted least cost edge point
    # Weights are 80% on cost and 20% on heading differences
    # the idea being that the next point on the road should
    # be on a similar heading
    tb_cost <- dplyr::tibble(idx = 1:nrow(pts_edge),
                                diff = abs(abs(headings) - abs(heading)),
                                cost = cost[1,],
                                order_diff = 0,
                                order_cost = 0,
                                order_sum = 0)

    tb_cost <- tb_cost |>
      dplyr::arrange(diff) |>
      dplyr::mutate(order_diff = 1:nrow(pts_edge)) |>
      dplyr::arrange(cost) |>
      dplyr::mutate(order_cost = 1:nrow(pts_edge)) |>
      dplyr::mutate(order_sum = order_diff * 0.2 + order_cost * 0.8) |>
      dplyr::arrange(idx)

    idx_min <- which.min(tb_cost[["order_sum"]])
    cost_selected <- tb_cost[idx_min,][["cost"]]


    # Check if maximal cost criterea is met
    # before saving the least-cost pixel coordinates
    # or stopping the drive
    if(cost_selected > cost_max)
    {
      still_good <- FALSE
    } else {
      k <- k + 1
      list_coords[[k]] <- c(sf::st_coordinates(pts_edge[idx_min,]),
                            tb_cost[idx_min,][["cost"]])
    }
  }

  m_coords <- do.call(rbind, list_coords)
  newline <- m_coords[,1:2] |>
    sf::st_linestring() |>
    sf::st_sfc(crs = crs)

  return(list(line = newline,
              cost = m_coords[,3]))
}