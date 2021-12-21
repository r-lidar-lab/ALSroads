#' Check amplitude of differences between corrected and uncorrected roads
#'
#' Check amplitude of differences between corrected and uncorrected roads
#'
#' @param roads  multiples lines (\code{sf} format). Corrected but unconnected roads.
#' @param roads_ori  multiples lines (\code{sf} format). Original non-corrected but connected roads.
#' @param field  character. Unique identifier field in both road datasets.
#'
#' @return data.frame with metrics about each roads.
#' @export
check_road_difference <- function(roads, roads_ori, field)
{
  IDs1 <- sort(unique(roads[[field]]))
  IDs2 <- sort(unique(roads_ori[[field]]))
  if (length(IDs1) != length(roads[[field]])) stop("Values in unique identifier field are not unique for 'roads'.", call. = FALSE)
  if (length(IDs2) != length(roads_ori[[field]])) stop("Values in unique identifier field are not unique for 'roads_ori'.", call. = FALSE)
  if (!all(IDs1 == IDs2)) stop("Values in unique identifier field are not the same in both road datasets.", call. = FALSE)

  # Arrange both road datasets to make sure
  # that line indices will match between them
  roads <- dplyr::arrange(roads, field)
  roads_ori <- dplyr::arrange(roads_ori, field)

  # Compute metrics
  ratios_area_perimeter <- mapply(diff_area_perimeter, sf::st_geometry(roads_ori), sf::st_geometry(roads))
  quantiles_along_road <- mapply(diff_along_road, sf::st_geometry(roads_ori), sf::st_geometry(roads), SIMPLIFY = FALSE) |> do.call(what = rbind)

  # Format results
  df_results <- data.frame(
    field = roads_ori[[field]],
    area_over_perimeter = ratios_area_perimeter) |>
    cbind(quantiles_along_road)
  
  names(df_results)[1] <- field

  return(df_results)
}


#' Compute a difference index between two roads by considering them as a polygon 
#'
#' Connect the two roads by their ends to construct a polygon from which
#' its area and perimeter will be computed. A ratio area/perimeter calculated.
#' The larger the value, the larger the differences between the two roads. It
#' must be noted that in some extreme cases, a low value doesn't mean low differences.
#'
#' @param road_ori  line (\code{sf} format). Original non-corrected road.
#' @param road_cor  line (\code{sf} format). Corrected road.
#' @param graph  boolean. Whether of not to display graphics.
#'
#' @return numeric. Ratio of the area over the perimeter of the constructed polygon.
#' @noRd
diff_area_perimeter <- function(road_ori, road_cor, graph = FALSE)
{
  if (nrow(road_ori) > 1 | nrow(road_cor) > 1) stop("'road_ori' and 'road_cor' must contain only one feature.", call. = FALSE)

  # Extract vertices
  coords_ori <- sf::st_coordinates(road_ori)[,-3]
  coords_cor <- sf::st_coordinates(road_cor)[,-3]
  
  
  # Adjust coordinates in order to create a valid sequence
  # of vertices for a polygon
  dist_to_start <- dist(rbind(coords_ori[1,], coords_cor[1,]))[1]
  dist_to_end <- dist(rbind(coords_ori[1,], coords_cor[nrow(coords_cor),]))[1]
  closest_vertex <- which.min(c(dist_to_start, dist_to_end))
  
  if (closest_vertex == 1) coords_cor <- apply(coords_cor, 2, rev)
  
  
  # Construct polygon
  poly <- rbind(coords_ori, coords_cor, coords_ori[1,]) |>
    list() |>
    sf::st_polygon() |>
    sf::st_sfc(crs = sf::st_crs(road_ori))
  
  
  # Compute area/perimeter ratio of polygon
  ratio <- as.numeric(sf::st_area(poly) / lwgeom::st_perimeter(poly))
  
  
  # Make graphical representation of the differences
  if (graph)
  {
    coords <- rbind(coords_ori, coords_cor)
    limits <- list(x = range(coords[,1]), y = range(coords[,2]))
    bigger <- as.numeric(diff(limits$x) < diff(limits$y)) + 1
    offset <- mean(limits[[bigger]]) - limits[[bigger]][1]
    limits$x <- mean(limits$x) + c(-offset, offset)
    limits$y <- mean(limits$y) + c(-offset, offset)

    plot_poly <- ggplot2::ggplot() +
      ggplot2::geom_sf(data = poly, fill = "orange") +
      ggplot2::geom_sf(data = road_ori, size = 1, color = "red") +
      ggplot2::geom_sf(data = road_cor, size = 1, color = "darkgreen") +
      ggplot2::coord_sf(xlim = limits$x, ylim = limits$y, datum = sf::st_crs(road_ori)) +
      ggplot2::ggtitle(sprintf("Ratio area/perimeter: %.1f m²/m", ratio))
    
    print(plot_poly)
  }
  return(ratio)
}


#' Compute differences quantiles based on are far apart are the roads
#'
#' Sample points at regular interval each road and measure distances between
#' each pair of points. The larger the quantiles are, the larger the differences
#' between the two roads.
#'
#' @param road_ori  line (\code{sf} format). Original non-corrected road.
#' @param road_cor  line (\code{sf} format). Corrected road.
#' @param step  numeric (distance unit). Interval on \code{road_ori} at which sample differences.
#' @param graph  boolean. Whether of not to display graphics.
#'
#' @return named numeric. Quantiles are P50, P90, P100.
#' @noRd
diff_along_road <- function(road_ori, road_cor, step = 10, graph = FALSE)
{
  if (nrow(road_ori) > 1 | nrow(road_cor) > 1) stop("'road_ori' and 'road_cor' must contain only one feature.", call. = FALSE)

  road_ori <- sf::st_geometry(road_ori)
  road_cor <- sf::st_geometry(road_cor)
  
  
  st_full_line_sample <- function(line, n_steps)
  {
    coords <- sf::st_coordinates(line)[,-3]
    pts_extrems <- c(1, nrow(coords)) |>
      lapply(function(x) { sf::st_point(coords[x,]) }) |>
      sf::st_sfc(crs = sf::st_crs(road_ori))
    
    pts_middle <- sf::st_cast(sf::st_line_sample(line, n_steps), "POINT")
    pts <- c(pts_extrems[1], pts_middle, pts_extrems[2])
  }
  
  
  # Sample points along the two roads
  n_steps <- ceiling(sf::st_length(road_ori) / step) |> as.numeric()
  ls_pts <- lapply(c(road_ori, road_cor), st_full_line_sample, n_steps)
  names(ls_pts) <- c("road_ori", "road_cor")
  
  
  # Compute differences at P50, P90, P100
  dist_step <- sf::st_distance(ls_pts[["road_ori"]], ls_pts[["road_cor"]]) |>
    diag() |>
    as.numeric()
  
  p <- quantile(dist_step, probs = c(0.5, 0.90, 1))
  
  
  # Display graphical representation of the differences
  if (graph)
  {
    coords_ori <- sf::st_coordinates(ls_pts[["road_ori"]])
    coords_cor <- sf::st_coordinates(ls_pts[["road_cor"]])
    
    lines <- 1:nrow(coords_ori) |>
      lapply(function(x) { sf::st_linestring(rbind(coords_ori[x,], coords_cor[x,])) }) |>
      sf::st_sfc(crs = sf::st_crs(road_ori))
    
    df_dist <- data.frame(idx = 1:length(dist_step),
                          distance = dist_step)
    
    coords <- rbind(coords_ori, coords_cor)
    limits <- list(x = range(coords[,1]), y = range(coords[,2]))
    bigger <- as.numeric(diff(limits$x) < diff(limits$y)) + 1
    offset <- mean(limits[[bigger]]) - limits[[bigger]][1]
    limits$x <- mean(limits$x) + c(-offset, offset)
    limits$y <- mean(limits$y) + c(-offset, offset)

    plot_lines <- ggplot2::ggplot() +
      ggplot2::geom_sf(data = road_ori, size = 1, color = "red") +
      ggplot2::geom_sf(data = road_cor, size = 1, color = "darkgreen") +
      ggplot2::geom_sf(data = lines, size = 0.5, color = "purple") +
      ggplot2::coord_sf(xlim = limits$x, ylim = limits$y, datum = sf::st_crs(road_ori)) +
      ggplot2::ggtitle("Pairwise points along roads")

    plot_dots <- ggplot2::ggplot(df_dist, ggplot2::aes(x=idx, y=distance)) +
      ggplot2::geom_point() +
      ggplot2::geom_smooth(method = "loess", formula = "y ~ x", se = FALSE) +
      ggplot2::geom_hline(yintercept = p) +
      ggplot2::ggtitle("Pairwise distances along roads",
                       subtitle = sprintf("P50: %.1f m | P90: %.1f m | P100: %.1f m", p[1], p[2], p[3]))
    
    print(plot_lines)
    print(plot_dots)
  }
  return(p)
}
