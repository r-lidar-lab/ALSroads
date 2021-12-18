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
#' @param fov  numeric. Field of view (degrees) ahead of the search vector.
#' @param radius  numeric (distance unit). Search radius used to find the next most probable point on the road.
#' @param cost_max  numeric. Maximal cost allowed in the conductivity raster for point candidate to continue the search.
#' If no value is provided \code{radius * 6} will be.
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
#' res <- drive_road(start_line, conductivity)
#' 
#' raster::plot(conductivity, col = viridis::viridis(50))
#' raster::plot(res$line, add = TRUE, col = "red", lwd = 2)
#' raster::plot(start_line, add = TRUE, col = "green", lwd = 3)
#' plot(res$cost, type = "l")
drive_road <- function(start_line, conductivity, fov = 45, radius = 10, cost_max = NULL)
{
  if (is.null(cost_max)) cost_max <- radius * 6
  
  resolution <- terra::res(conductivity)
  if (resolution[1] != resolution[2]) stop("'conductivity' raster must have the same resolution in both X and Y axis.", call. = FALSE)
  resolution <- resolution[1]

  resolution_min <- 2  # Could be set as a parameter
  if (resolution < resolution_min)
  {
    agg_factor <- ceiling(resolution_min / resolution)
    resolution <- resolution * agg_factor
    cat("Step 0: downsample conductivity to", resolution, "m\n")
    conductivity <- terra::aggregate(conductivity, fact = agg_factor, fun = mean, na.rm = TRUE)
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
    X <- p2[1] + radius * cos(heading + angles_rad)
    Y <- p2[2] + radius * sin(heading + angles_rad)

    M <- data.frame(X,Y)
    ends <- sf::st_as_sf(M, coords = c("X", "Y")) |>
      sf::as_Spatial()
    
    # Find cost of each edge point
    start <- sf::st_point(p2) |>
      sf::st_sfc() |>
      sf::as_Spatial()
    cost <- as.numeric(gdistance::costDistance(trans, start, ends))
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


#' Generate and save to file tile of conductivity
#'
#' Generate and save to file tile of conductivity based on the extent
#' of an area of interest.
#'
#' @param bbox  named numeric vector. Bounding box vector with names \code{xmin}, \code{xmax}, \code{ymin}, \code{ymax}
#' @param dtm  raster (\code{raster} format). Digital terrain model covering \code{bbox}.
#' @param ctg  \code{LAScatalog} covering \code{bbox}.
#' @param buffer  numeric (distance unit). Buffer to be added outside of \code{bbox} when computing metrics. Won't affect extent of the produced tile.
#' @param outdir  character. Directory in which the produced tile will be saved.
#'
#' @return character representing filenames of generated tiles.
#' @export
#' @examples
#' library(sf)
#' library(raster)
#' library(lidR)
#' library(future)
#' library(future.apply)
#'
#' rootdir = "Y:/Developpement/Chemins_forestiers/drive_road"
#' setwd(rootdir)
#' 
#' outdir = getwd()
#' size = 500
#' buffer = 5
#' aoi = sf::st_read("aoi.gpkg")
#' ctg = readLAScatalog("LAZ")
#' dtm = raster("dtm.tif")
#' dtm = crs(aoi)
#' 
#' 
#' # Grid that will be used to make tiles
#' bbox <- st_bbox(aoi)
#' bbox[c("xmin","ymin")] <- floor(bbox[c("xmin","ymin")])
#' bbox[c("xmax","ymax")] <- ceiling(bbox[c("xmax","ymax")])
#' grid <- st_make_grid(bbox, size)
#' 
#' # List of bboxes for parallelisation
#' bboxes <- lapply(grid, function(x) {
#'   coords <- st_coordinates(x)
#'   bbox <- c(range(coords[,"X"]), range(coords[,"Y"]))
#'   names(bbox) <- c("xmin","xmax","ymin","ymax")
#'   bbox})
#' 
#' # Generate tiles
#' future::plan(multisession)
#' filenames <- future.apply::future_sapply(bboxes, tile_conductivity, dtm, ctg, buffer, outdir)
#' future::plan(sequential)
#' 
#' # Create VRT from tiles
#' vrtfile <- file.path(outdir, "conductivity.vrt")
#' terra::vrt(file.path(outdir, filenames), vrtfile)
tile_conductivity <- function(bbox, dtm, ctg, buffer, outdir)
  {
    # Buffer added to bbox to ensure valid calculations at edges and thus
    # a smooth transition between tiles
    res <- raster::xres(dtm)
    bbox_buf <- c(bbox["xmin"], bbox["xmax"], bbox["ymin"], bbox["ymax"]) + c(-buffer*res, buffer*res, -buffer*res, buffer*res)

    # Crop DTM and point cloud
    cropped <- raster::crop(dtm, raster::extent(bbox_buf))
    las <- lidR::clip_rectangle(ctg, bbox_buf["xmin"], bbox_buf["ymin"], bbox_buf["xmax"], bbox_buf["ymax"])
    
    # Compute and write to file final conductivity layer
    conductivities <- grid_conductivity(las, cropped, road = NULL, water = NULL)
    conductivity <- raster::crop(conductivities$conductivity, raster::extent(bbox))

    filename <- paste0("conductivity_", bbox["xmin"], "_", bbox["ymin"], ".tif")
    raster::writeRaster(conductivity, file.path(outdir, filename))

    return(filename)
  }