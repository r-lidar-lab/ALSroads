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
#' @param conductivity  raster (\code{terra} format)
#' @param fov  numeric. Field of view (degrees) ahead of the search vector.
#' @param radius  numeric (distance unit). Search radius used to find the next most probable point on the road.
#' @param cost_max  numeric. Maximal cost allowed in the conductivity raster for point candidate to continue the search.
#' If no value is provided \code{radius * 6} will be used.
#'
#' @return list, \code{line} being the road found (as \code{sfc}) and \code{cost} a numeric vector
#' of cost at each invidual vertices of the line.
#' @export
#' @examples
#' library(terra)
#' library(sf)
#'
#' conductivity <- system.file("extdata", "drived_conductivity.tif", package = "MFFProads")
#' drived_road <- system.file("extdata", "drived_road.gpkg", package = "MFFProads")
#' conductivity <- rast(conductivity)
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


  resolution <- terra::res(conductivity)
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
    terra::vect()

  conductivity_crop <- terra::crop(conductivity, aoi)

  resolution_min <- 2  # Could be set as a parameter
  if (resolution < resolution_min)
  {
    agg_factor <- ceiling(resolution_min / resolution)
    resolution <- resolution * agg_factor
    cat("Step 0: downsample conductivity to", resolution, "m\n")
    conductivity_crop <- terra::aggregate(conductivity_crop, fact = agg_factor, fun = mean, na.rm = TRUE)
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
    raster::raster() |>
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
    val <- terra::extract(conductivity_crop, M)
    if (any(is.nan(val[,2])))
    {
      # Make a newly cropped conductivity raster further ahead
      line <- c(p1[1:2], p2[1:2]) |>
        matrix(ncol = 2, byrow = TRUE) |>
        sf::st_linestring() |>
        sf::st_sfc(crs = sf::st_crs(starting_road))

      aoi <- line |>
        st_extend_line(c(buf_dist["ahead"], buf_dist["behind"])) |>
        sf::st_buffer(buf_dist["side"], endCapStyle = "FLAT") |>
        terra::vect()

      conductivity_crop <- terra::crop(conductivity, aoi)

      resolution <- terra::res(conductivity)[1]
      resolution_min <- 2  # Could be set as a parameter
      if (resolution < resolution_min)
      {
        agg_factor <- ceiling(resolution_min / resolution)
        resolution <- resolution * agg_factor
        conductivity_crop <- terra::aggregate(conductivity_crop, fact = agg_factor, fun = mean, na.rm = TRUE)
      }

      # Check again if some ends fall outside of the newly cropped conductivity raster
      val <- terra::extract(conductivity_crop, M)
      if (any(is.nan(val[,2])))
      {
        warning("Drive stopped early. Edge of conductivity raster has been reached.", call. = FALSE)
        cost_max <- -Inf
      }

      # Generate a new transition matrix
      trans <- conductivity_crop |>
        raster::raster() |>
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


#' Generate and save to file tile of conductivity
#'
#' Generate and save to file tile of conductivity based on the extent
#' of an area of interest.
#'
#' @param bbox  named numeric vector. Bounding box vector with names \code{xmin}, \code{xmax}, \code{ymin}, \code{ymax}
#' @param dtm  raster (\code{terra} format) or path to raster file. Digital terrain model covering \code{bbox}.
#' @param ctg  \code{LAScatalog} covering \code{bbox}.
#' @param buffer  numeric (pixel unit). Buffer to be added outside of \code{bbox} when computing metrics. Won't affect extent of the produced tile.
#' @param outdir  character. Directory in which the produced tile will be saved.
#'
#' @return character representing filenames of generated tiles.
#' @export
#' @examples
#' library(sf)
#' library(terra)
#' library(lidR)
#'
#' dir <- system.file("extdata", "", package = "MFFProads")
#' path_dtm <- system.file("extdata", "dtm_1m.tif", package = "MFFProads")
#'
#' ctg <- readLAScatalog(dir)
#' dtm <- rast(path_dtm)
#' outdir <- getwd()
#'
#' # Define parameters for grid of tiles
#' # Internal buffer is only to mess things up a bit
#' aoi <- st_bbox(ctg) |>
#'   st_as_sfc() |>
#'   st_buffer(-14.7)
#' size <- 500  # (in distance unit)
#'
#' # Generate grid of tiles
#' bbox <- st_bbox(aoi)
#' bbox[c("xmin","ymin")] <- floor(bbox[c("xmin","ymin")])
#'
#' xlength <- bbox["xmax"] - bbox["xmin"]
#' bbox["xmax"] <- bbox["xmin"] + ceiling(xlength / xres(dtm)) * xres(dtm)
#' ylength <- bbox["ymax"] - bbox["ymin"]
#' bbox["ymax"] <- bbox["ymin"] + ceiling(ylength / yres(dtm)) * yres(dtm)
#'
#' grid <- st_make_grid(bbox, size)
#'
#' # Generate list of bounding boxes
#' bboxes <- lapply(grid, function(x)
#' {
#'   coords <- st_coordinates(x)
#'   bbox <- c(range(coords[,"X"]), range(coords[,"Y"]))
#'   names(bbox) <- c("xmin","xmax","ymin","ymax")
#'   bbox
#'  })
#'
#' # Generate tiles
#' buffer <- 5
#' filenames <- sapply(bboxes, tile_conductivity, path_dtm, ctg, outdir, buffer)
#'
#' # Create VRT from tiles for future use (allow a similar workflow as using a LAScatalog for LAS files)
#' path_filenames <- file.path(outdir, filenames)
#' path_vrtfile <- file.path(outdir, "conductivity.vrt")
#' vrt(path_filenames, path_vrtfile)
tile_conductivity <- function(bbox, dtm, ctg, outdir, buffer = 0, param = mffproads_default_parameters)
{
  # Check if dtm is a path or a SpatRaster
  if (class(dtm)[1] != "SpatRaster")
  {
    if (is.character(dtm))
      {
        dtm <- terra::rast(dtm)
      } else {
        stop("Invalid 'dtm'. Must be either object 'SpatRaster' or path to raster file", call. = FALSE)
      }
  }

  # Buffer added to bounding box to ensure valid calculations at edges and
  # thus a smooth transition between tiles
  resolution <- terra::xres(dtm)
  bbox_buf <- c(bbox["xmin"], bbox["xmax"], bbox["ymin"], bbox["ymax"]) + resolution * c(-buffer, buffer, -buffer, buffer)

  # Crop DTM and point cloud
  cropped <- terra::crop(dtm, bbox_buf) |> raster::raster()
  las <- lidR::clip_rectangle(ctg, bbox_buf["xmin"], bbox_buf["ymin"], bbox_buf["xmax"], bbox_buf["ymax"])

  # Compute and write to file final conductivity layer
  conductivities <- grid_conductivity(las, centerline = NULL, dtm = cropped, param = param)
  conductivity <- terra::crop(terra::rast(conductivities$conductivity), terra::ext(bbox))

  filename <- paste0("conductivity_", bbox["xmin"], "_", bbox["ymin"], ".tif")
  terra::writeRaster(conductivity, file.path(outdir, filename), overwrite = TRUE, gdal = "COMPRESS=DEFLATE")

  return(filename)
}


#' Find potential secondary roads from main road
#'
#' Tries to find secondary roads (branches) joining the main road using
#' a conductivity raster and the line geometry of the main road.
#'
#' @param road  line (\code{sf} format). Knowed road.
#' @param conductivity  raster (\code{terra} format). Conductivity raster covering the road.
#'
#' @return lines (\code{sf} format) representing potential starting points of branching roads.
#' @export
#' @examples
#' library(sf)
#' library(terra)
#'
#' road <- system.file("extdata", "drived_road.gpkg", package="MFFProads")
#' conductivity <- system.file("extdata", "drived_conductivity.tif", package="MFFProads")
#'
#' road <- st_read(road, "drived_road", quiet = TRUE)
#' conductivity <- rast(conductivity)
#'
#' road_branches <- find_road_branches(road, conductivity)
#' plot(road, col = "green", lwd = 3)
#' plot(road_branches, col = "red", add = TRUE)
find_road_branches <- function(road, conductivity)
{
  # Subset conductivity raster around the main road
  # Only consider conductivity values over 0.2 as good enough
  # to be considered as being potentially part of a road
  buffer40m <- sf::st_buffer(road, 40)
  conduc_crop <- terra::crop(conductivity, buffer40m)
  conductivity_threshold <- terra::mask(conduc_crop > 0.2, terra::vect(buffer40m))

  # Extract boolean raster which should mainly contain
  # pixels from the input road and the potential side jonctions
  # It's kind of doing a region growing operation with a seed on
  # one of the vertices of the main road
  coords_road <- sf::st_coordinates(road)
  m_seed <- coords_road[round(nrow(coords_road)/2), -3] |>
    matrix(nrow = 1)

  region_grow <- terra::patches(conductivity_threshold, directions = 8, zeroAsNA = TRUE)
  value_seed <- terra::extract(region_grow, m_seed)[1,1]
  conductivity_region <- terra::mask(conduc_crop, region_grow == value_seed, inverse = TRUE, maskvalues = 1)

  # Create (raster) search zone as a band around main road where potential jonction roads could be.
  # The search zone is between 20 m and 30 m on each side of the main road
  buffer20m <- sf::st_buffer(road, 20)
  buffer30m <- sf::st_buffer(road, 30)

  search_zone <- sf::st_difference(buffer30m, buffer20m) |>
    terra::vect() |>
    terra::rasterize(conductivity_region)

  # Inside the search zone raster, only keep pixels with high conductivity
  conductivity_search_zone <- terra::mask(conductivity_region, search_zone, inverse = TRUE, maskvalues = 1)
  conductivity_search_zone <- terra::mask(conductivity_region, conductivity_search_zone > 0.4, inverse = TRUE, maskvalues = 1)

  # Find clump of pixels that are big enough to be considered
  # as having a good potential to be part of a road
  # These pixels are then converted to points
  pts_clump <- conductivity_search_zone |>
    terra::patches(directions = 8, zeroAsNA = TRUE) |>
    terra::as.points() |>
    sf::st_as_sf()

  pts_clump_filtered <- pts_clump |>
    dplyr::group_by_at(1) |>
    dplyr::add_count() |>
    dplyr::filter(n > 10)

  # Find centroid of each filtered clump of pixels/points
  # to narrow a bit more the number of possible jonction
  pts_clump_centroid <- pts_clump_filtered |>
    dplyr::group_by_at(1) |>
    dplyr::group_split() |>
    lapply(function(x) { colMeans(sf::st_coordinates(x)) |> c(n = x[1,][["n"]]) }) |>
    do.call(what = rbind) |>
    as.data.frame() |>
    sf::st_as_sf(coords = c("X","Y"), crs = sf::st_crs(road))

  # Link each clump centroid to the closest point on main road
  # by taking the least-cost path
  starting_roads_pot <- sf::st_nearest_points(pts_clump_centroid, road)

  starts <- lwgeom::st_startpoint(starting_roads_pot) |> sf::as_Spatial()
  ends <- lwgeom::st_endpoint(starting_roads_pot) |> sf::as_Spatial()

  trans <- conductivity_region |>
    raster::raster() |>
    transition()

  branches_full <- lapply(seq_along(starts), function(i) { gdistance::shortestPath(trans, starts[i], ends[i], output = "SpatialLines") }) |>
    suppressWarnings() |>
    do.call(what = rbind) |>
    sf::st_as_sf() |>
    sf::st_set_crs(sf::st_crs(road)) |>
    suppressWarnings()

  # As the closest point on the main road might force a longer path
  # (in fact, in all cases where the main road and the jonction road
  # don't meet at 90°), a small buffer around the main road is used
  # to clip the least-cost path. The remaining part of the path should
  # be more representative of the real jonction point between the two roads.
  buffer2m <- sf::st_buffer(road, 2)
  branches_clip <- sf::st_difference(branches_full, buffer2m)

  # Only keep first LINESTRING (which will always be the one containing
  # the starting point) if the difference operation generated
  # many. Also reverse vertex order at the end so that branches
  # are going away from the main road
  branches_first <- branches_clip |>
    sf::st_geometry() |>
    lapply(function(x) { sf::st_cast(x, "LINESTRING") }) |>
    suppressWarnings() |>
    sf::st_sfc() |>
    sf::st_as_sf(crs = sf::st_crs(road)) |>
    sf::st_reverse()

  # Compute the final cost of each path.
  # The cost is divided by the length of each path in order
  # to exclude some paths that are inefficient and most likely
  # false positives
  starts <- lwgeom::st_startpoint(branches_first) |> sf::as_Spatial()
  ends <- lwgeom::st_endpoint(branches_first) |> sf::as_Spatial()

  costs <- sapply(seq_along(starts), function(i) { as.numeric(gdistance::costDistance(trans, starts[i], ends[i])) })
  lengths <- as.numeric(sf::st_length(branches_first))

  branches_cost <- branches_first |>
    dplyr::mutate(cost_m = costs/lengths) |>
    dplyr::filter(cost_m < 1.7) |>
    sf::st_simplify(dTolerance = 2)

  return(branches_cost)
}


#' Conductivity raster from boolean raster
#'
#' Conductivity raster from boolean raster. Can be useful if one has a road raster generated
#' by classification algorithm and whant to extract centerlines from it. The moving window
#' will enhance the center of an object. In cases of roads (narrow elongated shape), the centerline
#' will becone more apparent.
#'
#' @param x  raster (\code{terra} format). Boolean raster (1/0) where all 1 represent a pixel that might be categorized as a road.
#' @param w_max  numeric. Maximum window size to use (start from 3x3). Must be an odd number over 3.
#'
#' @return  raster (\code{terra} format) of conductivity with values ranging from 0 to 1.
#' @export
#' @examples
#' library(terra)
#'
#' segmented_road <- system.file("extdata", "segmented_road.tif", package="MFFProads")
#' segmented_road <- rast(segmented_road)
#'
#' conductivity <- conductivity_from_bool(segmented_road)
#' plot(conductivity)
conductivity_from_bool <- function(x, w_max = 15)
{
  if ((w_max %% 2 == 0) | (w_max <= 3)) stop("'w_max' must be an odd number over 3.", call. = FALSE)

  conductivity <- seq(3, w_max, 2) |>
    lapply(function(w) terra::focal(x, w)/w^2) |>
    terra::sds() |>
    terra::app(fun = mean)

  return(conductivity)
}
