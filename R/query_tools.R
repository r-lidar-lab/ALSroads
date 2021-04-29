clip_orectangle_with_index = function(las, xmin, ymin, xmax, ymax, angle)
{
  index = las@index$quadtree
  ids = filter_orectangle_with_index(index, xmin, xmax, ymin, ymax, angle)
  return(las[ids])
}

clip_longline = function(ctg, line, buffer, chunk = 300)
{
  poly <- sf::st_buffer(line, buffer)
  grid <- sf::st_make_grid(poly, cellsize = chunk, square = TRUE)
  u <- sf::st_intersection(sf::st_geometry(poly),grid)
  u <- sf::st_as_sf(u)
  las <- lidR::clip_roi(ctg, u)
  if (methods::is(las, 'LAS')) return(las)
  las <- do.call(rbind, las)
  return(las)
}

extract_road = function(las, road, param)
{
  # Automatically extract the points of the road if the input is a LAScatalog
  if (methods::is(las, "LAScatalog"))
  {
    cat("Extracting the point cloud...\n")
    width <- param[["extraction"]][["road_buffer"]]
    las <- clip_longline(las, road, (width + 5)/2)
  }

  # Create a spatial index for fast consecutive range queries
  if (methods::is(las, "LAS"))
  {
    if (is.null(las@index[["quadtree"]]))
    {
      cat("Creating a spatial index...\n")
      las@index[["quadtree"]] <- quadtree(las)
    }
  }

  return(las)
}
