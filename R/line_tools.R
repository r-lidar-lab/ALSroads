#' Get heading of both ends of a line
#'
#' Retreive heading of both ends of a line, pointing away from it.
#'
#' @param line  line (\code{sfc} or \code{sfg} format)
#'
#' @return numeric of length 2 expressing respectively the head and tail heading in degrees
#' (range [-180°, 180°] using \code{atan2()}.
#' @noRd
st_ends_heading <- function(line)
{
  M <- sf::st_coordinates(line)
  n = nrow(M)
  
  headings <- sapply(c(2, n), function(i) {
    if (i == 2) {j <- 1} else {j <- -1 ; i <- n-1}
    Ax = M[i-j,1]
    Ay = M[i-j,2]
    Bx = M[i,1]
    By = M[i,2]
    atan2(Ay-By, Ax-Bx)*180/pi
  })
  
  return(unname(headings))
}


#' Extend line by given distance
#'
#' Extend one or both ends of a line by a given distance.
#' No new vertices are added, instead, the first/last is moved
#' to its new location following the same heading as before.
#'
#' @param line  line (\code{sfc} or \code{sfg} format)
#' @param distance  numeric; distance by which the line will be extended. In case of \code{end = "BOTH"},
#' a vector of length 2 can be provided to extend respectively the head and tail end by different values.
#' @param end  character; (\code{BOTH}, \code{HEAD} or \code{TAIL}; define which end will be extended.
#'
#' @return Same line as input but now extended.
#' @noRd
st_extend_line <- function(line, distance, end = "BOTH")
{
  if (!(end %in% c("BOTH", "HEAD", "TAIL")) | length(end) != 1) stop("'end' must be 'BOTH', 'HEAD' or 'TAIL'")

  M <- sf::st_coordinates(line)[,-3]

  ends <- c(1, nrow(M))
  headings <- st_ends_heading(line)/180*pi
  distances <- distance[1:2]

  if (end == "HEAD") {
    ends <- ends[1]
    headings <- headings[1]
    distances[2] <- distance[1]
  } else if (end == "TAIL") {
    ends <- ends[2]
    headings <- headings[2]
    distances[2] <- distance[1]
  }
  
  M[ends, 1:2] <- M[ends, 1:2] + distances * c(cos(headings), sin(headings))
  newline <- sf::st_linestring(M)

  # If input is sfc_LINESTRING and not sfg_LINESTRING
  if (is.list(line)) newline <- sf::st_sfc(newline, crs = sf::st_crs(line))
  
  return(newline)
}


#' Get point from distance on a line
#'
#' This is essentially the reverse of rgeos::gProject(). It must be noted
#' that due floating point precision issue, the point returned won't be
#' exactly on the line and thus won't split it.
#'
#' @param distance  distance on the \code{line} at which the coordinates will be retreived.
#' @param line  line (\code{sfc} format)
#' @param from_end  logical; (\code{FALSE} by default); if \code{TRUE} start measuring from the end
#' of \code{line} instead of from the beginning.
#'
#' @return point (as \code{sfc_POINT}) as close as possible of being on the \code{line} at the
#' given distance
#' @noRd
st_point_on_line <- function(distance, line, from_end = FALSE)
{
  if (distance > as.numeric(sf::st_length(line))) stop("'distance' must be smaller than the total length of 'line'")

  coords <- sf::st_coordinates(line)

  if (from_end) coords <- apply(coords, 2, rev)
  
  dist_lagged <- sqrt(diff(coords[,1])^2 + diff(coords[,2])^2)
  dist_along_line <- data.frame(dist = dist_lagged,
                                cumsum = cumsum(dist_lagged))
  
  idx <- which(dist_along_line[,2] >= distance)[1]
  
  theta <- atan2(coords[idx+1,2] - coords[idx,2],
                 coords[idx+1,1] - coords[idx,1])
  hypo <- distance - dist_along_line[idx,2] + dist_along_line[idx,1]
  
  x_junction <- coords[idx,1] + hypo * cos(theta)
  y_junction <- coords[idx,2] + hypo * sin(theta)

  pt_junction <- c(x_junction, y_junction) |>
    sf::st_point() |>
    sf::st_sfc(crs = sf::st_crs(line))

  return(pt_junction)
}


#' Split line at a given point
#'
#' Split line at a specified point. The point must be
#' within the tolerance value of the line for a split to occur.
#'
#' @param split_point  point (\code{sfc} format) at which the \code{line} will be split.
#' @param line  line (\code{sfc} format)
#' @param tolerance  numeric; maximum distance allowed between the point and the line
#' for a split to occur.
#'
#' @return geometries (as \code{sfc_LINESTRING}) resulting from the split.
#' @noRd
st_split_at_point <- function(split_point, line, tolerance = 0.01)
{
  # Make a small perpendicular linestring for the splitting
  # operation as the split point provided won't split line
  # as it can't be exactly on the line, except at one of the
  # line's vertices. The reason is that the segment between two
  # vertices is define by a function and thus (almost) no point
  # on it can be described by a 64-bit float value.
  blade <- sf::st_nearest_points(split_point, line)

  if (as.numeric(sf::st_length(blade)) > tolerance) stop("'split_point' too far from 'line' to perform split")

  # Sharpen blade (by slightly extending it)!
  coords <- sf::st_coordinates(blade)
  theta <- atan2(diff(coords[,2]), diff(coords[,1]))

  xmax <- max(coords[,"X"]) + abs(0.01 * cos(theta))
  ymax <- max(coords[,"Y"]) + abs(0.01 * sin(theta))
  xmin <- min(coords[,"X"]) - abs(0.01 * cos(theta))
  ymin <- min(coords[,"Y"]) - abs(0.01 * sin(theta))
  
  # Split line and extract the two parts
  split_line <- rbind(c(xmax, ymax), c(xmin, ymin)) |>
    sf::st_linestring() |>
    lwgeom::st_split(x = line) |>
    sf::st_collection_extract("LINESTRING")
  
  return(split_line)
}


sinuosity <- function(x)
{
  UseMethod("sinuosity", x)
}

sinuosity.matrix = function(x)
{
  m  <- nrow(x)

  dx <- diff(x[,1])
  dy <- diff(x[,2])
  L  <- sum(sqrt(dx^2 + dy^2))

  dx <- diff(x[c(1,m),1])
  dy <- diff(x[c(1,m),2])
  D  <- sum(sqrt(dx^2 + dy^2))

  return(L/D)
}

sinuosity.sf <- function(x)
{
  lines <- x
  n <- nrow(lines)
  S <- numeric(n)
  for (i in 1:n)
  {
    line <- lines[i,]

    if (sf::st_geometry_type(line) != "LINESTRING")
      stop(paste0("Geometry ", i, " is not a LINESTRING"))

    S[i] <- sinuosity.sfc_LINESTRING(line)
  }

  round(S,2)
}

sinuosity.sfc_LINESTRING <- function(x)
{
  XY <- sf::st_coordinates(x)
  return(sinuosity.matrix(XY))
}

adjust_spline = function(points)
{
  # Adjust a spline to create a smooth line from points
  xroad <- sf::st_coordinates(points)[,1]
  yroad <- sf::st_coordinates(points)[,2]
  troad <- points$distance_to_start
  wroad <- points$find_score

  ux <- stats::smooth.spline(troad, xroad, spar = 0.4, all.knots = TRUE)
  uy <- stats::smooth.spline(troad, yroad, spar = 0.4, all.knots = TRUE)

  # Sometime one point is missing in the spline for an unknown reason. If the output
  # does not have the same length than the input we resize the input to fit the output
  # and prevent failures
  rm <- FALSE
  if (length(ux$x) < length(troad))
  {
    rm <- troad %in% ux$x
    xroad = xroad[rm]
    yroad = yroad[rm]
    troad = troad[rm]
    wroad = wroad[rm]
  }

  # Maybe it is possible to have an extra point. This case is handled here but never observed
  if (length(ux$x) > length(troad))
  {
    print(dput(points))
    stop("Different length")
  }

  # The adjusted spline can be far from some outlier points. This is the role of the spline
  # to smooth the road. But if some important outliers are detected far from the spline (more than
  # 10 m) we remove those points and fit again the spline
  rm2 <- FALSE
  rmx <- abs(ux$y - xroad) < 10
  rmy <- abs(uy$y - yroad) < 10

  if (any(!rmx) | any(!rmy))
  {
    rm2   <- rmx & rmy
    xroad <- xroad[rm2]
    yroad <- yroad[rm2]
    troad <- troad[rm2]
    wroad <- wroad[rm2]

    if (length(xroad) > 3)
    {
      ux <- stats::smooth.spline(troad, xroad, w = wroad)
      uy <- stats::smooth.spline(troad, yroad, w = wroad)
    }
    else
    {
      rm2 <- FALSE
    }
  }

  spline <- sf::st_drop_geometry(points)
  if (!isFALSE(rm))  spline <- spline[rm,]
  if (!isFALSE(rm2)) spline <- spline[rm2,]
  spline$x <- ux$y
  spline$y <- uy$y
  spline <- sf::st_as_sf(spline, coords = c("x", "y"))

  return(sf::st_sfc(sf::st_linestring(sf::st_coordinates(spline))))
}

st_merge_line = function(x)
{
  u <- sf::st_geometry(sf::st_linestring(sf::st_coordinates(sf::st_cast(sf::st_geometry(x), "POINT"))))
  u <- sf::st_set_crs(u, sf::st_crs(x))
  return(u)
}

st_is_loop = function(line)
{
  # Exit early for loops because does not work (yet?)
  p1 <- lwgeom::st_startpoint(line)
  p2 <- lwgeom::st_endpoint(line)
  d  <- as.numeric(sf::st_distance(p1,p2)[1,1])
  return(d < 2)
}

st_angles <- function(line)
{
  M <- sf::st_coordinates(line)
  M <- M[,-3]

  n = nrow(M)
  if (n < 3) return(0)

  angles <- numeric(n)
  angles[] <- NA_real_
  for (i in 2:(n-1))
  {
    Ax = M[i-1,1]
    Ay = M[i-1,2]
    Bx = M[i,1]
    By = M[i,2]
    Cx = M[i+1,1]
    Cy = M[i+1,2]
    u <- c(Bx-Ax, By-Ay)
    v <- c(Cx-Bx, Cy-By)
    angles[i] <- angle(u,v)
  }

   return(abs(angles[-c(1,n)]))
}

angle <- function(x,y)
{
  dot.prod <- x%*%y
  norm.x <- norm(x,type="2")
  norm.y <- norm(y,type="2")
  r <- as.numeric(dot.prod / (norm.x * norm.y))
  if (r > 0.999999999) return(0)
  theta <- acos(r)
  return(theta*180/pi)
}
