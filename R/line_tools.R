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
  lines = x
  n <- nrow(lines)
  S <- numeric(n)
  for (i in 1:n)
  {
    line <- lines[i,]

    if (sf::st_geometry_type(line) != "LINESTRING")
      stop(paste0("Geometry ", i, " is not a LINESTRING"))

    S[i] <- sinuosity(line)
  }

  round(S,2)
}

sinuosity.sfc_LINESTRING <- function(x)
{
  XY <- sf::st_coordinates(x)
  return(sinuosity(XY))
}

cut_line <- function(line, at, metrics)
{
  if (is.integer(at))
  {
    v <- rle(at)
    j <- c(1, cumsum(v$length))
    p <- cumsum(v$lengths/sum(v$lengths))
  }
  else
  {
    p <- at
  }

  n <- length(p)
  p <- c(0, p)
  o <- vector("list", n)
  for(i in 1:n)
  {
    if (!is.null(metrics))m = metrics[j[i]:j[i+1]]
    l = lwgeom::st_linesubstring(line, p[i], p[i+1])
    if (!is.null(metrics))w = 4-m$state
    l = l[1]
    if (!is.null(metrics))l$ROADWIDTH = round(sum(w*m$road_width)/sum(w),1)
    if (!is.null(metrics))l$DRIVABLEWIDTH = round(mean(m$drivable_width),1)
    if (!is.null(metrics))l$STATE = v$values[i]
    o[[i]] = l
  }

  lines = do.call(rbind, o)
  return(lines)
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
