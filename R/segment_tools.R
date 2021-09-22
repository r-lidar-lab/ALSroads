#' Computes size and state metrics of a given section
#'
#' Computes
#' @noRd
segment_road_metrics = function(nlas_segment, param)
{
  # We need both normalized and normalized versions
  las_segment <- lidR::unnormalize_height(nlas_segment)

  # In measurement mode the road is on the centre
  xc <- 0

  # Compute some profiles of vegetation to detect some open areas
  dd2 <- compute_veg_profiles(nlas_segment, res = param[["extraction"]][["profile_resolution"]])

  # Retrieve a more accurate road centre at +/- 5m around current centre
  # in case the map provided is no strictly exact
  Xr <- Y <- NULL
  D <- dd2[Xr >= xc-7 & Xr <= xc + 7]
  sdZ  <- ma(D$sdZ, 12)
  sdZ  <-  sdZ/max(sdZ, na.rm= T)
  ngnd <- ma(D$ngnd, 12)
  ngnd <- ngnd/max(ngnd, na.rm = T)
  ngnd <- 1 - (ngnd - min(ngnd, na.rm = T))
  cvi  <- ma(D$cvi, 12)
  cvi  <-  cvi/max(cvi, na.rm = T)
  all  <- sdZ + 2*ngnd + cvi
  lm1  <- stats::lm(all~poly(D$Xr,2))
  coe  <- stats::coefficients(lm1)
  u    <- -coe[2]/(2*coe[3])
  if (coe[3] > 0 & u < 5 & u > -5) xc <- u


  # Normalize but relatively to the road only to get a flat and horizontal road
  dtm <- lidR::grid_terrain(las_segment, param[["extraction"]][["profile_resolution"]], lidR::knnidw(), fast = TRUE, Wdegenerated = FALSE)
  road_dtm <- make_road_dtm(dtm, xc)
  road_norm_segment <- lidR::normalize_height(las_segment, road_dtm, na.rm = TRUE)

  # Compute some profiles from ground points
  dd <- compute_gnd_profiles(road_norm_segment, res = param[["extraction"]][["profile_resolution"]])

  # Find the embankments i.e. some slopes on each side of xc that delimit the roads
  accotement <- find_accotement(dd, xc, th_slope = param[["embankments"]][["min_slope"]])
  accotement_edges <- c(accotement$left$start, accotement$right$start)

  # We may not have found embankments in this case we need a rescue method
  rescue_edges = c(NA_real_, NA_real_)
  if (!accotement$left$embankment | !accotement$right$embankment)
    rescue_edges <- find_road_edges(dd, xc, accotement, thsd = param[["terrain"]][["max_sd_ground_points"]], thz = param[["terrain"]][["max_elevation_ground_points"]])

  # The road edges are computed using first the embankments, then the rescue edges if no
  # embankment
  left = if (!accotement$left$embankment) rescue_edges[1] else accotement$left$start
  right = if (!accotement$right$embankment) rescue_edges[2] else accotement$right$start
  if (left > right) stop("Internal error when measuring road size")
  road_edges <- c(left, right)

  # Using the road edges and vegetation profile we can estimate the drivable edges
  drivable_edges <- find_drivable_edges(dd2, road_edges, xc)

  # Compute different estimation of the width
  accotement_width = round(diff(accotement_edges), 1)
  rescue_width = round(diff(rescue_edges), 1)
  drivable_width = round(diff(drivable_edges$edges), 1)
  road_width = diff(road_edges)

  # Compute some metrics within the boundaries of the road
  road <- lidR::filter_poi(nlas_segment, Y > road_edges[1], Y < road_edges[2])
  n    <- lidR::npoints(road)

  if (n > 0)
  {
    n5            <- sum(road$Z < 5)
    nabove2       <- sum(road$Z > 2 & road$Z < 5)
    nabove05      <- sum(road$Z > 0.5 & road$Z < 5)
    nground       <- sum(road$Classification == lidR::LASGROUND)
    pzabove2      <- round(nabove2 / n5 * 100, 1)
    pzabove05     <- round(nabove05 / n5 * 100, 1)
    pground       <- round(nground / n * 100, 1)

    # Compute some metrics within the boundaries of the drivable road
    road          <- lidR::filter_poi(road, Y > drivable_edges$edges[1], Y < drivable_edges$edges[2])
    n             <- lidR::npoints(road)
    n5            <- sum(road$Z < 5)
    nabove2       <- sum(road$Z > 2 & road$Z < 5)
    nabove05      <- sum(road$Z > 0.5 & road$Z < 5)
    nground       <- sum(road$Classification == lidR::LASGROUND)
    pzabove2_drive  <- round(nabove2 / n5 * 100, 1)
    pzabove05_drive <- round(nabove05 / n5 * 100, 1)
    pground_drive   <- round(nground / n * 100, 1)
  }
  else
  {
    pground         <- NA
    pzabove05       <- 0
    pzabove2        <- 100
    pground_drive   <- NA
    pzabove05_drive <- 0
    pzabove2_drive  <- 100
  }

  # Compute the elevation of the road
  zroad = road_dtm[1]

  # Location of the road in the original CRS
  rot <- attributes(las_segment)$rotation
  off <- attributes(las_segment)$offset

  yc <- mean(range(las_segment$X))
  xyc <- cbind(yc, xc)%*%t(rot)
  xyc[,1] <- xyc[,1] + off[1]
  xyc[,2] <- xyc[,2] + off[2]

  # Store the metrics
  m = list(
    # Use to plot only (discarded when m is returned)
    accotement = accotement,
    accotement_edges = accotement_edges,
    rescue_edges = rescue_edges,
    drivable_edges = drivable_edges$edges,
    road_edges = road_edges,

    # Metrics
    accotement_width = accotement_width,
    rescue_width = rescue_width,
    road_width = road_width,
    drivable_width = drivable_width,
    pzabove05 = pzabove05,
    pzabove2 = pzabove2,
    pzabove05_drive = pzabove05_drive,
    pzabove2_drive = pzabove2_drive,
    xroad = xyc[,1],
    yroad = xyc[,2],
    zroad = zroad,
    number_accotements = accotement$left$embankment + accotement$right$embankment,
    xc = xc)

  if (getOption("MFFProads.debug.measuring")) plot_road_width(road_norm_segment, nlas_segment, m, dd, dd2)

  return(m[-c(1:5)])
}

compute_gnd_profiles = function(road_norm_pslice, res = 0.25)
{
  Z <- Y <-Classification <- ReturnNumber <- .N <- Xr <- NULL
  dd <- road_norm_pslice@data[Classification %in% c(2L, 9L), list(sdZ = fsd(Z), avgZ = mean(Z), pfgnd = sum(ReturnNumber == 1)/.N), by = list(Xr = lidR:::round_any(Y, res))]
  data.table::setkey(dd, Xr)
  dd$sdZ <- zoo::na.approx(dd$sdZ, na.rm = F)
  dd$ssdZ <- ma(dd$sdZ)
  dd$avgZ <- zoo::na.approx(dd$avgZ, na.rm = F)
  dd$avgZ <- ma(dd$avgZ)
  dd
}

compute_veg_profiles = function(nlas_segment, res = 0.5)
{
  right <- Classification <- Y <- Z <- Intensity <- .N <- Xr <- NULL

  dd <- nlas_segment@data[, list(ngnd = sum(Classification == 2),
                                 nabove2 = sum(Z >= 2 & Z < 5),
                                 nabove05 = sum(Z >= 0.5 & Z < 5),
                                 nabove01 = sum(Z >= 0.1 & Z < 5),
                                 npoint = .N,
                                 npoint2 = sum(Z < 5),
                                 cvi = diff(range(Intensity)),
                                 sdZ = fsd(Z)),
                          by = list(Xr = lidR:::round_any(Y, res))]

  data.table::setkey(dd, Xr)
  dd$ngnd[is.na(dd$ngnd)] = 0
  dd$pabove2 <- dd$nabove2/dd$npoint2
  dd$pabove2[is.nan(dd$pabove2)] <- 0
  dd$pabove2 <- ma(dd$pabove2)

  dd$pabove05 <- dd$nabove05/dd$npoint2
  dd$pabove05[is.nan(dd$pabove05)] <- 0
  dd$pabove05 <- ma(dd$pabove05)

  dd$pabove01 <- dd$nabove01/dd$npoint2
  dd$pabove01[is.nan(dd$pabove01)] <- 0
  dd$pabove01 <- ma(dd$pabove01)

  return(dd)
}

find_accotement = function(profiles, xc, th_slope = 10)
{
  dd <- profiles

  y <- dd$avgZ
  x <- dd$Xr
  s <- deriv(x,y)
  s <- abs(s)
  s <- round(atan(s)*180/pi,2)
  dd$s <- s

  # Adjust xc
  Xr <- NULL
  dd2 <- dd[Xr >= xc-2 & Xr <= xc+2]
  if (nrow(dd2) > 2)
  {
    i <- which.min(ma(dd2$s))
    if (length(i) > 0)
      xc = dd2$Xr[i]
  }

  right <- which(x > xc)
  left <- rev(which(x < xc))
  f <- function(s, x, idx, th_slope)
  {
    start = 0
    end = 0
    j = 1
    i = 0
    while(end == 0 & i != data.table::last(idx))
    {
      i = idx[j]
      j = j+1

      if (is.na(s[i]))
       next

      if (s[i] >= th_slope & start == 0 )
        start = i

      if (s[i] <= th_slope & start != 0 )
        end = i
    }

    if (start == 0)
      start = idx[length(idx)]

    if (end == 0)
      end = start

    pos = x[c(start, end)]
    dx = abs(diff(pos))
    dz = round(abs(diff(y[c(start,end)])), 2)
    slope = round(atan(dz/dx)*180/pi,2)
    if (is.nan(slope)) slope = 0
    width = abs(dx)
    embankment = width >= 1 & slope >= 10
    return(list(start = pos[1], end = pos[2], dz = dz, slope = slope, width = width, embankment = embankment))
  }

  default <- list(start = xc, end = xc, dz = 0, slope = 0, width = 0, embankment = FALSE)
  left_edges = if(length(left) > 0) f(dd$s, x, left, th_slope) else default
  right_edges = if(length(right) > 0) f(dd$s, x, right, th_slope) else default
  return(list(left = left_edges, right = right_edges, center = (left_edges$start + right_edges$start)/2))
}

find_road_edges = function(profiles, xc, accotement, thsd = 0.1, thz = 0.15)
{
  dd = profiles

  #xr = if (accotement$right$embankment) accotement$right$start else xc +5
  #xl = if (accotement$left$embankment) accotement$left$start else xc - 5

  #dd2 = dd[Xr >= xl & Xr <= xr]
  #i = which.min(dd2$ssdZ)
  #xc = dd2$Xr[i]

  sdZ = dd$ssdZ
  avgZ = dd$avgZ

  right = which(dd$Xr >= xc)
  left = rev(which(dd$Xr <= xc))

  if (length(right) > 0)
    pos1 = screen_profile(right, sdZ, avgZ, thsd, thz)
  else
    pos1 = nrow(dd)

  if (length(right) > 0)
    pos2 = screen_profile(left, sdZ, avgZ, thsd, thz)
  else
    pos2 = 1

  xi = dd$Xr[pos2]
  xf = dd$Xr[pos1]

  edges = c(xi,xf)
  return(edges)
}

find_drivable_edges = function(profiles, true_edges, xc)
{
  Xr <- NULL
  dd2 = profiles
  dd2 = dd2[Xr >= true_edges[1] & Xr <= true_edges[2]]

  idx = which(dd2$pabove05 < min(dd2$pabove05, na.rm = TRUE)+0.01)
  if (is.na(idx[1])) idx = which(dd2$pabove2 < min(dd2$pabove2, na.rm = TRUE)+0.01)
  if (is.na(idx[1])) return(list(center = xc, edges = c(xc,xc)))

  right = which(dd2$Xr >= xc)
  left = rev(which(dd2$Xr <= xc))

  if (length(right) > 0)
    pos1 = screen_profile(right, dd2$pabove05, NULL, 0.15, NULL)
  else
    pos1 = left[1]

  if (length(left) > 0)
    pos2 = screen_profile(left, dd2$pabove05, NULL, 0.15, NULL)
  else
    pos2 = right[1]

  xi = dd2$Xr[pos2]
  xf = dd2$Xr[pos1]

  edges = c(xi, xf)

  return(list(center = xc, edges = edges))
}

screen_profile = function(idx, y1, y2, thsd, thz)
{
  if (is.null(y2))
  {
    y2 = rep(0, length(y1))
    thz = 0
  }

  pos = 0
  for (i in idx)
  {
    if ((is.na(y1[i]) | is.na(y2[i])) & pos == 0) {
      pos = i
      next
    }

    if ((y1[i] > thsd | abs(y2[i]) > thz) & pos == 0)
      pos = i
  }

  if (i == data.table::last(idx) & pos == 0)
    pos = data.table::last(idx)

  return(pos)
}


