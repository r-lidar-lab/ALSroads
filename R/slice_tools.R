#' Computes size and state metrics of a given slice
#'
#' This function takes a properly oriented slice perpendicular to the road and (1) relocated more accurately
#' the centerline (2) measure the width of the road
#'
#' @param nlas_slice LAS. Normalized slice that will be analysed either as is or using raw value, i.e.
#' after unnormalization
#' @param list. parameters
#'
#' @noRd
slice_metrics = function(nlas_slice, param)
{
  pr <- param[["extraction"]][["profile_resolution"]]

  # We need both normalized and normalized versions
  las_slice <- lidR::unnormalize_height(nlas_slice)

  # Compute some profiles of vegetation to detect some open areas
  veg_profiles <- compute_veg_profiles(nlas_slice, res = pr)

  # Retrieve a more accurate road centre at +/- 5m around current centreline
  # because the least cost path is computed with a resolution of 2 m only and has a tendency
  # to go straight
  xc <- relocate_centerline(veg_profiles)

  # Normalize but relatively to the road only to get a flat and horizontal road
  # (Not described in the paper)
  dtm <- lidR::grid_terrain(las_slice, pr, lidR::knnidw(), fast = TRUE, Wdegenerated = FALSE, full_raster = TRUE)
  road_dtm <- make_road_dtm(dtm, xc)
  road_norm_slice <- lidR::normalize_height(las_slice, road_dtm, na.rm = TRUE)

  # Compute some profiles from ground points
  gnd_profiles <- compute_gnd_profiles(road_norm_slice, res = pr)

  # Find the shoulder i.e. some slopes on each side of the centerline that delimit the roads
  min_slope <- param[["embankments"]][["min_slope"]]
  shoulders <- find_shoulders(gnd_profiles, xc, th_slope = min_slope)
  shoulders_edges <- c(shoulders$left$start, shoulders$right$start)

  # We may not have found the shoulders in this case we need a rescue method to find the edges
  # of the road. It arises if one or two shoulders are missing
  rescue_edges = c(NA_real_, NA_real_)
  if (!shoulders$left$embankment | !shoulders$right$embankment)
    rescue_edges <- find_road_edges(gnd_profiles, xc, shoulders, thsd = param[["terrain"]][["max_sd_ground_points"]], thz = param[["terrain"]][["max_elevation_ground_points"]])

  # The road edges are computed using preliminary the shoulders, then the rescue edges if no
  # shoulder is found
  left  <- if (!shoulders$left$embankment)  rescue_edges[1] else shoulders$left$start
  right <- if (!shoulders$right$embankment) rescue_edges[2] else shoulders$right$start
  if (left > right) stop("Internal error when measuring road size")
  road_edges <- c(left, right)

  # Using the road edges and vegetation profile we can estimate the derivable edges
  drivable_edges <- find_drivable_edges(veg_profiles, road_edges, xc)

  # Compute different estimation of the width
  shoulders_width <- round(diff(shoulders_edges), 1)
  rescue_width    <- round(diff(rescue_edges), 1)
  drivable_width  <- round(diff(drivable_edges$edges), 1)
  road_width      <- diff(road_edges)

  # Compute some metrics within the boundaries of the road by extracting only the point within
  # the limits of the road
  Y    <- NULL
  road <- lidR::filter_poi(nlas_slice, Y > road_edges[1], Y < road_edges[2])
  n    <- lidR::npoints(road)

  # Coverage extra metrics we want to compute in addition to the width
  pground         <- NA   # Percentage of ground points (useless?)
  pzabove05       <- 0    # Percentage of point between 0.5 and below 5 (P)
  pzabove2        <- 100  # Percentage of point between 2 and below 5 (useless?)
  pground_drive   <- NA   # The same but within the drivable zone (useless?)
  pzabove05_drive <- 0
  pzabove2_drive  <- 100

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

  # Compute the elevation of the road
  zroad <- road_dtm[1]

  # Computation of the location of the road in the original CRS (currently we only have xc = +/- 5m)
  rot <- attributes(las_slice)$rotation
  off <- attributes(las_slice)$offset
  yc  <- mean(range(las_slice$X))
  xyc <- cbind(yc, xc) %*% t(rot)
  xyc[,1] <- xyc[,1] + off[1]
  xyc[,2] <- xyc[,2] + off[2]

  # Strore the metrics of interest. Some are only useful for display and debugging. Some seem
  # to no longer be useful and need cleaning
  m = list(
    # Use to plot only (discarded when returned)
    accotement = shoulders,
    accotement_edges = shoulders_edges,
    rescue_edges = rescue_edges,
    drivable_edges = drivable_edges$edges,
    road_edges = road_edges,

    # More or less interesting metrics
    accotement_width = shoulders_width,
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
    number_accotements = shoulders$left$embankment + shoulders$right$embankment,
    xc = xc)

  if (getOption("MFFProads.debug.measuring")) plot_road_width(road_norm_slice, nlas_slice, m, gnd_profiles, veg_profiles)

  return(m[-c(1:5)])
}

#' Compute profiles from road normalized oriented point-cloud
#'
#' The metrics useful are std. Z and avg Z. The profiles are smoothed. sd is computed
#' with fast custom fsd function.
#' @noRd
compute_gnd_profiles = function(road_norm_pslice, res)
{
  Z <- Y <-Classification <- ReturnNumber <- .N <- Xr <- NULL
  gnd_profiles <- road_norm_pslice@data[Classification %in% c(2L, 9L),
                                        list(sdZ = fsd(Z), avgZ = mean(Z)),
                                        by = list(Xr = lidR:::round_any(Y, res))]
  data.table::setkey(gnd_profiles, Xr)
  gnd_profiles$ssdZ <- ma(gnd_profiles$sdZ)
  gnd_profiles$avgZ <- ma(gnd_profiles$avgZ)
  return(gnd_profiles)
}

compute_veg_profiles = function(nlas_slice, res = 0.5)
{
  right <- Classification <- Y <- Z <- Intensity <- .N <- Xr <- NULL

  if (!"Intensity" %in% names(nlas_slice))
    nlas_slice@data[["Intensity"]] <- 0L

  veg_profiles <- nlas_slice@data[, list(ngnd = sum(Classification == 2),
                                 nabove2 = sum(Z >= 2 & Z < 5),
                                 nabove05 = sum(Z >= 0.5 & Z < 5),
                                 nabove01 = sum(Z >= 0.1 & Z < 5),
                                 npoint = .N,
                                 npoint2 = sum(Z < 5),
                                 cvi = diff(range(Intensity)),
                                 sdZ = fsd(Z)),
                          by = list(Xr = lidR:::round_any(Y, res))]

  data.table::setkey(veg_profiles, Xr)
  veg_profiles$ngnd[is.na(veg_profiles$ngnd)] = 0
  veg_profiles$pabove2 <- veg_profiles$nabove2/veg_profiles$npoint2
  veg_profiles$pabove2[is.nan(veg_profiles$pabove2)] <- 0
  veg_profiles$pabove2 <- ma(veg_profiles$pabove2)

  veg_profiles$pabove05 <- veg_profiles$nabove05/veg_profiles$npoint2
  veg_profiles$pabove05[is.nan(veg_profiles$pabove05)] <- 0
  veg_profiles$pabove05 <- ma(veg_profiles$pabove05)

  veg_profiles$pabove01 <- veg_profiles$nabove01/veg_profiles$npoint2
  veg_profiles$pabove01[is.nan(veg_profiles$pabove01)] <- 0
  veg_profiles$pabove01 <- ma(veg_profiles$pabove01)

  return(veg_profiles)
}

find_shoulders = function(profiles, xc, th_slope = 10)
{
  gnd_profiles <- profiles

  y <- gnd_profiles$avgZ
  x <- gnd_profiles$Xr
  s <- deriv(x,y)
  s <- abs(s)
  s <- round(atan(s)*180/pi,2)
  gnd_profiles$s <- s

  # Adjust xc
  Xr <- NULL
  veg_profiles <- gnd_profiles[Xr >= xc-2 & Xr <= xc+2]
  if (nrow(veg_profiles) > 2)
  {
    i <- which.min(ma(veg_profiles$s))
    if (length(i) > 0)
      xc = veg_profiles$Xr[i]
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
  left_edges = if(length(left) > 0) f(gnd_profiles$s, x, left, th_slope) else default
  right_edges = if(length(right) > 0) f(gnd_profiles$s, x, right, th_slope) else default
  return(list(left = left_edges, right = right_edges, center = (left_edges$start + right_edges$start)/2))
}

find_road_edges = function(profiles, xc, accotement, thsd = 0.1, thz = 0.15)
{
  gnd_profiles = profiles

  #xr = if (accotement$right$embankment) accotement$right$start else xc +5
  #xl = if (accotement$left$embankment) accotement$left$start else xc - 5

  #veg_profiles = gnd_profiles[Xr >= xl & Xr <= xr]
  #i = which.min(veg_profiles$ssdZ)
  #xc = veg_profiles$Xr[i]

  sdZ = gnd_profiles$ssdZ
  avgZ = gnd_profiles$avgZ

  right = which(gnd_profiles$Xr >= xc)
  left = rev(which(gnd_profiles$Xr <= xc))

  if (length(right) > 0)
    pos1 = screen_profile(right, sdZ, avgZ, thsd, thz)
  else
    pos1 = nrow(gnd_profiles)

  if (length(left) > 0)
    pos2 = screen_profile(left, sdZ, avgZ, thsd, thz)
  else
    pos2 = 1

  xi = gnd_profiles$Xr[pos2]
  xf = gnd_profiles$Xr[pos1]

  edges = c(xi,xf)
  return(edges)
}

find_drivable_edges = function(profiles, true_edges, xc)
{
  Xr <- NULL
  veg_profiles = profiles
  veg_profiles = veg_profiles[Xr >= true_edges[1] & Xr <= true_edges[2]]

  if (nrow(veg_profiles) == 1) return(list(center = xc, edges = c(xc,xc)))

  idx = which(veg_profiles$pabove05 < min(veg_profiles$pabove05, na.rm = TRUE)+0.01)
  if (is.na(idx[1])) idx = which(veg_profiles$pabove2 < min(veg_profiles$pabove2, na.rm = TRUE)+0.01)
  if (is.na(idx[1])) return(list(center = xc, edges = c(xc,xc)))

  right = which(veg_profiles$Xr >= xc)
  left = rev(which(veg_profiles$Xr <= xc))

  if (length(right) > 0)
    pos1 = screen_profile(right, veg_profiles$pabove05, NULL, 0.15, NULL)
  else
    pos1 = left[1]

  if (length(left) > 0)
    pos2 = screen_profile(left, veg_profiles$pabove05, NULL, 0.15, NULL)
  else
    pos2 = right[1]

  xi = veg_profiles$Xr[pos2]
  xf = veg_profiles$Xr[pos1]

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


relocate_centerline <- function(veg_profiles)
{
  xc <- 0
  Xr <- Y <- NULL
  D  <- veg_profiles[Xr >= xc-7 & Xr <= xc + 7]

  #
  sdZ  <- ma(D$sdZ, 12)
  if (!all(is.na(sdZ)))
    sdZ  <-  sdZ/max(sdZ, na.rm = T)

  ngnd_max <- max(D$ngnd, na.rm = TRUE)
  if (ngnd_max == 0) ngnd_max = 1

  ngnd <- ma(D$ngnd, 12)
  ngnd <- ngnd/ngnd_max
  ngnd <- 1 - (ngnd - min(ngnd, na.rm = T))

  cvi  <- ma(D$cvi, 12)
  cvi  <- cvi/max(cvi, na.rm = T)
  all  <- sdZ + 2*ngnd + cvi

  if (sum(!is.na(all)) > 3)
  {
    tryCatch(
    {
      lm1  <- stats::lm(all~ D$Xr + I(D$Xr^2))
    },
    error = function(e)
    {
      stop("Unexpected linear model failure")
    })

    coe  <- stats::coefficients(lm1)
    u    <- -coe[2]/(2*coe[3])
    if (coe[3] > 0 & u < 5 & u > -5) xc <- u
    if (xc < min(D$Xr)) xc <- min(D$Xr)
    if (xc > max(D$Xr)) xc <- max(D$Xr)
  }

  return(xc)
}
