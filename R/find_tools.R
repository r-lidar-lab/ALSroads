find_road = function(nlas_segment, param, prev_xc = NULL)
{
  . <- NULL

  # We use both normalized and non normalized point-cloud to assess the position of the
  # road in a given segment
  las_segment <- lidR::unnormalize_height(nlas_segment)

  # Adjust intensity for multiple flightline.
  if (data.table::uniqueN(nlas_segment$PointSourceID) > 1)
  {
    PointSourceID <- .N <- Intensity <- NULL
    imean <- nlas_segment@data[, .(I = mean(Intensity), n = .N), by = PointSourceID]
    iref <- imean$I[which.max(imean$n)]
    imean$corr <- iref / imean$I
    data <- nlas_segment@data[imean, on = 'PointSourceID']
    nlas_segment$Intensity <- as.integer(data$Intensity * data$corr)
  }

  # The method relies on the computation of 'probability' of the existence of the road
  # along the x direction. These 'probabilities' are computed at these value of x
  x <- seq(ceiling(min(las_segment$Y)), floor(max(las_segment$Y)), 1)

  # Compute the profiles that are used to estimate the location of the road.
  gnd <- compute_gnd_profiles(las_segment, res = param[["extraction"]][["profile_resolution"]])
  veg <- compute_veg_profiles(nlas_segment,res =  param[["extraction"]][["profile_resolution"]])

  # Compute the slope of the terrain from the average of ground points
  slope  <- deriv(gnd$Xr, gnd$avgZ)
  sslope <- ma(slope*180/pi, 5)
  slope  <- abs(slope)
  slope  <- round(atan(slope)*180/pi,2)

  # Perpendicularly to the road the terrain is expected to be flat.
  # We allocate probabilities to slope < 5 degrees
  xr <- gnd$Xr
  w_flat <- allocate_probabilities(x, xr, slope, param[["terrain"]][["max_flat_slope"]])

  # Water bodies are extremely flat and might be interpreted as a road
  # If the point cloud is water classified the following allows to multiplies
  # the allocated probabilities by 0 in water so water cannot be detected as road
  water <- veg$water
  dtx <- data.table::data.table(x)
  dtw <- data.table::data.table(x = veg$Xr, water)
  data.table::setkey(dtx, x)
  data.table::setkey(dtw, x)
  water <- dtw[dtx,nomatch = NA]$water
  water[is.na(water)] <- 1
  f_water <- ma(water, 5)
  f_water[is.na(f_water)] <- 1
  f_water <- 1 - f_water

  # From the slope we can estimate if there is likely something that look like a road.
  # We allocate probabilities to places that are likely between two embankments
  slope_th = param[["embankments"]][["min_slope"]]
  width_max = param[["embankments"]][["max_width"]]
  is_embankement <- is.embankement(sslope) | is.embankement(-sslope)
  w_bank = allocate_probabilities(x, xr, !is_embankement, 1)

  # Allocate probabilities to places with vegetation > 50 cm below 25%
  xr <- veg$Xr
  u <- veg$pabove05
  w_up50 <- allocate_probabilities(x, xr, u, param[["vegetation"]][["max_percentage_point_above_50cm"]])

  # Allocate probabilities to spikes in the profile vegetation > 50 cm
  u <- stats::approx(xr,u,x)$y
  u <- double_pass_robust_peak_detection(u, param[["peak"]][["lag"]], param[["peak"]][["threshold"]], param[["peak"]][["influence"]]) <= -1
  w_sp50 <- allocate_probabilities(x, x, 1-u, 1)

  # Allocate probabilities to places with vegetation > 10 cm below 25%
  xr <- veg$Xr
  u <- veg$pabove01
  w_up10 <- allocate_probabilities(x, xr, u, param[["vegetation"]][["max_percentage_point_above_10cm"]])

  # Allocate probabilities to spikes in the profile vegetation > 10 cm
  u <- stats::approx(xr,u,x)$y
  u <- double_pass_robust_peak_detection(u, param[["peak"]][["lag"]], param[["peak"]][["threshold"]], param[["peak"]][["influence"]]) <= -1
  w_sp10 <- allocate_probabilities(x, x, 1-u, 1)

  # Allocate probabilities to places with  stdZ std below 0.15
  xr <- veg$Xr
  u <- veg$sdZ
  w_stdz <- allocate_probabilities(x, xr, u, param[["vegetation"]][["max_stdz"]])

  # Allocate probabilities to spikes in the profile stdZ std
  u <- stats::approx(xr,u,x)$y
  u <- double_pass_robust_peak_detection(u, param[["peak"]][["lag"]], param[["peak"]][["threshold"]], param[["peak"]][["influence"]]) <= -1
  w_spstdz <- allocate_probabilities(x, x, 1-u, 1)

  # Allocate probabilities to places with CV of intensity below x
  xr <- veg$Xr
  u <- veg$cvi
  mnu <- min(u, na.rm = T)
  mxu <- max(u, na.rm = T)
  du <- mxu-mnu
  w_cvi <- allocate_probabilities(x, xr, u, mnu+0.1*du)

  # Allocate probabilities to spike of cvi
  u <- stats::approx(xr,u,x)$y
  u <- abs(double_pass_robust_peak_detection(u, param[["peak"]][["lag"]], param[["peak"]][["threshold"]], param[["peak"]][["influence"]])) <= -1
  w_spcvi <- allocate_probabilities(x, x, 1-u, 1)

  # Total of allocated probabilities. f_water comes in factor to clear possible high
  # probabilities in water. Some metrics are weighted because they are usually more robust.
  pf = f_water*(0.5*w_flat+w_bank+w_up50+2*w_up10+2*w_stdz+w_cvi+2*w_sp50+2*w_sp10+w_spcvi)
  pf = ma(pf, 3)

  # pf is expected to have high probability spots. We will find them and the highest is
  # the location of the road
  if (all(pf == 0))
  {
      xc <- 0
      score <- 0
  }
  else
  {
    # We must also consider the previous locations of the road that are likely
    # to be close to the current position. This way if the current location is not
    # the global maximum we can find it anyway. If there is no previous location
    # we use the global maximum. Otherwise we find the high spots and determine which one
    # match best.
    if (is.infinite(prev_xc[1]))
    {
      i <- which.max(pf)
      xc <- x[i]
      score <- pf[i]
    }
    else
    {
      # We find the high probability spots
      fsp <- find_maxima(x, pf, FALSE)

      h <- fsp$Height
      nn1 <- abs(fsp$X - prev_xc[1])
      i1 <- which.min(nn1)
      i2 <- NA

      if (!is.infinite(prev_xc[2]))
      {
        nn2 <- abs(fsp$X - prev_xc[2])
        i2 <- which.min(nn2)
      }

      # The high spot are multiplied if they are close to the previous location
      h[i1] <- h[i1]*2
      if (!is.na(i2)) h[i2] = h[i2]*1.5

      # We now search for the highest spot. Which is not necessarily the one
      # that is the closest to the previous location. Very high probabilities
      # can save previous errors
      i <- which.max(h)
      xc <- fsp$X[i]
      score <- h[i]
    }
  }

  # Display almost everything we did so far
  if (getOption("MFFProads.debug.finding"))
  {
    opar = graphics::par(mfrow=c(4,2))
    on.exit(graphics::par(opar))

    minz = min(las_segment$Z)

    plot(las_segment$Y,las_segment$Z, col = lidR:::set.colors(nlas_segment$Z, lidR::height.colors(50)), asp = 1, pch = 19, cex = 0.2, main = "LiDAR", xlab = "X", ylab = "Z")
    graphics::abline(v = xc)
    graphics::abline(v = prev_xc[1], lty = 2)
    graphics::abline(v = prev_xc[2], lty = 3)
    graphics::text(xc - 6, max(las_segment$Z) + 10, label = "Road")
    graphics::par(new=TRUE)
    plot(x, f_water, type = "l", col = "purple", lwd = 2, axes = FALSE, ylim = c(0, 3))
    graphics::par(new=TRUE)
    plot(x, water, type = "S", col = "deepskyblue", lwd = 1, axes = FALSE, ylim = c(0, 3))

    plot(gnd$Xr, gnd$avgZ-minz, col = "black", asp = 1, type = "l", main = "Slope", xlab = "X", ylab = "Z")
    graphics::par(new=TRUE)
    plot(gnd$Xr, slope, col = "red", type = "l")
    graphics::par(new=TRUE)
    plot(x,w_flat, type = "l", col = "purple", lwd = 2, axes = FALSE, ylim = c(0, 3))
    graphics::axis(4)

    plot(gnd$Xr, gnd$avgZ-minz, col = "black", asp = 1, type = "l", main = "Embankment", xlab = "X", ylab = "Slope")
    graphics::par(new=TRUE)
    plot(gnd$Xr, sslope, type = "l", col = "red")
    graphics::par(new=TRUE)
    plot(x, w_bank, type = "l", col = "purple", lwd = 2, axes = FALSE, ylim = c(0, 3))
    graphics::axis(4)

    plot(veg$Xr, veg$sdZ, col = "green", type = "l", main = "Std. Z", ylim = c(0, max(veg$sdZ, na.rm = T)), xlab = "X", ylab = "std. Z")
    graphics::par(new=TRUE)
    plot(x,w_stdz, type = "l", col = "purple", axes = FALSE, lwd = 2, ylim = c(0, 3))
    graphics::lines(x,w_spstdz, type = "l", col = "deeppink3", lwd = 2, ylim = c(0, 3))
    graphics::axis(4)

    plot(veg$Xr, veg$pabove05, type = "l", col = "blue", main = "% above 0.5 m", ylim = c(0, max(veg$pabove05, na.rm = T)))
    graphics::par(new=TRUE)
    plot(x,w_up50,type = "l", col = "purple", lwd = 2, axes = FALSE, ylim = c(0, 3))
    graphics::lines(x,w_sp50, col = "deeppink3", lwd = 2)
    graphics::axis(4)

    plot(veg$Xr, veg$pabove01, type = "l", col = "deepskyblue", main = "% above 0.1 m", ylim = c(0, max(veg$pabove01, na.rm = T)))
    graphics::par(new=TRUE)
    plot(x, w_up10, type = "l", col = "purple", lwd = 2, axes = FALSE, ylim = c(0, 3))
    graphics::lines(x, w_sp10, col = "deeppink3", lwd = 2)
    graphics::axis(4)

    plot(veg$Xr, (veg$cvi-min(veg$cvi, na.rm = T)), type = "l", col = "orange", main = "std. intensity", xlab = "X", ylab = "Std I")
    graphics::par(new=TRUE)
    plot(x,w_cvi, type = "l", col = "purple", lwd = 2, axes = FALSE, ylim = c(0, 3))
    graphics::lines(x,w_spcvi, col = "deeppink3", lwd = 2)
    graphics::axis(4)

    if (is.null(find_maxima(x, pf, TRUE)) ) plot(x,pf, type = "l", col = "purple", lwd = 2)
    graphics::abline(v = xc)
    graphics::abline(v = prev_xc, lty = 3)

    graphics::par(opar)
  }

  return(list(xc = xc, score = score))
}

allocate_probabilities = function(x, xr, u, th)
{
  stopifnot(length(xr) == length(u))

  p = numeric(length(x))
  for (i in seq_along(xr))
  {
    if (!is.na(u[i]) && u[i] < th)
    {
      # u[i] < th create a distrib with a std [1, 5] propotionnal
      # to of low is u[i]. If u[i] == 0 -> std = 2 , if u[i] == th -> std = 5
      perc = (1-(th - u[i])/th)
      stdmin = 2
      stdmax = 4
      std = (stdmax-stdmin)*perc + stdmin
      p = p + stats::dnorm(x, xr[i], std)
    }
  }

  return(p)
}

find_maxima = function(x, y, disp = TRUE)
{
  sy = sign(diff(y))
  i = which(sy == 0)
  sy[i] <- 1
  i = which(diff(sy) == -2) + 1

  if (length(i) == 0) return(NULL)

  spikes = vector("list", length(i))

  if (disp) plot(x,y, type = "l")

  for (k in 1:length(i))
  {
    ii = i[k]

    if (disp) graphics::lines(c(x[ii], x[ii]), c(0, y[ii]), lty = 3, col = "red")

    j = i[k]
    while (y[j-1] <= y[j] && j > 2)
      j = j-1

    j0 = j

    j = i[k]
    while(y[j+1] <= y[j] && j < length(y)-1)
      j = j+1

    j1 = j

    jj = c(j0, j1)


    Ampl1 =  y[ii] - y[j0]
    Ampl2 =  y[ii] - y[j1]
    Width = x[j1] - x[j0]
    Area = sum(y[j0:j1])

    if (disp)
    {
      graphics::points(x[jj], y[jj], col = "blue", pch = 19)
      graphics::lines(c(x[j0], x[j0]), c(0, y[j0]), lty = 3, col = "blue")
      graphics::lines(c(x[j1], x[j1]), c(0, y[j1]), lty = 3, col = "blue")
      graphics::lines(c(x[ii], x[ii]), c(0, (Ampl1 + Ampl2)/2), col = "gray", lwd = 3)
    }

    spikes[[k]] <- list(X = x[ii], left = y[j0], right = y[j1], Ampl1 = Ampl1, Ampl2 = Ampl2, Ampl = (Ampl1 + Ampl2)/2, Width = Width, Area = Area, Height = y[ii])
  }

  spikes <- data.table::rbindlist(spikes)
  spikes
}

is.embankement = function(sslope, slope_th = 10, width_max = 20)
{
  # Search embankment range
  start_left = NA
  start_right = NA
  n <- length(sslope)
  i <- 1
  lim <- slope_th
  sslope[is.na(sslope)] <- 0
  is_embankement = rep(FALSE, n)
  while (i <= n)
  {
    if (sslope[i] > lim)
    {
      while(sslope[i] > lim && i <= n) i = i + 1
      start_left = i
    }

    if (sslope[i] > lim && !is.na(start_left))
    {
      while(sslope[i] > lim && i <= n) i = i + 1
      start_left = i
    }

    if (sslope[i] < -lim && !is.na(start_left))
    {
      start_right = i
      while(sslope[i] < -lim && i <= n) i = i + 1
    }

    if (!is.na(start_left) && !is.na(start_right))
    {
      if (start_right - start_left < width_max)
        is_embankement[start_left:start_right] = TRUE

      start_right = NA
      start_left = NA
    }

    i = i + 1
  }

  return(is_embankement)
}
