plot_road_width = function(road_norm_segment, nlas_segment, m, profiles_gnd, profiles_veg)#, profiles_veg)
{
  Xr <- NULL

  dd = profiles_gnd
  dd2 = profiles_veg
  X = road_norm_segment$Y
  Z = road_norm_segment$Z
  col = lidR:::set.colors(Z, lidR::height.colors(25))
  col[road_norm_segment$Classification == 2] = "purple"
  cex = rep(0.3, lidR::npoints(road_norm_segment))
  cex[road_norm_segment$Classification == 2] = 0.5

  plot(X, Z, asp = 1, cex = cex, pch = 19, col = col, ylim = c(min(c(-7,Z)), max(15, max(Z))), xlab = "X (m)", ylab = "Z (m)")
  graphics::abline(h = 0)
  graphics::lines(c(0,0), c(0, 15), lty = 3)
  graphics::lines(c(m$xc, m$xc), c(0, -8), lty = 3, col = "black")
  graphics::text(m$xc, -9, "Adjusted\ncenterline", cex = 0.8)
  graphics::text(0, 16, "Pathfinder\ncenterline", cex = 0.8)

  #lines(dd$Xr,dd$sdZ*10, col = "darkgreen")
  #graphics::lines(dd$Xr, dd$ssdZ*10, col = "red")
  #graphics::lines(dd$Xr, dd$avgZ, col = "blue")
  #graphics::lines(dd2$Xr, dd2$pabove05*2, col = "darkgreen")

  if (!isFALSE(getOption("MFFProads.debug.measuring.size")))
  {
    zaccotementstart = -3
    zaccotementend = 13
    zdrivable = 9
    zrescue = 8
    zedge = 11

    #graphics::abline(h = c(1,2), lty = 3, col = "gray")
    graphics::lines(c(m$accotement$left$start, m$accotement$left$start), c(0, zaccotementstart), lty = 3, col = "blue")
    graphics::lines(c(m$accotement$right$start, m$accotement$right$start), c(0, zaccotementstart), lty = 3, col = "blue")
    graphics::arrows(m$accotement$left$start, zaccotementstart, m$accotement$right$start, zaccotementstart, code = 3, length = 0.1, col = "blue")
    graphics::text(mean(c(m$accotement$left$start, m$accotement$right$start)), zaccotementstart - 0.5, paste0(m$accotement$right$start-m$accotement$left$start, " m"), col = "blue")

    graphics::lines(c(m$accotement$left$end, m$accotement$left$end), c(0, zaccotementend), lty = 3, col = "purple")
    graphics::lines(c(m$accotement$right$end, m$accotement$right$end), c(0, zaccotementend), lty = 3, col = "purple")
    graphics::arrows(m$accotement$left$end, zaccotementend, m$accotement$right$end, zaccotementend, code = 3, length = 0.1, col = "purple")
    graphics::text(mean(c(m$accotement$left$end, m$accotement$right$end)), zaccotementend+0.5, paste0(m$accotement$right$end-m$accotement$left$end, " m"), col = "purple")

    graphics::lines(c(m$drivable_edges[1], m$drivable_edges[1]), c(0, zdrivable), lty = 3, col = "darkgreen")
    graphics::lines(c(m$drivable_edges[2], m$drivable_edges[2]), c(0, zdrivable), lty = 3, col = "darkgreen")
    if (m$drivable_width > 0)
      graphics::arrows(m$drivable_edges[1], zdrivable, m$drivable_edges[2], zdrivable, code = 3, length = 0.1, col = "darkgreen")
    graphics::text(mean(m$drivable_edges), zdrivable+0.5, paste0(m$drivable_width, " m"), col = "darkgreen")

    graphics::lines(c(m$rescue_edges[1], m$rescue_edges[1]), c(0,zrescue), lty = 3, col = "pink")
    graphics::lines(c(m$rescue_edges[2], m$rescue_edges[2]), c(0,zrescue), lty = 3, col = "pink")
    if (!is.na(m$rescue_width)) {
      graphics::arrows(m$rescue_edges[1], -3, m$rescue_edges[2], -3, code = 3, length = 0.1, col = "pink")
      graphics::text(mean(m$rescue_edges), -3.5, paste0(m$rescue_width, " m"), col = "pink")
    }

    graphics::lines(c(m$road_edges[1], m$road_edges[1]), c(0, zedge), lty = 3, col = "red", lwd = 1)
    graphics::lines(c(m$road_edges[2], m$road_edges[2]), c(0, zedge), lty = 3, col = "red", lwd = 1)
    if (m$road_width > 0)
      graphics::arrows(m$road_edges[1], zedge, m$road_edges[2], zedge, code = 3, length = 0.1, col = "red")
    graphics::text(mean(m$road_edges), zedge+0.5, paste0(m$road_width, " m"), col = "red")


    graphics::legend(x = "topright",
                     legend = c("Shoulders edges", "Rescue edges", "Right of way edges", "Drivable edges", "Road edges"),
                     col = c("blue", "pink", "purple","darkgreen", "red"),
                     lty = c(3,3,3,3,3),
                     lwd = c(1,1,1,1,1))
  }

  if (!isFALSE(getOption("MFFProads.debug.measuring.pos")))
  {
    D <- dd2[Xr >= -7 & Xr <= + 7]
    sdZ = ma(D$sdZ, 12)
    sdZ =  sdZ/max(sdZ, na.rm = T)
    ngnd = ma(D$ngnd, 12)
    ngnd = ngnd/max(ngnd, na.rm = T)
    ngnd = 1- (ngnd - min(ngnd, na.rm = T))
    cvi = ma(D$cvi, 12)
    cvi =  cvi/max(cvi, na.rm = T)
    all = sdZ + 2*ngnd + cvi
    lm1  <- stats::lm(all ~ D$Xr + I(D$Xr^2))
    coe  <- stats::coefficients(lm1)
    y = stats::predict(lm1, D)

    graphics::lines(D$Xr, sdZ*5, col = "red")
    graphics::lines(D$Xr, ngnd*5, col = "blue")
    graphics::lines(D$Xr, cvi*5, col = "orange")
    graphics::lines(D$Xr, all/4*5, col = "purple", lwd = 2)
    graphics::lines(D$X, y/4*5, col = "green", lwd = 3)

    graphics::legend(x = "topleft",
                     legend = c("Std. Z", "Ground points", "Intensity", "Average", "ax\u00B2+bx+c"),
                     col = c("red", "blue", "orange", "purple", "green"),
                     lty = c(1,1,1,1,1),
                     lwd = c(1,1,1,2,3))
  }
}

plot_road_metrics = function(chemin, road_metrics, segment_metrics)
{
  opar = graphics::par(mfrow=c(2,2))
  on.exit({graphics::par(opar)})

  #plot((seq_along(zroad)-1)*10, zroad, type = "l", asp = 1, xlab = "X", ylab = "Z")

  seuils = c(0, param$state$percentage_veg_thresholds)
  cols = c("darkgreen", "orange", "red")

  # col = cols[data.table::last(which(seuils < road_metrics$PABOVE2))]
  # plot(segment_metrics$pzabove2, type = "l", ylim = c(0, 100), main = "Points above 2 m", ylab = "%", xlab = NULL)
  # graphics::abline(h = road_metrics$PABOVE2, col = col)
  # graphics::abline(h = seuils[1], col = cols[1], lty = 3)
  # graphics::abline(h = seuils[2], col = cols[2], lty = 3)
  # graphics::abline(h = seuils[3], col = cols[3], lty = 3)

  col = cols[data.table::last(which(seuils < road_metrics$PABOVE05))]
  plot(segment_metrics$pzabove05, type = "l", ylim = c(0, 100), main = "Points in [0.5, 3] m", ylab = "%", xlab = NULL)
  graphics::abline(h = road_metrics$PABOVE05, col = col)
  graphics::abline(h = seuils[1], col = cols[1], lty = 3)
  graphics::abline(h = seuils[2], col = cols[2], lty = 3)
  graphics::abline(h = seuils[3], col = cols[3], lty = 3)


  plot(segment_metrics$road_width, type = "l", main = "Road width", ylim = c(0, max(segment_metrics$road_width)), ylab = "Width (m)", xlab = NULL)
  graphics::abline(h = road_metrics$ROADWIDTH, col = "blue", lwd = 2)

  cols = rev(cols)
  seuils = c(0, param$state$drivable_width_thresholds)
  col = cols[data.table::last(which(seuils < road_metrics$DRIVABLEWIDTH))]

  plot(segment_metrics$drivable_width, type = "l", main = "Drivable width", ylim = c(0, max(segment_metrics$drivable_width)), ylab = "Width (m)", xlab = NULL)
  graphics::abline(h = mean(road_metrics$DRIVABLEWIDTH), col = col, lwd = 2)
  graphics::abline(h = seuils[1], col = cols[1], lty = 3)
  graphics::abline(h = seuils[2], col = cols[2], lty = 3)
  graphics::abline(h = seuils[3], col = cols[3], lty = 3)

  plot(segment_metrics$number_accotements, type = "l", main = "Shoulders", ylim = c(0, 2), ylab = "n", xlab = NULL)
  graphics::abline(h = mean(segment_metrics$number_accotements), col = "red", lwd = 2)

  seuils = c(param$state$score_thresholds, 1)
  cols = c("red", "orange", "darkgreen")

  #score = chemin$SCORE
  #col = cols[data.table::first(which(seuils > score))]
  #plot(sf::st_geometry(chemin), col = col, main = paste0("Score = ", score))
  #graphics::axis(1)
  #graphics::axis(2)
  #graphics::box()
}

add_road_insert = function(chemin, i, p1, p2, col)
{
  graphics::title(paste("segment", i))
  opar = graphics::par( fig = c(.7, .95, .7, .95), mar=.1+c(0,0,0,0), new = TRUE )
  plot(sf::st_geometry(chemin))
  graphics::points(p1[1], p1[2], cex = 0.5, pch = 19, col = col)
  graphics::points(p2[1], p2[2], cex = 0.5, pch = 19, col = col)
  graphics::box()
  graphics::par(opar)
}

plot2D = function(las)
{
  plot(las$Y, las$Z, asp = 1, cex = 0.2, col = (las$Classification == 2) + 1, ylim = c(min(c(-7,las$Z)), max(15, max(las$Z))))
}

