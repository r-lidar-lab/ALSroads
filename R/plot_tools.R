plot_road_width = function(road_norm_segment, nlas_segment, m, profiles_gnd, profiles_veg)#, profiles_veg)
{
  Xr <- NULL

  dd = profiles_gnd
  dd2 = profiles_veg
  X = road_norm_segment$Y
  Z = road_norm_segment$Z

  plot(X, Z, asp = 1, cex = 0.3, pch = 19, col = (road_norm_segment$Classification == 2) + 1, ylim = c(min(c(-7,Z)), max(15, max(Z))))
  graphics::abline(h = 0)
  graphics::abline(v = 0, lty = 3)
  graphics::abline(v = m$xc, lty = 1, col = "black", lwd = 2)

  #lines(dd$Xr,dd$sdZ*10, col = "darkgreen")
  #graphics::lines(dd$Xr, dd$ssdZ*10, col = "red")
  #graphics::lines(dd$Xr, dd$avgZ, col = "blue")
  #graphics::lines(dd2$Xr, dd2$pabove05*2, col = "darkgreen")



  if (!isFALSE(getOption("MFFProads.debug.measuring.size")))
  {
  #graphics::abline(h = c(1,2), lty = 3, col = "gray")
  graphics::abline(v = m$accotement$left$start, lty = 3, col = "blue")
  graphics::abline(v = m$accotement$right$start, lty = 3, col = "blue")
  graphics::abline(v = m$accotement$left$end, lty = 3, col = "lightskyblue")
  graphics::abline(v = m$accotement$right$end, lty = 3, col = "lightskyblue")
  graphics::abline(v = m$drivable_edges, lty = 3, col = "darkgreen")
  graphics::abline(v = m$rescue_edges[1], lty = 3, col = "red")
  graphics::abline(v = m$rescue_edges[2], lty = 3, col = "red")
  graphics::abline(v = m$road_edges[1], lty = 2, col = "pink", lwd = 2)
  graphics::abline(v = m$road_edges[2], lty = 2, col = "pink", lwd = 2)


  graphics::legend(x = "topright",
         legend = c("Shoulders edge", "Rescue edges", "Right of way edges", "Road center", "Drivable edges", "Road edges"),
         col = c("blue", "red", "lightskyblue", "purple", "darkgreen", "pink"),
         lty = c(3,3,3,1,3,1),
         lwd = c(1,1,1,2,1,2))


  graphics::arrows(m$accotement$left$start, -2, m$accotement$right$start, -2, code = 3, length = 0.1, col = "blue")
  graphics::text(mean(c(m$accotement$left$start, m$accotement$right$start)), -2.5, paste0(m$accotement$right$start-m$accotement$left$start, " m"), col = "blue")

  graphics::arrows(m$accotement$left$end, -4, m$accotement$right$end, -4, code = 3, length = 0.1, col = "lightskyblue")
  graphics::text(mean(c(m$accotement$left$end, m$accotement$right$end)), -4.5, paste0(m$accotement$right$end-m$accotement$left$end, " m"), col = "lightskyblue")


  if (m$drivable_width > 0)
    graphics::arrows(m$drivable_edges[1], 12.5, m$drivable_edges[2], 12.5, code = 3, length = 0.1, col = "darkgreen")
  graphics::text(mean(m$drivable_edges), 13.2, paste0(m$drivable_width, " m"), col = "darkgreen")

  if (!is.na(m$rescue_width))
  {
    graphics::arrows(m$rescue_edges[1], -3, m$rescue_edges[2], -3, code = 3, length = 0.1, col = "red")
    graphics::text(mean(m$rescue_edges), -3.5, paste0(m$rescue_width, " m"), col = "red")
  }

  if (m$road_width > 0)
    graphics::arrows(m$road_edges[1], 11, m$road_edges[2], 11, code = 3, length = 0.1, col = "pink")
  graphics::text(mean(m$road_edges), 11.5, paste0(m$road_width, " m"), col = "pink")
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

  graphics::legend(x = "bottomright",
                   legend = c("Std. Z", "Ground points", "Intensity", "Average", "ax\u00B2+bx+c"),
                   col = c("red", "blue", "orange", "purple", "green"),
                   lty = 1,
                   lwd = c(1,1,1,2,3))
  }
}

plot_road_metrics = function(chemin, road_metrics, segment_metrics)
{
  opar = graphics::par(mfrow=c(2,3))
  on.exit({graphics::par(opar)})

  #plot((seq_along(zroad)-1)*10, zroad, type = "l", asp = 1, xlab = "X", ylab = "Z")

  seuils = c(0, param$state$percentage_veg_thresholds)
  cols = c("darkgreen", "orange", "red")

  col = cols[data.table::last(which(seuils < road_metrics$PABOVE2))]
  plot(segment_metrics$pzabove2, type = "l", ylim = c(0, 100), main = "% points above 2 m", ylab = "%")
  graphics::abline(h = road_metrics$PABOVE2, col = col)
  graphics::abline(h = seuils[1], col = cols[1], lty = 3)
  graphics::abline(h = seuils[2], col = cols[2], lty = 3)
  graphics::abline(h = seuils[3], col = cols[3], lty = 3)

  col = cols[data.table::last(which(seuils < road_metrics$PABOVE05))]
  plot(segment_metrics$pzabove05, type = "l", ylim = c(0, 100), main = "% points above 0.5 m", ylab = "%")
  graphics::abline(h = road_metrics$PABOVE05, col = col)
  graphics::abline(h = seuils[1], col = cols[1], lty = 3)
  graphics::abline(h = seuils[2], col = cols[2], lty = 3)
  graphics::abline(h = seuils[3], col = cols[3], lty = 3)


  plot(segment_metrics$road_width, type = "l", xaxt="n", main = "Road width profile", ylim = c(0, max(segment_metrics$road_width)), ylab = "Width (m)")
  graphics::axis(1, at = 1:nrow(segment_metrics), las = 2)
  graphics::abline(h = road_metrics$ROADWIDTH, col = "blue", lwd = 2)

  cols = rev(cols)
  seuils = c(0, param$state$drivable_width_thresholds)
  col = cols[data.table::last(which(seuils < road_metrics$DRIVABLEWIDTH))]

  plot(segment_metrics$drivable_width, type = "l", xaxt="n", main = "Drivable width profile", ylim = c(0, max(segment_metrics$drivable_width)), ylab = "Width (m)")
  graphics::axis(1, at = 1:nrow(segment_metrics), las = 2)
  graphics::abline(h = mean(road_metrics$DRIVABLEWIDTH), col = col, lwd = 2)
  graphics::abline(h = seuils[1], col = cols[1], lty = 3)
  graphics::abline(h = seuils[2], col = cols[2], lty = 3)
  graphics::abline(h = seuils[3], col = cols[3], lty = 3)

  plot(segment_metrics$number_accotements, type = "l", xaxt="n", main = "Number of embankments", ylim = c(0, 2), ylab = "n")
  graphics::axis(1, at = 1:nrow(segment_metrics), las = 2)
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

