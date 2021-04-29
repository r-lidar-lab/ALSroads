plot_road_width = function(road_norm_segment, nlas_segment, m, profiles_gnd, profiles_veg)#, profiles_veg)
{
  dd = profiles_gnd
  dd2 = profiles_veg
  X = road_norm_segment$Y
  Z = road_norm_segment$Z

  plot(X, Z, asp = 1, cex = 0.2, col = (road_norm_segment$Classification == 2) + 1, ylim = c(min(c(-7,Z)), max(15, max(Z))))

  #lines(dd$Xr,dd$sdZ*10, col = "darkgreen")
  graphics::lines(dd$Xr, dd$ssdZ*10, col = "red")
  graphics::lines(dd$Xr, dd$avgZ, col = "blue")
  graphics::lines(dd2$Xr, dd2$pabove05*2, col = "darkgreen")

  graphics::abline(h = 0)
  graphics::abline(h = c(1,2), lty = 3, col = "gray")
  graphics::abline(v = m$accotement$left$start, lty = 3, col = "blue")
  graphics::abline(v = m$accotement$right$start, lty = 3, col = "blue")
  graphics::abline(v = m$accotement$left$end, lty = 3, col = "lightskyblue")
  graphics::abline(v = m$accotement$right$end, lty = 3, col = "lightskyblue")
  graphics::abline(v = m$road_center, lty = 1, col = "black")
  graphics::abline(v = m$drivable_edges, lty = 3, col = "darkgreen")
  graphics::abline(v = m$rescue_edges, lty = 3, col = "red")

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
    graphics::arrows(m$road_edges[1], 11, m$road_edges[2], 11, code = 3, length = 0.1, col = "purple")
  graphics::text(mean(m$road_edges), 11.5, paste0(m$road_width, " m"), col = "purple")

  x0 = min(X)+5
  mz = max(Z)
  graphics::text(x0, mz+2, paste0("> 2 m: ", m$pzabove2, "%"), col = "purple")
  graphics::text(x0, mz+4, paste0("> 50 cm: ", m$pzabove05, "%"), col = "purple")

  graphics::text(x0, mz+8, paste0("> 2 m: ", m$pzabove2_drive, "%"), col = "darkgreen")
  graphics::text(x0, mz+10, paste0("> 50 cm: ", m$pzabove05_drive, "%"), col = "darkgreen")

  cols = c("darkgreen", "orange", "red", "black")
  col = cols[m$state]
  graphics::points(x0, min(Z)-3, col = col, pch = 19, cex = 5)

  #lines(dd2$Xr, dd2$pabove2/10, col = "purple")
  #abline(v = m$drivable_edge, lty = 3, col = "purple")
  #arrows(m$drivable_edge[1], 14, m$drivable_edge[2], 14, code = 3, length = 0.1)
  #text(mean(m$drivable_edge), 14.7, paste0(m$drivable_width, " m"))
}

plot_road_metrics = function(chemin, segment_metrics, find_scores = NULL)
{
  opar = graphics::par(mfrow=c(3,3))
  on.exit({graphics::par(opar)})

  #plot((seq_along(zroad)-1)*10, zroad, type = "l", asp = 1, xlab = "X", ylab = "Z")

  seuils = c(0,20,40,70)
  cols = c("darkgreen", "orange", "red", "black")

  col = cols[data.table::last(which(seuils < chemin$PABOVE2))]
  plot(segment_metrics$pzabove2, type = "l", ylim = c(0, 100), main = "% points above 2 m", ylab = "%")
  graphics::abline(h = chemin$PABOVE2, col = col)
  graphics::abline(h = seuils[1], col = cols[1], lty = 3)
  graphics::abline(h = seuils[2], col = cols[2], lty = 3)
  graphics::abline(h = seuils[3], col = cols[3], lty = 3)
  graphics::abline(h = seuils[4], col = cols[4], lty = 3)

  col = cols[data.table::last(which(seuils < chemin$PABOVE05))]
  plot(segment_metrics$pzabove05, type = "l", ylim = c(0, 100), main = "% points above 0.5 m", ylab = "%")
  graphics::abline(h = chemin$PABOVE05, col = col)
  graphics::abline(h = seuils[1], col = cols[1], lty = 3)
  graphics::abline(h = seuils[2], col = cols[2], lty = 3)
  graphics::abline(h = seuils[3], col = cols[3], lty = 3)
  graphics::abline(h = seuils[4], col = cols[4], lty = 3)

  plot(segment_metrics$road_width, type = "l", xaxt="n", main = "Road width profile", ylim = c(0, max(segment_metrics$road_width)))
  graphics::axis(1, at = 1:nrow(segment_metrics), las = 2)
  graphics::abline(h = chemin$ROADWIDTH, col = "blue", lwd = 2)

  col = rev(cols)[data.table::first(which(seuils > chemin$DRIVABLEWIDTH))]
  plot(segment_metrics$drivable_width, type = "l", xaxt="n", main = "Drivable width profile", ylim = c(0, max(segment_metrics$drivable_width)))
  graphics::axis(1, at = 1:nrow(segment_metrics), las = 2)
  graphics::abline(h = mean(chemin$DRIVABLEWIDTH), col = col, lwd = 2)

  plot(segment_metrics$number_accotements, type = "l", xaxt="n", main = "Number of embankments", ylim = c(0, 2))
  graphics::axis(1, at = 1:nrow(segment_metrics), las = 2)
  graphics::abline(h = mean(segment_metrics$number_accotements), col = "red", lwd = 2)

  if (!is.null(find_scores))
  {
    plot(find_scores, type = "l", main = "Finding scores", ylim = c(0,50))
    abline(h = mean(find_scores, na.rm = TRUE))
  }
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

