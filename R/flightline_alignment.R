flightlines_z_misalignment_matrix <- function(las, res = 5)
{
  PointSourceID <- NULL

  lay <- lidR:::raster_layout(las, res)
  lay <- lidR:::raster_materialize(lay)
  ids <- unique(las$PointSourceID)
  res <- vector("list", 0)
  data.table::setindex(las@data, PointSourceID)

  # If not populated
  if(length(ids) == 1 && ids == 0)
  {
    las <- lidR::retrieve_flightlines(las)
    las$PointSourceID <- las$flightlineID
    ids <- unique(las$PointSourceID)
  }

  if (length(ids) == 1) {
    M <- matrix(0,1,1)
    rownames(M) = ids
    colnames(M) = ids
    return(list(offsets = M, count = M))
  }

  for (i in ids) {
    psi <- las@data[PointSourceID == i, .(X,Y,Z)]
    psi <- lidR::LAS(psi, las@header, check = FALSE)
    res[[as.character(i)]] <- lidR:::rasterize_fast(psi, lay, method = "min")
  }

  #plot(terra::rast(res))

  M <- matrix(0, length(ids), length(ids))
  rownames(M) = ids
  colnames(M) = ids
  N <- M

  for (i in 1:length(ids)) {
    for (j in 1:length(ids)) {
      if (i != j) {
        R <- res[[i]] - res[[j]]
        M[i,j] <- round(mean(R[], na.rm = TRUE), 3)
        N[i,j] <- sum(!is.na(R[]))
      }
    }
  }

  M[is.nan(M)] <- NA
  N[is.nan(N)] <- NA
  return(list(offsets = M, count = N))
}


flightlines_z_realignment <- function(las, M)
{
  .N <- PointSourceID <- NULL

  offset_tot = apply(M$offsets,1, \(x) sum(abs(x)))
  k = which(!is.na(offset_tot))

  if (length(k) == 0)
  {
    M <- fill_misalignment_matrix(M)
    island <- anyNA(M)
    M$offsets[is.nan(M$offsets)] = 0

    if (!island)
    {
      offset_tot = apply(M$offsets, 1, \(x) sum(abs(x)))
      k = which(!is.na(offset_tot))
    }
    else
    {
      offset_tot = -apply(M$count, 1, \(x) sum(abs(x)))
      k = which(!is.na(offset_tot))
    }
  }

  offset_tot[is.na(offset_tot)] = Inf
  k <- which.min(offset_tot)

  u <- las@data[, .N, by = PointSourceID]
  idx <- match(las$PointSourceID, u$PointSourceID)
  #k <- which.max(u$N)

  offsets <- as.numeric(M$offsets[,k])
  newZ <- las$Z - offsets[idx]
  lidR:::quantize(newZ, las[["Z scale factor"]], las[["Z offset"]], by_reference = TRUE)
  las@data[["Z"]] <- newZ
  return(las)
}

fill_misalignment_matrix <- function(M)
{
  N = M$count
  O = M$offsets
  Tr <- O
  Tr[is.na(Tr)] <- 0
  Tr[Tr != 0] <- 1
  G <- igraph::graph_from_adjacency_matrix(Tr)

  E <- O
  E[] <- 0

  for (i in 1:nrow(E)) {
    for (j in 1:ncol(E)) {
      if (i != j) {
        paths <- igraph::all_simple_paths(G, i, j, cutoff = 3)
        delta <- vector("list", length(paths))
        err <- get_all_path_errors(M, paths)
        E[i,j] <- round(sum(err$deltas*err$weight)/sum(err$weight), 3)
      }
    }
  }

  return(list(offsets = E, count = M$count))
}

get_all_path_errors = function(M, paths)
{
  deltas <- vector("numeric", length(paths))
  weight <- deltas

  for(i in seq_along(paths)) {
    idx <- as.numeric(paths[[i]])
    e <- 0
    w <- 0
    for (j in 1:(length(idx)-1))
    {
      e <- e + M$offsets[idx[j],idx[j+1]]
      w <- w + M$count[idx[j],idx[j+1]]
    }

    deltas[i] <- e
    weight[i] <- w
  }

  return(list(deltas = deltas, weight = weight))
}
