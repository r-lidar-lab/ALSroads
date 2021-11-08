make_road_dtm = function(dtm, xc)
{
  row <- raster::rowFromY(dtm, xc)
  roadz <- dtm[raster::cellFromRow(dtm, row)]
  m = matrix(roadz, nrow(dtm), ncol(dtm), byrow = TRUE)
  road_dtm = dtm
  road_dtm[] = m
  return(road_dtm)
}
