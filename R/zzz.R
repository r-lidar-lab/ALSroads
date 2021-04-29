.onLoad <- function(libname, pkgname)
{
  op <- options()
  op.MFFProads<- list(
    MFFProads.debug.finding = FALSE,
    MFFProads.debug.measuring = FALSE,
    MFFProads.debug.metrics = FALSE,
    MFFProads.debug = FALSE)

  toset <- !(names(op.MFFProads) %in% names(op))
  if (any(toset)) options(op.MFFProads[toset])
}

.onUnload <- function(libpath)
{
  library.dynam.unload("MFFProads", libpath)
}

.datatable.aware = TRUE
