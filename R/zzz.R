.onLoad <- function(libname, pkgname)
{
  op <- options()
  op.MFFProads<- list(
    MFFProads.debug.finding = FALSE,
    MFFProads.debug.measuring = FALSE,
    MFFProads.debug.metrics = FALSE,
    MFFProads.debug.verbose = FALSE,
    MFFProads.debug.progress = TRUE,
    MFFProads.debug = FALSE)

  toset <- !(names(op.MFFProads) %in% names(op))
  if (any(toset)) options(op.MFFProads[toset])
}

.onUnload <- function(libpath)
{
  library.dynam.unload("MFFProads", libpath)
}

.datatable.aware = TRUE

verbose = function(...)
{
  if (getOption("MFFProads.debug.verbose") & getOption("MFFProads.debug.progress"))
  {
    message("verbose and progress cannot be both TRUE. progress set to FALSE")
    options(MFFProads.debug.progress = FALSE)
  }

  if (getOption("MFFProads.debug.verbose")) cat(...)
  if (getOption("MFFProads.debug.progress")) cat(".")
}
