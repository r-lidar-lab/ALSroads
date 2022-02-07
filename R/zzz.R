.onLoad <- function(libname, pkgname)
{
  op <- options()
  op.ALSroads<- list(
    ALSroads.debug.finding = FALSE,
    ALSroads.debug.measuring = FALSE,
    ALSroads.debug.metrics = FALSE,
    ALSroads.debug.verbose = FALSE,
    ALSroads.debug.progress = TRUE,
    ALSroads.debug = FALSE)

  toset <- !(names(op.ALSroads) %in% names(op))
  if (any(toset)) options(op.ALSroads[toset])
}

.onUnload <- function(libpath)
{
  library.dynam.unload("ALSroads", libpath)
}

.datatable.aware = TRUE

verbose = function(...)
{
  if (getOption("ALSroads.debug.verbose") & getOption("ALSroads.debug.progress"))
  {
    message("verbose and progress cannot be both TRUE. progress set to FALSE")
    options(ALSroads.debug.progress = FALSE)
  }

  if (getOption("ALSroads.debug.verbose")) cat(...)
  if (getOption("ALSroads.debug.progress")) cat(".")
}
