#' Moving average
#' @noRd
ma <- function(x, n = 3)
{
  if (anyNA(x))x <- zoo::na.approx(x, na.rm = FALSE)
  as.numeric(stats::filter(x, rep(1 / n, n), sides = 2))
}

percent <- function(x)
{
  if (length(x) == 0L) return(0)
  sum(x)/length(x)
}

#' First order derivative
#' @noRd
deriv = function(x,y)
{
  n = length(x)
  dx = x[-(1:2)] - x[-((n-1):n)]
  dy = y[-(1:2)] - y[-((n-1):n)]

  d1 = (y[2]-y[1])/(x[2]-x[1])
  dn = (y[n]-y[n-1])/(x[n]-x[n-1])

  c(d1,dy/dx,dn)
}

#' Generic activation function
#'
#' @param x input vector
#' @param th thresold where the activation function changes
#' @param mode two modes of activation function are necessary in this project.  Linear and stairway
#' @param asc bool. The activation goes from 0 to 1 or from 1 to 0.
#' @noRd
activation <- function(x, th, mode = c("thresholds", "piecewise-linear"), asc = TRUE)
{
  o = x
  x = x[]
  mode <- match.arg(mode)
  y <- numeric(length(x))
  start <- if (asc) 0 else 1
  delta <- if (asc) 1 else -1

  if (mode == "piecewise-linear")
  {
    stopifnot(length(th) == 2)

    a <- delta/(th[2]-th[1])
    b <- start -a*th[1]
    y <- a*x + b
    y[x > th[2]] <- 1 - start
    y[x < th[1]] <- start
  }
  else
  {
    val <- seq(from = 0, to = 1, length.out = length(th) + 1)
    val <- val[-c(1, length(val))]
    if (!asc) val <- rev(val)

    for (i in seq_along(th))
      y[x >= th[i]] <- val[i]

    y[x < th[1]] <- start
    y[x >= th[length(th)]] <- 1-start
  }

  o[] <- round(y,3)
  return(o)
}

activation2 <- function(r)
{
  u <- quantile(r[], probs = seq(0,1, length.out = 12), na.rm = TRUE)
  y <- r
  y[] <- findInterval(y[], u)
  ymin = min(y[], na.rm = TRUE)
  ymax = max(y[], na.rm = TRUE)
  y <- 1-(y-ymin)/(ymax-ymin)
  y
}
