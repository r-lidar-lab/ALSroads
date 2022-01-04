#' Moving average
#' @noRd
ma <- function(x, n = 3)
{
  if (anyNA(x))x <- zoo::na.approx(x, na.rm = FALSE)
  as.numeric(stats::filter(x, rep(1 / n, n), sides = 2))
}

#' Mode of a distribution
#'
#' The mode is the value that appears most often
#'
#' @param x numeric
#' @noRd
stat_mode <- function(x, na.rm = FALSE)
{
  if(na.rm){
    x = x[!is.na(x)]
  }

  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}


#' Moving mode of a distribution
#'
#' Moving mode of a distribution
#'
#' @examples
#' x = c(1,1,1,2,1,1,2,2,1,2,2,4,2,2,4,4,4,2,1,4,4,4,1,4)
#' moving_mode(x, size = 6)
#' @noRd
moving_mode <- function(x, size = 5)
{
  nas <- rep(NA, size)
  y <- integer(length(x))
  x <- c(nas, x, nas)
  for (i in (size+1):(length(x)-size))
  {
    xx <- x[(i-size):(i+size)]
    y[i-size] = stat_mode(xx, na.rm = TRUE)
  }

  return(y)
}

#' Moving mode until stability
#'
#' Moving mode repetitively until it is stable
#'
#' @examples
#' x = c(1,1,2,2,2,1,2,1,1,2,2,1,2,2,4,2,2,3,4,2,4,2,4,2,1,4,3,4,1,4)
#' moving_mode(x, size = 6)
#' stable_moving_mode(x, size = 6)
#' @noRd
stable_moving_mode = function(x, size = 5)
{
  if (length(x) <= 2) return(rep(x[1], length(x)))

  y <- moving_mode(x)
  i <- 1
  while(!all(x == y, na.rm = T) & i < 10)
  {
    x <- y
    y <- moving_mode(x, size = )
    i <- i+1
  }
  return(y)
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

#' Robust peak detection algorithm
#'
#' Detects statistical peaks in a time serie
#'
#' @param lag: the lag parameter determines how much your data will be smoothed and how
#' adaptive the algorithm is to changes in the long-term average of the data. The more
#' stationary your data is, the more lags you should include (this should improve the
#' robustness of the algorithm). If your data contains time-varying trends, you should
#' consider how quickly you want the algorithm to adapt to these trends. I.e., if you put
#' lag at 10, it takes 10 'periods' before the algorithm's treshold is adjusted to any
#' systematic changes in the long-term average. So choose the lag parameter based on the
#' trending behavior of your data and how adaptive you want the algorithm to be.
#' @param influence: this parameter determines the influence of signals on the algorithm's
#' detection threshold. If put at 0, signals have no influence on the threshold, such that
#' future signals are detected based on a threshold that is calculated with a mean and
#' standard deviation that is not influenced by past signals. If put at 0.5, signals have
#' half the influence of normal data points. Another way to think about this is that if you
#' put the influence at 0, you implicitly assume stationarity (i.e. no matter how many signals
#' there are, you always expect the time series to return to the same average over the long term).
#' If this is not the case, you should put the influence parameter somewhere between 0 and 1,
#' depending on the extent to which signals can systematically influence the time-varying trend
#' of the data. E.g., if signals lead to a structural break of the long-term average of the time
#' series, the influence parameter should be put high (close to 1) so the threshold can react to
#' structural breaks quickly.
#' @param threshold: the threshold parameter is the number of standard deviations from the
#' moving mean above which the algorithm will classify a new datapoint as being a signal.
#' For example, if a new datapoint is 4.0 standard deviations above the moving mean and
#' the threshold parameter is set as 3.5, the algorithm will identify the datapoint as a
#' signal. This parameter should be set based on how many signals you expect. For example,
#' if your data is normally distributed, a threshold (or: z-score) of 3.5 corresponds to a
#' signaling probability of 0.00047 (from this table), which implies that you expect a signal
#' once every 2128 datapoints (1/0.00047). The threshold therefore directly influences how
#' sensitive the algorithm is and thereby also determines how often the algorithm signals.
#' Examine your own data and choose a sensible threshold that makes the algorithm signal
#' when you want it to (some trial-and-error might be needed here to get to a good threshold
#' for your purpose).
#'
#' @references
#' https://stackoverflow.com/a/22640362/8442410
#' @examples
#' y <- runif(200, 0, 1)
#' y[55:60] <- runif(7,1,4)
#' y[110:115] <- runif(7,2,3)
#' lag       <- 30
#' threshold <- 5
#' influence <- 0
#' robust_peak_detection(y, lag, threshold, influence, TRUE)
#' @noRd
robust_peak_detection <- function(y, lag, threshold, influence, display = FALSE)
{
  signals <- rep(0,length(y))
  filteredY <- y[0:lag]
  avgFilter <- NULL
  stdFilter <- NULL
  avgFilter[lag] <- mean(y[0:lag], na.rm=TRUE)
  stdFilter[lag] <- stats::sd(y[0:lag], na.rm=TRUE)
  for (i in (lag+1):length(y))
  {
    if (abs(y[i]-avgFilter[i-1]) > threshold*stdFilter[i-1])
    {
      if (y[i] > avgFilter[i-1])
      {
        signals[i] <- 1;
      }
      else
      {
        signals[i] <- -1;
      }
      filteredY[i] <- influence*y[i]+(1-influence)*filteredY[i-1]
    }
    else
    {
      signals[i] <- 0
      filteredY[i] <- y[i]
    }
    avgFilter[i] <- mean(filteredY[(i-lag):i], na.rm=TRUE)
    stdFilter[i] <- stats::sd(filteredY[(i-lag):i], na.rm=TRUE)
  }

  result = list("signals"=signals,"avgFilter"=avgFilter,"stdFilter"=stdFilter)

  if (display)
  {
    ymin <- min(result$avgFilter-threshold*result$stdFilter, y, na.rm = T)
    ymax <- max(result$avgFilter+threshold*result$stdFilter, y, na.rm = T)
    plot(1:length(y),y,type="l",ylab="",xlab="", ylim = c(ymin, ymax))
    graphics::lines(1:length(y),result$avgFilter,type="l",col="cyan",lwd=2)
    graphics::lines(1:length(y),result$avgFilter+threshold*result$stdFilter,type="l",col="green",lwd=2)
    graphics::lines(1:length(y),result$avgFilter-threshold*result$stdFilter,type="l",col="green",lwd=2)
    graphics::lines(result$signals+ymin,type="S",col="red",ylab="",xlab="",ylim=c(-1.5,1.5),lwd=2)
  }

  return(result)
}

#' Robust peak detection algorithm
#'
#' Detects statistical peaks in a time serie like robust_peak_detection() but makes 2 passes
#' in two direction to also be able to detect peaks before the lag
#'
#' @examples
#' y <- runif(200, 0, 1)
#' y[15:20] <- runif(7,1,4)
#' y[110:115] <- runif(7,2,3)
#' lag       <- 30
#' threshold <- 5
#' influence <- 0
#' double_pass_robust_peak_detection(y, lag, threshold, influence, TRUE)
#' @noRd
double_pass_robust_peak_detection = function(y, lag, threshold, influence, display = FALSE)
{
  nas = is.na(y)

  y[nas] = mean(y, na.rm = T)
  ry <- rev(y)

  if (display)
  {
    opar = graphics::par(mfrow = c(2,2),oma = c(2,2,0,0) + 0.1,mar = c(0,0,2,1) + 0.2)
    on.exit(graphics::par(opar))
  }

  result <- robust_peak_detection(y, lag, threshold, influence, display)
  rresult <- robust_peak_detection(ry, lag, threshold, influence, display)

  s1 = result$signals
  s2 = rev(rresult$signals)
  s = s1+s2
  s[nas] <- 0

  if (display)
  {
    plot(1:length(y),y,type="l",ylab="",xlab="", ylim = c(min(-2, y), max(2,y)))
    graphics::lines(s - min(y, na.rm = T), col = "red", lwd = 2)
  }

  return(s)
}

percent <- function(x)
{
  if (length(x) == 0L) return(0)
  sum(x)/length(x)
}

activation <- function(x, th, mode = c("thresholds", "piecewise-linear"), asc = TRUE)
{
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

  return(y)
}
