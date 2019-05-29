#' Calculates atmospheric pressure from pressure at sea-level
#'
#' @description
#' `pressure.height` applies a simple height adjustment to sea-level pressure data
#'
#' @param psl an optional single numeric value, raster, two-dimensional array or matrix of
#' sea-level atmospheric pressures (Pa)
#' @param dem  a single numeric value, raster, two-dimensional array or matrix of elevations (m)

#' @return a single numeric value, two-dimensional array or matrix of atmospheric pressure (Pa)
#'
#' @details Height adjustments are performed using a simplification of the ideal gas law,
#' assuming 20°C for a standard atmosphere. If psl is not supplied, a typical mean value of sea-level
#' pressure is assumed.
#'
#' @export
#' @examples
#' plot(pressure.height(dem = dtm100m))
#'
pressure.height <- function(psl = 101300, dem) {
  psl <- is_raster(psl)
  psl <- psl / 1000
  z <- is_raster(dem)
  p <- psl * ((293 - 0.0065 * z) / 293)
  p <- p * 1000
  if_raster(p, dem)
}
#' Calculates crop reference evapotranspiration from hourly data
#'
#' @description
#' `cre.hourly` applies the Penman-Monteith equation to derive hourly crop reference evapotranspiration
#'
#' @param Rn a single numeric value, vector, two-dimensional array or matrix of net radiation (MJ / m^2 / hr)
#' @param tc a single numeric value, vector, two-dimensional array or matrix of temperature (deg C)
#' @param u1 a single numeric value, vector, two-dimensional array or matrix of wind speed at one metre above the ground (m /s)
#' @param rh a single numeric value, vector, two-dimensional array or matrix of relative humidity (percentage)
#' @param u1 a single numeric value, vector, two-dimensional array or matrix of wind speed at one metre above the ground (m /s)
#' @param p  a single numeric value, vector, two-dimensional array or matrix of pressure (Pa)
#' @param dn binary variable indicating whether day or night (day - 1, night = 0). Must be
#' a single numeric value or have the same dimensions as Rn
#'
#' @return a single numeric value, vector, two-dimensional array or matrix of crop reference evapotranspiration (mm / hr)
#'
#' @details Applies the method detailed in Allen et al (1998). If relatively fine-resolution spatial data are needed,
#' radiation, temperature, wind and relative humidity can be downscaled using the microclima package. The variable `dn`
#' and humidity conversions can also be obtained using this package. If only sea-level pressure is known, then
#' `p` can be derived using function [pressure.height()]. To convert radiation from watts m^-2 to MJ m^-2 hr^-1 multiply
#' by 0.0036.
#'
#' @export
#' @seealso cre.daily for deriving crop reference evapotranspiration from daily data
#'
#' @examples
#' cre.hourly(1.1, 20, 1.3, 70, 101300, 1) # day
#' cre.hourly(-0.3, 5, 1.3, 100, 101300, 0) # night (cold and wet)
#' cre.hourly(-0.1, 30, 1.3, 60, 101300, 0) # night (hot and dry)
cre.hourly <- function(Rn, tc, u1, rh, p, dn) {
  dn <- Rn * 0 + dn
  pk <- p / 1000
  delta <- 4098 * (0.6108 * exp(17.27 * tc / (tc + 237.3)) / (tc + 237.3)^2)
  G <- 0.1 * Rn
  G0 <- 0.5 * Rn
  G[dn == 0] <- G0[dn == 0]
  Gamma <- 0.000665 * pk
  es <- 0.6108 * exp(17.27 * tc / (tc + 237.3))
  ea <- es * rh / 100
  u2 <- u1 * 4.87 / 4.133245
  cre <- (0.408 * delta * (Rn - G) + Gamma * (37 / (tc + 273)) * u2 * (es - ea))  /
    (delta + Gamma * (1 + 0.34 * u2))
  ifelse(cre < 0, 0, cre)
}
#' Calculates crop reference evapotranspiration from daily data
#'
#' @description
#' `cre.daily` applies the Penman-Monteith equation to derive hourly crop reference evapotranspiration
#'
#' @param Rn a single numeric value, vector, two-dimensional array or matrix of net radiation (MJ / m^2 / day)
#' @param tmin a single numeric value, vector, two-dimensional array or matrix of daily minimum temperature (deg C)
#' @param tmax a single numeric value, vector, two-dimensional array or matrix of daily maximum temperature (deg C)
#' @param u1 a single numeric value, vector, two-dimensional array or matrix of wind speed at one metre above the ground (m /s)
#' @param rhmin a  single numeric value, vector, two-dimensional array or matrix of daily minimum relative humidity (percentage)
#' @param rhmax a  single numeric value, vector, two-dimensional array or matrix of daily maximum relative humidity (percentage)
#' @param p  a single numeric value, vector, two-dimensional array or matrix of pressure (Pa)

#' @return a single numeric value, vector, two-dimensional array or matrix of crop reference evapotranspiration (mm / day)
#'
#' @details Applies the method detailed in Allen et al (1998). If relatively fine-resolution spatial data are needed,
#' radiation, temperature, wind and relative humidity can be downscaled using the microclima package. Humidity conversions
#' can also be performed using this package. If only sea-level pressure is known, then
#' `p` can be derived using function [pressure.height()]. To convert radiation from watts m^-2 to MJ m^-2 day^-1 multiply
#' by 0.0864.
#'
#' @export
#' @seealso cre.hourly for deriving crop reference evapotranspiration from hourly data
#'
#' @examples
#' cre.daily(13.2, 10, 20, 1.3, 60, 100, 101300)
cre.daily <- function(Rn, tmin, tmax, u1, rhmin, rhmax, p) {
  pk <- p / 1000
  tc <- (tmax + tmin) / 2
  delta <- 4098 * (0.6108 * exp(17.27 * tc / (tc + 237.3)) / (tc + 237.3)^2)
  Gamma <- 0.000665 * pk
  esmax <- 0.6108 * exp(17.27 * tmax / (tmax + 237.3))
  esmin <- 0.6108 * exp(17.27 * tmin / (tmin + 237.3))
  es <- (esmax + esmin) / 2
  ea <- (esmin * rhmax / 100 + esmax * rhmin / 100) / 2
  u2 <- u1 * 4.87 / 4.133245
  cre <- (0.408 * delta * Rn + Gamma * (900 / (tc + 273)) * u2 * (es - ea))  /
    (delta + Gamma * (1 + 0.34 * u2))
  cre
}
#' Calculates elevation correction factor to apply to fine-resolution rainfall data
#'
#' @description
#' `raincorrect` Calculates an altitudinal correction factor for application to fine-resolution e.g.
#' daily rainfall derived using bilinear interpolation.
#'
#' @param demc a coarse-resolution raster of digital elevation data
#' @param demf a fine-resolution raster of digital elevation data for the region for which correction factors are required
#' @param rainc a coarse-resolution raster of e.g. monthly or annual rainfall covering the same extent
#' as `demc`
#' @param rainf a fine-resolution raster of e.g. monthly or annual rainfall derived by resampling
#' coarse-resolution rainfall data, covering the same extent as `demf`.

#' @return a raster of correction factors
#'
#' @details Fine-scale e.g. daily or hourly ainfall data derived by resampling or interpolating
#' coarse-resolution data fails to adequately capture elevation effects on rainfall. This
#' function fits a thin plate spline model fitted to coarse-resolution rainfall data with
#' elevation as a predictor. The model is then applied to fine-resolution data and the results
#' compared to those obtained using simple raster resampling using bilinear interpolation.
#' The correction factor can then be applied to e.g. daily rainfall.
#'
#' @importFrom rgcvpack fitTps
#' @importFrom rgcvpack predict.Tps
#' @export
#'
#' @examples
#' rainf <- resample(cornwallrain, dtm100m)
#' plot(mask(rainf, dtm100m)) # rainfall derived using bilinear interpolation
#' cf <- raincorrect(dtm1km, dtm100m, cornwallrain, rainf) # takes ~ 20 seconds to run
#' plot(rainf * cf) # rainfall with correction factor applied
#'
raincorrect <- function(demc, demf, rainc, rainf) {
  xy <- data.frame(xyFromCell(demc, 1:ncell(demc)))
  z <- extract(demc, xy)
  xyz <- cbind(xy, z)
  v <- extract(rainc, xy)
  sel <- which(is.na(v) == F)
  v <- v[is.na(v) == F]
  xyz <- xyz[sel, ]
  tps <- fitTps(xyz, v, m = 2)
  xy <- data.frame(xyFromCell(demf, 1:ncell(demf)))
  z <- extract(demf, xy)
  xyz <- cbind(xy, z)
  v <- extract(demf, xy)
  sel <- which(is.na(v) == F)
  v <- v[is.na(v) == F]
  xyz <- xyz[sel, ]
  xy$z <- NA
  xy$z[sel] <- predict.Tps(tps, xyz)
  rf2 <- rasterFromXYZ(xy)
  cf <- rf2 / rainf
  cf
}
# sequence of rainfalls
# Z higher level (e.g. storm) sum
.propadjust <- function(rainseq, Z) {
  Xs <- rainseq * (Z / sum(rainseq))
  Xs
}
# runs Bartlett-Lewis until sequence of L wet days is generated
# BLest = Bartlett lewis paramaters
# dailysteps = number of values per day (i.e. 24 for hourly)
# dn = number of days
.level0 <- function(BLest, dailyvals, dn) {
  wet <- function(rainseq) {
    wd <- ifelse(sum(rainseq) > 0, 1, 0)
    wd
  }
  sr <- NA
  w <- NULL
  iter <- 0
  while (length(w) < 1) {
    sim<-BLRPM.sim(BLest[1] / dailyvals, BLest[2] / dailyvals, BLest[3] / dailyvals,
                   BLest[4] / dailyvals, BLest[5] / dailyvals, dn * dailyvals * 100,
                   dn * dailyvals * 100, 1, 0)
    hr <- t(matrix(sim$RR, nrow = dailyvals))
    dr <- apply(hr, 1, wet)
    x <- rle(dr)
    w <- which(x$lengths == dn & x$values == 1)
    if (length(w) > 0) {
      v <- length(rep.int(x$values[1:w[1]], x$lengths[1:w[1]])) - dn + 1
      sr <- hr[v:(v + dn -1),]
    }
    iter <- iter + 1
  }
  return(sr)
}
.level1 <- function(rainseq, BLest, dailyvals, dlim, maxiter) {
  d <- dlim * 2
  iter <- 0
  while (d > dlim & iter < maxiter) {
    dn <- length(rainseq)
    l0 <- .level0(BLest, dailyvals, dn)
    if (length(rainseq) > 1) {
      dr <- apply(l0, 1, sum)
    } else dr <- sum(l0)
    d <- sum(log((rainseq + 0.1) / (dr + 0.1))^2)^0.5
    iter <- iter + 1
  }
  if (d > dlim) l0<-NA
  return(l0)
}
.oneday <- function(dayrain, BLest, dailyvals, dlim, maxiter) {
  l1 <- NA
  while (is.na(sum(l1))) {
    l1 <- .level1(dayrain, BLest, dailyvals, dlim, maxiter)
  }
  l1
}
.level3 <- function(rainseq, BLest, dailyvals, dlim, maxiter) {
  l1 <- matrix(NA, nrow = length(rainseq), ncol = 24)
  r <- matrix(rainseq, nrow = 1)
  while (is.na(sum(l1))) {
    for (i in 1:dim(r)[1]) {
      r1 <- r[i,]
      r1 <- r1[is.na(r1) == F]
      st <- 1
      if (i > 1) {
        xx <- r[1:(i-1),]
        xx <- xx[is.na(xx) == F]
        st <- st + length(xx)
      }
      xx <- r[i,]
      xx <- xx[is.na(xx) == F]
      ed <- st + length(xx) - 1
      if (length(r1) == 1) {
        if (is.na(sum(l1[st:ed,]))) l1[st:ed,] <- .oneday(r1, BLest, dailyvals, dlim, maxiter)
      }  else {
        if (is.na(sum(l1[st:ed,]))) {
          l1[st:ed,] <- .level1(r1, BLest, dailyvals, dlim, maxiter)
          if (is.na(sum(l1[st:ed,]))) {
            ii <- which(r1 == min(r1[1:(length(r1) - 1)]))
            r11 <- r1[1:ii]
            r12 <- r1[(ii + 1):length(r1)]
            rr <- matrix(r[-i,], ncol = ncol(r))
            nc <- max(length(r11),length(r12))
            dm <- dim(rr)[1]
            mxc <- apply(rr, 1, function(x) length(x[is.na(x) ==F]))
            if (dm > 0) nc <- max(nc, max(mxc))
            r <- matrix(NA, nrow = dm + 2, ncol = nc)
            if (dm > 0) r[1:dm,1:mxc] <- rr[1:dm,1:mxc]
            r[dm + 1, 1:length(r11)] <- r11
            r[dm + 2, 1:length(r12)] <- r12
            l1a <- .level1(r11, BLest, dailyvals, dlim, maxiter)
            l1b <- .level1(r12, BLest, dailyvals, dlim, maxiter)
            st2 <- st + length(r11)
            l1[st:(st2 -1), ] <- l1a
            l1[st2:(st2 + length(r12) - 1),] <- l1b
          }
        }
      }
    }
  }
  l1
}
#' Estimate sub-daily rainfall from daily rainfall
#'
#' @description
#' `subdailyrain` estimate sub-daily rainfall using Bartlett-Lewis rectangular pulse rainfall model.
#'
#' @param rain a vector time-series of rainfall
#' @param BLpar a data.frame of Bartlett-Lewis parameters as returned by [findBLpar()].
#' @param dailyvals the number of rainfall values required for each day (i.e. 24 for hourly).
#'
#' @return A matrix with `length(rain)` rows and `dailyvals` columns of sub-daily rainfall.
#'
#' @export
#'
#' @details The function is based on the Bartlett-Lewis Rectangular Pulse model described by
#' Rodriguez-Iturbe (1987 & 1988). The model has six parameters (see [findBLpar()]) and is
#' characterized as a particular form of clustering process in which each cluster of rainfall events
#' (hereafter storms) consists of one or more rainfall cells being generated in the start of the
#' process. The parameters of `BLpar` governs the frequency of storms, the start and end of rainfall
#' events associated with each storms, the intensity of rainfall associated with storms variation
#' in the duration of storms, and can be used to generate data for any time-interval. Since these
#' vary seasonally, or by month, it is wise to generate sb-daily data seperately for each month using
#' different parameter estimates.
#'
#' Singificant element sof the coding have been borrowed from from the HyetosMinute package, and
#' the library must be loaded and attached, i.e. `library(HyetosMinute)' as the function calls C++ code
#' included with the package. The package is not available on CRAN and must be obtained or installed
#' directly from here: http://www.itia.ntua.gr/en/softinfo/3/.
#'
#' @references
#' Rodriguez-Iturbe I, Cox DR & Isham V (1987) Some models for rainfall based on stochastic point
#' processes. Proc. R. Soc. Lond., A 410: 269-288.
#'
#' Rodriguez-Iturbe I, Cox DR & Isham V (1988) A point process model for rainfall: Further
#' developments. Proc. R. Soc. Lond., A 417: 283-298.
#'
#' @examples
#' # =========================================== #
#' # ~~~ Generate hourly data for March 2015 ~~~ #
#' # =========================================== #
#' # ~~~~ Get paramaters for March
#' tme <- as.POSIXlt(dailyrain$obs_time)
#' marchrain <- dailyrain$precipitation[which(tme$mon + 1 == 3)]
#' BLpar <- findBLpar(marchrain) # Takes ~ 30 seconds
#' # ~~~~ Generate hourly data for March 2015
#' sel <- which(tme$mon + 1 == 3 & tme$year + 1900 == 2015)
#' march2015 <- dailyrain$precipitation[sel]
#' hourly <- subdailyrain(march2015, BLpar)
#' # ~~~~ Plots comparing hourly and daily / 24 data
#' o <- as.vector(t(matrix(rep(c(1:31), 24), nrow = 31, ncol = 24)))
#' marchhfd <- march2015[o] / 24
#' hourlyv <-as.vector(t(hourly))
#' dd <- c(1:(31 * 24)) / 24
#' plot(hourlyv ~ dd, type = "l", ylim = c(0, max(hourlyv)),
#'      xlab = "Decimal day", ylab = "Rain (mm / hr)", col = "red")
#' par(new = T)
#' plot(marchhfd ~ dd, type = "l", ylim = c(0, max(hourlyv)),
#'      xlab = "", ylab = "", col = "blue", lwd = 2)

#'
subdailyrain <- function(rain, BLest, dailyvals = 24, dlim = 0.2, maxiter = 1000, splitthreshold = 0.2, trace = TRUE) {
  rain[is.na(rain)] <- 0
  srain <- matrix(0, ncol = dailyvals, nrow = length(rain))
  wd <- ifelse(rain > splitthreshold, 1, 0)
  st <- which(diff(c(0, wd)) == 1)
  ed <- which(diff(c(wd, 0)) == -1)
  if (trace) cat(paste("Number of rainfall events:",length(ed),"\n"))
  for (i in 1:length(ed)) {
    r <- rain[st[i]:ed[i]]
    rseq <- .level3(r, BLest, dailyvals, dlim, maxiter)
    for (j in 1:length(r)) rseq[j,] <- .propadjust(rseq[j,], r[j])
    srain[st[i]:ed[i],] <- rseq
    if (trace) cat(paste("Completed rainfall event:",i,"\n"))
  }
  if (trace) cat("Processing days with rain < splitthreshold \n")
  sel <- which(rain <= splitthreshold & rain > 0)
  if (length(sel) > 0) {
    for (i in 1:length(sel)) {
      r <- rain[sel[i]]
      rseq <- .level3(r, BLest, dailyvals, dlim, maxiter)
      srain[sel[i],] <-  .propadjust(rseq, r)
    }
  }
  srain
}
plotrain <- function(daily, subdaily) {
  dailyvals <- length(subdaily) / length(daily)
  d24 <- as.vector(t(matrix(rep(daily, dailyvals), ncol = dailyvals)))
  d24 <- d24 / dailyvals
  day <- c(0:(length(d24) - 1)) / dailyvals
  sday <-  as.vector(t(subdaily))
  xs <- c(day, max(day), 0)
  ys1 <- c(d24, 0, 0)
  ys2 <- c(sday, 0, 0)
  par(mar=c(5,5,5,5))
  plot(NA, xlim = c(0,max(day)), ylim = c(0,max(sday)),
       ylab = "Rainfall (mm)", xlab = "Day", cex.axis = 2, cex.lab = 2)
  polygon(xs, ys1, col = rgb(1,0,0,0.5))
  polygon(xs, ys2, col = rgb(0,0,1,0.5))
}
