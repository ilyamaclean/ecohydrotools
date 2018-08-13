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
  psl <- is.raster(psl)
  psl <- psl / 1000
  z <- is.raster(dem)
  p <- psl * ((293 - 0.0065 * z) / 293)
  p <- p * 1000
  if.raster(p, dem)
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
#' #' @importFrom rgcvpack fitTps
#' #' @importFrom rgcvpack predict.Tps
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


