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
#' Internal funcion used to calculate mean rainfall
#'
meanMBLRPM<-function(a,l,v,k,f,mx,h=1) {
  x<-(h*l*mx*v*(1+k/f))/(a-1)
  return(x)
}
#' Internal funcion used to calculate variance of rainfall
#'
varMBLRPM<-function(a,l,v,k,f,mx,h=1) {
  A<-(2*l*(1+k/f)*(mx^2)*(v^a))/((f^2)*((f^2)-1)*(a-1)*(a-2)*(a-3))
  B<-(2*(f^2)-2+k*f)*(f^2)*((a-3)*h*(v^(2-a))-(v^(3-a))+((v+h)^(3-a)))
  C<-k*(f*(a-3)*h*(v^(2-a))-(v^(3-a))+((v+f*h)^(3-a)))
  D<-A*(B-C)
  return(D)
}
#' Internal funcion used to calculate covariance of rainfall
#'
covarMBLRPM<-function(a,l,v,k,f,mx,h=1,lag=1) {
  A<-(l*(1+k/f)*(mx^2)*(v^a))/((f^2)*((f^2)-1)*(a-1)*(a-2)*(a-3))
  B<-(2*(f^2)-2+k*f)*(f^2)*(((v+(lag+1)*h)^(3-a))-2*((v+lag*h)^(3-a))+((v+(lag-1)*h)^(3-a)))
  C<-k*(((v+(lag+1)*h*f)^(3-a))-(2*((v+h*lag*f)^(3-a)))+((v+(lag-1)*h*f)^(3-a)))
  D<-A*(B-C)
  return(D)
}
#' Internal funcion used to calculate probability of being dry
#'
pdrMBLRPM<-function(a,l,v,k,f,h=1) {
  mt<-((1+(f*(k+f))-(0.25*f*(k+f)*(k+4*f))+((f/72)*(k+f)*(4*(k^2)+27*k*f+72*(f^2))))*v)/(f*(a-1))
  G00<-((1-k-f+1.5*k*f+(f^2)+0.5*(k^2))*v)/(f*(a-1))
  A<-(f+(k*(v/(v+(k+f)*h))^(a-1)))/(f+k)
  D<-exp(l*(-h-mt+G00*A))
  return(D)
}
#' Internal objective funcion used to derive Bartlett-Lewis model parameters
#'
objf1 <- function(x, means, vars, covs, pdrs, period = "daily") {
  S1 <- 0
  S6 <- 0
  S12 <- 0
  S48 <- ((meanMBLRPM(x[1], x[2], x[3], x[4], x[5], x[6], 48) / means[5]) - 1)^2 +
    ((varMBLRPM(x[1], x[2], x[3], x[4], x[5], x[6], 48) / vars[5]) - 1)^2 +
    ((covarMBLRPM(x[1], x[2], x[3], x[4], x[5], x[6], 48)/ covs[5]) - 1)^2 +
    ((pdrMBLRPM(x[1], x[2], x[3], x[4], x[5], 48) / pdrs[5]) - 1)^2
  S24 <- ((meanMBLRPM(x[1], x[2], x[3], x[4], x[5], x[6], 24) / means[4]) - 1)^2 +
    ((varMBLRPM(x[1], x[2], x[3], x[4], x[5], x[6], 24) / vars[4]) - 1)^2 +
    ((covarMBLRPM(x[1], x[2], x[3], x[4], x[5], x[6], 24)/ covs[4]) - 1)^2 +
    ((pdrMBLRPM(x[1], x[2], x[3], x[4], x[5], 24) / pdrs[4]) - 1)^2
  if (period == "sixhourly" | period == "hourly") {
    S12 <- ((meanMBLRPM(x[1], x[2], x[3], x[4], x[5], x[6], 12) / means[3]) - 1)^2 +
      ((varMBLRPM(x[1], x[2], x[3], x[4], x[5], x[6], 12) / vars[3]) - 1)^2 +
      ((covarMBLRPM(x[1], x[2], x[3], x[4], x[5], x[6], 12)/ covs[3]) - 1)^2 +
      ((pdrMBLRPM(x[1], x[2], x[3], x[4], x[5], 12) / pdrs[3]) - 1)^2
    S6 <- ((meanMBLRPM(x[1], x[2], x[3], x[4], x[5], x[6], 6) / means[2]) - 1)^2 +
      ((varMBLRPM(x[1], x[2], x[3], x[4], x[5], x[6], 6) / vars[2]) - 1)^2 +
      ((covarMBLRPM(x[1], x[2], x[3], x[4], x[5], x[6], 6)/ covs[2]) - 1)^2 +
      ((pdrMBLRPM(x[1], x[2], x[3], x[4], x[5], 6) / pdrs[2]) - 1)^2
  }
  if (period == "hourly") {
    S1 <- ((meanMBLRPM(x[1], x[2], x[3], x[4], x[5], x[6], 1) / means[1]) - 1)^2 +
      ((varMBLRPM(x[1], x[2], x[3], x[4], x[5], x[6], 1) / vars[1]) - 1)^2 +
      ((covarMBLRPM(x[1], x[2], x[3], x[4], x[5], x[6], 1)/ covs[1]) - 1)^2 +
      ((pdrMBLRPM(x[1], x[2], x[3], x[4], x[5], 1) / pdrs[1]) - 1)^2
  }
  S <- S1 + S6 + S12 + S24 + S48
  if(is.infinite(S)) {S<-10^40}
  if(is.na(S)) {S<-10^40}
  return(S)
}
#' Estimate Bartlett-Lewis model parameters for deriving sub-daily rainfall
#'
#' @description
#' `findBLpar` is used to estimate Bartlett-Lewis model parameters needed by function [subdailyrain()] to derive sub-daily rainfall.
#'
#' @param raindseq a vector time-series of rainfall
#' @param period a desciption of the time-interval between successive values of `rainseq`.
#' One of `hourly`, `sixhourly` or `daily`. Other time intervals are not supported.
#'
#' @return A data.frame with 10 elements. The first six are the Bartlett-Lewis model parameters
#' (see details). Additionally, `mean.rain`: the modelled mean intensity of rainfall in `rainseq`
#' as derived using the parameters of Bartlett-Lewis model; `var.rain`: the modelled variance in
#' rainfall;  `cvar.rain`, the modelled covariance between rainfall and rainfall lagged by one
#' time-step; and  `prop.dry`, the modelled proportion of elements in `rainseq` that have no rain. These are not used for generating sub-daily rainfall, but are used to check that parameter estimates are sensible (see details).
#' @export
#'
#' @details The function is based on the Bartlett-Lewis Rectangular Pulse model described by
#' Rodriguez-Iturbe (1987 & 1988). The model has six parameters and is characterized as a
#' particular form of clustering process in which each cluster of rainfall events (hereafter storms)
#' consists of one or more rainfall cells being generated in the start of the process. The
#' parameter `l` governs the frequency of storms; `k` is used to control to start of rainfall
#' events associated with each storms and `f` is used to is used to control the time at which
#' rainfall ceases; `mx` represents intensity of rainfall associated with storms and `a` and `v`
#' are scale and shape parameters used to control variation in the duration of storms.
#'
#' Once know, the parameters can be used to calculate, for any time-interval, the mean and variance
#' in the intensity of rain, the autocorrelation rainfall intensities for any lag period, and the
#' proportion of dry days or hours. Since these vary seasonally, or by month, normally the model
#' parameters are derived separately for each month using time-series of data covering multiple
#' years (see example 2).
#'
#' These variables can also be calculated analytically from the sequence of rainfall data, but at
#' only at time intervals ≥ that of the data. This function derives the parameters through
#' iteration. In each iteration, 1000 parameters are generated randomly from uniform distributions
#' within finite limits and the mean, variance, covariance and dry proportion estimated and
#' compared to those obtained analytically using an objective function. This comparison is used to
#' narrow the limits from which parameters can be generated and algorithm iterated until
#' convergence.
#'
#' It remains challenging to generate six parameters using a limited amount of information against
#' which to ensure a unique optimum of the objective, and hence convergence of during optimization.
#' Initial parameter limits of chosen carefully to assist with optimisation, but the function is
#' crude, not particularly computationally efficient, and not widely tested. It does, however, work
#' in most circumstances, and the estimated mean, variance, covariance and dry proportion of rainfall
#' events are also given to enable sensible cross-checking. Some element of the coding from the
#' HyetosMinute package, available here under a GNU General Public License:
#' http://www.itia.ntua.gr/en/softinfo/3/, have been used.
#'
#' @references
#' Rodriguez-Iturbe I, Cox DR & Isham V (1987) Some models for rainfall based on stochastic point
#' processes. Proc. R. Soc. Lond., A 410: 269-288.
#'
#' Rodriguez-Iturbe I, Cox DR & Isham V (1988) A point process model for rainfall: Further
#' developments. Proc. R. Soc. Lond., A 417: 283-298.
#'
#' @examples
#' # ========= Example 1: ============ #
#' # Find Bartlett-Lewis model parameters for March rainfall
#' tme <- as.POSIXlt(dailyrain$obs_time)
#' marchrain <- dailyrain$precipitation[which(tme$mon + 1 == 3)]
#' findBLpar(marchrain) # Takes ~ 30 seconds

#' # ========= Example 2: ============ #
#' # Find Bartlett-Lewis model parameters for each month
#' # Warning: Takes a few minutes to run
#' tme <- as.POSIXlt(dailyrain$obs_time)
#' pars <- list()
#' for (mon in 1:12) {
#'   monthrain <- dailyrain$precipitation[which((tme$mon + 1) == mon)]
#'   pars[[mon]] <- findBLpar(marchrain)
#'  print(pars[[mon]])
#' }
#'
findBLpar <- function(rainseq, period = "daily") {
  rainseq[is.na(rainseq)] <- 0
  xmin <- c(1,0,0,0,0,0)
  xmax <- c(20,0.1,20,20,1,15)
  means <- c(0:4) * 0
  vars <- c(0:4) * 0
  covs <- c(0:4) * 0
  pdrs <- c(0:4) * 0
  rain <- list()
  if (period == "hourly") {
    if (length(rainseq)%%48 != 0) {docuemnt
      rainseq <- rainseq[1:(floor(length(rainseq) / 48) * 48)]
    }
    rain[[1]] <- rainseq
    rain6 <- matrix(rainseq, nc = length(rainseq) / 6)
    rain[[2]] <- apply(rain6, 2, sum)
    rain12 <- matrix(rainseq, nc = length(rainseq) / 12)
    rain[[3]] <- apply(rain12, 2, sum)
    rain24 <- matrix(rainseq, nc = length(rainseq) / 24)
    rain[[4]] <- apply(rain24, 2, sum)
    rain48 <- matrix(rainseq, nc = length(rainseq) / 48)
    rain[[5]] <- apply(rain48, 2, sum)
    for (i in 1:5) {
      r <- rain[[i]]
      means[i] <- mean(r)
      vars[i]  <- var(r)
      covs[i] <- cov(r[2:length(r)], r[1:(length(r) -1)])
      pdrs[i] <-length(which(r > 0)) / length(r)
    }
  }
  if (period == "sixhourly") {
    if (length(rainseq)%%8 != 0) {
      warning("Incomplete two-days. Truncating data")
      rainseq <- rainseq[1:(floor(length(rainseq) / 8) * 8)]
    }
    rain[[2]] <- rainseq
    rain12 <- matrix(rainseq, nc = length(rainseq) / 2)
    rain[[3]] <- apply(rain12, 2, sum)
    rain24 <- matrix(rainseq, nc = length(rainseq) / 4)
    rain[[4]] <- apply(rain24, 2, sum)
    rain48 <- matrix(rainseq, nc = length(rainseq) / 8)
    rain[[5]] <- apply(rain48, 2, sum)
    for (i in 2:5) {
      r <- rain[[i]]
      means[i] <- mean(r)
      vars[i]  <- var(r)
      covs[i] <- cov(r[2:length(r)], r[1:(length(r) -1)])
      pdrs[i] <-length(which(r > 0)) / length(r)
    }
  }
  if (period == "daily") {
    if (length(rainseq)%%2 != 0) {
      warning("Incomplete two-days. Truncating data")
      rainseq <- rainseq[1:(floor(length(rainseq) / 2) * 2)]
    }
    rain[[4]] <- rainseq
    rain48 <- matrix(rainseq, nc = length(rainseq) / 2)
    rain[[5]] <- apply(rain48, 2, sum)
    for (i in 4:5) {
      r <- rain[[i]]
      means[i] <- mean(r)
      vars[i]  <- var(r)
      covs[i] <- cov(r[2:length(r)], r[1:(length(r) -1)])
      pdrs[i] <-length(which(r > 0)) / length(r)
    }
  }
  max.iter <- 2000
  iter <- 1
  repeat {
    pop <- mat.or.vec(nr = 1000, nc = 6)
    for (i in 1:6)
      pop[, i] <- runif(n = 1000, min = xmin[i], max = xmax[i])
    fpop <- apply(pop, 1, objf1, means, vars, covs, pdrs, period)
    o<-order(fpop)
    pop<-pop[o[1:200],]
    xmin <- apply(pop,2,min)
    xmax <- apply(pop,2,max)
    dif <- xmax - xmin
    if (iter > max.iter  | max(dif) < 0.001) {break}
    iter <- iter + 1
  }
  x <- (xmin + xmax) / 2
  x <- data.frame(matrix(x, nr = 1, nc = 6))
  names(x) <- c("a", "l", "v", "k", "f", "mx")
  x$mean.rain <- as.numeric(meanMBLRPM(x[1], x[2], x[3], x[4], x[5], x[6], h = 24))
  x$var.rain <- as.numeric(varMBLRPM(x[1], x[2], x[3], x[4], x[5], x[6], 24))
  x$cvar.rain <- as.numeric(covarMBLRPM(x[1], x[2], x[3], x[4], x[5], x[6], 24))
  x$prop.dry <- as.numeric(pdrMBLRPM(x[1], x[2], x[3], x[4], x[5], 24))
  x$l <- x$l * 24
  x$mx <- x$mx * 24
  x$v <- x$v / 24
  return(x)
}
#' Internal funcion used by [subdailyrain()]
#'
#' @import HyetosMinute
#'
IMlevel2 <- function (L, l, f, k, a, v, mx, Zhist, subtimescale) {
  sx <- mx
  TotalNumRep <- 1000
  TotalNumOfRepeat <- 0
  if (L != 1) {
    repeat {
      q <- level0(L = L, l = l, f = f, k = k, a = a, v = v,
                  MaxTotalRep = TotalNumRep)
      if (q$NumOfLevel0Rep >= 10 * TotalNumRep) {
        #warning("Rainday cluster split")
        i <- which(Zhist == min(Zhist[1:(length(Zhist) - 1)]))[1]
        L1 <- i
        L2 <- L - i
        Zhist1 <- Zhist[1:i]
        Zhist2 <- Zhist[-(1:i)]
        r <- list(c(L1, Zhist1), c(L2, Zhist2))
        return(r)
        break
      }
      else {
        b <- level1(Zhist = Zhist, q = q, mx = mx, sx = sx,
                    TFweib = FALSE, iot = NA, da = 0.1, L = L,
                    F = 20, MinNumOfLevel1Rep = 50,
                    MaxNumOfTotalRep = TotalNumRep - TotalNumOfRepeat)
        TotalNumOfRepeat <- TotalNumOfRepeat + b$NumOfLevel1Rep
        if (b$LogDist <= 0.1) {
          if (subtimescale == 1) {
            Synth.Depth <- discr.hourly(cellorigin = q[["Sequence"]]$CellOrigins,
                                        cellend = q[["Sequence"]]$CellEnds, cellint = b$CellIntensities,
                                        L = L)
          }
          else {
            Synth.Depth <- discr.subhourly(cellorigin = q[["Sequence"]]$CellOrigins,
                                           cellend = q[["Sequence"]]$CellEnds, cellint = b$CellIntensities,
                                           L = L, timescale = 60 * subtimescale)
          }
          AdjustedHourlyDepth <- adjusting.proc(Zhist = Zhist,
                                                X = Synth.Depth[1:(24/subtimescale * L)],
                                                Z = b$DailyDepth, tms = subtimescale)
          r <- list(AdjustedHourlyDepth)
          return(r)
          break
        }
        if (TotalNumOfRepeat >= TotalNumRep) {
          x <- 1:L
          i <- sample(x[x < L], size = 1)
          L1 <- i
          L2 <- L - i
          Zhist1 <- Zhist[1:i]
          Zhist2 <- Zhist[-(1:i)]
          r <- list(c(L1, Zhist1), c(L2, Zhist2))
          return(r)
          break
        }
      }
    }
  }
  else {
    repeat {
      q <- level0(L = L, l = l, f = f, k = k, a = a, v = v,
                  MaxTotalRep = TotalNumRep)
      if (q$NumOfLevel0Rep >= 10 * TotalNumRep) {
        cat("Impossible to form a wet day \n")
        cat("NA values will be returned for this day \n")
        AdjustedHourlyDepth <- rep(NA, L * (24/subtimescale))
        r <- list(AdjustedHourlyDepth)
        return(r)
        break
      }
      else {
        b <- level1(Zhist = Zhist, q = q, mx = mx, sx = sx,
                    TFweib = FALSE, iot = NA, da = 0.1, L = L,
                    F = 20, MinNumOfLevel1Rep = 50,
                    MaxNumOfTotalRep = TotalNumRep - TotalNumOfRepeat)
        TotalNumOfRepeat <- TotalNumOfRepeat + b$NumOfLevel1Rep
        if (b$LogDist <= 0.1) {
          if (subtimescale == 1) {
            Synth.Depth <- discr.hourly(cellorigin = q[["Sequence"]]$CellOrigins,
                                        cellend = q[["Sequence"]]$CellEnds, cellint = b$CellIntensities,
                                        L = L)
          }
          else {
            Synth.Depth <- discr.subhourly(cellorigin = q[["Sequence"]]$CellOrigins,
                                           cellend = q[["Sequence"]]$CellEnds, cellint = b$CellIntensities,
                                           L = L, timescale = 60 * subtimescale)
          }
          AdjustedHourlyDepth <- adjusting.proc(Zhist = Zhist,
                                                X = Synth.Depth[1:(24/subtimescale * L)],
                                                Z = b$DailyDepth, tms = subtimescale)
          r <- list(AdjustedHourlyDepth)
          return(r)
          break
        }
        if (TotalNumOfRepeat >= 10 * TotalNumRep) {
          cat("Impossible to obtain a sequence of cell intensities that resembles the real one \n")
          cat("Possibly too high daily rainfall depth \n")
          cat("NA values will be returned for this day \n")
          AdjustedHourlyDepth <- rep(NA, L * (24/subtimescale))
          r <- list(AdjustedHourlyDepth)
          return(r)
          break
        }
      }
    }
  }
}
#' Internal funcion used by [subdailyrain()]
#'
#' @import HyetosMinute
#'
IMHyetosRepScheme <- function (L, l, f, k, a, v, mx, Zhist, subtimescale) {
  i <- IMlevel2(L = L, l = l, f = f, k = k, a = a, v = v, Zhist = Zhist, mx = mx,
                subtimescale = subtimescale)
  if (length(i) != 1) {
    repeat {
      i <- lapply(i, function(o) if (!is.list(o)) {
        list(o)
      } else o)
      w <- numeric()
      for (j in 1:length(i)) {w <- c(w, i[[j]])}
      i <- w
      i <- lapply(i, function(r) if (r[1] != length(r) - 1 || is.na(r[1])) {
        r
      } else {
        IMlevel2(L = r[1], l = l, f = f, k = k, a = a,
                 v = v, Zhist = r[2:(r[1] + 1)], mx = mx, subtimescale = subtimescale)
      })
      if (length(unlist(i)) == L * 24/subtimescale) {break}
    }
  }
  i <- as.vector(unlist(i))
  return(i)
}
#' Estimate sub-daily rainfall from daily rainfall
#'
#' @description
#' `subdailyrain` estimate sub-daily rainfall using Bartlett-Lewis rectangular pulse rainfall model.
#'
#' @param rain a vector time-series of rainfall
#' @param BLpar a data.frame of Bartlett-Lewis as returned by [findBLpar()].
#' @param dailyvals the number of rainfall values required for each day (i.e. 24 for hourly).
#'
#' @return A matrix with `length(rain)` rows and `dailyvals` columns of sub-daily rainfall.
#'
#' @import HyetosMinute
#' @export
#'
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
subdailyrain <- function (rain, BLpar, dailyvals = 24) {
  rain[is.na(rain)] <- 0
  rd <- ifelse(rain > 0 , 1, 0)
  st <- c(1, which(rd[2:length(rd)] == 1 &
                     rd[1:(length(rd) - 1)] == 0) + 1)
  if (rd[1] == 0) st <- st[-1]
  ed <- c(which(rd[1:(length(rd) - 1)] == 1 &
                  rd[2:length(rd)] == 0), length(rd))
  if (rd[length(rd)] == 0) ed <- ed[-length(ed)]
  allclusters <-  vector("list", length(st))
  for (i in 1:length(st))
    allclusters[[i]] <- rain[st[i]:ed[i]]
  N <- length(allclusters)
  disaglist <- vector("list", N)
  subdaily <- array(0, dim = c(length(rain), dailyvals))
  for (i in 1:N) {
    clustN <- allclusters[[i]]
    disagdata <- IMHyetosRepScheme(L = length(clustN), l = BLpar$l, f = BLpar$f, k = BLpar$k,
                                   a = BLpar$a, v = BLpar$v, mx = BLpar$mx, Zhist = clustN, subtimescale = 1 / dailyvals * 24)
    subdaily[st[i]:ed[i],] <- t(array(disagdata, dim = c(dailyvals, length(clustN))))
  }
  if (is.na(sum(subdaily)))
    warning("NAs in dataset. Some wet days NOT disaggregated successfully")
  return(subdaily)
}

