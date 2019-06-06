#' Gets van Genuchten soil paramaters
#'
#' @description
#' `getsoilparams` Obtains van Genuchten soil paramaters for a given soil type from an
#' internal dataset
#'
#' @param soiltype a single character, matrix or raster of characters aslisted in [soilparams()]
#' @return a list of van Genuchten parameters (see details)
#'
#' @export
#' @seealso [soilparams()]
#' @details  This function derives the following van Genuchten paramaters used to estimate
#' hydrualic conductivity and diffusivity (see van Genuchten 1980 Soil Science Society of
#' America Journal 44:892-898). The returned list contains the following elements:
#' `n` the pore size distribution parameter (dimensionless, > 1), `alpha` van Genuchten
#' shape paramater (cm^-1), `Ksat` saturated hydraulic conductivity (cm / day), `Smin`
#' Residual water content (cm^3 / cm^3), `Smax` Water content at saturation (cm^3 / cm^3).
#' Each element of the list is of the same class as `soiltype`.
#'
#' @examples
#' getsoilparams("Loam")
#'
getsoilparams <- function(soiltype) {
  data(soilparams)
  r <- soiltype
  soiltype <- is_raster(soiltype)
  u <- unique(as.vector(soiltype))
  u <- u[is.na(u) == F]
  n <- soiltype
  alpha <- soiltype
  Ksat <- soiltype
  Smin <- soiltype
  Smax <- soiltype
  for (i in 1:length(u)) {
    sel <- which(soilparams$Soil.type == u[i])
    s1 <- which(soiltype == u[i])
    n[s1] <- soilparams$n[sel]
    alpha[s1] <- soilparams$alpha[sel]
    Ksat[s1] <- soilparams$Ksat[sel]
    Smin[s1] <- soilparams$Smin[sel]
    Smax[s1] <- soilparams$Smax[sel]
  }
  if (length(as.vector(soiltype)) == 1) {
    lo <- list(n = as.numeric(n),
               alpha = as.numeric(alpha),
               Ksat = as.numeric(Ksat),
               Smin = as.numeric(Smin),
               Smax = as.numeric(Smax))
  } else {
    lo <- list(n = if_raster(array(as.numeric(n), dim = dim(soiltype)), r),
               alpha = if_raster(array(as.numeric(alpha), dim = dim(soiltype)), r),
               Ksat = if_raster(array(as.numeric(Ksat), dim = dim(soiltype)), r),
               Smin = if_raster(array(as.numeric(Smin), dim = dim(soiltype)), r),
               Smax = if_raster(array(as.numeric(Smax), dim = dim(soiltype)), r))
  }
  return(lo)
}

#' Hydraulic conductivity
#'
#' @description
#' `Ktheta` calculates hydrualic conductivity
#'
#' @param sm Soil water fraction
#' @param sp a list of van Genuchten soil parameters as returned by [getsoilparams()]
#' @param timestep number of seconds in each time step of model run (default one day)
#' @return Hydraulic conductivity (cm / time-step)
#'
#' @seealso [soilparams()]
#' @details  Applies the van Genuchten equation to estimate the hydrualic conductivity of
#' unsaturated soil (van Genuchten 1980 Soil Science Society of America Journal 44:892-898).
#' Paramaters of the van Genuchten equation can be found in the dataset [soilparams()]. Ksat
#' and the returned value are adjusted so that returned units are in cm per `timestep`.
#'
#' @export
#' @examples
#' sm <- c(74:422) / 1000
#' plot(Ktheta(sm, getsoilparams("Loam")) ~ sm)
#'
Ktheta <- function(sm, sp, timestep = 86400) {
  n <- sp$n
  Ksat <- sp$Ksat
  Smin <- sp$Smin
  Smax <- sp$Smax
  Ksats <- (Ksat * timestep) / (24 * 3600)
  theta <- (sm - Smin) / (Smax - Smin)
  m <- 1 - 1/n
  Btheta <- theta^0.5 * (1 - (1 - theta^(1/m))^m)^2 * Ksats
  Btheta
}
#' Hydraulic diffusivity
#'
#' @description
#' `Dtheta` calculates hydrualic diffusivity
#'
#' @param sm Soil water fraction
#' @param sp a list of van Genuchten soil parameters as returned by [getsoilparams()]
#' @param timestep number of seconds in each time step of model run (default one day)
#' @return Hydraulic conductivity (cm / time-step)
#'
#' @seealso [soilparams()]
#' @details  Applies the van Genuchten equation to estimate the hydrualic diffusivity
#' (van Genuchten 1980 Soil Science Society of America Journal 44:892-898).
#' Paramaters of the van Genuchten equation can be found in the dataset [soilparams()].
#' Ksat and the returned value are adjusted so that returned units are in cm per `timestep`.
#' `Dtheta` is contrained to be finite by ensuring that `sm` marginally less than `Smax`
#' when saturated and is assumed to be zero, when `sm` = `Smin`.
#'
#'
#' @export
#' @examples
#' sm <- c(74:422) / 1000
#' plot(log(Dtheta(sm, getsoilparams("Loam")), 10) ~ sm)
#'
Dtheta <- function(sm, sp, timestep = 86400) {
  n <- sp$n
  alpha <- sp$alpha
  Ksat <- sp$Ksat
  Smin <- sp$Smin
  Smax <- sp$Smax
  Ksats <- (Ksat * timestep) / (24 * 3600)
  theta <- (sm - Smin) / (Smax - Smin)
  theta <- ifelse(theta > 0.999, 0.999, theta)
  m <- 1 - 1/n
  Btheta <- ((1 - m) * Ksats) / (alpha * m * (Smax - Smin)) *
    theta^(0.5 - 1 / m) * ((1 - theta^(1/m)) ^ -m +
                             (1 - theta^(1/m))^m - 2)
  Btheta[is.na(Btheta)] <- 0
  Btheta
}
#' Adjusts run-off curve number by antecedent soil moisture
.SCN_adjust <- function(scn, sm, Smin, Smax) {
  theta <- (sm - Smin) / (Smax - Smin)
  scn2 <- scn / 100
  scn2 <- ifelse(theta <= 0.25,  0.242536 * scn2 + 0.678166 * scn2^2 + 0.011486, scn2)
  scn2 <- ifelse(theta >= 0.75,  1.61044 * scn2 - 0.69114 * scn2^2 + 0.06730, scn2)
  scn2 <- scn2 * 100
  scn2 <- ifelse(scn == 100, scn, scn2)
  scn2
}
#' Calculates runoff
#'
#' @description
#' `runoff` calculates direct runoff based on USDA Run-off curve number
#'
#' @param rain rainfall in mm
#' @param cn run-off curve number (see details)
#' @param sm Soil water fraction
#' @param Smin Residual water content (cm^3 / cm^3)
#' @param Smax Volumetric water content at saturation (cm^3 / cm^3)
#' @return direct runoff in mm
#'
#' @details  Applies the United States Department of Agriculture method for
#' calculating direct runoff. Runoff curve numbers can be obtained from here:
#' https://www.wcc.nrcs.usda.gov/ftpref/wntsc/H&H/training/runoff-curve-numbers1.pdf
#' Adjustemnts for antecedent moisture condition are automatically applied, so the
#' supplied value of `cn` should be for average soil moisture (AMC II). `Smin` and `Smax`
#' values for soil types are given in [soilparams()].
#' @export
#'
#' @examples
#' rain <- c(1:150)
#' ro1 <- runoff(rain, cn = 50, sm = 0.4, Smin = 0.074, Smax = 0.422)
#' ro2 <- runoff(rain, cn = 90, sm = 0.4, Smin = 0.074, Smax = 0.422)
#' plot(ro1 ~ rain, type = "l")
#' plot(ro2 ~ rain, type = "l")
runoff <- function(rain, cn, sm, Smin, Smax) {
  cn2 <- .SCN_adjust(cn, sm, Smin, Smax)
  P <- 0.0393701 * rain
  S <- (1000 / cn2) - 10
  Ia <- 0.2 * S
  Q <- ifelse(P > Ia, (P - Ia)^2 / (P - Ia + S), 0)
  Qmm <- Q / 0.0393701
  Qmm
}
#' Soil evapotranspiration
#' @details
#' `soilevap` calculates soil evapotranspiration from crop reference evapotranspiration.
#'
#' @param cre Crop reference evapotranspiration as returned by [cre.hourly()] or [cre.daily()]
#' @param sm Soil water fraction
#' @param cover Fraction of vegetation cover
#' @param topratio Ratio of root water uptake from top layer relative to bottom layer
#' (1 = all from top layer, 0 = all form bottom layer)
#' @return a list of with two elements: `top` evaporation and evapotranspiration from
#' the top soil layer, `btm` evapotranspiration from the bottom layer. Units are as for `CRE`
#' @export
#'
#' @seealso [cre.hourly()] [cre.daily()]
#'
#' @details As the soil water fraction decreases below the saturation value, evaporation from
#' bare soil is usually assumed to continue at the potential rate until the water content
#' decreases below a critical value. Thereafter, evaporation is assumed to decrease linearly
#' with decreasing water content, vanishing to zero when water content is low. The paramater
#' `cover` aportions the ratio of evapotranspiration to bare soil evaporation. The paramater
#' `topratio` aportions evapotranspiration between the top and bottom layer and depends on
#' the chosen soil layer depths and root depths.
#'
#' @examples
#' sm <- c(0:100) / 100
#' plot(soilevap(1, sm, cover = 1, surface = 1)$top ~ sm, type = "l")
#' plot(soilevap(1, sm, cover = 0, surface = 1)$top ~ sm, type = "l")
soilevap <- function(cre, sm, cover  = 0.5, topratio = 0.5) {
  rat <- 1.08723 * sm ^ 0.56007 + 0.1
  rat <- ifelse(rat > 1, 1, rat)
  se <- rat * cre
  top <- (1 - cover) * se + cover * topratio * cre
  btm <- cover * (1 - topratio) * cre
  return(list(top = top, btm = btm))
}
.onehydro <- function(rain, evap, surface, theta1, theta2, z1, z2, n, alpha, Ksat,
                      Smin, Smax, cn, cover, topratio, timestep = 86400, n2, Ksat2) {
  sp <- list(n = n, alpha = alpha, Ksat = Ksat, Smin = Smin, Smax = Smax)
  sp2 <- list(n = n2, alpha = alpha, Ksat = Ksat2, Smin = Smin, Smax = Smax)
  # add surplus to rainfall
  su2 <- ifelse(theta2 > Smax, z2 * (theta2 - Smax), 0)
  theta1 <- theta1 + su2 / z1
  su1 <- ifelse(theta1 > Smax, z1 * 10 * (theta1 - Smax), 0)
  surface <- surface + su1
  rain <- rain + surface
  # Run model
  thetab <- ifelse(theta2 > theta1, theta2, theta1)
  Dth <- Dtheta(thetab, sp, timestep)
  Kth <- Ktheta(thetab, sp, timestep)
  c1 <- which(theta2 > theta1)
  Kth[c1] <- -Kth[c1]
  c1 <- which(theta2 == theta1)
  Kth[c1] <- 0
  evp <- soilevap(evap, theta1, cover, topratio)
  evp1 <- evp$top / 10
  evp2 <- evp$btm / 10
  ro <- runoff(rain, cn, theta1, Smin, Smax)
  In <- (rain - ro) / 10
  # Delta due to hydrualic conductivity and diffusivity (cm)
  xx <- 2 * Dth * ((theta1 - theta2) / (z1 + z2)) + Kth
  # Constrain xx to do no more than equalise water content of layers
  xmx <- abs(theta1 - theta2) * z1 * (z2 / (z1 + z2))
  xx <- ifelse(xx > xmx, xmx, xx)
  xx <- ifelse(xx < -xmx, -xmx, xx)
  delta.sm1 <- In - evp1 - xx
  delta.sm2 <- Ktheta(theta2, sp2, timestep) - evp2  + xx
  delta.sm1 <- delta.sm1 / z1
  delta.sm2 <- delta.sm2 / z2
  theta1 <- theta1 + delta.sm1
  theta2 <- theta2 + delta.sm2
  # surplus
  su2 <- ifelse(theta2 > Smax, z2 * (theta2 - Smax), 0)
  theta1 <- theta1 + su2 / z1
  su1 <- ifelse(theta1 > Smax, z1 * 10 * (theta1 - Smax), 0)
  ro <- ro + su1
  theta1 <- ifelse(theta1 > Smax, Smax, theta1)
  theta2 <- ifelse(theta2 > Smax, Smax, theta2)
  theta1 <- ifelse(theta1 < Smin, Smin, theta1)
  theta2 <- ifelse(theta2 < Smin, Smin, theta2)
  return(list(theta1 = theta1, theta2 = theta2, surplus = ro))
}
#' Implementation of Mahrt and Pan two-layer model of soil hydrology
#'
#' @description `hydromodel_time` is a time-series implementation of the Mahrt and Pan two-layer
#' model of soil hydrology for a point location.
#'
#' @param rain a vector of  rainfall (mm per time step)
#' @param cre a vector of crop reference evapotranspirations (mm per timestep) as returned by
#' [cre.daily()] or [cre.hourly()]
#' @param sp a list of van Genuchten soil parameters as returned by [getsoilparams()]
#' @param theta1 initial soil water fraction of top layer
#' @param theta2 initial soil water fraction of bottom layer
#' @param surface surface water depth (mm)
#' @param add.surplus logical indicating whether to add surface water to rainfall in each time-step (see details)
#' @param z1 assumed depth of top soil layer 1 (cm)
#' @param z2 assumed depth of bottom soil layer 2 (cm)
#' @param cn Runoff curve number (see details)
#' @param cover Fraction of vegetation cover
#' @param topratio ratio of root water uptake from top layer relative to bottom layer
#' @param timestep number of seconds in each time step of model run (default one day)
#' @param n2 optional pore size distribution parameter for controlling
#' ground-water seepage (see details)
#' @param Ksat2 optional saturated hydraulic conductivity parameter for controlling
#' ground-water seepage (see details)
#'
#' @return a list of three vectors: `sm1` a vector of soil water fractions in the top soil layer,
#' `sm2` a vector of soil water fractions in the bottom soil layer, `surplus` a vector of
#' direct run-off (mm)
#' @export
#' @seealso [soilparams()] [Dtheta()] [Ktheta()] [runoff()] [soilevap()]
#' @details Implementation of Mahrt and Pan (1984) two-layer model of soil hydrology for applications
#'  where only computer time and complexity are restricted. Volumetric soil water is
#'  computed in a thin upper layer for use in calculation of surface evaporation. Water
#'  storage is computed for an underlying deeper layer. Precipitation enters the top soil layer,
#'  but any precipitation that cannot infiltrate or re-evaporate is specified to be runoff. The rate and
#'  direction of exchange of water between the soil layers is determined by soil moisture-dependent
#'  hydraulic diffusivity and conductivity and by the difference in soil moisture between the two layers.
#'  Evapotranspiration is assumed to occur from the vegetated portion of the basin, and can be partioned
#'  between both layers.
#'
#'  Soil moisture-dependent hydraulic diffusivity and conductivity are
#'  calculated using the van Genuchten equation (van Genuchten 1980 Soil Science Society of America Journal
#'  44:892-898), and the parameters `n`, `alpha`, `Ksat`, `Smax` and `Smin` (see [getsoilparams()]) control the relationship
#'  with soil moisture. The van Genuchten parameters for major soil types are shown in soilparams.
#'  The paramaters `n2` and `Ksat2` are equivelent to `n` and `Ksat`, but control ground water
#'  seepage rates. Setting `Ksat2` to zero assumes ground water seepage. Paramaters for the
#'  van Genuchten equation should be in units of centimeters and days. Correction for shorter time-scales
#'  is automatically applied.
#'
#'  For simplicity, the diffusion of soil water across the bottom of the lower layer is neglected since gradients
#'  and resulting fluxes at this depth are generally unimportant, except over longer time-scales.
#'  Exceptions include a high water table and the advance of a deep ‘wetting front’. In the
#'  latter case, the thickness of the model should be increased.
#'
#'  Surface runoff is dependent on the porosity of the soil and is controlled by a run-off
#'  curve number, itself dependent on soil type, land cover and soil condition. Runoff curve numbers can be
#'  obtained from here: https://www.wcc.nrcs.usda.gov/ftpref/wntsc/H&H/training/runoff-curve-numbers1.pdf.
#'  If add.surplus is set to false (the default), runoff is assumed to leave the basin. If set to true,
#'  it is added to rainfall at each time step.
#'
#'
#' @examples
#' # Rainfall and evapotranspiration time series
#' data(rainevap)
#' rain <- rainevap$rain
#' cre <- rainevap$evap
#' # Get van Genuchten soil paramaters
#' sp <- getsoilparams("Clay loam")
#' # Run model for clay-loam with run-off curve number for fair condition pasture
#' hm <- hydromodel_time(rain, cre, sp = sp, cn = 84, cover = 1)
#' # ============== Rainfall and run-off ==================== #
#' plot(rain, type = "l", lwd = 1, col = "blue", ylim = c(0,40))
#' par(new=T)
#' plot(hm$sur, type = "l", lwd = 1, xlab = "", ylab = "", ylim = c(0,40))
#' # Soil moisture: top layer
#' plot (hm$sm1, type = "l", ylim = c(0, 0.5))
#' # Soil moisture: bottom layer
#' plot (hm$sm2, type = "l", ylim = c(0, 0.5))
#'
#'
hydromodel_time <- function(rain, cre, sp, theta1 = 0.35, theta2 = 0.35, surface = 0, add.surplus = F,
                            z1 = 5, z2 = 95, cn = 82, cover = 0.8, topratio = 0.5, timestep = 86400,
                            n2 = 1.1, Ksat2 = 0) {
  if (class(rain) != "integer" & class(rain) != "numeric") {
    stop("rain must a vector")
  }
  if (class(cre) != "integer" & class(cre) != "numeric") {
    stop("evap must a vector")
  }
  if(length(rain) != length(cre)) {
    stop("rain and cre miust have same length")
  }
  n <- sp$n
  alpha <- sp$alpha
  Ksat <- sp$Ksat
  Smin <- sp$Smin
  Smax <- sp$Smax
  sm1 <- 0
  sm2 <- 0
  sur <- surface
  for (i in 1:length(rain)) {
    if (add.surplus) rain[i] <- rain[i] + sur[i]
    oh <- .onehydro(rain[i], cre[i], surface, theta1, theta2, z1, z2, n, alpha, Ksat,
                    Smin, Smax, cn, cover, topratio, timestep = 86400, n2, Ksat2)
    sm1[i] <- oh$theta1
    sm2[i] <- oh$theta2
    theta1 <- sm1[i]
    theta2 <- sm2[i]
    sur[i + 1] <- oh$surplus
  }
  return(list(sm1= sm1, sm2 = sm2, surplus = sur))
}
#' Implementation of Mahrt and Pan two-layer model of soil hydrology: spatial
#'
#' @description `hydromodel_time` is a spatial implementation of the Mahrt and Pan two-layer
#' model of soil hydrology for a point location.
#'
#' @param dem a raster or matrix of digital elevation data
#' @param rain a vector or 3-dimensional array ofrainfall (mm per time step)
#' @param cre a vector or 3-dimensional array of crop reference evapotranspirations (mm per timestep) as returned by
#' [cre.daily()] or [cre.hourly()].
#' @param sp van Genuchten soil paramaters as returned by [getsoilparams()]
#' @param basins optional raster or matrix of hydrological basins as returned by `basindelin`.
#' Calculated if not supplied.
#' @param theta1 initial soil water fraction of top layer
#' @param theta2 initial soil water fraction of bottom layer
#' @param surface initial surface water depth (mm)
#' @param p1 a single numeric value of the power adjustment to apply to topographic wetness index values  for surface layer (see [topdist()])
#' @param p2 a single numeric value of the power adjustment to apply to topographic wetness index values  for sub-surface layer (see [topdist()])
#' @param p3 a single numeric value of the power adjustment to apply to topographic wetness index values  for surface water (see [topdistw()])
#' @param z1 a single numeric value, raster or matrix of assumed depth(s) of the top soil layer 1 (cm)
#' @param z2 a single numeric value, raster or matrix of assumed depth(s) of bottom soil layer 2 (cm)
#' @param cn a single numeric value, raster or matrix of runoff curve number(s) (see details)
#' @param cover a single numeric value, raster or matrix of fractional of vegetation cover
#' @param topratio a single numeric value, raster or matrix of ratio(s) of root water uptake from top layer relative to bottom layer
#' @param timestep number of seconds in each time step of model run (default one day)
#' @param n2 optional single numeric value, raster or matrix of pore size distribution parameter(s) for controlling
#' ground-water seepage (see details)
#' @param Ksat2 optional single numeric value, raster or matrix of saturated hydraulic conductivity parameter(s) for controlling
#' ground-water seepage (see details)
#' @param Trace logical indicating whether to produce plots tracking progress
#' @return a list of three arrays: `sm1` soil water fractions in the top soil layer,
#' `sm2` soil water fractions in the bottom soil layer, `surface` surface water depth (mm)
#' @import dplyr raster stringr
#' @export
#' @seealso [getsoilparams()] [Dtheta()] [Ktheta()] [runoff()] [soilevap()] [topdist()] [topdistw()] [hydromodel_time()]
#'
#' @details Modified spatial implementation of Mahrt and Pan (1984) two-layer model of soil hydrology. FRactional Soil water
#' content is computed in a thin upper layer for use in calculation of surface evaporation. Water
#' storage is computed for an underlying deeper layer. Precipitation enters the top soil layer,
#' but any precipitation that cannot infiltrate or re-evaporate is specified to be runoff. The rate and
#' direction of exchange of water between the soil layers is determined by soil moisture-dependent
#' hydraulic diffusivity and conductivity and by the difference in soil moisture between the two layers.
#' Evapotranspiration is assumed to occur from the vegetated portion of the basin, and can be partioned
#' between both layers.
#'
#' In this variant of the model, the spatial distributions of soil and surface water are explicitely accounted
#' for. Hydrological basins are delineated from the dem using [basindelin()] and the basin volumes, pour points and
#' basins to which overflow would accrue are calculated. Starting with basin with the highest elevation pour point, the model
#' is run iteratively for each basin. Within each time-step soil and surface water are distributed across the basin
#' by topographic wetness index values calculated using [topidx()]. As such, run-off, hydraulic conductivity and diffusivity
#' are permitted to vary within the basin. Surplus surface water remains within the basin unless the basin volume is exceeded,
#' in which case it is accrued to adjoining basin at the pour point.
#'
#' Soil moisture-dependent hydraulic diffusivity and conductivity are calculated using the van Genuchten equation
#' (van Genuchten 1980 Soil Science Society of America Journal 44:892-898), and the parameters `n`, `alpha`, `Ksat`, `Smax` and `Smin`, obtained using [getsoilparams()]
#' control the relationship with soil moisture. The paramaters `n2` and `Ksat2` are equivelent to `n` and `Ksat`, but control
#' ground water seepage rates. Setting `Ksat2` to zero assumes ground water seepage. Paramaters for the
#' van Genuchten equation should be in units of centimeters and days. Correction for shorter time-scales
#' is automatically applied.
#'
#' For simplicity, the diffusion of soil water across the bottom of the lower layer is neglected since gradients
#' and resulting fluxes at this depth are generally unimportant, except over longer time-scales.
#' Exceptions include a high water table and the advance of a deep ‘wetting front’. In the
#' latter case, the thickness of the model should be increased.
#'
#' Surface runoff is dependent on the porosity of the soil and is controlled by a run-off
#' curve number, itself dependent on soil type, land cover and soil condition. Runoff curve numbers can be
#' obtained from here: https://www.wcc.nrcs.usda.gov/ftpref/wntsc/H&H/training/runoff-curve-numbers1.pdf.
#' If add.surplus is set to false (the default), runoff is assumed to leave the basin. If set to true,
#' it is added to rainfall at each time step.
#'
#' The paramaters `rain` and `cre`  can be provided as vectors or 3-D arrays. If provides as vectors,
#' spatially uniform rainfall and evapotranspiration is assumed, and number of time-steps over which
#' the model is run is determined the length of these vectors. If provided as 3D arrays, the the first
#' and second dimensions indicated spatial variability, and the third the number of timesteps over which
#' the model is run. The remaining paramaters, with the exception `timestep`, and `Trace` can all be
#' provided as single numeric values (no spatial vaiability) or as matrices or rasters (spatial variability).
#'
#' The model outputs three arrays, each covering the same extent as dem and the third the surface and soil water
#' at each time-step. Note however, that the model relies on the computation of slope angles using the raster
#' package, which sets edge cells to NA. In consequence, the edge cells of the model output are also NA. For the
#' model to work `dem` must have a planer projection such that x and y values are equidistant.
#'
#' @examples
#' # Rainfall and evapotranspiration time series
#' data(rainevap)
#' rain <- rainevap$rain
#' cre <- rainevap$evap
#' # Get van Genuchten soil paramaters
#' sp <- getsoilparams("Clay loam")
#' # ================================================================== #
#' # Run model for Lizard Peninsula over 365 days assuming clay-loam soils
#' # and with default settings. Trace is `TRUE` so progress is tracked
#' # MOdel takes ~ 5 mins to run
#' # ================================================================== #
#' hm <- hydromodel_spatial(dtm100m, rain, cre, sp)
#' # ============== Calculate and plot means  ==================== #
#' soil1 <- apply(hm$sm1, c(1, 2), mean) # soil moisture in layer 1
#' soil2 <- apply(hm$sm2, c(1, 2), mean) # soil moisture in layer 2
#' sfcew <- apply(hm$surface, c(1, 2), mean) # standing water (mm)
#' plot(if_raster(soil1, dtm100m), main = "Layer 1")
#' plot(if_raster(soil2, dtm100m), main = "Layer 2")
#'
hydromodel_spatial <- function(dem, rain, cre, sp, basins = NA, theta1 = 0.35, theta2 = 0.35, surface = 0,
                               p1 = 0.2, p2 = 0.1, p3 = 0.25, z1 = 5, z2 = 95, cn = 82, cover = 0.8,
                               topratio = 0.5, timestep = 86400, n2 = 1.1, Ksat2 = 0.5, Trace = TRUE) {
  n <- sp$n
  alpha <- sp$alpha
  Ksat <- sp$Ksat
  Smin <- sp$Smin
  Smax <- sp$Smax
  cat("Delineating basins and calculating basin attributes\n")
  if (class(basins) == "logical") {
    basins <- basindelin(dem)
  }
  mb <- is_raster(basins)
  md <- is_raster(dem)
  tpidx <- topidx(dem)
  tpidxm <- is_raster(tpidx)
  bc <- .basinchars(md, mb, sea = TRUE)
  sel <- which(bc$pourpointbasin == -999)
  bc$pourpointhgt[sel] <- 0
  bc <- bc[order(bc$pourpointhgt, decreasing = T), ]
  bc$basindepth <- bc$pourpointhgt - bc$basinmindepth
  bc$basindepth <- ifelse(bc$basindepth < 0, 0, bc$basindepth)
  # # Calculate basin volume
  bc$basinvolume <- 0
  for (i in 1:length(bc$basin)) {
    sel <- which(mb == bc$basin[i])
    bdem <- md[sel]
    bdem <- bdem[bdem < bc$pourpointhgt[i]]
    vols <- (bc$pourpointhgt[i] - bdem) / 1000
    vols <- vols * res(basins)[1] * res(basins)[2]
    bc$basinvolume[i] <- sum(vols)
  }
  u <- unique(bc$basin)
  u <- u[is.na(u) == F]
  # # # Convert rainfaill to array
  if (class(rain) == "array") {
    if (dim(rain)[1] != dim(cre)[1]) stop ("Rain and cre must have same dimensions")
    if (dim(rain)[2] != dim(cre)[2]) stop ("Rain and cre must have same dimensions")
  }
  if (class(rain) != "array") {
    if (length(rain) != length(cre)) stop ("Rain and cre vectors must have identical lengths\n")
  }
  if (class(rain) != "array") {
    cat("Converting rainfall and evaporation to arrays\n")
    raina <- rep(rain, each = dim(mb)[1] * dim(mb)[2])
    raina <- array(raina, dim = c(dim(mb)[1], dim(mb)[2], length(rain)))
    crea <- rep(cre, each = dim(mb)[1] * dim(mb)[2])
    crea <- array(crea, dim = c(dim(mb)[1], dim(mb)[2], length(rain)))
  } else {
    raina <- rain
    crea <- cre
  }
  # Convert initial values to matrix
  if (length(as.vector(theta1)) == 1) {
    theta1 <- array(theta1, dim = dim(md))
  }
  if (length(as.vector(theta2)) == 1) {
    theta2 <- array(theta2, dim = dim(md))
  }
  if (length(as.vector(surface)) == 1) {
    surface <- array(theta1, dim = dim(md))
  }
  if (length(as.vector(n)) == 1) {
    n <- array(n, dim = dim(md))
  }
  if (length(as.vector(alpha)) == 1) {
    alpha <- array(alpha, dim = dim(md))
  }
  if (length(as.vector(Ksat)) == 1) {
    Ksat <- array(Ksat, dim = dim(md))
  }
  if (length(as.vector(Smin)) == 1) {
    Smin <- array(Smin, dim = dim(md))
  }
  if (length(as.vector(Smax)) == 1) {
    Smax <- array(Smax, dim = dim(md))
  }
  if (length(as.vector(cn)) == 1) {
    cn <- array(cn, dim = dim(md))
  }
  if (length(as.vector(cover)) == 1) {
    cover <- array(cover, dim = dim(md))
  }
  if (length(as.vector(topratio)) == 1) {
    topratio <- array(topratio, dim = dim(md))
  }
  if (length(as.vector(n2)) == 1) {
    n2 <- array(n2, dim = dim(md))
  }
  if (length(as.vector(Ksat2)) == 1) {
    Ksat2 <- array(Ksat2, dim = dim(md))
  }
  if (length(as.vector(z1)) == 1) {
    z1 <- array(z1, dim = dim(md))
  }
  if (length(as.vector(z2)) == 1) {
    z2 <- array(z2, dim = dim(md))
  }
  if (length(as.vector(p1)) == 1) {
    p1 <- array(p1, dim = dim(md))
  }
  if (length(as.vector(p2)) == 1) {
    p2 <- array(p2, dim = dim(md))
  }
  if (length(as.vector(p3)) == 1) {
    p3 <- array(p3, dim = dim(md))
  }
  # Convert from raster
  surface <- is_raster(surface)
  theta1 <- is_raster(theta1)
  theta2 <- is_raster(theta2)
  z1 <- is_raster(z1)
  z2 <- is_raster(z2)
  n <- is_raster(n)
  alpha <- is_raster(alpha)
  Ksat <- is_raster(Ksat)
  Smin <- is_raster(Smin)
  Smax <- is_raster(Smax)
  cn <- is_raster(cn)
  cover <- is_raster(cover)
  topratio <- is_raster(topratio)
  n2 <- is_raster(n2)
  Ksat2 <- is_raster(Ksat2)
  # Probably need to create arrays for storing data here
  runff <- array(0, dim = c(length(u), dim(raina)[3]))
  sm1a <- raina * NA
  sm2a <- sm1a
  sm0a <- sm1a
  cat(paste("Running model over", length(u), "basins\n"))
  mrn <- mean(raina, na.rm = T)
  mcr <- mean(crea, na.rm = T)
  for (bsn in 1:length(u)) {
    sel <- which(mb == u[bsn])
    tx <- tpidxm[sel]
    if (is.na(mean(tx, na.rm = T)) == F) {
      bSmin <- Smin[sel]
      bSmax <- Smax[sel]
      bp1 <- p1[sel]
      bp2 <- p2[sel]
      bp3 <- p3[sel]
      th1 <- mean(theta1[sel], na.rm = T)
      th2 <- mean(theta2[sel], na.rm = T)
      th1 <- topdist(tx, th1, mean(bp1, na.rm = T), mean(bSmin, na.rm = T), mean(bSmax, na.rm = T))
      th2 <- topdist(tx, th2, mean(bp2, na.rm = T), mean(bSmin, na.rm = T), mean(bSmax, na.rm = T))
      swv <- sum((surface[sel] / 1000) * xres(basins) * yres(basins), na.rm = T)
      sfce <- topdistw(tx, swv, mean(bp3, na.rm = T), xres(basins), yres(basins))
      bz1 <- mean(z1[sel], na.rm = T)
      bz2 <- mean(z2[sel], na.rm = T)
      bn <- mean(n[sel], na.rm = T)
      balpha <- mean(alpha[sel], na.rm = T)
      bKsat <- mean(Ksat[sel], na.rm = T)
      bcn <- cn[sel]
      bcover <- mean(cover[sel], na.rm = T)
      btopratio <- mean(topratio[sel], na.rm = T)
      bn2 <- mean(n2[sel], na.rm = T)
      bKsat2 <- mean(Ksat2[sel], na.rm = T)
      for (day in 1:dim(raina)[3]) {
        rd <- raina[,,day][sel]
        cd <- crea[,,day][sel]
        rd[is.na(rd)] <- mrn
        cd[is.na(cd)] <- mcr
        rof <-  topdistw(tx, runff[bsn, day], 0.35, res(basins)[1], res(basins)[2])
        rd <- rd + rof
        oh <- .onehydro(rd, cd, sfce, th1, th2, bz1, bz2, bn, balpha, bKsat,
                        bSmin, bSmax, bcn, bcover, btopratio, timestep, bn2, bKsat2)
        # partiton between standing water and run-off  (values in m3)
        swv <- sum((oh$surplus / 1000) * xres(basins) * yres(basins), na.rm = T)
        msel <- which(bc$basin == bsn)
        ro <- swv - bc$basinvolume[msel]
        ro <- ifelse(ro < 0, 0, swv)
        sfce <- topdistw(tx, swv - ro, mean(bp3, na.rm = T), xres(basins), yres(basins))
        # Assign surpluss to adjoining basin in mm per cell
        tb <- bc$pourpointbasin[bsn]
        if (tb != -999) {
          runff[u[tb], day] <- ro
        }
        # Recycle & save
        th1 <- topdist(tx, mean(oh$theta1, na.rm = T), mean(bp1, na.rm = T),
                       mean(bSmin, na.rm = T), mean(bSmax, na.rm = T))
        th2 <- topdist(tx, mean(oh$theta2, na.rm = T), mean(bp2, na.rm = T),
                       mean(bSmin, na.rm = T), mean(bSmax, na.rm = T))
        sel2 <- sel + (day - 1) * dim(md)[1] * dim(md)[2]
        sm0a[sel2] <- sfce
        sm1a[sel2] <- th1
        sm2a[sel2] <- th2
      } # end day
    } # if check
    if (Trace) {
      iter <- ceiling(length(u) / 20)
      if (bsn%%iter == 0 | bsn == length(u)) {
        i <- ceiling(dim(raina)[3] / 2)
        r <- if_raster(sm1a[,,i], dem)
        plot(r, main = paste(bsn, "of", length(u), "done"))
      }
    }
  } # end bsn
  if (Trace) {
    s1 <- apply(sm1a, 3, mean, na.rm = T)
    s2 <- apply(sm2a, 3, mean, na.rm = T)
    plot(s1, type = "l", col = "red", xlab = "Time-steps", ylab = "Moisture",
         ylim = c(min(Smin, na.rm = T), max(Smax, na.rm = T)))
    par(new = T)
    plot(s2, type = "l", col = "blue", xlab = "", ylab = "",
         ylim = c(min(Smin, na.rm = T), max(Smax, na.rm = T)))
  }
  return(list(sm1 = sm1a, sm2 = sm2a, surface = sm0a))
}
