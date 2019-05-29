#' A raster of basins numbered as integers
#'
#' A raster object of basins numbered as integers for the area bounded by 160000, 181400,
#' 11300, 30000  (xmin, xmax, ymin, ymax) using the Ordance Survey GB Grid Reference system
#' (CRS: +init=epsg:27700) as returned by [basindelin()].
#'
#' @format A raster object with 187 rows and 214 columns.
"basins100m"

#' A 1 km resolution raster object of rainfall in west Cornwall, UK in January 2015.
#'
#' A raster object containing rainfall (mm) with sea coded as NA for the area bounded by
#' 82500, 198500, 4500, 81500  (xmin, xmax, ymin, ymax) using the Ordance Survey GB Grid
#' Reference system (CRS: +init=epsg:27700).
#'
#' @format A raster object with 77 rows and 116 columns.
#' @source \url{https://eip.ceh.ac.uk/rainfall/}
"cornwallrain"
#' Daily precipitation
#'
#' A dataset containing daily precipitation on the Lizard Peninsula for the period 1982-2015
#'
#' @format A data.frame with 12053 rows and two variables.
#' \describe{
#'   \item{obs_time}{Date of observation}
#'   \item{precipitation}{Daily precipitation total (mm)}
#' }
#' @source \url{https://eip.ceh.ac.uk/rainfall/}
#'
"dailyrain"

#' A 1 m resolution raster object of elevation for part of the Lizard Peninsula, Cornwall, UK.
#'
#' A raster object containing elevation in metres with sea coded as NA for the area bounded by
#' 169000, 170000, 12000, 13000  (xmin, xmax, ymin, ymax) using the Ordance
#' Survey GB Grid Reference system (CRS: +init=epsg:27700).
#'
#' @format A raster object with 1000 rows and 1000 columns.

#' @source \url{http://www.tellusgb.ac.uk/}
"dtm1m"

#' A 100 m resolution raster object of elevation for the Lizard Peninsula, Cornwall, UK.
#'
#' A raster object containing elevation in metres with sea coded as NA for the area bounded by
#' 160000, 181400, 11300, 30000  (xmin, xmax, ymin, ymax) using the Ordance Survey GB Grid
#' Reference system (CRS: +init=epsg:27700).
#'
#' @format A raster object with 187 rows and 214 columns.

#' @source \url{http://www.tellusgb.ac.uk/}
"dtm100m"

#' A 1 km resolution raster object of elevations in west Cornwall, UK.
#'
#' A raster object containing elevation in metres with sea coded as NA for the area bounded by
#' 169000, 170000, 12000, 13000   (xmin, xmax, ymin, ymax) using the Ordance Survey GB Grid
#' Reference system (CRS: +init=epsg:27700).
#'
#' @format A raster object with 77 rows and 116 columns.

#' @source \url{http://www.tellusgb.ac.uk/}
"dtm1km"

#' Demo time-series of rainfall and evapotranspiration
#'
#' A dataset of daily rainfall and evapotranspiration in mm for the Lizard Peninsula in
#' Cornwall in 2017.
#'
#'  @format A data frame with 365 rows and 2 variables
#'  \describe{
#'   \item {rain} {daily rainfall (mm)}
#'   \item {evap} {daily crop reference evapotranspiration (mm)}
#'  }
#' "rainevap"

#' Van Genuchten model parameters
#'
#' A dataset containing Van Genuchten model parameters and Hydrologic Soil Groupings for
#' major soil types derived by combining data from Carsel and Parrish (1988) Water Resources
#' Research 24:755–769 and Schaap et al. (1998) Soil Science Society of America Journal
#' 62:847–855.
#'
#' @format A data frame with 12 rows and 7 variables:
#' \describe{
#'  \item {soil.type} {description of soil type}
#'  \item {Smax} {Volumetric water content at saturation (cm^3 / cm^3)}
#'  \item {Smin} {Residual water content (cm^3 / cm^3)}
#'  \item {alpha} {Shape parameter of the van Genuchten mode (cm^-1)}
#'  \item {n} {Pore size distribution parameter (dimensionless, > 1)}
#'  \item {Ksat} {Saturated hydraulic conductivity (cm / day)}
#'  \item {HSG} {USDA Hydrological soil group}
#' }
#' @source: \url{https://onlinelibrary.wiley.com/doi/full/10.1002/ird.1751}
#' "soilparams"
