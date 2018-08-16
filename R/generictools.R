#' Flexible conversion to raster object
#'
#' @description
#' `if.raster` is used to permit flexibility in the use of rasters, matrices or arrays in many functions.
#'
#' @param x an R object
#' @param r an R object
#'
#' @return if `r` is a raster, `x` is converted to a raster with the same attributes as `r`, otherwise returns `x`
#' @import raster
#'
#' @examples
#' r <- is.raster(dtm100m)
#' r1 <- if.raster(r, dtm100m)
#' r2 <- if.raster(r, r)
#' class(r1) # is a RasterLayer
#' class(r2) # is a matrix
if.raster <- function(x, r) {
  if (class(r) == "RasterLayer")
    x <- raster(x, template = r)
  x
}
#' Checks whether object is a raster and returns a matrix if yes

#' @description
#' `is.raster` is used to permit flexibility in the use of rasters, matrices or arrays in many functions.
#'
#' @param r an R object
#'
#' @return if `r` is a raster, returns a matrix containing all values of `r`, otherwise returns `r`
#' @import raster
#'
#' @examples
#' r <- is.raster(dtm100m)
#' class(dtm100m) # is a RasterLayer
#' class(r) # is a matrix
#' plot(r) # not a raster
#' plot(raster(r)) # converts to raster
is.raster <- function(r) {
  if (class(r) == "RasterLayer")
    r <- getValues(r, format = "matrix")
  r
}
