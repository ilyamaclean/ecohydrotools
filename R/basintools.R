#' Orders drainage basins by elevation
#'
#' @description
#' `basinsort` is an internal function used by [basindelin()] and [basinmerge()] to number
#' drainage basins sequentially from lowest elevation (of lowest point) to highest.
#'
#' @param dem a raster object, two-dimensional array or matrix of elevations.
#' @param basins a raster object, two-dimensional array or matrix of basins numbered as integers.
#'
#' @return a raster object, two-dimensional array or matrix of sequentially numbered basins
#' @importFrom dplyr left_join
#'
#' @examples
#' basins <- matrix(c(1:4), nrow = 2, ncol = 2)
#' dem <- matrix(c(4:1), nrow = 2, ncol = 2)
#' basinsort(dem, basins)
basinsort <- function(dem, basins) {
  mdf <- function(u) {
    sel <- which(basins == u)
    min(dem[sel], na.rm = TRUE)
  }
  r <- dem
  dem <- is_raster(dem)
  basins <- is_raster(basins)
  u <- unique(as.vector(basins))
  u <- u[is.na(u) == F]
  u2 <- c(1:length(u))
  df1 <- data.frame(old = as.vector(basins))
  df2 <- data.frame(old = c(NA, u), new = c(NA, u2))
  df3 <- left_join(df1, df2, by = "old")
  basins <- array(df3$new, dim = dim(basins))
  u <- unique(as.vector(basins))
  u <- u[is.na(u) == F]
  mnd <- sapply(u, mdf)
  o <- order(mnd)
  u2 <- u[o]
  df1 <- data.frame(old = as.vector(basins))
  df2 <- data.frame(old = c(NA, u), new = c(NA, u2))
  df3 <- left_join(df1, df2, by = "old")
  bm2 <- array(df3$new, dim = dim(basins))
  if_raster(bm2, r)
}
#' Delineates hydrological basins
#'
#' @description `basindelin` uses digital elevation data to delineate hydrological basins.
#'
#' @param dem a raster object, two-dimensional array or matrix of elevations.
#'
#' @return a raster object, two-dimensional array or matrix  with individual basins numbered sequentially as integers.
#' @import raster
#' @export
#'
#' @seealso [basindelin_big()] for working with large datasets.
#'
#' @details
#' If `dem` is a raster object, a raster onbject is returned.
#' This function is used to delineate hydrological basins.
#' It iteratively identifies the lowest elevation pixel of `dem`, and
#' assigns any of the eight adjoining pixels to the same basin if higher.
#' The process is repeated on all assigned adjoining pixels until no further
#' higher pixels are found. The next lowest unassigned pixel is then
#' identified, a new basin identity assigned and the processes repeated until
#' all pixels are assigned to a basin. Relative to heuristic algorithms, it is
#' slow and run time increases exponentially with size of `dtm`. However, in
#' contrast to many such algorithms, all basins are correctly seperated by
#' boundaries >0. With large datasets, with e.g. > 160,000 pixels, the calculations
#' will be slow and [basindelin_big()] should be used instead.
#'
#' @examples
#' dem <- aggregate(dtm1m, 20)
#' basins <- basindelin (dem)
#' plot(basins, main = "Basins")
basindelin <- function(dem) {
  hgttongbr <- function(m) {
    hgt <- rep(m[2, 2], 8)
    neighbours <- c(m[1, 2], m[1, 3], m[2, 3], m[3, 3], m[3, 2], m[3, 1],
                    m[2, 1], m[1, 1])
    htn <- ifelse(hgt > neighbours, 1, 0)
    intcode <- sum(2 ^ (which(rev(unlist(strsplit(as.character(htn), "")) ==
                                    1)) - 1))
    intcode
  }
  integertobinary8 <- function(i) {
    a <- 2 ^ (0:9)
    b <- 2 * a
    binc <- format(sapply(i, function(x) sum(10 ^ (0:9)[(x %% b) >= a])),
                   scientific = FALSE)
    if (nchar(binc) > 8) warning("Integer > 8 bit binary")
    binc <- str_pad(binc, 8, pad = "0")
    binc
  }
  updown <- function(dem) {
    m <- dem
    m2 <- array(9999, dim = c(dim(dem)[1] + 2, dim(dem)[2] + 2))
    m2[2:(dim(dem)[1] + 1), 2:(dim(dem)[2] + 1)] <- m
    updownm <- matrix(rep(NA, length(dem)), nrow = dim(dem)[1])
    for (y in 2:(dim(m2)[2] - 1)) {
      for (x in 2:(dim(m2)[1] - 1)) {
        focalcell <- m2[x, y]
        if (is.na(focalcell) == FALSE) {
          m9 <- matrix(c(m2[x - 1, y - 1], m2[x, y - 1], m2[x + 1, y - 1],
                         m2[x - 1, y], focalcell, m2[x + 1, y],
                         m2[x - 1, y + 1], m2[x, y + 1], m2[x + 1, y + 1]),
                       nrow = 3)
          updownm[(x - 1), (y - 1)] <- hgttongbr(m9)
        }
      }
    }
    updownm
  }
  r <- dem
  dem <- is_raster(dem)
  ngbrow <- c(-1, -1, 0, 1, 1, 1, 0, -1)
  ngbcol <- c(0, 1, 1, 1, 0, -1, -1, -1)
  updownm <- updown(dem)
  basinsm <- ifelse(is.na(dem), NA, 0)
  donem <- dem
  basincell <- order(dem)[1]
  basin <- 1
  while (basincell != -999) {
    while (basincell != -999) {
      i <- basincell
      irow <- arrayInd(i, dim(dem))[1]
      icol <- arrayInd(i, dim(dem))[2]
      basinsm[irow, icol] <- basin
      neighbours <- integertobinary8(updownm[i])
      for (n in 1:8) {
        if ((irow + ngbrow[n]) > 0 & (irow + ngbrow[n]) <= dim(basinsm)[1] &
            (icol + ngbcol[n]) > 0 & (icol + ngbcol[n]) <= dim(basinsm)[2]) {
          if (!is.na(basinsm[(irow + ngbrow[n]), (icol + ngbcol[n])])) {
            if (substr(neighbours, n, n) == "0" &
                basinsm[(irow + ngbrow[n]), (icol + ngbcol[n])] == 0) {
              basinsm[(irow + ngbrow[n]), (icol + ngbcol[n])] <- basin
            }
          }
        }
      }
      donem[i] <- NA
      if (length(which(basinsm == basin & !is.na(donem))) > 0) {
        basincell <- min(which(basinsm == basin & !is.na(donem) ))
      }
      else basincell <- -999
    }
    if (length(which(!is.na(donem)) > 0)) {
      basincell <- which(donem == min(donem, na.rm = TRUE))[1]
    }
    else basincell <- -999
    basin <- basin + 1
  }
  basinsm <- if_raster(basinsm, r)
  basinsort(r, basinsm)
  }
#' Delineates hydrological basins for large datasets
#'
#' @description
#' `basindelin_big` is for use with large digital elevation datasets, to
#' delineate hydrological basins.
#'
#' @param dem a raster object of elevations.
#' @param dirout an optional character vector containing a single path directory for temporarily storing tiles. Deleted after use. Tilde expansion (see [path.expand()]) is done.
#' @param trace a logical value indicating whether to plot and report on progress.
#'
#' @return a raster object with individual basins numbered sequentially as integers.
#' @import raster
#' @export
#' @seealso [basindelin()] for working with smaller datasets.
#'
#' @details `basindelin_big` divides the large dataset into tiles and then
#' uses [basindelin()] to delineate basins for each tile before mosaicing back
#' together and merging basins along tile edges if not seperated by a boundary
#' > 0. If `dirout` is unspecified, then a directory `basinsout` is
#' temporarily created within the working directory. If `trace` is TRUE (the
#' default) then progress is tracked during three stages: (1) the basins
#' of each tile are plotted, (2) basins after mosaicing, but prior
#' to merging are plotted and (3) on each merge iteration, the number of basins
#' to merge is printed and processed basin is plotted.
#'
#' @examples
#' basins <- basindelin_big(dtm1m)
#' plot(basins, main = "Basins")
basindelin_big <- function(dem, dirout = NA, trace = TRUE) {
  onebasin_merge <- function(md, mb, b) {
    mb2 <- array(NA, dim = c(dim(mb)[1] + 2, dim(mb)[2] + 2))
    mb2[2:(dim(mb)[1] + 1), 2:(dim(mb)[2] + 1)] <- mb
    md2 <- array(NA, dim = c(dim(md)[1] + 2, dim(md)[2] + 2))
    md2[2:(dim(md)[1] + 1), 2:(dim(md)[2] + 1)] <- md
    bm <-  b
    sel <- which(mb == b)
    x <- arrayInd(sel, dim(md))[, 1]
    y <- arrayInd(sel, dim(md))[, 2]
    for (i in 1:length(x)) {
      mb9 <- mb2[x[i]:(x[i] + 2), y[i]:(y[i] + 2)]
      md9 <- md2[x[i]:(x[i] + 2), y[i]:(y[i] + 2)]
      sel2 <- which(mb9 > b)
      if (length(sel2) > 0) {
        for (j in 1:length(sel2)) {
          xi <- arrayInd(sel2[j], dim(mb9))[, 1]
          yi <- arrayInd(sel2[j], dim(mb9))[, 2]
          x2 <- c(xi - 1, xi, xi + 1)
          y2 <- c(yi - 1, yi, yi + 1)
          x2 <- x2[x2 > 0 & x2 < 4]
          y2 <- y2[y2 > 0 & y2 < 4]
          mb3 <- mb9[x2, y2]
          md3 <- md9[x2, y2]
          md3[mb3 == b] <- NA
          u <- unique(mb3[mb3 > b])
          u <- u[is.na(u) == F]
          for (k in 1:length(u)) {
            vrs <- md3[mb3  == u[k]]
            if (max(vrs, na.rm = T) > md9[2, 2]) bm <- c(bm, u[k])
          }
        }
      }
    }
    unique(bm)
  }
  dem <- trim(dem)
  dmsx <- ceiling(dim(dem)[2] / 200) - 1
  dmsy <- ceiling(dim(dem)[1] / 200) - 1
  if (dmsx < 1 & dmsy < 1) {
    basins <- basindelin(dem)
  }
  else {
    xres <- xres(dem)
    yres <- yres(dem)
    ed <- extent(dem)
    if (is.na(dirout)) dirout <- "basinsout/"
    dir.create(dirout)
    fol <- ""
    ii <- 1
    for (i in 0:dmsx) {
      for (j in 0:dmsy) {
        xmn <- ed@xmin + i * 200 * xres
        xmx <- min((xmn + 200 * xres), ed@xmax)
        ymn <- ed@ymin + j * 200 * yres
        ymx <- min((ymn + 200 * yres), ed@ymax)
        e <- extent(c(xmn, xmx, ymn, ymx))
        r <- crop(dem, e)
        v <- getValues(r)
        if (is.na(mean(v, na.rm = TRUE)) == FALSE) {
          b <- basindelin(r)
          if (trace) {
            progress <- paste0(round(ii / ((dmsx + 1) * (dmsy + 1)) * 100, 1),
                               "%")
            plot(b, main = paste0("basin slice progress: ", progress))
          }
          fo <- paste0(dirout, "b", i, "_", j, ".tif")
          fol <- c(fol, fo)
          writeRaster(b, file = fo, overwrite = TRUE)
        }
        ii <- ii + 1
      }
    }
    fol <- fol[fol != ""]
    basins <- raster(fol[1])
    ta <- max(getValues(basins), na.rm = TRUE)
    if (length(fol) > 1) {
      for (i in 2:length(fol)) {
        r <- raster(fol[i]) + ta
        basins <- mosaic(basins, r, fun = mean)
        ta <- max(getValues(basins), na.rm = TRUE)
      }
    }
    if (trace) {
      plot(basins, main = "Basin mosaic complete")
    }
    unlink(dirout)
    iter <- 1
    test <- 0
    basins[is.na(basins)] <- 9999
    dem[is.na(dem)] <- 9999
    while (test != 1) {
      basins <- basinsort(dem, basins)
      bmm <- getValues(basins, format = "matrix")
      lst <- as.list(c(1:max(getValues(basins), na.rm = TRUE)))
      ii <- 1
      if (dmsx > 0) {
        for (i in 1:dmsx) {
          xmn <- ed@xmin + i * 200 * xres - 1
          xmx <- min((xmn + 2 * xres), ed@xmax)
          e <- extent(c(xmn, xmx, ed@ymin, ed@ymax))
          r <- crop(basins, e)
          ds <- crop(dem, e)
          md <- getValues(ds, format = "matrix")
          mb <- getValues(r, format = "matrix")
          u <- unique(getValues(r))
          u <- u[is.na(u) == F]
          u <- u[order(u)]
          if (length(u) > 0) {
            for (k in 1:length(u)) {
              lst[[ii]] <- onebasin_merge(md, mb, u[k])
              ii <- ii + 1
            }
          }
        }
      }
      if (dmsy > 0) {
        for (j in 1:dmsy) {
          ymn <- ed@ymin + j * 200 * yres - 1
          ymx <- ymn + 2 * yres
          e <- extent(c(ed@xmin, ed@xmax, ymn, ymx))
          r <- crop(basins, e)
          ds <- crop(dem, e)
          md <- getValues(ds, format = "matrix")
          mb <- getValues(r, format = "matrix")
          u <- unique(getValues(r))
          u <- u[is.na(u) == F]
          u <- u[order(u)]
          if (length(u) > 0) {
            for (k in 1:length(u)) {
              lst[[ii]] <- onebasin_merge(md, mb, u[k])
              ii <- ii + 1
            }
          }
        }
      }
      lst <- lst[1:ii]
      ub <- unique(unlist(lst))
      ub <- ub[order(ub)]
      ub2 <- ub
      for (i in 1:length(ub)) {
        for (j in 1:ii) {
          v <- lst[[j]]
          tst <- which(v == ub[i])
          if (length(tst) > 0) ub2[i] <- min(ub2[i], min(v))
        }
      }
      sel <- which(ub != ub2)
      if (trace) {
        print(paste0("Basins to merge: ", length(sel)))
      }
      if (length(sel) == 0) {
        test <- 1
        if (trace) plot(basins, main = "Merge complete")
      }
      else {
        mbout <- bmm
        for (i in 1:length(sel)) {
          mbout[bmm == ub[sel[i]]] <- ub2[sel[i]]
        }
        basins <- raster(mbout, template = basins)
        if (trace) {
          plot(basins, main = paste0("Basin merge iteration: ", iter))
        }
        iter <- iter + 1
      }
    }
  }
  basins[dem == 9999] <- NA
  if (trace) plot (basins, main = "Merge complete")
  basins <- basinsort(dem, basins)
  basins
}
#' Internal function for calculating basin characteristics
#' @export
.basinchars <- function(md, mb, sea = F) {
  mb2 <- array(NA, dim = c(dim(mb)[1] + 2, dim(mb)[2] + 2))
  mb2[2:(dim(mb)[1] + 1), 2:(dim(mb)[2] + 1)] <- mb
  sel <- which(is.na(mb2) == T)
  mb2[sel] <- (-999)
  md2 <- array(NA, dim = c(dim(md)[1] + 2, dim(md)[2] + 2))
  md2[2:(dim(md)[1] + 1), 2:(dim(md)[2] + 1)] <- md
  if (sea) {
    md2[sel] <- -9999
  } else md2[sel] <- 9999
  dfm <- data.frame(basin = NA, basinmindepth = NA,
                    pourpointbasin = NA, pourpointhgt = NA)
  for (b in 1:max(mb, na.rm = T)) {
    sel <- which(mb == b)
    x <- arrayInd(sel, dim(md))[, 1]
    y <- arrayInd(sel, dim(md))[, 2]
    df1 <- data.frame(basin = b,
                      basinmindepth = min(md[sel], na.rm = T),
                      pourpointbasin = NA, pourpointhgt = 0)
    b2 <- 0
    pph <- 0
    for (i in 1:length(x)) {
      mb9 <- mb2[x[i]:(x[i] + 2), y[i]:(y[i] + 2)]
      md9 <- md2[x[i]:(x[i] + 2), y[i]:(y[i] + 2)]
      sel <- which(as.vector(mb9) != b)
      b2[i] <- NA
      pph[i] <- NA
      if (length(sel) > 0) {
        mb9 <- mb9[sel]
        md9 <- md9[sel]
        sel2 <- which(as.vector(md9) == min(as.vector(md9)))
        b2[i] <- mb9[sel2][1]
        pph[i] <- md9[sel2][1]
      }
    }
    sel <- which(pph == min(pph, na.rm = T))
    df1$pourpointbasin <- b2[sel[1]]
    df1$pourpointhgt <- min(pph, na.rm = T)
    dfm <- rbind(dfm, df1)
  }
  dfm <- dfm[which(is.na(dfm$basin) == F), ]
  dfm
}
#' Merges adjoining basins
#'
#' @description
#' `basinmerge` merges adjoining basins if the height differences between the
#' bottom of the basin and the pour point is less than than that
#' specified by `boundary`.
#'
#' @param dem a raster object, two-dimensional array or matrix of elevations.
#' @param basins a raster object, two-dimensional array or matrix with basins numbered as integers as returned by [basindelin()].
#' @param boundary a single numeric value. Basins seperated by boundaries below this height are merged (should have same units as `dtm`).
#'
#' @return a raster object, two-dimensional array or matrix with basins numbered as integers.
#' @import raster
#' @export
#'
#' @details
#' If `dem` is a raster object, then a raster object is returned.
#' If the differences in height between the pour-point and bottom of the basin is
#' less than that specified by `boundary` the basin is merged with basin to which
#' water or air would pour.
#'
#' @examples
#' basins2 <- basinmerge(dtm100m, basins100m, 1)
#' par(mfrow=c(1, 2))
#' plot(basins100m, main = "Basins")
#' plot(basins2, main = "Merged basins")
basinmerge <- function(dem, basins, boundary) {
  bvarsrem <- function(bvars, boundary) {
    sel <- which(bvars$pourpointbasin != -999)
    if (length(sel) == 0) warning("all basins flow into sea")
    bvars <- bvars[sel, ]
    bvars$basindepth <- bvars$pourpointhgt - bvars$basinmindepth
    bvars <- bvars[bvars$basindepth < boundary, ]
    o <- order(bvars$pourpointhgt, decreasing = TRUE)
    bvars <- bvars[o, ]
    bvars
  }
  if (all.equal(dim(basins)[1:2], dim(dem)[1:2]) == FALSE)
    stop ("basins and dem have different dimensions")
  r <- basins
  dem <- is_raster(dem)
  basins <- is_raster(basins)
  test <- F
  while (test == F) {
    mb2 <- basins
    bvars <- .basinchars(dem, basins)
    bkeep <- bvarsrem(bvars, boundary)
    if (dim(bkeep)[1] == 0) test <- T
    for (b in 1:(dim(bkeep)[1])) {
      sel <- which(bkeep$pourpointbasin == bkeep$basin[b])
      if (length(sel) > 0) {
        bkeep$pourpointbasin[sel] <- bkeep$pourpointbasin[b]
      }
      sel <- which(basins == bkeep$basin[b])
      mb2[sel] <- bkeep$pourpointbasin[b]
    }
    u <- unique(as.vector(mb2))
    sel <- which(is.na(u) == F)
    u <- u[sel]
    for (i in 1:length(u)) {
      sel <- which(mb2 == u[i])
      basins[sel] <- i
    }
  }
  basins <- basinsort(dem, basins)
  if_raster(basins, r)
}
#' Calculates flow direction
#'
#' @description
#' `flowdir` is used to calculate the direction of flow from every cell of a digital
#' elevation dataset
#'
#' @param dem a raster object, two-dimensional array or matrix of elevations.
#'
#' @return a raster object, two-dimensional array or matrix of flow directions (1-9).
#' @details Flow direction is expressed as given by `array(c(1:9), dim = c(3,3))`. I.e.
#' 1 = NW, 2 = W, 3 = SW, 4 = N, 5 = to self, 6 = S, 7 = NE, 8 = E and 9 = SE. Flow direction is determined
#' by which adjacent cell has the lowest elevation.
#' @export
#'
#' @examples
#' dem <- aggregate(dtm1m, 20)
#' plot(flowdir(dem), main = "Flow direction")
#'
flowdir <- function(dem) {
  md <- is_raster(dem)
  fd <- md * 0
  md2 <- array(NA, dim = c(dim(md)[1] + 2, dim(md)[2] + 2))
  md2[2:(dim(md)[1] + 1), 2:(dim(md)[2] + 1)] <- md
  v <- c(1:length(md))
  v <- v[is.na(md) == F]
  x <- arrayInd(v, dim(md))[, 1]
  y <- arrayInd(v, dim(md))[, 2]
  for (i in 1:length(x)) {
    md9 <- md2[x[i]:(x[i] + 2), y[i]:(y[i] + 2)]
    fd[x[i], y[i]] <- round(mean(which(md9 == min(md9, na.rm = T))), 0)
  }
  if_raster(fd, dem)
}
#' Calculates accumulated flow
#'
#' @description
#' `flowacc` is used by [pcad()] to calculate accumulated flow to each cold air drainage
#' basin
#'
#' @param dem a raster object, two-dimensional array or matrix of elevations.
#' @param basins a raster object, two-dimensional array or matrix with basins numbered as integers as returned by [basindelin()].
#'
#' @return a raster object, two-dimensional array or matrix of accumulated flow.
#' @details Accumulated flow is expressed in terms of number of cells.
#' @export
#'
#' @examples
#' dem <- aggregate(dtm1m, 20)
#' plot(flowacc(dem), main = "Accumulated flow")
#'
flowacc <- function(dem) {
  dm <- is_raster(dem)
  fd <- flowdir(dm)
  fa <- fd * 0 + 1
  o <- order(dm, decreasing = T, na.last = NA)
  for (i in 1:length(o)) {
    x <- arrayInd(o[i], dim(dm))[1]
    y <- arrayInd(o[i], dim(dm))[2]
    f <- fd[x, y]
    x2 <- x + (f - 1) %% 3 - 1
    y2 <- y + (f - 1) %/% 3 - 1
    if (x2 > 0 & x2 < dim(dm)[1] & y2 > 0 & y2 < dim(dm)[2])
      fa[x2, y2] <- fa[x, y] + 1
  }
  if_raster(fa, dem)
}
#' Calculates Bevan and Kirkby's topographic wetness index
#'
#' @description
#' `topidx` is used to calculate the topographic wetness index described in Bevan and Kirkby (1979)
#' basin
#'
#' @param dem a raster object of elevations.
#' @param minslope an optional positive value specifying the minimum slope allowed (see details).
#'
#' @return a raster object of topographic wetness.
#'
#' @details If slope is zero, infinite topographic wetness index values are returned. Thus,
#' flat areas are assigned a value specifed by `minslope`. The default value, given by
#' atan(0.005 / mean(res(dem))), assumes that height differences between adjacent cells exceed
#' 0.005. Accumulated flow is multiplied by the cell dimensions, thus`dem` should have an
#' equal area projection.
#'
#' @export
#' @references Bevan KJ & Kirkby MJ (1979) A physically based, variable contributing area
#' model of basin hydrology. Hydrological Sciences Journal 24: 43-69.
#'
#' @examples
#' dem <- aggregate(dtm1m, 20)
#' plot(log(topidx(dem)), main = "Topographic wetness")
#'
topidx <- function(dem, minslope = atan(0.005 / mean(res(dem)))) {
  slope <- terrain(dem)
  B <- is_raster(slope)
  B[B < minslope] <- minslope
  a <- is_raster(flowacc(dem))
  a <- a * res(dem)[1] * res(dem)[2]
  tpx <-  a / tan(B)
  if_raster(tpx, dem)
}
#' Calculates slope perpendicular to direction of flow
#'
#' @description
#' `perpslope` is an internal function used by [flowvelocity()] to calculate the slope
#' perpendicular to the direction of flow.
#'
#' @param dem a raster object of elevations.
#'
#' @return a raster object of slope angle in radians perpendicular to the direction of flow
#' @importFrom dplyr left_join
#' @export
#'
#' @examples
#' dem <- aggregate(dtm1m, 20)
#' plot(perpslope(dem) * (180 / pi), main = "slope")
#'
perpslope <- function(dem) {
  md <- is_raster(dem)
  md2 <- array(NA, dim = c(dim(md)[1] + 2, dim(md)[2] + 2))
  md2[2:(dim(md)[1] + 1), 2:(dim(md)[2] + 1)] <- md
  m <- array(c(1:9), dim = c(3, 3))
  df1 <- data.frame(old = as.vector(flowdir(md)))
  df2 <- data.frame(old = as.vector(m), new = as.vector(t(m[3:1, ])))
  df3 <- left_join(df1, df2, by = "old")
  fdl <- array(df3$new, dim = dim(md))
  df2 <- data.frame(old = as.vector(m), new = as.vector(t(m[, 3:1])))
  df3 <- left_join(df1, df2, by = "old")
  fdr <- array(df3$new, dim = dim(md))
  ad <- array(NA, dim = c(dim(md), 9))
  xs <- c(-1, 0, 1, -1, 0, 1, -1, 0, 1)
  ys <- c(-1, -1, -1, 0, 0, 0, 1, 1, 1)
  # check this bit
  for (i in 1:9)
    ad[,,i] <- md2[(2 - ys[i]):(dim(md)[1] - ys[i] + 1), (2 - xs[i]):(dim(md)[2] - xs[i] + 1)]
  hl <- ad
  hr <- ad
  for (i in 1:9) hl[,,i] <- ifelse(fdl == i, ad[,,i], 0)
  for (i in 1:9) hr[,,i] <- ifelse(fdr == i, ad[,,i], 0)
  hl <- apply(hl, c(1,2), sum)
  hr <- apply(hr, c(1,2), sum)
  hd <- pmax(sqrt((hl - md)^2), sqrt((hr - md)^2))
  theta <- atan(hd / mean(res(dem)))
  if_raster(theta, dem)
}
#' Finds coefficient for scaling topographic wetness to flow velocity
#'
#' @description
#' `findm` is an internal function used by [flowvelocity()] to scale topographic
#' wetness to flow velocity.
#'
#' @param tpx a raster object, two-dimensional array or matrix of topographic wetness
#' as returned by [topidx()]
#'
#' @return a numeric scaling ceofficient
#' @export
#'
#' @examples
#' findm(topidx(dtm100m))
#'
findm <- function(tpx) {
  tpx <- log(is_raster(tpx))
  mn <- 0
  for (i in 0:100) {
    a <- i/5 - 10
    ltpx <- 1 / (1 + exp(-log(tpx) - a))
    mn[i+1] <- mean(ltpx, na.rm = T)
  }
  i <- c(0:100)
  a <- i/5 - 10
  mn2 <- log(mn / (1 - mn))
  m1 <- lm(a ~ mn2)
  as.numeric(m1$coef[1])
}
#' Calculates flow velocity
#'
#' @description
#' `flowvelocity` applies Manning's equation to calculate flow velocity across all cells
#' of a digital elevation model. It is used to calculate which cells contribute to peak
#' flow.
#'
#' @param dem a raster object of elevations
#' @param nr an optional single numeric value, or a raster object, two-dimensional array or matrix
#' of Manning's roughness coefficients (see details)
#' @param wf n optional single numeric value, or a raster object, two-dimensional array or matrix of
#' mean water volumes (see details)
#' @param minslope an optional positive value specifying the minimum slope allowed (see details).
#'
#' @return a raster object of flow velocities (m / s)
#'
#' @details If `wf` is a single value it should approximate the mean proportion of land wetted.
#' If `wf` is an array, matrix or raster object, then values should approximate the mean proportion of
#' land wetted seperately for each basin, but should not be variable within basins.
#' If slope is zero, infinite topographic wetness index values are returned. Thus,
#' flat areas are assigned a value specifed by `minslope`. The default value, given by
#' atan(0.005 / mean(res(dem))), assumes that height differences between adjacent cells exceed
#' 0.005. Accumulated flow is multiplied by the cell dimensions, thus`dem` should have an
#' equal area projection. Values of `nr` fora variaty of substrate types can be found
#' here: \code{\link{http://www.fsl.orst.edu/geowater/FX3/help/8_Hydraulic_Reference/Mannings_n_Tables.htm}}.
#'
#'
#' @export
#'
#' @examples
#' plot(flowvelocity(dtm100m), main = "Flow velocity")
#'
flowvelocity <- function(dem, nr = 0.025, wf = 0.0001, minslope = atan(0.005 / mean(res(dem)))) {
  res <- mean(res(dem))
  nr <- is_raster(nr)
  wf <- is_raster(wf)
  pslope <- is_raster(perpslope(dem))
  slope <- is_raster(terrain(dem))
  pslope[pslope < minslope] <- minslope
  slope[slope < minslope] <- minslope
  tpx <- is_raster(topidx(dem))
  a <- log(wf / (1 - wf)) + findm(tpx)
  W <- 1 / (1 + exp(-log(tpx) - a))
  hd <- res * tan(pslope)
  A <- 0.5 * hd * res
  P <- sqrt(hd ^ 2 + res^2)
  H <- (A / P) * W
  V <- (1 / nr) * H^(2/3) * sqrt(atan(slope))
  if_raster(V, dem)
}
#' Calculates flow time to basin bottom
#'
#' @description
#' `flowtimes` calculates flow paths and applies the flow velocity function to calculate the
#' the time taken for water from each cell to reach the lowest point in the basin. It is used
#' to calculate which cells contribute to peak flow.
#'
#' @param dem a raster object of elevations
#' @param nr an optional single numeric value, or a raster object, two-dimensional array or matrix
#' of Manning's roughness coefficients (see details)
#' @param wf n optional single numeric value, or a raster object, two-dimensional array or matrix of
#' mean water volumes (see details)
#' @param maxiter a positive integer specifying the maximum number of iterations allowed when flow
#' path algorithm gets stuck in a loop.
#' @param minslope an optional positive value specifying the minimum slope allowed (see details).
#'
#' @return a raster object of flow times (s)
#'
#' @details If `wf` is a single value it should approximate the mean proportion of land wetted.
#' If `wf` is an array, matrix or raster object, then values should approximate the mean proportion of
#' land wetted seperately for each basin, but should not be variable within basins.
#' If slope is zero, infinite topographic wetness index values are returned. Thus,
#' flat areas are assigned a value specifed by `minslope`. The default value, given by
#' atan(0.005 / mean(res(dem))), assumes that height differences between adjacent cells exceed
#' 0.005. Accumulated flow is multiplied by the cell dimensions, thus`dem` should have an
#' equal area projection. Values of `nr` fora variaty of substrate types can be found
#' here: \code{\link{http://www.fsl.orst.edu/geowater/FX3/help/8_Hydraulic_Reference/Mannings_n_Tables.htm}}.
#'
#'
#' @export
#'
#' @examples
#' plot(log(flowtimes(dtm100m)), main = "Log flow time (s)")
#'
flowtimes <- function(dem, nr = 0.025, wf = 0.0005, maxiter = 100,
                      minslope = atan(0.005 / mean(res(dem)))) {
  res <- mean(res(dem))
  fd <- is_raster(flowdir(dem))
  fv <- is_raster(flowvelocity(dem, nr, wf, minslope))
  fd2 <- array(5, dim = dim(fd) + 2)
  fd2[2:(dim(fd)[1] + 1), 2:(dim(fd)[2] + 1)] <- fd
  fv2 <- array(100, dim = dim(fd) + 2)
  fv2[2:(dim(fv)[1] + 1), 2:(dim(fv)[2] + 1)] <- fv
  fv2[is.na(fv2)] <- 100
  fd2[is.na(fd2)] <- 5
  md <- is_raster(mask(dem, if_raster(fv, dem)))
  mdo <- order(md, decreasing = T, na.last = NA)
  # For each pixel in turn
  tms <- array(NA, dim = dim(fd))
  for (cell in 1:length(mdo)) {
    dne <- T
    x <- arrayInd(mdo[cell], dim(md))[, 1]
    y <- arrayInd(mdo[cell], dim(md))[, 2]
    tme <- 1 / fv2[x + 1, y +1]
    iter <- 0
    while (dne & iter < maxiter) {
      fdc <- fd2[x + 1, y + 1]
      if (is.na(fdc) | fdc == 5) dne <- F
      x <- x + (fdc - 1) %% 3 - 1
      y <- y + (fdc - 1) %/% 3 - 1
      tme <- tme + 1 / fv2[x + 1, y + 1]
      iter <- iter + 1
    }
    tms[mdo[cell]] <- tme
  }
  tms <- tms * res
  tms <- if_raster(tms, dem)
  mask(tms, dem)
}
#' Calculates the contribution of each cell to peak flow
#'
#' @description
#' `contributiontopeak` calculates the contribution of each cell of a digital elevation dataset to
#' peak flow
#'
#' @param ft a raster object, two-dimensional array or matrix of flow times as returned by [flowtimes()]
#' @param basins a raster object, two-dimensional array or matrix with basins numbered as integers as
#' returned by [basindelin()].
#' @param bysize an optional logical indicating whether or not to weight values by basin size
#'
#' @return a raster object of the contribution of each cell to peak flow
#'
#' @details contribution to peak flow is calculated by fitting a log-normal distribution
#' to flow times, estimating the probability density function (pdf) of this distribution, weighted
#' by the maximum value of of the pdf. If `bysize` is TRUE, output values are multiplied
#' by the size of the basin thereby making the assumption that larger basins have higher
#' peak flow. Alternatively real peak flows can be estimated from meterological
#' data.
#'
#' @importFrom MASS fitdistr
#' @export
#'
#' @examples
#' basins <- basindelin(dtm100m) # Takes ~20 seconds to run
#' ft <- flowtimes(dtm100m)
#' plot(contributiontopeak(ft, basins))
#' plot(contributiontopeak(ft, basins, bysize = F))
#'
contributiontopeak <- function(ft, basins, bysize = TRUE) {
  fm <- is_raster(ft)
  bm <-  is_raster(basins)
  zm <- fm * 0
  u <- unique(as.vector(bm))
  u <- u[is.na(u) == F]
  for (i in 1:length(u)) {
    sel <- which(bm == u[i])
    ftb <- fm[sel]
    sel2 <- which(is.na(ftb) == F)
    z <- ftb
    ftb <- ftb[sel2]
    if (length(ftb) > 2) {
      f1 <- fitdistr(log(ftb), densfun = "normal")
      z[sel2] <- dnorm(log(ftb), mean=f1$estimate[1], sd=f1$estimate[2])
      z <- z  / max(z, na.rm=T)
      if (bysize) z <- z * length(sel)
      zm[sel] <- z
    }
  }

  if_raster(zm, ft)
}
#' Distributes soil moisture by topographic wetness
#'
#' @description
#' `topdist` distributes soil moisture, typically within a basin, by the topograhic wetness index of each basin
#' grid cell
#'
#' @param tx a vector or marix of topographic wetness index values
#' @param sm a single numeric value of the mean soil water fraction
#' @param p a single numeric value of to power adjustment to apply to topographic wetness
#' index values (see details)
#' @param Smin a single numeric value of the residual soil water fraction (see [soilparams()])
#' @param Smax a single numeric value of the saturated  soil water fraction (see [soilparams()])
#'
#' @return a vector or matrix of soil moisture values of each grid cell of tx
#' @seealso [topidx()], [soilparams()]
#'
#' @details the coefficient p scales the relationship between topographic wetness and the
#' distribution of soil moisture values. If p = 1, soil moisture is strongly concentrated in
#' cells with a high topographic wetness. If p = 0, soil moisture is evenly distributed across
#' the basin. The correct choice of p depends on sub-grid scale variation in topographic wetness, and on
#' complex lateral flows, and thus best determined empirically, by measuring soil moisture across
#' a catchment. A choice of 0.25 for surface soil layers and 0.15 for sub-surface soil layers
#' is likely to provide a reasonable approximation.
#'
#' @export
#'
#' @examples
#' library(raster)
#' tpx <- topidx(dtm100m)
#' tp1 <- topdist(is_raster(tpx), 0.3, p = 0.25)
#' tp2 <- topdist(is_raster(tpx), 0.3, p = 0.15)
#' # NB - not distributed by basin
#' plot(if_raster(tp1, dtm100m))
#' plot(if_raster(tp2, dtm100m))
topdist <- function(tx, sm, p = 0.25, Smin = 0.091, Smax = 0.419) {
  sm <- ifelse(sm > Smax, Smax, sm)
  sm <- ifelse(sm < Smin, Smin, sm)
  rscl <- (sm - Smin) / (Smax - Smin)
  tx2 <- tx ^ p
  tx2 <- tx2[is.na(tx2) == F]
  av <- log(rscl / (1 - rscl))
  sout <- tx2 - mean(tx2) + av
  sout <- 1 / (1 + exp(-1 * sout))
  sout <- sout * (Smax - Smin) + Smin
  sm2 <- mean(sout)
  rat <- sm / sm2
  sout2 <- rat * sout
  if (rat > 1) {
    x1 <- sum(sout2) - length(sout2[sout2 > Smax]) * Smax
    x2 <- sum(sout2[sout2 <= Smax])
    sout2 <- sout2 * (x1 / x2)
    sout2[sout2 > Smax] <- Smax
  } else {
    x1 <- sum(sout2) - length(sout2[sout2 < Smin]) * Smin
    x2 <- sum(sout2[sout2 >= Smin])
    sout2 <- sout2 * (x1 / x2)
    sout2[sout2 < Smin] <- Smin
  }
  if (round(rat, 3) == 1) sout2 <- sout
  so <- tx
  sel <- which(is.na(tx) == F)
  so[sel] <- sout2
  so
}

#' Distributes surface water by topographic wetness
#'
#' @description
#' `topdistw` surface water, typically within a basin, by the topograhic wetness index of each basin
#' grid cell
#'
#' @param tx a vector or matrix of topographic wetness index values
#' @param wvol a single numeric value of basin  water volume (m^3)
#' @param p a single numeric value of to power adjustment to apply to topographic wetness
#' index values (see details)
#' @param xres x grid cell resolution (m)
#' @param yres y grid cell resolution (m)
#'
#' @return a vector or matrix of surface water depth (mm) for each grid cell of tx
#' @seealso [topidx()]
#'
#' @details the coefficient p scales the relationship between topographic wetness and the
#' distribution of surface water. If p = 1, soil moisture is strongly concentrated in
#' cells with a high topographic wetness. If p = 0, soil moisture is evenly distributed across
#' the basin. The correct choice of p depends on sub-grid scale variation in basin depressions, and on
#' complex lateral flows, and thus best determined empirically. A choice of 0.3 for is likely to provide a reasonable approximation.
#'
#' @export
#'
#' @examples
#' library(raster)
#' tpx <- topidx(dtm100m)
#' # average depth on 10 mm
#' sel <- which(is.na(is_raster(tpx) == F))
#' wvol <- 10 * 100 * 100 * length(sel) / 1000
#' # NB - not distributed by basin
#' swd1 <- topdistw(is_raster(tpx), wvol, p = 0.3, xres = 100, yres = 100)
#' swd2 <- topdistw(is_raster(tpx), wvol, p = 0.1, xres = 100, yres = 100)
#' plot(if_raster(swd1, dtm100m))
#' plot(if_raster(swd2, dtm100m))
topdistw <- function(tx, wvol, p = 0.3, xres, yres) {
  if (wvol > 0) {
    tx2 <- tx ^ p
    swtr <- tx2 * wvol
    swtr <- swtr * wvol / sum(swtr, na.rm = T)
    # Convert to mm per grid cell
    swtr <- (swtr * 1000) / (xres * yres)
  } else swtr <- rep(0, length(tx))
  swtr <-ifelse(swtr < 0, 0, swtr)
}

