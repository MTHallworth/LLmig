#' Construct different fishnet types over a polygon surface.
#'
#' @param \code{areaOfInterest} a \code{SpatialPolygon} defining the boundaries of the fishnet area.
#' @param \code{type} the shape of the fishnet grid - must be "square" or "hexagon"
#' @param \code{cell_width} width of the cell in meters
#' @param \code{cell_area} area of the cell in meters
#' @param \code{clipToArea} logical indicating whether to clip fishnet to \code{areaOfInterest}
#'
#' @export

fishnet <- function(areaOfInterest, type, cell_width, cell_area, clipToArea = FALSE) {
if (!type %in% c("square", "hexagon")) {
  stop("Type must be either 'square' or 'hexagon'")
}

if (missing(cell_width)) {
  if (missing(cell_area)) {
    stop("Provide either cell_width or cell_area")
  } else {
    if (type == "square") {
      cell_width <- sqrt(cell_area)
    } else if (type == "hexagonal") {
      cell_width <- sqrt(2 * cell_area / sqrt(3))
    }
  }
}
# buffered extent of study area to define cells over
ext <- as(raster::extent(x) + cell_width, "SpatialPolygons")
raster::projection(ext) <- raster::projection(x)
# generate grid
if (type == "square") {
  g <- raster::raster(ext, resolution = cell_width)
  g <- as(g, "SpatialPolygons")
} else if (type == "hexagonal") {
  # generate array of hexagon centers
  g <- sp::spsample(ext, type = "hexagonal", cellsize = cell_width, offset = c(0, 0))
  # convert center points to hexagons
  g <- sp::HexPoints2SpatialPolygons(g, dx = cell_width)
}

# clip to boundary of study area
if (clip) {
  g <- rgeos::gIntersection(g, x, byid = TRUE)
} else {
  g <- g[x, ]
}
# clean up feature IDs
row.names(g) <- as.character(1:length(g))
return(g)
}
