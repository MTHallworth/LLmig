#' Construct different fishnet types over a polygon surface.
#'
#' @param \code{areaOfInterest} a \code{SpatialPolygon} defining the boundaries of the fishnet area.
#' @param \code{type} the shape of the fishnet grid - must be "square" or "hexagon"
#' @param \code{cell_width} width of the cell in meters
#' @param \code{cell_area} area of the cell in meters
#' @param \code{clipToArea} logical indicating whether to clip fishnet to \code{areaOfInterest}
#'
#' @export

fishnet <- function(areaOfInterest, type, cell_width = NA, cell_area = NA, clipToArea = FALSE) {
if (!type %in% c("square", "hexagon")) {
  stop("type must be either 'square' or 'hexagon'")
}

if (is.na(cell_width)) {
  if (is.na(cell_area)) {
    stop("Provide either cell_width or cell_area")
}
}

if(!is.na(cell_area) & is.na(cell_width)){
if (type == "square") {
      cell_area <- cell_area*1e6
      cell_width <- sqrt(cell_area)
    }
}

if(!is.na(cell_area) & is.na(cell_width)){
   if (type == "hexagon") {
    cell_area <- cell_area*1e6
    cell_width <- sqrt(2 * cell_area / sqrt(3))
  }
}

# buffered extent of study area to define cells over
ext <- as(raster::extent(areaOfInterest) + cell_width, "SpatialPolygons")
raster::projection(ext) <- raster::projection(areaOfInterest)
# generate grid
if (type == "square") {
  g <- raster::raster(ext, resolution = cell_width)
  g <- as(g, "SpatialPolygons")

} else if (type == "hexagon") {
  # generate array of hexagon centers
  g <- sp::spsample(ext, type = "hexagonal", cellsize = cell_width, offset = c(0, 0))
  # convert center points to hexagons
  g <- sp::HexPoints2SpatialPolygons(g, dx = cell_width)
}
raster::crs(g) <- raster::crs(NH)

# clip to boundary of study area
if (clipToArea) {
  g <- rgeos::gIntersection(g, areaOfInterest, byid = TRUE)
} else {
  g <- g[areaOfInterest, ]
}
# clean up feature IDs
row.names(g) <- as.character(1:length(g))
return(g)
}
