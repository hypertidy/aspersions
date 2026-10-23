#' aspersions: cast shade and sight lines on terrain
#'
#' Two kinds of casting over a terrain grid.
#'
#' * [viewshed()]: what can be seen from a point. Exact, on the terrain as triangles:
#'   the visible part of every triangle as a polygon, a visible fraction per cell, and the
#'   skyline from the point. [viewshed_mesh()] does the same for any triangulation, in
#'   the front-to-back order from [mesh_order()].
#' * [cast_shade()]: where the terrain shades itself for one sun position, and
#'   [sun_hours()] for the hours of direct sun over a day, with [sun_position()] for the
#'   sun.
#'
#' Grids are a matrix of heights with row 1 at the north edge, square cells, and an
#' extent `c(xmin, xmax, ymin, ymax)` in a projected coordinate system in metres; or
#' directly the result of `gdalraster::read_ds()`, which carries its own extent. Results
#' are matrices in the same layout, ready for `ximage::ximage(m, extent = ...)`.
#'
#' @keywords internal
"_PACKAGE"
