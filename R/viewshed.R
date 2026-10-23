#' What can be seen from a point
#'
#' The ground visible from an eye at height `eye` above the terrain at `xy`, computed
#' exactly on the terrain as triangles (vertices at cell centres, each square of four
#' split on the diagonal from top-left to bottom-right).
#'
#' The method is an object-space visibility map. Around the eye, four 90 degree wedges
#' of azimuth each use the central (gnomonic) projection, in which straight edges stay
#' straight. Triangles are visited front to back, and a horizon is kept per narrow
#' sector of azimuth as an upper envelope. A triangle is visible exactly where it rises
#' above the horizon, which it then raises. The visible pieces map back to the ground as
#' convex polygons, and the horizon left at the end is the skyline.
#'
#' @param z matrix of heights, row 1 at the north edge, square cells
#' @param extent `c(xmin, xmax, ymin, ymax)` of `z`
#' @param xy the eye, `c(x, y)`, inside `z` (inside the hull of its cell centres)
#' @param eye height of the eye above the ground (m)
#' @param far,far_extent optional coarser grid surrounding `z` (and its extent), for
#'   terrain further out; its cells that touch `z` are skipped
#' @param max_dist ignore terrain further than this (m)
#' @param curvature lower the terrain by `d^2 / 2R` (earth curvature with standard
#'   refraction, k = 0.13)
#' @param supersample visible fraction per cell is sampled at `supersample^2` points per
#'   cell (0 to skip)
#' @param pieces keep the visible polygons
#' @param sectors sectors per 90 degree wedge (more: less bookkeeping per sector)
#' @return a list with `visible` (matrix like `z`, fraction of each cell visible), `area`
#'   (visible area, m^2, including `far`), `skyline` (data frame of `azimuth`, degrees
#'   clockwise from grid north, and `elevation` of the horizon, degrees), `pieces` (data
#'   frame of polygon vertices `x`, `y` with `piece`, `triangle` and `ring` ids),
#'   `eye` (x, y, z of the eye) and counts. `NULL` if `xy` is outside `z`.
#' @export
#' @examples
#' z <- t(volcano)[ncol(volcano):1, ]
#' v <- viewshed(z, c(0, 870, 0, 610), c(300, 200))
#' v
#' image(v$visible[nrow(z):1, ] |> t())
viewshed <- function(z, extent, xy, eye = 1.5, far = NULL, far_extent = NULL, max_dist = Inf,
                     curvature = TRUE, supersample = 2L, pieces = TRUE, sectors = 64L) {
  g <- grid_args(z, extent)
  f <- far_args(far, far_extent)
  xy <- eye_xy(xy)
  r <- rs_viewshed(g$z, g$nrow, g$ncol, g$x0, g$ytop, g$cell,
                   f$z, f$nrow, f$ncol, f$x0, f$ytop, f$cell,
                   xy[1L], xy[2L], as.double(eye), as.double(max_dist), isTRUE(curvature),
                   as.integer(sectors), as.integer(supersample), isTRUE(pieces))
  if (is.null(r)) return(NULL)
  vis <- if (length(r$visible)) matrix(r$visible, g$nrow, g$ncol, byrow = TRUE) else NULL
  structure(list(visible = vis, area = r$area,
                 skyline = data.frame(azimuth = 0:359, elevation = r$skyline),
                 pieces = pieces_df(r$pieces),
                 eye = c(x = xy[1L], y = xy[2L], z = r$z_eye),
                 n_triangles = r$n_tri, n_visible_triangles = r$n_vis_tri),
            class = "aspersions_viewshed")
}

#' @export
print.aspersions_viewshed <- function(x, ...) {
  cat("<aspersions viewshed>\n")
  cat(sprintf("  eye at (%g, %g), %g m\n", x$eye[["x"]], x$eye[["y"]], x$eye[["z"]]))
  cat(sprintf("  visible: %.4g m^2 from %d of %d triangles\n", x$area, as.integer(x$n_visible_triangles), as.integer(x$n_triangles)))
  if (!is.null(x$visible)) cat(sprintf("  %.1f%% of the grid visible\n", 100 * mean(x$visible)))
  if (!is.null(x$pieces)) cat(sprintf("  %d visible pieces\n", length(unique(x$pieces$piece))))
  cat(sprintf("  skyline: %.1f to %.1f degrees\n", min(x$skyline$elevation), max(x$skyline$elevation)))
  invisible(x)
}

#' Viewshed over any triangulation
#'
#' As [viewshed()], for an explicit mesh: a constrained or adaptive triangulation, or a
#' grid with other diagonals. Triangles are fed in the front-to-back order from
#' [mesh_order()] unless `order` is given. The mesh must cover a convex region.
#'
#' @param x,y,z vertex coordinates and heights
#' @param triangles integer matrix of vertex indices (1-based), one row per triangle
#' @param order front-to-back order of triangles (1-based), or `NULL` to compute it
#' @inheritParams viewshed
#' @return a list with `area`, `skyline`, `pieces` (as [viewshed()]), `eye`, and
#'   `cycles` (triangles that had to be forced into order; 0 for Delaunay meshes)
#' @export
#' @examples
#' z <- t(volcano)[ncol(volcano):1, ]
#' xy <- expand.grid(x = seq(5, 865, by = 10), y = seq(605, 5, by = -10))
#' nc <- ncol(z); nr <- nrow(z)
#' q <- expand.grid(c = 1:(nc - 1), r = 1:(nr - 1))
#' tl <- (q$r - 1) * nc + q$c
#' tri <- rbind(cbind(tl, tl + 1, tl + nc), cbind(tl + 1, tl + nc + 1, tl + nc))
#' v <- viewshed_mesh(xy$x, xy$y, as.vector(t(z)), tri, c(300, 200))
#' v$area
viewshed_mesh <- function(x, y, z, triangles, xy, eye = 1.5, order = NULL, curvature = TRUE, sectors = 64L) {
  tri <- mesh_tri(triangles, length(x))
  xy <- eye_xy(xy)
  ord <- if (is.null(order)) integer(0) else as.integer(order) - 1L
  r <- rs_viewshed_mesh(as.double(x), as.double(y), as.double(z), tri, xy[1L], xy[2L],
                        as.double(eye), isTRUE(curvature), as.integer(sectors), ord)
  if (is.null(r)) return(NULL)
  structure(list(visible = NULL, area = r$area,
                 skyline = data.frame(azimuth = 0:359, elevation = r$skyline),
                 pieces = pieces_df(r$pieces), eye = c(x = xy[1L], y = xy[2L], z = r$z_eye),
                 n_triangles = r$n_tri, n_visible_triangles = r$n_vis_tri, cycles = r$cycles),
            class = "aspersions_viewshed")
}

#' Front-to-back order of a mesh's triangles
#'
#' Along any plan ray from the eye, the triangles crossed are neighbours across shared
#' edges, each crossed from the side facing the eye. So "the triangle on the eye's side
#' of a shared edge goes first", over all shared edges, is the whole in-front relation;
#' this is a topological sort of it. Delaunay meshes never have a cycle; other meshes
#' can, and triangles left in a cycle are appended nearest first and counted.
#'
#' @param x,y vertex coordinates
#' @inheritParams viewshed_mesh
#' @return a list with `order` (1-based triangle indices) and `cycles`
#' @export
mesh_order <- function(x, y, triangles, xy) {
  tri <- mesh_tri(triangles, length(x))
  xy <- eye_xy(xy)
  r <- rs_mesh_order(as.double(x), as.double(y), tri, xy[1L], xy[2L])
  list(order = as.integer(r$order) + 1L, cycles = as.integer(r$cycles))
}

mesh_tri <- function(triangles, nv) {
  tri <- as.integer(t(as.matrix(triangles)))
  if (length(tri) %% 3L != 0L) stop("'triangles' needs 3 vertex indices per triangle", call. = FALSE)
  if (anyNA(tri) || min(tri) < 1L || max(tri) > nv) stop("'triangles' indices must be in 1..length(x)", call. = FALSE)
  tri - 1L
}
