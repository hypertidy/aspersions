# a matrix and extent to the arguments the Rust side wants (row-major, row 0 north)
# z is a matrix (row 1 north) with extent c(xmin, xmax, ymin, ymax), or a vector from
# gdalraster::read_ds() whose "gis" attribute carries dim and bbox
as_grid <- function(z, extent = NULL, label = "z") {
  gis <- attr(z, "gis")
  if (!is.matrix(z) && !is.null(gis)) {
    z <- matrix(as.vector(z)[seq_len(gis$dim[1] * gis$dim[2])], gis$dim[2], gis$dim[1], byrow = TRUE)
    if (is.null(extent)) extent <- gis$bbox[c(1, 3, 2, 4)]
  }
  if (!is.matrix(z)) stop(sprintf("'%s' must be a matrix of heights (row 1 at the north edge) or a gdalraster::read_ds() result", label), call. = FALSE)
  if (is.null(extent)) stop(sprintf("'%s' needs an extent, c(xmin, xmax, ymin, ymax)", label), call. = FALSE)
  list(z = z, extent = as.double(extent))
}
grid_args <- function(z, extent, label = "z") {
  zz <- as_grid(z, extent, label); z <- zz$z; extent <- zz$extent
  if (length(extent) != 4L || !all(is.finite(extent))) stop("an extent is c(xmin, xmax, ymin, ymax)", call. = FALSE)
  nr <- nrow(z); nc <- ncol(z)
  if (nr < 2L || nc < 2L) stop(sprintf("'%s' needs at least 2 rows and 2 columns", label), call. = FALSE)
  dx <- (extent[2L] - extent[1L]) / nc; dy <- (extent[4L] - extent[3L]) / nr
  if (!(dx > 0 && dy > 0)) stop("extent must have xmax > xmin and ymax > ymin", call. = FALSE)
  if (abs(dx - dy) > 1e-6 * max(dx, dy)) stop(sprintf("cells must be square (%s: %g x %g)", label, dx, dy), call. = FALSE)
  v <- as.vector(t(z))
  v[!is.finite(v)] <- 0
  list(z = as.double(v), nrow = as.integer(nr), ncol = as.integer(nc), x0 = extent[1L], ytop = extent[4L], cell = dx,
       extent = extent)
}
no_far <- list(z = double(0), nrow = 0L, ncol = 0L, x0 = 0, ytop = 0, cell = 1)
far_args <- function(far, far_extent) {
  if (is.null(far)) return(no_far)
  if (is.null(far_extent)) stop("'far' needs 'far_extent'", call. = FALSE)
  grid_args(far, far_extent, "far")
}
eye_xy <- function(xy) {
  xy <- as.double(unlist(xy))
  if (length(xy) != 2L || !all(is.finite(xy))) stop("'xy' is one point, c(x, y)", call. = FALSE)
  xy
}
pieces_df <- function(p) {
  if (is.null(p)) return(NULL)
  np <- length(p$tri)
  nv <- diff(p$off)
  m <- matrix(p$xy, ncol = 2L, byrow = TRUE)
  data.frame(x = m[, 1L], y = m[, 2L], piece = rep(seq_len(np), nv),
             triangle = rep(p$tri + 1, nv), ring = rep(p$tag + 1L, nv))
}
