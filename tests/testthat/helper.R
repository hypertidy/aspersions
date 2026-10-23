volcano_grid <- function() {
  list(z = t(datasets::volcano)[ncol(datasets::volcano):1, ], extent = c(0, 870, 0, 610))
}
# height on the grid's TIN (vertices at cell centres, quads split top-left to bottom-right)
tin_height <- function(z, extent, x, y) {
  cell <- (extent[2] - extent[1]) / ncol(z)
  u <- (x - extent[1]) / cell - 0.5; v <- (extent[4] - y) / cell - 0.5
  c0 <- pmin(ncol(z) - 1, floor(u)); r0 <- pmin(nrow(z) - 1, floor(v))
  fx <- u - c0; fy <- v - r0
  tl <- z[cbind(r0 + 1, c0 + 1)]; tr <- z[cbind(r0 + 1, c0 + 2)]
  bl <- z[cbind(r0 + 2, c0 + 1)]; br <- z[cbind(r0 + 2, c0 + 2)]
  ifelse(fx >= fy, tl + fx * (tr - tl) + fy * (br - tr), tl + fy * (bl - tl) + fx * (br - bl))
}
# brute-force line of sight on the TIN: the most the terrain rises above the sight line
los_margin <- function(z, extent, eye_xyz, x, y) {
  cell <- (extent[2] - extent[1]) / ncol(z)
  z1 <- tin_height(z, extent, x, y)
  d <- sqrt((x - eye_xyz[1])^2 + (y - eye_xyz[2])^2)
  n <- ceiling(d / (cell / 16))
  t <- seq_len(n - 1) / n
  t <- t[d * (1 - t) > cell * 0.05]
  if (!length(t)) return(-Inf)
  zz <- tin_height(z, extent, eye_xyz[1] + t * (x - eye_xyz[1]), eye_xyz[2] + t * (y - eye_xyz[2]))
  max(zz - (eye_xyz[3] + t * (z1 - eye_xyz[3])))
}
point_in_pieces <- function(pc, x, y) {
  for (p in split(pc, pc$piece)) {
    n <- nrow(p); j <- c(2:n, 1)
    s <- sum(p$x * p$y[j] - p$x[j] * p$y)
    cr <- (p$x[j] - p$x) * (y - p$y) - (p$y[j] - p$y) * (x - p$x)
    if (all(cr * sign(s) >= -1e-7)) return(TRUE)
  }
  FALSE
}
tri_id <- function(z, extent, x, y) {
  cell <- (extent[2] - extent[1]) / ncol(z)
  u <- (x - extent[1]) / cell - 0.5; v <- (extent[4] - y) / cell - 0.5
  c0 <- floor(u); r0 <- floor(v)
  (r0 * (ncol(z) - 1) + c0) * 2 + ifelse(u - c0 >= v - r0, 0, 1) + 1
}
