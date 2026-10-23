#' Where the terrain shades itself
#'
#' Cast shadows for one sun position, by a single directional sweep over the grid: a
#' "shadow surface" is carried along the light's direction, dropping by
#' `step * tan(elevation)` each cell, and a cell is in shadow where the surface arriving
#' at it is above the ground. With `far`, a coarse surrounding grid casts first (with
#' earth curvature) and its shadow surface comes in over `z` from distant ridges.
#'
#' @inheritParams viewshed
#' @param azimuth sun azimuth, degrees clockwise from true north
#' @param elevation sun elevation, degrees above the horizon
#' @param convergence grid azimuth of true north (meridian convergence, degrees), added
#'   to `azimuth`; 0 for a grid aligned with true north. In polar stereographic
#'   (EPSG:3031, EPSG:3413 style) it is the longitude, or minus it in the north
#' @return logical matrix like `z`, `TRUE` in shadow (all `TRUE` when the sun is down),
#'   with attribute `extent`. Missing heights cast as 0 m and are missing in the result.
#' @export
#' @examples
#' z <- volcano[nrow(volcano):1, ncol(volcano):1]
#' ex <- c(2667400, 2668010, 6478700, 6479570)
#' s <- cast_shade(z, ex, azimuth = 300, elevation = 15)
#' mean(s)
#' if (requireNamespace("ximage", quietly = TRUE)) {
#'   ximage::ximage(z, extent = ex, col = grey.colors(64), asp = 1)
#'   ximage::ximage(s * 1, extent = ex, add = TRUE, col = c(NA, adjustcolor("navy", 0.5)))
#' }
cast_shade <- function(z, extent = NULL, azimuth, elevation, far = NULL, far_extent = NULL, convergence = 0) {
  if (length(azimuth) != 1L || length(elevation) != 1L) stop("one sun position at a time; see sun_hours() for many", call. = FALSE)
  g <- as_grid(z, extent)
  w <- sun_weight(g$z, g$extent, azimuth, elevation, 1, far, far_extent, convergence)
  m <- matrix(w == 0, nrow(g$z), ncol(g$z), byrow = TRUE)
  m[!is.finite(g$z)] <- NA
  structure(m, extent = g$extent)
}

sun_weight <- function(z, extent, azimuth, elevation, weight, far, far_extent, convergence) {
  g <- grid_args(z, extent)
  f <- far_args(far, far_extent)
  rs_sun_weight(g$z, g$nrow, g$ncol, g$x0, g$ytop, g$cell, f$z, f$nrow, f$ncol, f$x0, f$ytop, f$cell,
                as.double(azimuth), as.double(elevation), as.double(weight), as.double(convergence))
}

#' Hours of direct sun over a day
#'
#' Sun positions every `step` minutes over the local day (at the middle of each step),
#' one shadow sweep each, summed.
#'
#' @inheritParams cast_shade
#' @param date the day, a `Date` or `"YYYY-MM-DD"`, local to `tz`
#' @param lon,lat where the sun is computed from (the middle of the grid, say), degrees
#' @param tz time zone of `date`
#' @param step minutes between sun positions
#' @return matrix like `z` of hours of direct sun (missing where `z` is), with attributes
#'   `daylight` (hours the sun is above a flat horizon) and `extent`
#' @export
#' @examples
#' ## Maungawhau, Auckland
#' z <- volcano[nrow(volcano):1, ncol(volcano):1]
#' ex <- c(2667400, 2668010, 6478700, 6479570)
#' h <- sun_hours(z, ex, "2026-06-21", lon = 174.764, lat = -36.876,
#'                tz = "Pacific/Auckland", step = 20)
#' range(h)
sun_hours <- function(z, extent = NULL, date, lon, lat, tz = "UTC", step = 10, far = NULL, far_extent = NULL, convergence = 0) {
  t0 <- as.POSIXct(paste(format(as.Date(date)), "00:00:00"), tz = tz)
  mins <- seq(step / 2, 1440 - step / 2, by = step)
  sp <- sun_position(t0 + mins * 60, lon, lat)
  g <- as_grid(z, extent)
  w <- sun_weight(g$z, g$extent, sp$azimuth, sp$elevation, rep(step / 60, length(mins)), far, far_extent, convergence)
  m <- matrix(w, nrow(g$z), ncol(g$z), byrow = TRUE)
  m[!is.finite(g$z)] <- NA
  structure(m, daylight = sum(sp$elevation > 0) * step / 60, extent = g$extent)
}

#' Where the sun is
#'
#' Sun azimuth and elevation by the NOAA solar calculator equations, with standard
#' atmospheric refraction.
#'
#' @param time `POSIXct` times
#' @param lon,lat degrees
#' @return data frame of `azimuth` (degrees clockwise from true north) and `elevation`
#'   (degrees above a flat horizon)
#' @export
#' @examples
#' sun_position(as.POSIXct("2026-06-21 12:00", tz = "Australia/Hobart"), 147.3, -42.9)
sun_position <- function(time, lon, lat) {
  rad <- pi / 180
  ms <- as.numeric(as.POSIXct(time)) * 1000
  jd <- ms / 86400000 + 2440587.5; T <- (jd - 2451545) / 36525
  L0 <- (280.46646 + T * (36000.76983 + T * 0.0003032)) %% 360
  M <- 357.52911 + T * (35999.05029 - 0.0001537 * T)
  e <- 0.016708634 - T * (0.000042037 + 0.0000001267 * T)
  Mr <- M * rad
  C <- sin(Mr) * (1.914602 - T * (0.004817 + 0.000014 * T)) + sin(2 * Mr) * (0.019993 - 0.000101 * T) + sin(3 * Mr) * 0.000289
  omega <- (125.04 - 1934.136 * T) * rad
  lambda <- (L0 + C - 0.00569 - 0.00478 * sin(omega)) * rad
  eps0 <- 23 + (26 + (21.448 - T * (46.815 + T * (0.00059 - T * 0.001813))) / 60) / 60
  eps <- (eps0 + 0.00256 * cos(omega)) * rad
  decl <- asin(sin(eps) * sin(lambda))
  y <- tan(eps / 2)^2; L0r <- L0 * rad
  eot <- 4 / rad * (y * sin(2 * L0r) - 2 * e * sin(Mr) + 4 * e * y * sin(Mr) * cos(2 * L0r) -
                      0.5 * y * y * sin(4 * L0r) - 1.25 * e * e * sin(2 * Mr))
  utc_min <- (ms %% 86400000) / 60000
  ha <- ((utc_min + eot + 4 * lon) / 4 - 180) * rad; phi <- lat * rad
  cz <- pmax(-1, pmin(1, sin(phi) * sin(decl) + cos(phi) * cos(decl) * cos(ha)))
  elev <- 90 - acos(cz) / rad
  az <- atan2(sin(ha), cos(ha) * sin(phi) - tan(decl) * cos(phi)) / rad + 180
  up <- elev > -1
  elev[up] <- elev[up] + 1.02 / tan((elev[up] + 10.3 / (elev[up] + 5.11)) * rad) / 60
  data.frame(azimuth = (az + 360) %% 360, elevation = elev)
}
