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
#'   to `azimuth`; 0 for a grid aligned with true north
#' @return logical matrix like `z`, `TRUE` in shadow (all `TRUE` when the sun is down)
#' @export
#' @examples
#' z <- t(volcano)[ncol(volcano):1, ]
#' s <- cast_shade(z, c(0, 870, 0, 610), azimuth = 300, elevation = 15)
#' mean(s)
cast_shade <- function(z, extent, azimuth, elevation, far = NULL, far_extent = NULL, convergence = 0) {
  if (length(azimuth) != 1L || length(elevation) != 1L) stop("one sun position at a time; see sun_hours() for many", call. = FALSE)
  w <- sun_weight(z, extent, azimuth, elevation, 1, far, far_extent, convergence)
  matrix(w == 0, nrow(z), ncol(z), byrow = TRUE)
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
#' @return matrix like `z` of hours of direct sun, with attribute `daylight` (hours the
#'   sun is above a flat horizon)
#' @export
#' @examples
#' z <- t(volcano)[ncol(volcano):1, ]
#' h <- sun_hours(z, c(0, 870, 0, 610), "2026-06-21", lon = 174.77, lat = -36.88,
#'                tz = "Pacific/Auckland", step = 20)
#' range(h)
sun_hours <- function(z, extent, date, lon, lat, tz = "UTC", step = 10, far = NULL, far_extent = NULL, convergence = 0) {
  t0 <- as.POSIXct(paste(format(as.Date(date)), "00:00:00"), tz = tz)
  mins <- seq(step / 2, 1440 - step / 2, by = step)
  sp <- sun_position(t0 + mins * 60, lon, lat)
  w <- sun_weight(z, extent, sp$azimuth, sp$elevation, rep(step / 60, length(mins)), far, far_extent, convergence)
  structure(matrix(w, nrow(z), ncol(z), byrow = TRUE), daylight = sum(sp$elevation > 0) * step / 60)
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
