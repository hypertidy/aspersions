# brute-force shadow by ray marching on the TIN
brute_shade <- function(z, extent, azimuth, elevation, x, y) {
  cell <- (extent[2] - extent[1]) / ncol(z)
  a <- azimuth * pi / 180; tn <- tan(elevation * pi / 180)
  z0 <- tin_height(z, extent, x, y)
  d <- seq(cell / 2, 2000, by = cell / 4)
  xs <- x + sin(a) * d; ys <- y + cos(a) * d
  ok <- xs > extent[1] + cell / 2 & xs < extent[2] - cell / 2 & ys > extent[3] + cell / 2 & ys < extent[4] - cell / 2
  if (!any(ok)) return(FALSE)
  any(tin_height(z, extent, xs[ok], ys[ok]) > z0 + d[ok] * tn + 1e-3)
}

test_that("cast_shade agrees with brute-force ray marching", {
  g <- volcano_grid(); cell <- 10
  set.seed(2)
  for (sun in list(c(300, 15), c(45, 8), c(180, 25))) {
    s <- cast_shade(g$z, g$extent, sun[1], sun[2])
    expect_equal(dim(s), dim(g$z))
    rc <- cbind(sample(2:(nrow(g$z) - 1), 300, TRUE), sample(2:(ncol(g$z) - 1), 300, TRUE))
    b <- apply(rc, 1, function(k) brute_shade(g$z, g$extent, sun[1], sun[2],
                                              g$extent[1] + (k[2] - 0.5) * cell, g$extent[4] - (k[1] - 0.5) * cell))
    expect_gt(mean(s[rc] == b), 0.95)
  }
})

test_that("no sun, all shade; high sun, little", {
  g <- volcano_grid()
  expect_true(all(cast_shade(g$z, g$extent, 0, -5)))
  expect_lt(mean(cast_shade(g$z, g$extent, 0, 80)), 0.05)
  expect_error(cast_shade(g$z, g$extent, c(1, 2), 10), "one sun position")
})

test_that("sun_position matches known values", {
  # Hobart, midwinter noon: the sun is low in the north, about 23.6 degrees up
  p <- sun_position(as.POSIXct("2026-06-21 12:00", tz = "Australia/Hobart"), 147.3, -42.9)
  expect_lt(abs(p$elevation - 23.6), 0.5)
  expect_true(p$azimuth < 10 || p$azimuth > 350)
  # midnight: below the horizon
  expect_lt(sun_position(as.POSIXct("2026-06-21 00:00", tz = "Australia/Hobart"), 147.3, -42.9)$elevation, 0)
})

test_that("sun_hours is bounded by daylight", {
  g <- volcano_grid()
  h <- sun_hours(g$z, g$extent, "2026-06-21", lon = 174.764, lat = -36.876, tz = "Pacific/Auckland", step = 20)
  expect_equal(dim(h), dim(g$z))
  expect_true(all(h >= 0 & h <= attr(h, "daylight") + 1e-9))
  expect_gt(attr(h, "daylight"), 9)
  expect_lt(attr(h, "daylight"), 11)
})
