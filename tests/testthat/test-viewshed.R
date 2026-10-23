test_that("viewshed agrees with brute-force line of sight on the same TIN", {
  g <- volcano_grid()
  set.seed(1)
  for (eye in list(c(300, 300, 1.5), c(401.7, 612.3, 10))) {
    v <- viewshed(g$z, g$extent, at(g, eye[1], eye[2]), eye = eye[3], curvature = FALSE)
    expect_s3_class(v, "aspersions_viewshed")
    by_tri <- split(v$pieces, v$pieces$triangle)
    n <- 0; agree <- 0
    for (k in 1:400) {
      x <- g$extent[1] + runif(1, 6, 604); y <- g$extent[3] + runif(1, 6, 864)
      m <- los_margin(g$z, g$extent, v$eye, x, y)
      if (abs(m) < 0.05) next
      pc <- by_tri[[as.character(tri_id(g$z, g$extent, x, y))]]
      ov <- !is.null(pc) && point_in_pieces(pc, x, y)
      n <- n + 1; agree <- agree + ((m < 0) == ov)
    }
    expect_gt(agree / n, 0.99)
  }
})

test_that("viewshed result has the documented shape", {
  g <- volcano_grid()
  v <- viewshed(g$z, g$extent, at(g, 300, 300))
  expect_equal(dim(v$visible), dim(g$z))
  expect_equal(v$extent, g$extent)
  expect_true(all(v$visible >= 0 & v$visible <= 1))
  expect_named(v$skyline, c("azimuth", "elevation"))
  expect_equal(nrow(v$skyline), 360)
  expect_true(all(v$skyline$elevation > -90))
  expect_named(v$pieces, c("x", "y", "piece", "triangle", "ring"))
  expect_gt(v$area, 0)
  expect_lt(v$area, 870 * 610)
  expect_null(viewshed(g$z, g$extent, at(g, -100, 200)))
  expect_output(print(v), "aspersions viewshed")
})

test_that("a higher eye sees at least as much", {
  g <- volcano_grid()
  a <- viewshed(g$z, g$extent, at(g, 300, 300), eye = 1.5, pieces = FALSE)$area
  b <- viewshed(g$z, g$extent, at(g, 300, 300), eye = 30, pieces = FALSE)$area
  expect_gt(b, a)
})

test_that("a far ring adds terrain outside the near grid", {
  g <- volcano_grid()
  ## rows 31:60 and columns 21:40 of the 87 x 61 grid
  near <- g$z[31:60, 21:40]
  near_ext <- c(g$extent[1] + 200, g$extent[1] + 400, g$extent[4] - 600, g$extent[4] - 300)
  eye <- c(g$extent[1] + 300, g$extent[4] - 450)
  v1 <- viewshed(near, near_ext, eye, pieces = FALSE)
  v2 <- viewshed(near, near_ext, eye, far = g$z, far_extent = g$extent, pieces = FALSE)
  expect_gt(v2$area, v1$area)
})

test_that("the mesh path matches the grid path", {
  g <- volcano_grid(); z <- g$z
  xy <- expand.grid(x = g$extent[1] + seq(5, 605, by = 10), y = g$extent[4] - seq(5, 865, by = 10))
  nc <- ncol(z); nr <- nrow(z)
  q <- expand.grid(c = 1:(nc - 1), r = 1:(nr - 1)); tl <- (q$r - 1) * nc + q$c
  tri <- rbind(cbind(tl, tl + 1, tl + nc + 1), cbind(tl, tl + nc + 1, tl + nc))
  e <- at(g, 300.3, 200.7)
  vg <- viewshed(z, g$extent, e)
  vm <- viewshed_mesh(xy$x, xy$y, as.vector(t(z)), tri, e)
  expect_equal(vm$area, vg$area, tolerance = 1e-9)
  expect_equal(vm$cycles, 0)
  o <- mesh_order(xy$x, xy$y, tri, e)
  expect_equal(sort(o$order), seq_len(nrow(tri)))
  expect_equal(o$cycles, 0L)
  expect_error(viewshed_mesh(xy$x, xy$y, as.vector(t(z)), tri + 1e6, e), "indices")
})

test_that("inputs are checked", {
  g <- volcano_grid(); e <- at(g, 300, 300)
  expect_error(viewshed(as.vector(g$z), g$extent, e), "matrix")
  expect_error(viewshed(g$z, xy = e), "extent")
  expect_error(viewshed(g$z, g$extent + c(0, 0, 0, -300), e), "square")
  expect_error(viewshed(g$z, g$extent, c(e, 1)), "one point")
  expect_error(viewshed(g$z, g$extent, e, far = g$z), "far_extent")
})

test_that("a gdalraster::read_ds() style vector works in place of matrix + extent", {
  g <- volcano_grid()
  v <- as.vector(t(g$z))
  attr(v, "gis") <- list(type = "raster", bbox = g$extent[c(1, 3, 2, 4)], dim = c(ncol(g$z), nrow(g$z), 1))
  e <- at(g, 300, 300)
  a <- viewshed(v, xy = e, pieces = FALSE)
  b <- viewshed(g$z, g$extent, e, pieces = FALSE)
  expect_equal(a$area, b$area)
  expect_equal(a$visible, b$visible)
  expect_equal(unclass(cast_shade(v, azimuth = 300, elevation = 15)), unclass(cast_shade(g$z, g$extent, 300, 15)))
})

test_that("missing heights are missing in the results", {
  g <- volcano_grid(); z <- g$z
  z[1:10, 1:10] <- NA
  v <- viewshed(z, g$extent, at(g, 300, 300), pieces = FALSE)
  expect_true(all(is.na(v$visible[1:10, 1:10])))
  expect_false(anyNA(v$visible[-(1:10), ]))
  s <- cast_shade(z, g$extent, 300, 15)
  expect_true(all(is.na(s[1:10, 1:10])))
  expect_false(anyNA(s[-(1:10), ]))
})
