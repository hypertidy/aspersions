test_that("viewshed agrees with brute-force line of sight on the same TIN", {
  g <- volcano_grid()
  set.seed(1)
  for (eye in list(c(300, 200, 1.5), c(612.3, 401.7, 10))) {
    v <- viewshed(g$z, g$extent, eye[1:2], eye = eye[3], curvature = FALSE)
    expect_s3_class(v, "aspersions_viewshed")
    by_tri <- split(v$pieces, v$pieces$triangle)
    n <- 0; agree <- 0
    for (k in 1:400) {
      x <- runif(1, 6, 864); y <- runif(1, 6, 604)
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
  v <- viewshed(g$z, g$extent, c(300, 200))
  expect_equal(dim(v$visible), dim(g$z))
  expect_true(all(v$visible >= 0 & v$visible <= 1))
  expect_named(v$skyline, c("azimuth", "elevation"))
  expect_equal(nrow(v$skyline), 360)
  expect_true(all(v$skyline$elevation > -90))
  expect_named(v$pieces, c("x", "y", "piece", "triangle", "ring"))
  expect_gt(v$area, 0)
  expect_lt(v$area, 870 * 610)
  expect_null(viewshed(g$z, g$extent, c(-100, 200)))
  expect_output(print(v), "aspersions viewshed")
})

test_that("a higher eye sees at least as much", {
  g <- volcano_grid()
  a <- viewshed(g$z, g$extent, c(300, 200), eye = 1.5, pieces = FALSE)$area
  b <- viewshed(g$z, g$extent, c(300, 200), eye = 30, pieces = FALSE)$area
  expect_gt(b, a)
})

test_that("a far ring adds terrain outside the near grid", {
  g <- volcano_grid()
  near <- g$z[21:40, 31:60]; near_ext <- c(300, 600, 210, 410)
  v1 <- viewshed(near, near_ext, c(450, 300), pieces = FALSE)
  v2 <- viewshed(near, near_ext, c(450, 300), far = g$z, far_extent = g$extent, pieces = FALSE)
  expect_gt(v2$area, v1$area)
})

test_that("the mesh path matches the grid path", {
  g <- volcano_grid(); z <- g$z
  xy <- expand.grid(x = seq(5, 865, by = 10), y = seq(605, 5, by = -10))
  nc <- ncol(z); nr <- nrow(z)
  q <- expand.grid(c = 1:(nc - 1), r = 1:(nr - 1)); tl <- (q$r - 1) * nc + q$c
  tri <- rbind(cbind(tl, tl + 1, tl + nc + 1), cbind(tl, tl + nc + 1, tl + nc))
  vg <- viewshed(z, g$extent, c(300.3, 200.7))
  vm <- viewshed_mesh(xy$x, xy$y, as.vector(t(z)), tri, c(300.3, 200.7))
  expect_equal(vm$area, vg$area, tolerance = 1e-9)
  expect_equal(vm$cycles, 0)
  o <- mesh_order(xy$x, xy$y, tri, c(300.3, 200.7))
  expect_equal(sort(o$order), seq_len(nrow(tri)))
  expect_equal(o$cycles, 0L)
  expect_error(viewshed_mesh(xy$x, xy$y, as.vector(t(z)), tri + 1e6, c(300, 200)), "indices")
})

test_that("inputs are checked", {
  g <- volcano_grid()
  expect_error(viewshed(as.vector(g$z), g$extent, c(300, 200)), "matrix")
  expect_error(viewshed(g$z, c(0, 870, 0, 300), c(300, 200)), "square")
  expect_error(viewshed(g$z, g$extent, c(300, 200, 1)), "one point")
  expect_error(viewshed(g$z, g$extent, c(300, 200), far = g$z), "far_extent")
})
