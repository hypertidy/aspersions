# Viewshed over any triangulation

As
[`viewshed()`](https://hypertidy.github.io/aspersions/reference/viewshed.md),
for an explicit mesh: a constrained or adaptive triangulation, or a grid
with other diagonals. Triangles are fed in the front-to-back order from
[`mesh_order()`](https://hypertidy.github.io/aspersions/reference/mesh_order.md)
unless `order` is given. The mesh must cover a convex region.

## Usage

``` r
viewshed_mesh(
  x,
  y,
  z,
  triangles,
  xy,
  eye = 1.5,
  order = NULL,
  curvature = TRUE,
  sectors = 64L
)
```

## Arguments

- x, y, z:

  vertex coordinates and heights

- triangles:

  integer matrix of vertex indices (1-based), one row per triangle

- xy:

  the eye, `c(x, y)`, inside `z` (inside the hull of its cell centres)

- eye:

  height of the eye above the ground (m)

- order:

  front-to-back order of triangles (1-based), or `NULL` to compute it

- curvature:

  lower the terrain by `d^2 / 2R` (earth curvature with standard
  refraction, k = 0.13)

- sectors:

  sectors per 90 degree wedge (more: less bookkeeping per sector)

## Value

a list with `area`, `skyline`, `pieces` (as
[`viewshed()`](https://hypertidy.github.io/aspersions/reference/viewshed.md)),
`eye`, and `cycles` (triangles that had to be forced into order; 0 for
Delaunay meshes)

## Examples

``` r
z <- volcano[nrow(volcano):1, ncol(volcano):1]
## vertices at cell centres, row by row from the north-west corner
xy <- expand.grid(x = 2667400 + seq(5, 605, by = 10), y = 6479570 - seq(5, 865, by = 10))
nc <- ncol(z); nr <- nrow(z)
q <- expand.grid(c = 1:(nc - 1), r = 1:(nr - 1))
tl <- (q$r - 1) * nc + q$c
tri <- rbind(cbind(tl, tl + 1, tl + nc), cbind(tl + 1, tl + nc + 1, tl + nc))
v <- viewshed_mesh(xy$x, xy$y, as.vector(t(z)), tri, c(2667700, 6478900))
v$area
#> [1] 39870.68
```
