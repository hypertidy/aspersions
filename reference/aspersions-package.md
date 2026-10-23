# aspersions: cast shade and sight lines on terrain

Two kinds of casting over a terrain grid.

## Details

- [`viewshed()`](https://hypertidy.github.io/aspersions/reference/viewshed.md):
  what can be seen from a point. Exact, on the terrain as triangles: the
  visible part of every triangle as a polygon, a visible fraction per
  cell, and the skyline from the point.
  [`viewshed_mesh()`](https://hypertidy.github.io/aspersions/reference/viewshed_mesh.md)
  does the same for any triangulation, in the front-to-back order from
  [`mesh_order()`](https://hypertidy.github.io/aspersions/reference/mesh_order.md).

- [`cast_shade()`](https://hypertidy.github.io/aspersions/reference/cast_shade.md):
  where the terrain shades itself for one sun position, and
  [`sun_hours()`](https://hypertidy.github.io/aspersions/reference/sun_hours.md)
  for the hours of direct sun over a day, with
  [`sun_position()`](https://hypertidy.github.io/aspersions/reference/sun_position.md)
  for the sun.

Grids are a matrix of heights with row 1 at the north edge, square
cells, and an extent `c(xmin, xmax, ymin, ymax)` in a projected
coordinate system in metres; or directly the result of
[`gdalraster::read_ds()`](https://firelab.github.io/gdalraster/reference/read_ds.html),
which carries its own extent. Results are matrices in the same layout,
ready for `ximage::ximage(m, extent = ...)`.

## See also

Useful links:

- <https://github.com/hypertidy/aspersions>

- Report bugs at <https://github.com/hypertidy/aspersions/issues>

## Author

**Maintainer**: Michael Sumner <mdsumner@gmail.com>
