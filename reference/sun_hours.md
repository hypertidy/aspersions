# Hours of direct sun over a day

Sun positions every `step` minutes over the local day (at the middle of
each step), one shadow sweep each, summed.

## Usage

``` r
sun_hours(
  z,
  extent = NULL,
  date,
  lon,
  lat,
  tz = "UTC",
  step = 10,
  far = NULL,
  far_extent = NULL,
  convergence = 0
)
```

## Arguments

- z:

  matrix of heights, row 1 at the north edge, square cells; or the
  result of
  [`gdalraster::read_ds()`](https://firelab.github.io/gdalraster/reference/read_ds.html)
  (a vector with a `gis` attribute), in which case `extent` is taken
  from it

- extent:

  `c(xmin, xmax, ymin, ymax)` of `z`

- date:

  the day, a `Date` or `"YYYY-MM-DD"`, local to `tz`

- lon, lat:

  where the sun is computed from (the middle of the grid, say), degrees

- tz:

  time zone of `date`

- step:

  minutes between sun positions

- far, far_extent:

  optional coarser grid surrounding `z` (and its extent), for terrain
  further out; its cells that touch `z` are skipped

- convergence:

  grid azimuth of true north (meridian convergence, degrees), added to
  `azimuth`; 0 for a grid aligned with true north. In polar
  stereographic (EPSG:3031, EPSG:3413 style) it is the longitude, or
  minus it in the north

## Value

matrix like `z` of hours of direct sun (missing where `z` is), with
attributes `daylight` (hours the sun is above a flat horizon) and
`extent`

## Examples

``` r
## Maungawhau, Auckland
z <- volcano[nrow(volcano):1, ncol(volcano):1]
ex <- c(2667400, 2668010, 6478700, 6479570)
h <- sun_hours(z, ex, "2026-06-21", lon = 174.764, lat = -36.876,
               tz = "Pacific/Auckland", step = 20)
range(h)
#> [1] 0.000000 9.666667
```
