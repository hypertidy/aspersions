# Where the sun is

Sun azimuth and elevation by the NOAA solar calculator equations, with
standard atmospheric refraction.

## Usage

``` r
sun_position(time, lon, lat)
```

## Arguments

- time:

  `POSIXct` times

- lon, lat:

  degrees

## Value

data frame of `azimuth` (degrees clockwise from true north) and
`elevation` (degrees above a flat horizon)

## Examples

``` r
sun_position(as.POSIXct("2026-06-21 12:00", tz = "Australia/Hobart"), 147.3, -42.9)
#>    azimuth elevation
#> 1 3.137254  23.63753
```
