# Front-to-back order of a mesh's triangles

Along any plan ray from the eye, the triangles crossed are neighbours
across shared edges, each crossed from the side facing the eye. So "the
triangle on the eye's side of a shared edge goes first", over all shared
edges, is the whole in-front relation; this is a topological sort of it.
Delaunay meshes never have a cycle; other meshes can, and triangles left
in a cycle are appended nearest first and counted.

## Usage

``` r
mesh_order(x, y, triangles, xy)
```

## Arguments

- x, y:

  vertex coordinates

- triangles:

  integer matrix of vertex indices (1-based), one row per triangle

- xy:

  the eye, `c(x, y)`, inside `z` (inside the hull of its cell centres)

## Value

a list with `order` (1-based triangle indices) and `cycles`
