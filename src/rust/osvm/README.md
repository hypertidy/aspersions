# osvm

*Vendored into aspersions from [mdsumner/lunashade](https://github.com/mdsumner/lunashade) `rust/osvm`, where the
parity tests against the JavaScript live. Relative links below point into that repository.*

Viewsheds as an object-space visibility map: the exact visible part of every triangle of a terrain mesh, seen
from a point, as convex polygons in plan, plus the skyline from that point. No dependencies.

A port of [`js/osvm.js`](../../js/osvm.js) from [lunashade](https://github.com/mdsumner/lunashade); the method
is written up in [docs/viewshed.md](../../docs/viewshed.md). The two implementations are kept in step and
checked against each other: `node test/make-parity.js` (from the repository root) writes
`tests/fixtures/parity.txt`, and `cargo test` checks the crate against it. At the time of writing they agree to
the last digit of the visible area, and on every count.

## Layers

- `Horizon`, the kernel: triangles in, nearest first; area, pieces and skyline out. It knows nothing about
  where the triangles come from.
- Sources: `grid_source` (a north-up grid, vertices at cell centres, walked by Manhattan offset from the eye:
  a valid front-to-back order with no sort) and `mesh_source` (any triangulation as vertex and index arrays).
- Orders: `order_topo` for any mesh of a convex region. It does a topological sort of "the triangle on the
  eye's side of each shared edge goes first", and counts triangles caught in cycles. Delaunay meshes have none.
- `viewshed` (nested grids, near to far) and `viewshed_mesh` put them together; `coverage` samples the pieces
  back onto a grid.

```rust
use osvm::{viewshed, Grid, ViewshedOptions};
let g = Grid { z: &heights, w, h, cell: 4.0, x0, ytop };
let res = viewshed(&[g], e, n, &ViewshedOptions { eye: 1.5, ..Default::default() }).unwrap();
println!("{} m2 visible, skyline due north {} deg", res.area, res.skyline[0]);
```

`cargo run --release --example bench` times the app-sized case (5.1 M triangles): about 1.4 s, against
about 2.5 s for the JavaScript in node.

## Where it is going

- Bindings: R (extendr), and WebAssembly for lunashade in the browser. The R side takes meshes as plain
  vertex and triangle tables, so mesh producers (constrained Delaunay, CGAL, incremental Delaunay) and this
  kernel stay independent.
- Level of detail by angle, not metres: an incremental Delaunay mesh refined where the error, seen from
  the eye, is largest.
