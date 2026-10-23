//! # osvm: viewsheds as an object-space visibility map
//!
//! The ground visible from a point, computed exactly on a triangulated terrain: the
//! answer is the visible part of every triangle, as convex polygons in plan, plus the
//! skyline seen from the point.
//!
//! Three layers, kept apart so that mesh producers and the visibility kernel do not
//! depend on each other:
//!
//! * [`Horizon`], the kernel. Triangles go in, nearest first; it knows nothing about
//!   where they come from.
//! * Sources that feed it: [`grid_source`] (a north-up grid, vertices at cell centres)
//!   and [`mesh_source`] (any triangulation as vertex and index arrays).
//! * Orders: the Manhattan walk built into [`grid_source`], and [`order_topo`] for any
//!   mesh (a topological sort of "the triangle on the eye's side of each shared edge
//!   goes first", with cycles detected).
//!
//! [`viewshed`] and [`viewshed_mesh`] put them together.
//!
//! Cast shadows (a viewshed from the sun, where the rays are parallel) have their own,
//! much faster grid algorithm: [`sweep`], one directional pass per sun position, and
//! [`ShadeDomain`] for a fine near grid with a coarse far field.
//!
//! ## How the kernel works
//!
//! 1. *Gnomonic wedges.* Around the eye, four 90 degree wedges of azimuth (N, E, S, W),
//!    each with the central projection onto a vertical plane one unit ahead:
//!    `u = right.q / forward.q`, `v = z / forward.q`. Straight 3D edges stay straight, so
//!    a projected triangle is still a triangle; at fixed `u`, `v` increases with
//!    elevation angle. Each wedge is cut into narrow sectors of equal width in `u`; all
//!    cuts are vertical planes through the eye, so clipping is done in plan.
//! 2. *Front to back.* For a height field seen from above it, anything that can hide a
//!    point lies nearer along the same plan ray, so a front-to-back order exists for
//!    Delaunay meshes (Edelsbrunner's acyclicity theorem) and is cheap on grids.
//! 3. *The horizon is an upper envelope.* Everything drawn so far covers the view from
//!    the bottom up to a piecewise-linear horizon `v = E(u)` per sector. A triangle is
//!    visible exactly where it rises above `E`; then `E` becomes `max(E, its top edge)`.
//! 4. *Back to the ground.* Central projection between planes is projective, so each
//!    visible piece maps back to a convex polygon on its triangle.
//!
//! Heights can be lowered by `d^2 / 2R` for earth curvature (with standard refraction,
//! k = 0.13) before projection.
//!
//! This crate is a port of `js/osvm.js` in lunashade; the two are kept in step and
//! checked against each other (`tests/parity.rs`).

mod coverage;
mod grid;
mod horizon;
mod mesh;
mod shade;

pub use coverage::coverage;
pub use grid::{grid_source, sample_tin, Grid, GridSourceOptions};
pub use horizon::{Horizon, HorizonOptions, Pieces, ViewshedResult, R_EFF};
pub use mesh::{mesh_source, order_topo, sample_mesh, Mesh, TopoOrder};
pub use shade::{sweep, ShadeDomain, SweepOptions};

/// Options for [`viewshed`] and [`viewshed_mesh`].
#[derive(Clone, Debug)]
pub struct ViewshedOptions {
    /// Eye height above the ground, metres.
    pub eye: f64,
    /// Ignore terrain further than this (metres, plan distance). Grids only.
    pub max_dist: f64,
    /// Sectors per 90 degree wedge.
    pub sectors: usize,
    /// Lower heights by d^2 / 2R (with refraction) before projection.
    pub curvature: bool,
    /// Keep the visible polygons (otherwise only area, counts and skyline).
    pub keep_pieces: bool,
}
impl Default for ViewshedOptions {
    fn default() -> Self {
        ViewshedOptions {
            eye: 1.5,
            max_dist: f64::INFINITY,
            sectors: 64,
            curvature: true,
            keep_pieces: true,
        }
    }
}

/// Viewshed from (e, n) over nested grids: `rings[0]` is the finest and contains the eye,
/// each next ring is coarser and contains the previous one; quads of ring k > 0 touching
/// ring k - 1 are skipped. A plan ray leaves a convex ring once, so everything in a ring
/// is in front of everything outside it. Returns `None` if the eye is outside `rings[0]`.
/// Piece tags are the ring index.
pub fn viewshed(rings: &[Grid], e: f64, n: f64, opts: &ViewshedOptions) -> Option<ViewshedResult> {
    let zg = sample_tin(rings.first()?, e, n)?;
    let mut h = Horizon::new(
        zg + opts.eye,
        HorizonOptions {
            sectors: opts.sectors,
            curvature: opts.curvature,
            keep_pieces: opts.keep_pieces,
            e,
            n,
        },
    );
    for (ri, g) in rings.iter().enumerate() {
        let hole = if ri > 0 { Some(&rings[ri - 1]) } else { None };
        grid_source(
            &mut h,
            g,
            e,
            n,
            &GridSourceOptions {
                hole,
                max_dist: opts.max_dist,
                tag: ri as u8,
            },
        );
    }
    Some(h.result())
}

/// Viewshed from (e, n) over an explicit mesh, in the order from [`order_topo`] unless
/// `order` is given. Returns the result and the number of triangles that had to be
/// forced into order because of cycles (0 for Delaunay meshes).
pub fn viewshed_mesh(mesh: &Mesh, e: f64, n: f64, order: Option<&[u32]>, opts: &ViewshedOptions) -> Option<(ViewshedResult, usize)> {
    let zg = sample_mesh(mesh, e, n)?;
    let (ord, cycles) = match order {
        Some(o) => (o.to_vec(), 0),
        None => {
            let t = order_topo(mesh, e, n);
            (t.order, t.cycles)
        }
    };
    let mut h = Horizon::new(
        zg + opts.eye,
        HorizonOptions {
            sectors: opts.sectors,
            curvature: opts.curvature,
            keep_pieces: opts.keep_pieces,
            e,
            n,
        },
    );
    mesh_source(&mut h, mesh, &ord, e, n, 0);
    Some((h.result(), cycles))
}
