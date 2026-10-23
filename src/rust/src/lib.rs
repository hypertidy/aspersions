// aspersions: R bindings to the osvm crate (src/rust/osvm). Thin: grids arrive from R as
// row-major heights (row 0 north) with their origin and cell size; R does the checking
// and shaping (R/*.R), this file only moves numbers.
use extendr_api::prelude::*;
use osvm::{coverage, order_topo, viewshed, viewshed_mesh, Grid, Mesh, ShadeDomain, ViewshedOptions};

fn as_f32(z: &[f64]) -> Vec<f32> {
    z.iter().map(|&v| v as f32).collect()
}

fn pieces_list(res: &osvm::ViewshedResult) -> List {
    let p = &res.pieces;
    list!(
        xy = p.xy.clone(),
        off = p.off.iter().map(|&o| o as f64).collect::<Vec<f64>>(),
        tri = p.tri.iter().map(|&t| t as f64).collect::<Vec<f64>>(),
        tag = p.tag.iter().map(|&t| t as i32).collect::<Vec<i32>>()
    )
}

/// Viewshed over a grid and an optional coarser far grid (empty `far_z` for none).
/// @noRd
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_viewshed(
    z: &[f64], nrow: i32, ncol: i32, x0: f64, ytop: f64, cell: f64,
    far_z: &[f64], far_nrow: i32, far_ncol: i32, far_x0: f64, far_ytop: f64, far_cell: f64,
    e: f64, n: f64, eye: f64, max_dist: f64, curvature: bool, sectors: i32, ss: i32, keep_pieces: bool,
) -> Robj {
    let zn = as_f32(z);
    let zf = as_f32(far_z);
    let near = Grid { z: &zn, w: ncol as usize, h: nrow as usize, cell, x0, ytop };
    let mut rings = vec![near];
    if !zf.is_empty() {
        rings.push(Grid { z: &zf, w: far_ncol as usize, h: far_nrow as usize, cell: far_cell, x0: far_x0, ytop: far_ytop });
    }
    let opts = ViewshedOptions { eye, max_dist, sectors: sectors.max(1) as usize, curvature, keep_pieces: keep_pieces || ss > 0 };
    let res = match viewshed(&rings, e, n, &opts) {
        Some(r) => r,
        None => return Robj::from(()),
    };
    let cov: Vec<f64> = if ss > 0 {
        let k = (ss * ss) as f64;
        coverage(&res, &near, ss as usize, 0).into_iter().map(|c| c as f64 / k).collect()
    } else {
        Vec::new()
    };
    let pieces = if keep_pieces { pieces_list(&res).into_robj() } else { Robj::from(()) };
    list!(
        area = res.area,
        z_eye = res.z_eye,
        n_tri = res.n_tri as f64,
        n_vis_tri = res.n_vis_tri as f64,
        n_quick_hidden = res.n_quick_hidden as f64,
        skyline = res.skyline.clone(),
        visible = cov,
        pieces = pieces
    )
    .into_robj()
}

/// Viewshed over an explicit mesh; `tri` and `order` 0-based; empty `order` = order_topo.
/// @noRd
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_viewshed_mesh(x: &[f64], y: &[f64], z: &[f64], tri: &[i32], e: f64, n: f64, eye: f64, curvature: bool, sectors: i32, order: &[i32]) -> Robj {
    let t: Vec<u32> = tri.iter().map(|&v| v as u32).collect();
    let mesh = Mesh { x, y, z, tri: &t };
    let ord: Vec<u32> = order.iter().map(|&v| v as u32).collect();
    let opts = ViewshedOptions { eye, sectors: sectors.max(1) as usize, curvature, ..Default::default() };
    match viewshed_mesh(&mesh, e, n, if ord.is_empty() { None } else { Some(&ord) }, &opts) {
        None => Robj::from(()),
        Some((res, cycles)) => list!(
            area = res.area,
            z_eye = res.z_eye,
            n_tri = res.n_tri as f64,
            n_vis_tri = res.n_vis_tri as f64,
            cycles = cycles as f64,
            skyline = res.skyline.clone(),
            pieces = pieces_list(&res)
        )
        .into_robj(),
    }
}

/// Front-to-back order of a mesh's triangles from (e, n), 0-based.
/// @noRd
#[extendr]
fn rs_mesh_order(x: &[f64], y: &[f64], tri: &[i32], e: f64, n: f64) -> Robj {
    let t: Vec<u32> = tri.iter().map(|&v| v as u32).collect();
    let z = vec![0.0; x.len()];
    let o = order_topo(&Mesh { x, y, z: &z, tri: &t }, e, n);
    list!(order = o.order.iter().map(|&v| v as f64).collect::<Vec<f64>>(), cycles = o.cycles as f64).into_robj()
}

/// Minutes (or any weight) of direct sun per cell, summed over sun positions; one shadow
/// sweep per position. A single position with weight 1 gives the lit mask.
/// @noRd
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_sun_weight(
    z: &[f64], nrow: i32, ncol: i32, x0: f64, ytop: f64, cell: f64,
    far_z: &[f64], far_nrow: i32, far_ncol: i32, far_x0: f64, far_ytop: f64, far_cell: f64,
    azimuth: &[f64], elevation: &[f64], weight: &[f64], theta: f64,
) -> Vec<f64> {
    let zn = as_f32(z);
    let zf = as_f32(far_z);
    let near = Grid { z: &zn, w: ncol as usize, h: nrow as usize, cell, x0, ytop };
    let far = if zf.is_empty() {
        None
    } else {
        Some(Grid { z: &zf, w: far_ncol as usize, h: far_nrow as usize, cell: far_cell, x0: far_x0, ytop: far_ytop })
    };
    let dom = ShadeDomain::new(near, far);
    let mut acc = vec![0.0f64; zn.len()];
    for k in 0..azimuth.len() {
        if elevation[k] <= 0.0 {
            continue;
        }
        let sh = dom.shade_at(azimuth[k], elevation[k], theta);
        for (a, s) in acc.iter_mut().zip(sh) {
            if s == 0 {
                *a += weight[k];
            }
        }
    }
    acc
}

extendr_module! {
    mod aspersions;
    fn rs_viewshed;
    fn rs_viewshed_mesh;
    fn rs_mesh_order;
    fn rs_sun_weight;
}
