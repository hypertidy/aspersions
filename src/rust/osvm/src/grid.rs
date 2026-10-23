//! A north-up grid as a source of triangles.

use crate::horizon::Horizon;

/// A north-up grid: `z` has `w * h` heights, row 0 at the top (north), square cells of
/// `cell` metres, (`x0`, `ytop`) the outer corner of the top-left cell. As a TIN, the
/// vertices are the cell centres and each quad of four neighbouring centres is split on
/// the diagonal from top-left to bottom-right.
#[derive(Clone, Copy, Debug)]
pub struct Grid<'a> {
    pub z: &'a [f32],
    pub w: usize,
    pub h: usize,
    pub cell: f64,
    pub x0: f64,
    pub ytop: f64,
}

/// Height on the grid's TIN at (e, n), or `None` outside the hull of the cell centres.
pub fn sample_tin(g: &Grid, e: f64, n: f64) -> Option<f64> {
    let u = (e - g.x0) / g.cell - 0.5;
    let v = (g.ytop - n) / g.cell - 0.5;
    if u < 0.0 || v < 0.0 || u > (g.w - 1) as f64 || v > (g.h - 1) as f64 || g.w < 2 || g.h < 2 {
        return None;
    }
    let c = (u.floor() as usize).min(g.w - 2);
    let r = (v.floor() as usize).min(g.h - 2);
    let (fx, fy, w) = (u - c as f64, v - r as f64, g.w);
    let z = |i: usize| g.z[i] as f64;
    let (tl, tr, bl, br) = (z(r * w + c), z(r * w + c + 1), z((r + 1) * w + c), z((r + 1) * w + c + 1));
    Some(if fx >= fy {
        tl + fx * (tr - tl) + fy * (br - tr)
    } else {
        tl + fy * (bl - tl) + fx * (br - bl)
    })
}

/// Options for [`grid_source`].
#[derive(Clone, Copy, Debug)]
pub struct GridSourceOptions<'a> {
    /// A finer grid whose vertex hull is skipped (already fed, and nearer).
    pub hole: Option<&'a Grid<'a>>,
    /// Skip quads further than this (metres, plan).
    pub max_dist: f64,
    /// Carried to the pieces.
    pub tag: u8,
}
impl Default for GridSourceOptions<'_> {
    fn default() -> Self {
        GridSourceOptions {
            hole: None,
            max_dist: f64::INFINITY,
            tag: 0,
        }
    }
}

/// Feed every triangle of `g` to `h` front to back from the eye at (e, n).
///
/// The order: the quads a plan ray from the eye crosses on its way to quad Q have row and
/// column offsets from the eye's quad no larger than Q's (on the same side), so quads in
/// increasing `|dcol| + |drow|` are front to back, with no sort. Inside a quad the triangle
/// on the eye's side of the diagonal goes first. Triangle ids are `quad * 2 + half` with
/// `quad = row * (w - 1) + col`, half 0 = (TL, TR, BR), 1 = (TL, BR, BL).
pub fn grid_source(h: &mut Horizon, g: &Grid, e: f64, n: f64, opts: &GridSourceOptions) {
    let (w, hh, cell) = (g.w, g.h, g.cell);
    if w < 2 || hh < 2 {
        return;
    }
    let hole = opts.hole.map(|i| {
        (
            i.x0 + 0.5 * i.cell,
            i.x0 + (i.w as f64 - 0.5) * i.cell,
            i.ytop - (i.h as f64 - 0.5) * i.cell,
            i.ytop - 0.5 * i.cell,
        )
    });
    let pc = (e - g.x0) / cell - 0.5;
    let pr = (g.ytop - n) / cell - 0.5;
    let cp = (pc.floor() as isize).clamp(0, w as isize - 2);
    let rp = (pr.floor() as isize).clamp(0, hh as isize - 2);
    let (wm, hm) = (w as isize - 2, hh as isize - 2);
    let kmax = cp.max(wm - cp) + rp.max(hm - rp);
    let max_q = opts.max_dist / cell + 1.0;
    let md2 = opts.max_dist * opts.max_dist;
    let z = |i: usize| g.z[i] as f64;
    for kk in 0..=kmax {
        if (kk as f64) > max_q * 1.4143 + 2.0 {
            break;
        }
        let (dlo, dhi) = ((-kk).max(-cp), kk.min(wm - cp));
        for di in dlo..=dhi {
            let c = cp + di;
            let rem = kk - di.abs();
            let sides = if rem != 0 { 2 } else { 1 };
            for side in 0..sides {
                let r = if side == 1 { rp + rem } else { rp - rem };
                if r < 0 || r > hm {
                    continue;
                }
                let xl = g.x0 + (c as f64 + 0.5) * cell - e;
                let xr = xl + cell;
                let yt = g.ytop - (r as f64 + 0.5) * cell - n;
                let yb = yt - cell;
                if let Some((hx0, hx1, hy0, hy1)) = hole {
                    if xr + e >= hx0 && xl + e <= hx1 && yt + n >= hy0 && yb + n <= hy1 {
                        continue;
                    }
                }
                if opts.max_dist.is_finite() {
                    let mx = ((xl + 0.5 * cell).abs() - 0.5 * cell).max(0.0);
                    let my = ((yt - 0.5 * cell).abs() - 0.5 * cell).max(0.0);
                    if mx * mx + my * my > md2 {
                        continue;
                    }
                }
                let (c, r) = (c as usize, r as usize);
                let i0 = r * w + c;
                let (ztl, ztr, zbl, zbr) = (z(i0), z(i0 + 1), z(i0 + w), z(i0 + w + 1));
                let qi = ((r * (w - 1) + c) * 2) as u32;
                if (pc - c as f64) >= (pr - r as f64) {
                    h.tri(qi, opts.tag, xl, yt, ztl, xr, yt, ztr, xr, yb, zbr);
                    h.tri(qi + 1, opts.tag, xl, yt, ztl, xr, yb, zbr, xl, yb, zbl);
                } else {
                    h.tri(qi + 1, opts.tag, xl, yt, ztl, xr, yb, zbr, xl, yb, zbl);
                    h.tri(qi, opts.tag, xl, yt, ztl, xr, yt, ztr, xr, yb, zbr);
                }
            }
        }
    }
}
