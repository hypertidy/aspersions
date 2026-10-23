//! Visible pieces back onto a grid, as a sampled fraction per cell.

use crate::grid::Grid;
use crate::horizon::ViewshedResult;

/// Visible fraction of each cell of `g` from the pieces tagged `tag`, sampled on an
/// `ss` x `ss` grid of points per cell: counts 0 ..= ss * ss, row-major like `g.z`.
pub fn coverage(res: &ViewshedResult, g: &Grid, ss: usize, tag: u8) -> Vec<u8> {
    let ss = ss.max(1);
    let cap = (ss * ss).min(255) as u8;
    let step = g.cell / ss as f64;
    let mut out = vec![0u8; g.w * g.h];
    let p = &res.pieces;
    for k in 0..p.len() {
        if p.tag[k] != tag {
            continue;
        }
        let poly = p.polygon(k);
        let n = poly.len() / 2;
        let (mut xmin, mut xmax, mut ymin, mut ymax, mut s2) = (f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY, 0.0);
        for i in 0..n {
            let (x, y) = (poly[2 * i], poly[2 * i + 1]);
            xmin = xmin.min(x);
            xmax = xmax.max(x);
            ymin = ymin.min(y);
            ymax = ymax.max(y);
            let j = if i + 1 == n { 0 } else { i + 1 };
            s2 += x * poly[2 * j + 1] - poly[2 * j] * y;
        }
        let sgn = if s2 >= 0.0 { 1.0 } else { -1.0 };
        let k0 = ((xmin - g.x0) / step - 0.5).ceil().max(0.0) as isize;
        let k1 = (((xmax - g.x0) / step - 0.5).floor() as isize).min((g.w * ss) as isize - 1);
        let m0 = ((g.ytop - ymax) / step - 0.5).ceil().max(0.0) as isize;
        let m1 = (((g.ytop - ymin) / step - 0.5).floor() as isize).min((g.h * ss) as isize - 1);
        for m in m0..=m1 {
            let ny = g.ytop - (m as f64 + 0.5) * step;
            for kx in k0..=k1 {
                let ex = g.x0 + (kx as f64 + 0.5) * step;
                let inside = (0..n).all(|i| {
                    let j = if i + 1 == n { 0 } else { i + 1 };
                    let cr =
                        (poly[2 * j] - poly[2 * i]) * (ny - poly[2 * i + 1]) - (poly[2 * j + 1] - poly[2 * i + 1]) * (ex - poly[2 * i]);
                    cr * sgn >= -1e-9
                });
                if inside {
                    let o = (m as usize / ss) * g.w + kx as usize / ss;
                    if out[o] < cap {
                        out[o] += 1;
                    }
                }
            }
        }
    }
    out
}
