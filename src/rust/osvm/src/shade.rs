//! Cast shadows on a grid for one sun position: a single directional sweep, O(N), with an
//! optional far field. A port of `sweep()` and `Domain` in lunashade's `js/engine.js`.
//!
//! Light travels away from the sun. Visit the grid line by line along the light's
//! direction (columns or rows, whichever the light crosses faster); each cell's upstream
//! neighbour lies one line back at a fractional offset. Carry H, the height of the
//! "shadow surface": the lowest height a ray from the sun can pass at that cell without
//! having touched terrain upstream. One line downstream the surface drops by
//! `step * tan(elevation)` (plus the earth-curvature drop, if asked); a cell is in shadow
//! if the surface arriving at it is above the ground, and the surface leaving it is
//! `max(ground, arriving surface)`.
//!
//! Storage is f32 like the JavaScript (Float32Array), so the two agree cell for cell.

use crate::grid::Grid;
use crate::horizon::R_EFF;

const NEG: f64 = -1e9;
const RAD: f64 = std::f64::consts::PI / 180.0;

/// Options for [`sweep`].
#[derive(Clone, Copy, Debug, Default)]
pub struct SweepOptions<'a> {
    /// An incoming far-field shadow surface (w * h), or `None`.
    pub hfar: Option<&'a [f32]>,
    /// Carry distance-from-caster and add the earth-curvature drop.
    pub curvature: bool,
    /// Return the arriving surface (upstream casters only) rather than the leaving one.
    pub arrive: bool,
}

/// Shadow mask (1 = in shadow) and shadow surface for a grid of heights `z` (w x h,
/// row 0 north, square cells of `cell` m), for the sun at grid azimuth `az_grid`
/// (degrees clockwise from grid north) and elevation `elev` (degrees).
pub fn sweep(z: &[f32], w: usize, h: usize, cell: f64, az_grid: f64, elev: f64, opts: &SweepOptions) -> (Vec<u8>, Vec<f32>) {
    let n = w * h;
    let mut shade = vec![0u8; n];
    let mut hs = vec![0f32; n];
    let mut arrive = if opts.arrive { vec![0f32; n] } else { Vec::new() };
    let curv = opts.curvature;
    let mut dd = if curv { vec![0f32; n] } else { Vec::new() };
    if elev <= 0.0 {
        shade.fill(1);
        if opts.arrive {
            arrive.fill(1e9);
            return (shade, arrive);
        }
        return (shade, hs);
    }
    let a = az_grid * RAD;
    let (px, py) = (-a.sin(), a.cos()); // light direction in (col, row); rows run south
    let col_major = px.abs() >= py.abs();
    let (n1, n2, s1, d2, st1, st2): (usize, usize, isize, f64, usize, usize) = if col_major {
        (w, h, if px > 0.0 { 1 } else { -1 }, py / px.abs(), 1, w)
    } else {
        (h, w, if py > 0.0 { 1 } else { -1 }, px / py.abs(), w, 1)
    };
    let step = cell * (1.0 + d2 * d2).sqrt();
    let drop = step * (elev * RAD).tan();
    let mut i1: isize = if s1 > 0 { 0 } else { n1 as isize - 1 };
    // first line: nothing upstream but the far field
    for i2 in 0..n2 {
        let idx = i1 as usize * st1 + i2 * st2;
        let zi = z[idx] as f64;
        let cand = opts.hfar.map_or(NEG, |f| f[idx] as f64);
        if opts.arrive {
            arrive[idx] = cand as f32;
        }
        let sh = cand > zi;
        shade[idx] = sh as u8;
        hs[idx] = (if sh { cand } else { zi }) as f32;
        if curv {
            dd[idx] = 0.0;
        }
    }
    let j0 = (-d2).floor();
    let f = -d2 - j0;
    let j0 = j0 as isize;
    for _ in 1..n1 {
        let prev = i1 as usize * st1;
        i1 += s1;
        let cur = i1 as usize * st1;
        for i2 in 0..n2 {
            let idx = cur + i2 * st2;
            let zi = z[idx] as f64;
            let ja = i2 as isize + j0;
            let jb = ja + 1;
            let mut dp = 0.0f64;
            let hp: f64;
            let at = |j: isize| prev + j as usize * st2;
            if ja >= 0 && (jb as usize) < n2 {
                // take the higher of the two upstream surfaces where one is far lower
                // (avoids smearing a real surface into "no surface"); otherwise interpolate
                let (pa, pb) = (at(ja), at(jb));
                let (ha, hb) = (hs[pa] as f64, hs[pb] as f64);
                if ha < -1e8 || hb < -1e8 {
                    hp = ha.max(hb);
                    if curv {
                        dp = if ha > hb { dd[pa] as f64 } else { dd[pb] as f64 };
                    }
                } else {
                    hp = ha + f * (hb - ha);
                    if curv {
                        dp = dd[pa] as f64 + f * (dd[pb] as f64 - dd[pa] as f64);
                    }
                }
            } else if ja == n2 as isize - 1 && f < 1e-9 {
                hp = hs[at(ja)] as f64;
                if curv {
                    dp = dd[at(ja)] as f64;
                }
            } else if ja == -1 && jb == 0 {
                hp = hs[prev] as f64;
                if curv {
                    dp = dd[prev] as f64;
                }
            } else if ja == n2 as isize - 1 {
                hp = hs[at(ja)] as f64;
                if curv {
                    dp = dd[at(ja)] as f64;
                }
            } else {
                hp = NEG;
            }
            let mut cand = hp - drop;
            let mut dn = 0.0;
            if curv {
                dn = dp + step;
                cand -= (dn * dn - dp * dp) / (2.0 * R_EFF);
            }
            if let Some(fr) = opts.hfar {
                if fr[idx] as f64 > cand {
                    cand = fr[idx] as f64;
                    dn = 0.0;
                }
            }
            if opts.arrive {
                arrive[idx] = cand as f32;
            }
            if cand > zi {
                shade[idx] = 1;
                hs[idx] = cand as f32;
                if curv {
                    dd[idx] = dn as f32;
                }
            } else {
                shade[idx] = 0;
                hs[idx] = zi as f32;
                if curv {
                    dd[idx] = 0.0;
                }
            }
        }
    }
    if opts.arrive {
        (shade, arrive)
    } else {
        (shade, hs)
    }
}

/// A two-level domain for shadows: `near` (fine, the area of interest plus a margin) and
/// an optional `far` (coarse, wide). Far cells inside the near domain do not cast (the
/// near grid casts there at full detail); the far sweep's arriving surface, with earth
/// curvature, is interpolated onto the near grid as its incoming shadow surface.
pub struct ShadeDomain<'a> {
    near: Grid<'a>,
    far: Option<Grid<'a>>,
    zcast: Vec<f32>,
    cx: Vec<usize>,
    cf: Vec<f64>,
    rx: Vec<usize>,
    rf: Vec<f64>,
}

impl<'a> ShadeDomain<'a> {
    pub fn new(near: Grid<'a>, far: Option<Grid<'a>>) -> ShadeDomain<'a> {
        let mut d = ShadeDomain {
            near,
            far,
            zcast: Vec::new(),
            cx: Vec::new(),
            cf: Vec::new(),
            rx: Vec::new(),
            rf: Vec::new(),
        };
        if let Some(fg) = far {
            let mut zc: Vec<f32> = fg.z.to_vec();
            let nx1 = near.x0 + near.w as f64 * near.cell;
            let ny1 = near.ytop - near.h as f64 * near.cell;
            for r in 0..fg.h {
                let nn = fg.ytop - (r as f64 + 0.5) * fg.cell;
                if nn > near.ytop || nn < ny1 {
                    continue;
                }
                for c in 0..fg.w {
                    let e = fg.x0 + (c as f64 + 0.5) * fg.cell;
                    if e > near.x0 && e < nx1 {
                        zc[r * fg.w + c] = NEG as f32;
                    }
                }
            }
            d.zcast = zc;
            for i in 0..near.w {
                let u = ((near.x0 + (i as f64 + 0.5) * near.cell - fg.x0) / fg.cell - 0.5).clamp(0.0, fg.w as f64 - 1.000001);
                d.cx.push(u.floor() as usize);
                d.cf.push(u - u.floor());
            }
            for j in 0..near.h {
                let v = ((fg.ytop - (near.ytop - (j as f64 + 0.5) * near.cell)) / fg.cell - 0.5).clamp(0.0, fg.h as f64 - 1.000001);
                d.rx.push(v.floor() as usize);
                d.rf.push(v - v.floor());
            }
        }
        d
    }

    /// Shadow mask on the near grid (1 = shadow) for the sun at true azimuth `az` and
    /// elevation `elev` (degrees); `theta` is the grid azimuth of true north (meridian
    /// convergence), added to `az`.
    pub fn shade_at(&self, az: f64, elev: f64, theta: f64) -> Vec<u8> {
        let near = &self.near;
        let azg = az + theta;
        let mut hfar: Option<Vec<f32>> = None;
        if let Some(fg) = &self.far {
            let (_, fs) = sweep(
                &self.zcast,
                fg.w,
                fg.h,
                fg.cell,
                azg,
                elev,
                &SweepOptions {
                    hfar: None,
                    curvature: true,
                    arrive: true,
                },
            );
            let ww = fg.w;
            let mut hf = vec![0f32; near.w * near.h];
            for j in 0..near.h {
                let (r0, fy) = (self.rx[j], self.rf[j]);
                for i in 0..near.w {
                    let (c0, fx) = (self.cx[i], self.cf[i]);
                    let p = r0 * ww + c0;
                    let (a, b, cc, d) = (fs[p] as f64, fs[p + 1] as f64, fs[p + ww] as f64, fs[p + ww + 1] as f64);
                    hf[j * near.w + i] = if a < -1e8 || b < -1e8 || cc < -1e8 || d < -1e8 {
                        a.max(b).max(cc).max(d) as f32
                    } else {
                        ((a * (1.0 - fx) + b * fx) * (1.0 - fy) + (cc * (1.0 - fx) + d * fx) * fy) as f32
                    };
                }
            }
            hfar = Some(hf);
        }
        let opts = SweepOptions {
            hfar: hfar.as_deref(),
            curvature: false,
            arrive: false,
        };
        sweep(near.z, near.w, near.h, near.cell, azg, elev, &opts).0
    }
}
