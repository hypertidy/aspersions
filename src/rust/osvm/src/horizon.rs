//! The kernel: a horizon per sector, raised triangle by triangle, front to back.

/// Earth radius with standard refraction (k = 0.13), metres.
pub const R_EFF: f64 = 6371000.0 / (1.0 - 0.13);
const NEGV: f64 = -1e30;
// degrees <-> radians exactly as js/osvm.js does it (x * RAD, x / RAD), for parity
const RAD: f64 = std::f64::consts::PI / 180.0;
const TINY: f64 = 1e-13;
const NBIN: usize = 16;
const EPS_FORWARD: f64 = 1e-3;

/// Options for [`Horizon::new`].
#[derive(Clone, Debug)]
pub struct HorizonOptions {
    /// Sectors per 90 degree wedge.
    pub sectors: usize,
    /// Lower heights by d^2 / 2R before projection.
    pub curvature: bool,
    /// Keep the visible polygons.
    pub keep_pieces: bool,
    /// Plan position of the eye, added back to piece coordinates.
    pub e: f64,
    pub n: f64,
}
impl Default for HorizonOptions {
    fn default() -> Self {
        HorizonOptions {
            sectors: 64,
            curvature: true,
            keep_pieces: true,
            e: 0.0,
            n: 0.0,
        }
    }
}

/// Visible pieces: convex polygons in plan (absolute coordinates).
/// Polygon `p` has vertices `xy[2 * off[p] ..  2 * off[p + 1]]` as (x, y) pairs, came
/// from triangle `tri[p]` and carries `tag[p]` (ring or source).
#[derive(Clone, Debug, Default)]
pub struct Pieces {
    pub xy: Vec<f64>,
    pub off: Vec<u32>,
    pub tri: Vec<u32>,
    pub tag: Vec<u8>,
}
impl Pieces {
    pub fn len(&self) -> usize {
        self.tri.len()
    }
    pub fn is_empty(&self) -> bool {
        self.tri.is_empty()
    }
    /// Vertices of polygon p as (x, y) pairs.
    pub fn polygon(&self, p: usize) -> &[f64] {
        &self.xy[2 * self.off[p] as usize..2 * self.off[p + 1] as usize]
    }
}

/// What a viewshed returns.
#[derive(Clone, Debug)]
pub struct ViewshedResult {
    pub e: f64,
    pub n: f64,
    /// Height of the eye.
    pub z_eye: f64,
    /// Visible area, square metres (plan).
    pub area: f64,
    /// Triangles seen, triangles with a visible part, triangles rejected by the quick test.
    pub n_tri: usize,
    pub n_vis_tri: usize,
    pub n_quick_hidden: usize,
    /// Breakpoints in all the sector horizons at the end.
    pub envelope_points: usize,
    pub pieces: Pieces,
    /// Elevation angle of the horizon (degrees) at whole azimuths 0..359 clockwise from
    /// plan north; -90 where nothing was seen.
    pub skyline: Vec<f64>,
}

struct Sector {
    fx: f64,
    fy: f64,
    rx: f64,
    ry: f64,
    u0: f64,
    u1: f64,
    eu: Vec<f64>,
    ev: Vec<f64>,
    bmin: [f64; NBIN],
}

fn lower_bound(a: &[f64], x: f64) -> usize {
    a.partition_point(|&v| v < x)
}
fn upper_bound(a: &[f64], x: f64) -> usize {
    a.partition_point(|&v| v <= x)
}
fn seg_val(eu: &[f64], ev: &[f64], j: usize, x: f64) -> f64 {
    let (a, b) = (ev[j], ev[j + 1]);
    if a <= NEGV || b <= NEGV {
        return NEGV;
    }
    a + (b - a) * (x - eu[j]) / (eu[j + 1] - eu[j])
}

impl Sector {
    // envelope value just right of x (or at the right end, just left)
    fn eval_right(&self, x: f64) -> f64 {
        let n = self.eu.len();
        let j = upper_bound(&self.eu, x) as isize - 1;
        if j >= n as isize - 1 {
            return self.ev[n - 1];
        }
        seg_val(&self.eu, &self.ev, j.max(0) as usize, x)
    }
    fn bin_range(&self, umin: f64, umax: f64) -> (usize, usize) {
        let w = (self.u1 - self.u0) / NBIN as f64;
        let b0 = ((umin - self.u0) / w).floor().max(0.0) as usize;
        let b1 = (((umax - self.u0) / w).floor() as isize).clamp(0, NBIN as isize - 1) as usize;
        (b0.min(NBIN - 1), b1)
    }
    // conservative lower bound of the envelope in NBIN equal bins of u
    fn update_bins(&mut self, umin: f64, umax: f64) {
        let w = (self.u1 - self.u0) / NBIN as f64;
        let (b0, b1) = self.bin_range(umin, umax);
        let (eu, ev) = (&self.eu, &self.ev);
        for b in b0..=b1 {
            let x0 = self.u0 + b as f64 * w;
            let x1 = x0 + w;
            let mut j = (upper_bound(eu, x0) as isize - 1).max(0) as usize;
            let mut mn = f64::INFINITY;
            if j + 1 < eu.len() {
                mn = seg_val(eu, ev, j, x0);
            }
            j += 1;
            while j < eu.len() && eu[j] <= x1 {
                if ev[j] < mn {
                    mn = ev[j];
                }
                j += 1;
            }
            if j < eu.len() && j > 0 {
                let v1 = seg_val(eu, ev, j - 1, x1);
                if v1 < mn {
                    mn = v1;
                }
            }
            self.bmin[b] = mn;
        }
    }
    fn hidden_quick(&self, umin: f64, umax: f64, vmax: f64) -> bool {
        let (b0, b1) = self.bin_range(umin, umax);
        (b0..=b1).all(|b| vmax < self.bmin[b])
    }
}

// reusable buffers for the per-polygon step
#[derive(Default)]
struct Scratch {
    bp: Vec<f64>,
    nu: Vec<f64>,
    nv: Vec<f64>,
    new_u: Vec<f64>,
    new_v: Vec<f64>,
    su: Vec<f64>,
    sv: Vec<f64>,
    strips: Vec<[f64; 8]>,
}

enum Visible {
    None,
    Full,
    Strips,
}

// Visible part of convex polygon (pu, pv) above the sector's horizon (into scr.strips, or
// Full), then raise the horizon over it.
fn process_poly(sec: &mut Sector, pu: &mut [f64], pv: &[f64], scr: &mut Scratch) -> Visible {
    let n = pu.len();
    let (u0, u1) = (sec.u0, sec.u1);
    let (mut umin, mut umax) = (f64::INFINITY, f64::NEG_INFINITY);
    for u in pu.iter_mut() {
        // clamp, and snap vertices made by clipping on the sector's edges (they land a
        // rounding error inside, which would leave a sliver of no horizon at the edge)
        *u = if *u < u0 + 1e-12 {
            u0
        } else if *u > u1 - 1e-12 {
            u1
        } else {
            *u
        };
        umin = umin.min(*u);
        umax = umax.max(*u);
    }
    if umax - umin < TINY {
        return Visible::None;
    }
    // non-vertical edges (ua, va, ub, vb), ua < ub
    let mut edg = [[0.0f64; 4]; 16];
    let mut ne = 0;
    for i in 0..n {
        let j = if i + 1 == n { 0 } else { i + 1 };
        let (ua, ub) = (pu[i], pu[j]);
        if (ub - ua).abs() < TINY {
            continue;
        }
        edg[ne] = if ua < ub { [ua, pv[i], ub, pv[j]] } else { [ub, pv[j], ua, pv[i]] };
        ne += 1;
    }
    // breakpoints: polygon vertices and envelope vertices inside (umin, umax)
    let lo = lower_bound(&sec.eu, umin);
    let hi = upper_bound(&sec.eu, umax);
    scr.bp.clear();
    scr.bp.extend_from_slice(pu);
    for &x in &sec.eu[lo..hi] {
        if x > umin && x < umax {
            scr.bp.push(x);
        }
    }
    let bp = &mut scr.bp;
    for i in 1..bp.len() {
        let x = bp[i];
        let mut k = i;
        while k > 0 && bp[k - 1] > x {
            bp[k] = bp[k - 1];
            k -= 1;
        }
        bp[k] = x;
    }
    scr.nu.clear();
    scr.nv.clear();
    scr.strips.clear();
    let (eu, ev) = (&sec.eu, &sec.ev);
    let (mut any, mut all_full) = (false, true);
    let mut j_e = lo.saturating_sub(1);
    let mut prev = bp[0];
    for &b in bp.iter().skip(1) {
        let a = prev;
        if b - a < TINY {
            continue;
        }
        prev = b;
        let m = 0.5 * (a + b);
        // the two polygon edges over (a, b)
        let mut got = 0;
        let (mut xa, mut xb, mut ya, mut yb) = (0.0, 0.0, 0.0, 0.0);
        for e in edg.iter().take(ne) {
            if got >= 2 {
                break;
            }
            if e[0] <= m && e[2] >= m {
                let s = (e[3] - e[1]) / (e[2] - e[0]);
                let (va, vb) = (e[1] + s * (a - e[0]), e[1] + s * (b - e[0]));
                if got == 0 {
                    xa = va;
                    xb = vb;
                } else {
                    ya = va;
                    yb = vb;
                }
                got += 1;
            }
        }
        if got < 2 {
            continue;
        }
        let (la, lb, ua, ub) = if xa + xb <= ya + yb { (xa, xb, ya, yb) } else { (ya, yb, xa, xb) };
        // the envelope segment over (a, b)
        while j_e + 2 < eu.len() && eu[j_e + 1] <= m {
            j_e += 1;
        }
        while j_e + 1 < eu.len() && eu[j_e + 1] - eu[j_e] < TINY && eu[j_e + 1] <= m {
            j_e += 1;
        }
        let (ea, eb) = if ev[j_e] <= NEGV || ev[j_e + 1] <= NEGV {
            (NEGV, NEGV)
        } else {
            (seg_val(eu, ev, j_e, a), seg_val(eu, ev, j_e, b))
        };
        // split (a, b) where E crosses L or U
        let mut cut = [0.0f64; 4];
        let mut ncut = 1;
        if ea > NEGV {
            let (mut x_a, mut x_b) = (-1.0, -1.0);
            let (d0, d1) = (ea - la, eb - lb);
            if d0 * d1 < 0.0 {
                x_a = d0 / (d0 - d1);
            }
            let (d0, d1) = (ea - ua, eb - ub);
            if d0 * d1 < 0.0 {
                x_b = d0 / (d0 - d1);
            }
            if x_a > 0.0 && x_b > 0.0 && x_b < x_a {
                std::mem::swap(&mut x_a, &mut x_b);
            }
            if x_a > 0.0 {
                cut[ncut] = x_a;
                ncut += 1;
            }
            if x_b > 0.0 {
                cut[ncut] = x_b;
                ncut += 1;
            }
        }
        cut[ncut] = 1.0;
        ncut += 1;
        for c in 0..ncut - 1 {
            let (t0, t1) = (cut[c], cut[c + 1]);
            if t1 - t0 < 1e-12 {
                continue;
            }
            let (s0, s1, tm) = (a + (b - a) * t0, a + (b - a) * t1, 0.5 * (t0 + t1));
            let (l0, l1) = (la + (lb - la) * t0, la + (lb - la) * t1);
            let (u0v, u1v) = (ua + (ub - ua) * t0, ua + (ub - ua) * t1);
            let lerp_e = |t: f64| if ea <= NEGV { NEGV } else { ea + (eb - ea) * t };
            let (e0, e1, em) = (lerp_e(t0), lerp_e(t1), lerp_e(tm));
            let (lm, um) = (la + (lb - la) * tm, ua + (ub - ua) * tm);
            if em < lm {
                any = true;
                scr.strips.push([s0, l0, s1, l1, s1, u1v, s0, u0v]);
            } else if em < um {
                any = true;
                all_full = false;
                scr.strips.push([s0, e0, s1, e1, s1, u1v, s0, u0v]);
            } else {
                all_full = false;
            }
            scr.nu.push(s0);
            scr.nu.push(s1);
            if em < um {
                scr.nv.push(u0v);
                scr.nv.push(u1v);
            } else {
                scr.nv.push(e0);
                scr.nv.push(e1);
            }
        }
    }
    if !any {
        return Visible::None;
    }
    // splice the new envelope over [umin, umax]
    scr.new_u.clear();
    scr.new_v.clear();
    if lo > 0 {
        scr.new_u.push(umin);
        scr.new_v.push(seg_val(eu, ev, lo - 1, umin));
    }
    scr.new_u.extend_from_slice(&scr.nu);
    scr.new_v.extend_from_slice(&scr.nv);
    if hi < eu.len() {
        scr.new_u.push(umax);
        scr.new_v.push(seg_val(eu, ev, hi - 1, umax));
    }
    // tidy: drop repeated points and interior collinear points
    scr.su.clear();
    scr.sv.clear();
    for i in 0..scr.new_u.len() {
        let (nu_i, nv_i) = (scr.new_u[i], scr.new_v[i]);
        let k = scr.su.len();
        if k > 0 && (scr.su[k - 1] - nu_i).abs() < TINY && (scr.sv[k - 1] - nv_i).abs() <= 1e-12 * (1.0 + nv_i.abs()) {
            continue;
        }
        if k >= 2
            && nu_i - scr.su[k - 1] > TINY
            && scr.su[k - 1] - scr.su[k - 2] > TINY
            && scr.sv[k - 1] > NEGV
            && scr.sv[k - 2] > NEGV
            && nv_i > NEGV
        {
            let pr = scr.sv[k - 2] + (nv_i - scr.sv[k - 2]) * (scr.su[k - 1] - scr.su[k - 2]) / (nu_i - scr.su[k - 2]);
            if (pr - scr.sv[k - 1]).abs() <= 1e-12 * (1.0 + pr.abs()) {
                scr.su[k - 1] = nu_i;
                scr.sv[k - 1] = nv_i;
                continue;
            }
        }
        scr.su.push(nu_i);
        scr.sv.push(nv_i);
    }
    sec.eu.splice(lo..hi, scr.su.iter().copied());
    sec.ev.splice(lo..hi, scr.sv.iter().copied());
    sec.update_bins(umin, umax);
    if all_full {
        Visible::Full
    } else {
        Visible::Strips
    }
}

// Sutherland-Hodgman against A x + B y + C >= 0
#[allow(clippy::too_many_arguments)]
fn clip(ix: &[f64], iy: &[f64], iz: &[f64], n: usize, a: f64, b: f64, c: f64, ox: &mut [f64], oy: &mut [f64], oz: &mut [f64]) -> usize {
    let mut m = 0;
    for i in 0..n {
        let j = if i + 1 == n { 0 } else { i + 1 };
        let di = a * ix[i] + b * iy[i] + c;
        let dj = a * ix[j] + b * iy[j] + c;
        if di >= 0.0 {
            ox[m] = ix[i];
            oy[m] = iy[i];
            oz[m] = iz[i];
            m += 1;
        }
        if (di >= 0.0) != (dj >= 0.0) {
            let t = di / (di - dj);
            ox[m] = ix[i] + t * (ix[j] - ix[i]);
            oy[m] = iy[i] + t * (iy[j] - iy[i]);
            oz[m] = iz[i] + t * (iz[j] - iz[i]);
            m += 1;
        }
    }
    m
}

/// A viewshed under construction, for an eye at plan (0, 0) and height `z_eye`. Feed
/// triangles with [`Horizon::tri`] in a front-to-back order; read the answer with
/// [`Horizon::result`].
pub struct Horizon {
    z_eye: f64,
    k2: f64,
    s: usize,
    sec: Vec<Sector>,
    keep: bool,
    e0: f64,
    n0: f64,
    pieces: Pieces,
    area: f64,
    n_tri: usize,
    n_vis: usize,
    n_hid: usize,
    scr: Scratch,
}

impl Horizon {
    pub fn new(z_eye: f64, opts: HorizonOptions) -> Horizon {
        let s = opts.sectors.max(1);
        let mut sec = Vec::with_capacity(4 * s);
        for wq in 0..4 {
            let a = wq as f64 * 90.0 * RAD;
            let (fx, fy, rx, ry) = (a.sin(), a.cos(), a.cos(), -a.sin());
            for k in 0..s {
                let u0 = -1.0 + 2.0 * k as f64 / s as f64;
                let u1 = -1.0 + 2.0 * (k + 1) as f64 / s as f64;
                sec.push(Sector {
                    fx,
                    fy,
                    rx,
                    ry,
                    u0,
                    u1,
                    eu: vec![u0, u1],
                    ev: vec![NEGV, NEGV],
                    bmin: [NEGV; NBIN],
                });
            }
        }
        Horizon {
            z_eye,
            k2: if opts.curvature { 1.0 / (2.0 * R_EFF) } else { 0.0 },
            s,
            sec,
            keep: opts.keep_pieces,
            e0: opts.e,
            n0: opts.n,
            pieces: Pieces {
                off: vec![0],
                ..Default::default()
            },
            area: 0.0,
            n_tri: 0,
            n_vis: 0,
            n_hid: 0,
            scr: Scratch::default(),
        }
    }

    fn sector_of(&self, az: f64) -> usize {
        // az in [-45, 315)
        let s = self.s;
        let wq = (((az + 45.0) / 90.0).floor() as isize).clamp(0, 3) as usize;
        let u = ((az - wq as f64 * 90.0) * RAD).tan();
        wq * s + ((((u + 1.0) / 2.0 * s as f64).floor() as isize).clamp(0, s as isize - 1) as usize)
    }

    /// Add one triangle. `id` is the caller's triangle id and `tag` a small integer carried
    /// to the pieces. Plan coordinates are relative to the eye; heights are absolute.
    #[allow(clippy::too_many_arguments)]
    pub fn tri(&mut self, id: u32, tag: u8, ax: f64, ay: f64, za: f64, bx: f64, by: f64, zb: f64, cx: f64, cy: f64, zc: f64) {
        self.n_tri += 1;
        let k2 = self.k2;
        let az = za - self.z_eye - (ax * ax + ay * ay) * k2;
        let bz = zb - self.z_eye - (bx * bx + by * by) * k2;
        let cz = zc - self.z_eye - (cx * cx + cy * cy) * k2;
        // plane z = A x + B y + C
        let (ux, uy, uz, vx, vy, vz) = (bx - ax, by - ay, bz - az, cx - ax, cy - ay, cz - az);
        let (nx, ny, nz) = (uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx);
        if nz.abs() < 1e-12 {
            return;
        }
        let (pa, pb) = (-nx / nz, -ny / nz);
        let pc = az - pa * ax - pb * ay;
        // which sectors: the azimuth span of the vertices
        let ns = 4 * self.s;
        let (o1, o2, o3) = (ax * by - ay * bx, bx * cy - by * cx, cx * ay - cy * ax);
        let all = (o1 >= 0.0 && o2 >= 0.0 && o3 >= 0.0) || (o1 <= 0.0 && o2 <= 0.0 && o3 <= 0.0);
        let (s0, s1);
        if all {
            s0 = 0;
            s1 = ns - 1;
        } else {
            let wrap = |x: f64| (x % 360.0 + 360.0) % 360.0;
            let aa = ax.atan2(ay) / RAD;
            let ab = aa + (wrap(bx.atan2(by) / RAD - aa + 540.0) - 180.0);
            let ac = aa + (wrap(cx.atan2(cy) / RAD - aa + 540.0) - 180.0);
            let (mut amin, mut amax) = (aa.min(ab).min(ac), aa.max(ab).max(ac));
            while amin < -45.0 {
                amin += 360.0;
                amax += 360.0;
            }
            while amin >= 315.0 {
                amin -= 360.0;
                amax -= 360.0;
            }
            // widen by a hair: a vertex exactly on a sector boundary must not leave a gap there
            amin -= 1e-9;
            amax += 1e-9;
            if amin < -45.0 {
                amin += 360.0;
                amax += 360.0;
            }
            s0 = self.sector_of(amin);
            s1 = if amax >= 315.0 {
                ns + self.sector_of(amax - 360.0)
            } else {
                self.sector_of(amax)
            };
            if s0 == s1 {
                let sq = &self.sec[s0 % ns];
                let (ta, tb, tc) = (sq.fx * ax + sq.fy * ay, sq.fx * bx + sq.fy * by, sq.fx * cx + sq.fy * cy);
                if ta > EPS_FORWARD && tb > EPS_FORWARD && tc > EPS_FORWARD {
                    let qa = (sq.rx * ax + sq.ry * ay) / ta;
                    let qb = (sq.rx * bx + sq.ry * by) / tb;
                    let qc = (sq.rx * cx + sq.ry * cy) / tc;
                    let vmax = (az / ta).max(bz / tb).max(cz / tc);
                    if sq.hidden_quick(qa.min(qb).min(qc), qa.max(qb).max(qc), vmax) {
                        self.n_hid += 1;
                        return;
                    }
                }
            }
        }
        let mut visible = false;
        let (mut bx0, mut by0, mut bz0) = ([0.0f64; 16], [0.0f64; 16], [0.0f64; 16]);
        let (mut bx1, mut by1, mut bz1) = ([0.0f64; 16], [0.0f64; 16], [0.0f64; 16]);
        let mut pu = [0.0f64; 16];
        let mut pv = [0.0f64; 16];
        for si in s0..=s1 {
            let sc = &self.sec[si % ns];
            let (fx, fy, rx, ry, su0, su1) = (sc.fx, sc.fy, sc.rx, sc.ry, sc.u0, sc.u1);
            bx0[..3].copy_from_slice(&[ax, bx, cx]);
            by0[..3].copy_from_slice(&[ay, by, cy]);
            bz0[..3].copy_from_slice(&[az, bz, cz]);
            let n = clip(&bx0, &by0, &bz0, 3, rx - su0 * fx, ry - su0 * fy, 0.0, &mut bx1, &mut by1, &mut bz1);
            if n < 3 {
                continue;
            }
            let n = clip(&bx1, &by1, &bz1, n, su1 * fx - rx, su1 * fy - ry, 0.0, &mut bx0, &mut by0, &mut bz0);
            if n < 3 {
                continue;
            }
            let n = clip(&bx0, &by0, &bz0, n, fx, fy, -EPS_FORWARD, &mut bx1, &mut by1, &mut bz1);
            if n < 3 {
                continue;
            }
            for k in 0..n {
                let t = fx * bx1[k] + fy * by1[k];
                pu[k] = (rx * bx1[k] + ry * by1[k]) / t;
                pv[k] = bz1[k] / t;
            }
            let vis = process_poly(&mut self.sec[si % ns], &mut pu[..n], &pv[..n], &mut self.scr);
            let plane = (pa, pb, pc, fx, fy, rx, ry);
            match vis {
                Visible::None => {}
                Visible::Full => {
                    let mut poly = [0.0f64; 32];
                    for k in 0..n {
                        poly[2 * k] = pu[k];
                        poly[2 * k + 1] = pv[k];
                    }
                    visible |= self.emit(&poly[..2 * n], plane, id, tag);
                }
                Visible::Strips => {
                    let strips = std::mem::take(&mut self.scr.strips);
                    for st in &strips {
                        visible |= self.emit(st, plane, id, tag);
                    }
                    self.scr.strips = strips;
                }
            }
        }
        if visible {
            self.n_vis += 1;
        }
    }

    // image polygon (u, v pairs) -> plan, recorded if it has area
    fn emit(&mut self, uv: &[f64], plane: (f64, f64, f64, f64, f64, f64, f64), id: u32, tag: u8) -> bool {
        let (pa, pb, pc, fx, fy, rx, ry) = plane;
        let n = uv.len() / 2;
        let mut q = [0.0f64; 32];
        for k in 0..n {
            let (u, v) = (uv[2 * k], uv[2 * k + 1]);
            let (dx, dy) = (fx + u * rx, fy + u * ry);
            let den = v - pa * dx - pb * dy;
            let mut t = if den.abs() < 1e-15 { 0.0 } else { pc / den };
            if t.is_nan() || t <= 0.0 || t > 1e7 {
                t = 0.0;
            }
            q[2 * k] = t * dx;
            q[2 * k + 1] = t * dy;
        }
        let mut a2 = 0.0;
        for k in 0..n {
            let j = if k + 1 == n { 0 } else { k + 1 };
            a2 += q[2 * k] * q[2 * j + 1] - q[2 * j] * q[2 * k + 1];
        }
        let ar = a2.abs() / 2.0;
        if ar < 1e-6 {
            return false;
        }
        self.area += ar;
        if self.keep {
            for k in 0..n {
                self.pieces.xy.push(q[2 * k] + self.e0);
                self.pieces.xy.push(q[2 * k + 1] + self.n0);
            }
            self.pieces.off.push((self.pieces.xy.len() / 2) as u32);
            self.pieces.tri.push(id);
            self.pieces.tag.push(tag);
        }
        true
    }

    /// The horizon so far: elevation angle (degrees) at whole azimuths 0..359 clockwise
    /// from plan north, -90 where nothing has been seen.
    pub fn skyline(&self) -> Vec<f64> {
        (0..360)
            .map(|az| {
                let a2 = if az >= 315 { az as f64 - 360.0 } else { az as f64 };
                let si = self.sector_of(a2);
                let wq = si / self.s;
                let u = ((a2 - wq as f64 * 90.0) * RAD).tan();
                let v = self.sec[si].eval_right(u);
                if v <= NEGV {
                    -90.0
                } else {
                    (v / (1.0 + u * u).sqrt()).atan() / RAD
                }
            })
            .collect()
    }

    /// Finish: area, counts, pieces and skyline.
    pub fn result(self) -> ViewshedResult {
        let skyline = self.skyline();
        ViewshedResult {
            e: self.e0,
            n: self.n0,
            z_eye: self.z_eye,
            area: self.area,
            n_tri: self.n_tri,
            n_vis_tri: self.n_vis,
            n_quick_hidden: self.n_hid,
            envelope_points: self.sec.iter().map(|s| s.eu.len()).sum(),
            pieces: self.pieces,
            skyline,
        }
    }
}
