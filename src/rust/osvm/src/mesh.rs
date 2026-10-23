//! An explicit triangulation as a source of triangles, and a front-to-back order for it.

use crate::horizon::Horizon;
use std::collections::HashMap;

/// A triangulation: vertex coordinates (absolute) and three vertex indices per triangle.
#[derive(Clone, Copy, Debug)]
pub struct Mesh<'a> {
    pub x: &'a [f64],
    pub y: &'a [f64],
    pub z: &'a [f64],
    pub tri: &'a [u32],
}
impl Mesh<'_> {
    pub fn n_tri(&self) -> usize {
        self.tri.len() / 3
    }
}

/// A front-to-back order from [`order_topo`].
#[derive(Clone, Debug)]
pub struct TopoOrder {
    pub order: Vec<u32>,
    /// Triangles caught in cycles of the in-front relation, appended nearest first.
    pub cycles: usize,
    /// Number of shared-edge constraints.
    pub constraints: usize,
}

/// A front-to-back order for any triangulation of a convex region, seen from (e, n).
///
/// Along any plan ray from the eye, the triangles crossed are neighbours across shared
/// edges, each crossed from the side facing the eye. So "the triangle on the eye's side of
/// a shared edge goes first", over all shared edges, is the whole in-front relation, and a
/// topological sort of it is a valid order. Delaunay meshes never have a cycle
/// (Edelsbrunner's acyclicity theorem); other meshes can, and any triangles left in a
/// cycle are appended nearest first and counted. The region must be convex: a plan ray
/// must not leave the mesh and come back.
pub fn order_topo(mesh: &Mesh, e: f64, n: f64) -> TopoOrder {
    let (t_ix, xs, ys) = (mesh.tri, mesh.x, mesh.y);
    let nt = t_ix.len() / 3;
    let nv = xs.len() as u64;
    let mut edge: HashMap<u64, u32> = HashMap::with_capacity(nt * 2);
    let (mut from, mut to): (Vec<u32>, Vec<u32>) = (Vec::new(), Vec::new());
    for t in 0..nt {
        for k in 0..3 {
            let a = t_ix[3 * t + k] as usize;
            let b = t_ix[3 * t + (k + 1) % 3] as usize;
            let key = if a < b {
                a as u64 * nv + b as u64
            } else {
                b as u64 * nv + a as u64
            };
            let o = match edge.remove(&key) {
                None => {
                    edge.insert(key, t as u32);
                    continue;
                }
                Some(o) => o,
            };
            let opp = t_ix[3 * t + (k + 2) % 3] as usize;
            let (ex, ey) = (xs[b] - xs[a], ys[b] - ys[a]);
            let se = ex * (n - ys[a]) - ey * (e - xs[a]);
            let st = ex * (ys[opp] - ys[a]) - ey * (xs[opp] - xs[a]);
            if se.abs() <= 1e-9 * (ex.abs() + ey.abs()) * ((e - xs[a]).abs() + (n - ys[a]).abs() + 1.0) {
                continue; // the eye is on this edge's line
            }
            if (se > 0.0) == (st > 0.0) {
                from.push(t as u32);
                to.push(o);
            } else {
                from.push(o);
                to.push(t as u32);
            }
        }
    }
    let mut start = vec![0u32; nt + 1];
    let mut indeg = vec![0u32; nt];
    for i in 0..from.len() {
        start[from[i] as usize + 1] += 1;
        indeg[to[i] as usize] += 1;
    }
    for i in 0..nt {
        start[i + 1] += start[i];
    }
    let mut fill = start[..nt].to_vec();
    let mut adj = vec![0u32; from.len()];
    for i in 0..from.len() {
        let f = from[i] as usize;
        adj[fill[f] as usize] = to[i];
        fill[f] += 1;
    }
    let d2: Vec<f64> = (0..nt)
        .map(|t| {
            let (a, b, c) = (t_ix[3 * t] as usize, t_ix[3 * t + 1] as usize, t_ix[3 * t + 2] as usize);
            let cx = (xs[a] + xs[b] + xs[c]) / 3.0 - e;
            let cy = (ys[a] + ys[b] + ys[c]) / 3.0 - n;
            cx * cx + cy * cy
        })
        .collect();
    let by_dist = |p: &u32, q: &u32| d2[*p as usize].partial_cmp(&d2[*q as usize]).unwrap_or(std::cmp::Ordering::Equal);
    let mut seeds: Vec<u32> = (0..nt as u32).filter(|&t| indeg[t as usize] == 0).collect();
    seeds.sort_by(by_dist);
    let mut done = vec![false; nt];
    let mut order: Vec<u32> = Vec::with_capacity(nt);
    for &s in &seeds {
        done[s as usize] = true;
        order.push(s);
    }
    let mut head = 0;
    while head < order.len() {
        let u = order[head] as usize;
        head += 1;
        for k in start[u]..start[u + 1] {
            let v = adj[k as usize] as usize;
            indeg[v] -= 1;
            if indeg[v] == 0 && !done[v] {
                done[v] = true;
                order.push(v as u32);
            }
        }
    }
    let cycles = nt - order.len();
    if cycles > 0 {
        let mut rest: Vec<u32> = (0..nt as u32).filter(|&t| !done[t as usize]).collect();
        rest.sort_by(by_dist);
        order.extend(rest);
    }
    TopoOrder {
        order,
        cycles,
        constraints: from.len(),
    }
}

/// Feed the mesh's triangles to `h` in `order`, relative to the eye at (e, n).
pub fn mesh_source(h: &mut Horizon, mesh: &Mesh, order: &[u32], e: f64, n: f64, tag: u8) {
    let (t_ix, x, y, z) = (mesh.tri, mesh.x, mesh.y, mesh.z);
    for &t in order {
        let t = t as usize;
        let (a, b, c) = (t_ix[3 * t] as usize, t_ix[3 * t + 1] as usize, t_ix[3 * t + 2] as usize);
        h.tri(
            t as u32,
            tag,
            x[a] - e,
            y[a] - n,
            z[a],
            x[b] - e,
            y[b] - n,
            z[b],
            x[c] - e,
            y[c] - n,
            z[c],
        );
    }
}

/// Height on the mesh at (e, n) (linear scan), or `None` outside it.
pub fn sample_mesh(mesh: &Mesh, e: f64, n: f64) -> Option<f64> {
    let (t_ix, x, y, z) = (mesh.tri, mesh.x, mesh.y, mesh.z);
    for t in 0..t_ix.len() / 3 {
        let (a, b, c) = (t_ix[3 * t] as usize, t_ix[3 * t + 1] as usize, t_ix[3 * t + 2] as usize);
        let d = (y[b] - y[c]) * (x[a] - x[c]) + (x[c] - x[b]) * (y[a] - y[c]);
        let l1 = ((y[b] - y[c]) * (e - x[c]) + (x[c] - x[b]) * (n - y[c])) / d;
        let l2 = ((y[c] - y[a]) * (e - x[c]) + (x[a] - x[c]) * (n - y[c])) / d;
        let l3 = 1.0 - l1 - l2;
        if l1 >= -1e-12 && l2 >= -1e-12 && l3 >= -1e-12 {
            return Some(l1 * z[a] + l2 * z[b] + l3 * z[c]);
        }
    }
    None
}
