// Lyapunov exponent computation via tangent-space
// integration + QR (Gram-Schmidt) reorthogonalization.

use crate::ode::{
    self, Coeffs, State,
};

/// Data from Lyapunov computation.
#[derive(Debug, Clone)]
pub struct LyapData {
    pub spectrum: [f64; 3],
    pub ky_dim: f64,
}

/// 3-vector ops (inline, no alloc).
type V3 = [f64; 3];

#[inline(always)]
fn dot(a: &V3, b: &V3) -> f64 {
    a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
}

#[inline(always)]
fn norm(a: &V3) -> f64 {
    dot(a, a).sqrt()
}

#[inline(always)]
fn scale(a: &V3, s: f64) -> V3 {
    [a[0]*s, a[1]*s, a[2]*s]
}

#[inline(always)]
fn sub(a: &V3, b: &V3) -> V3 {
    [a[0]-b[0], a[1]-b[1], a[2]-b[2]]
}

#[inline(always)]
fn add(a: &V3, b: &V3) -> V3 {
    [a[0]+b[0], a[1]+b[1], a[2]+b[2]]
}

/// Advance tangent vector w by one RK4 step.
/// dw/dt = J(s) * w, where J is Jacobian at s.
/// We need J at intermediate RK4 points too,
/// so we co-integrate state and tangent.
#[inline(always)]
fn mat_vec(m: &[[f64; 3]; 3], v: &V3) -> V3 {
    [
        m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2],
        m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2],
        m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2],
    ]
}

/// RK4 co-integration of state + 3 tangent vectors.
/// Returns (new_state, new_tangent_vecs).
#[inline]
fn rk4_tangent(
    c: &Coeffs,
    s: &State,
    w: &[V3; 3],
    dt: f64,
) -> (State, [V3; 3]) {
    // State k1
    let f1 = ode::rhs(c, s);
    let j1 = ode::jacobian(c, s);
    let g1: [V3; 3] = [
        mat_vec(&j1, &w[0]),
        mat_vec(&j1, &w[1]),
        mat_vec(&j1, &w[2]),
    ];

    // Midpoint 1
    let s2 = [
        s[0] + 0.5*dt*f1[0],
        s[1] + 0.5*dt*f1[1],
        s[2] + 0.5*dt*f1[2],
    ];
    let w2 = [
        add(&w[0], &scale(&g1[0], 0.5*dt)),
        add(&w[1], &scale(&g1[1], 0.5*dt)),
        add(&w[2], &scale(&g1[2], 0.5*dt)),
    ];
    let f2 = ode::rhs(c, &s2);
    let j2 = ode::jacobian(c, &s2);
    let g2: [V3; 3] = [
        mat_vec(&j2, &w2[0]),
        mat_vec(&j2, &w2[1]),
        mat_vec(&j2, &w2[2]),
    ];

    // Midpoint 2
    let s3 = [
        s[0] + 0.5*dt*f2[0],
        s[1] + 0.5*dt*f2[1],
        s[2] + 0.5*dt*f2[2],
    ];
    let w3 = [
        add(&w[0], &scale(&g2[0], 0.5*dt)),
        add(&w[1], &scale(&g2[1], 0.5*dt)),
        add(&w[2], &scale(&g2[2], 0.5*dt)),
    ];
    let f3 = ode::rhs(c, &s3);
    let j3 = ode::jacobian(c, &s3);
    let g3: [V3; 3] = [
        mat_vec(&j3, &w3[0]),
        mat_vec(&j3, &w3[1]),
        mat_vec(&j3, &w3[2]),
    ];

    // Endpoint
    let s4 = [
        s[0] + dt*f3[0],
        s[1] + dt*f3[1],
        s[2] + dt*f3[2],
    ];
    let w4 = [
        add(&w[0], &scale(&g3[0], dt)),
        add(&w[1], &scale(&g3[1], dt)),
        add(&w[2], &scale(&g3[2], dt)),
    ];
    let f4 = ode::rhs(c, &s4);
    let j4 = ode::jacobian(c, &s4);
    let g4: [V3; 3] = [
        mat_vec(&j4, &w4[0]),
        mat_vec(&j4, &w4[1]),
        mat_vec(&j4, &w4[2]),
    ];

    // Combine
    let d6 = dt / 6.0;
    let sn = [
        s[0] + d6*(f1[0]+2.0*f2[0]+2.0*f3[0]+f4[0]),
        s[1] + d6*(f1[1]+2.0*f2[1]+2.0*f3[1]+f4[1]),
        s[2] + d6*(f1[2]+2.0*f2[2]+2.0*f3[2]+f4[2]),
    ];
    let mut wn = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            wn[i][j] = w[i][j] + d6 * (
                g1[i][j]
                + 2.0*g2[i][j]
                + 2.0*g3[i][j]
                + g4[i][j]
            );
        }
    }
    (sn, wn)
}

/// Gram-Schmidt orthonormalization of 3 vectors.
/// Returns norms before normalization (for exponents).
#[inline]
fn gram_schmidt(w: &mut [V3; 3]) -> [f64; 3] {
    // v0
    let n0 = norm(&w[0]);
    if n0 > 0.0 {
        w[0] = scale(&w[0], 1.0 / n0);
    }
    // v1 -= (v1·v0)*v0
    let d10 = dot(&w[1], &w[0]);
    w[1] = sub(&w[1], &scale(&w[0], d10));
    let n1 = norm(&w[1]);
    if n1 > 0.0 {
        w[1] = scale(&w[1], 1.0 / n1);
    }
    // v2 -= (v2·v0)*v0 + (v2·v1)*v1
    let d20 = dot(&w[2], &w[0]);
    let d21 = dot(&w[2], &w[1]);
    w[2] = sub(
        &sub(&w[2], &scale(&w[0], d20)),
        &scale(&w[1], d21),
    );
    let n2 = norm(&w[2]);
    if n2 > 0.0 {
        w[2] = scale(&w[2], 1.0 / n2);
    }
    [n0, n1, n2]
}

/// Compute maximal Lyapunov exponent only (fast).
pub fn max_lyapunov(
    c: &Coeffs,
    s0: &State,
    dt: f64,
    transient: usize,
    steps: usize,
) -> (f64, bool) {
    let mut s = *s0;
    // Transient
    for _ in 0..transient {
        s = ode::rk4_step(c, &s, dt);
        if s[0]*s[0]+s[1]*s[1]+s[2]*s[2] > 1e8 {
            return (f64::INFINITY, false);
        }
    }
    let mut w: V3 = [1.0, 0.0, 0.0];
    let mut sum = 0.0;
    let bound = 1e6;
    for _ in 0..steps {
        // Advance state
        let sn = ode::rk4_step(c, &s, dt);
        // Advance tangent: dw/dt = J*w
        let j = ode::jacobian(c, &s);
        // Simple Euler for tangent (good enough for
        // max exponent with frequent renorm)
        let dw = mat_vec(&j, &w);
        w = [
            w[0] + dt * dw[0],
            w[1] + dt * dw[1],
            w[2] + dt * dw[2],
        ];
        s = sn;
        let r2 = s[0]*s[0]+s[1]*s[1]+s[2]*s[2];
        if r2 > bound {
            return (f64::INFINITY, false);
        }
        // Renormalize every step
        let n = norm(&w);
        if n > 0.0 && n.is_finite() {
            sum += n.ln();
            w = scale(&w, 1.0 / n);
        } else {
            return (f64::NAN, false);
        }
    }
    let lyap = sum / (steps as f64 * dt);
    (lyap, true)
}

/// Full Lyapunov spectrum (3 exponents) via QR.
pub fn full_spectrum(
    c: &Coeffs,
    s0: &State,
    dt: f64,
    transient: usize,
    steps: usize,
    renorm_interval: usize,
) -> Option<LyapData> {
    let mut s = *s0;
    for _ in 0..transient {
        s = ode::rk4_step(c, &s, dt);
        let r2 = s[0]*s[0]+s[1]*s[1]+s[2]*s[2];
        if r2 > 1e8 {
            return None;
        }
    }
    let mut w: [V3; 3] = [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let mut sums = [0.0f64; 3];
    let mut count = 0u64;
    let bound = 1e6;
    for step in 0..steps {
        let (sn, wn) = rk4_tangent(c, &s, &w, dt);
        s = sn;
        w = wn;
        let r2 = s[0]*s[0]+s[1]*s[1]+s[2]*s[2];
        if r2 > bound || !r2.is_finite() {
            return None;
        }
        if (step + 1) % renorm_interval == 0 {
            let norms = gram_schmidt(&mut w);
            for i in 0..3 {
                if norms[i] > 0.0
                    && norms[i].is_finite()
                {
                    sums[i] += norms[i].ln();
                } else {
                    return None;
                }
            }
            count += 1;
        }
    }
    if count == 0 {
        return None;
    }
    let t = count as f64
        * renorm_interval as f64
        * dt;
    let spectrum = [
        sums[0] / t,
        sums[1] / t,
        sums[2] / t,
    ];
    let ky = kaplan_yorke(&spectrum);
    Some(LyapData { spectrum, ky_dim: ky })
}

/// Kaplan-Yorke dimension from sorted spectrum.
fn kaplan_yorke(spec: &[f64; 3]) -> f64 {
    let mut sorted = *spec;
    sorted.sort_by(|a, b| {
        b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal)
    });
    let mut sum = 0.0;
    for (j, &l) in sorted.iter().enumerate() {
        sum += l;
        if sum < 0.0 {
            if l.abs() < 1e-12 {
                return j as f64;
            }
            return j as f64
                + (sum - l) / l.abs();
        }
    }
    3.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ode::lorenz_coeffs;

    #[test]
    fn lorenz_max_lyapunov() {
        let c = lorenz_coeffs(10.0, 28.0, 8.0/3.0);
        let s0: State = [1.0, 1.0, 1.0];
        let (l, ok) = max_lyapunov(
            &c, &s0, 0.005, 2000, 20000,
        );
        assert!(ok, "integration failed");
        // Lorenz max LE ~ 0.9, accept 0.5..1.5
        assert!(
            l > 0.5 && l < 1.5,
            "Lorenz LE={l}, expected ~0.9"
        );
    }

    #[test]
    fn lorenz_full_spectrum() {
        let c = lorenz_coeffs(10.0, 28.0, 8.0/3.0);
        let s0: State = [1.0, 1.0, 1.0];
        let data = full_spectrum(
            &c, &s0, 0.005, 2000, 40000, 10,
        );
        let d = data.expect("spectrum failed");
        // λ1 > 0, λ2 ≈ 0, λ3 < 0
        assert!(
            d.spectrum[0] > 0.3,
            "λ1={} too small", d.spectrum[0]
        );
        // KY dimension ~ 2.06
        assert!(
            d.ky_dim > 1.5 && d.ky_dim < 2.8,
            "KY dim={}", d.ky_dim
        );
    }
}
