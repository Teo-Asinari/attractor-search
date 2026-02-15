// 3D quadratic ODE system + RK4 integrator.
// 10 terms per equation: 1,x,y,z,x²,y²,z²,xy,xz,yz
// 30 coefficients total. Zero allocations in hot path.

pub const NTERMS: usize = 10;
pub const NCOEFFS: usize = 30;

/// Coefficients for 3 equations, 10 terms each.
/// Layout: [eq0_term0..eq0_term9,
///          eq1_term0..eq1_term9,
///          eq2_term0..eq2_term9]
pub type Coeffs = [f64; NCOEFFS];
pub type State = [f64; 3];

/// Evaluate the 10 basis monomials at (x,y,z).
#[inline(always)]
pub fn basis(s: &State) -> [f64; NTERMS] {
    let (x, y, z) = (s[0], s[1], s[2]);
    [
        1.0,
        x,
        y,
        z,
        x * x,
        y * y,
        z * z,
        x * y,
        x * z,
        y * z,
    ]
}

/// Evaluate one equation: dot product of coeffs[off..off+10]
/// with the basis monomials.
#[inline(always)]
fn eval_eq(
    c: &Coeffs,
    off: usize,
    b: &[f64; NTERMS],
) -> f64 {
    let mut v = 0.0;
    let mut i = 0;
    while i < NTERMS {
        v += c[off + i] * b[i];
        i += 1;
    }
    v
}

/// Evaluate full RHS: ds/dt = f(s).
#[inline(always)]
pub fn rhs(c: &Coeffs, s: &State) -> State {
    let b = basis(s);
    [
        eval_eq(c, 0, &b),
        eval_eq(c, NTERMS, &b),
        eval_eq(c, 2 * NTERMS, &b),
    ]
}

/// Compute Jacobian df_i/dx_j at state s.
/// Returns 3x3 matrix as [[f64;3];3].
/// Row i = partials of equation i w.r.t. x,y,z.
#[inline(always)]
pub fn jacobian(
    c: &Coeffs,
    s: &State,
) -> [[f64; 3]; 3] {
    let (x, y, z) = (s[0], s[1], s[2]);
    let mut j = [[0.0f64; 3]; 3];
    for eq in 0..3 {
        let o = eq * NTERMS;
        // d/dx: terms x(1), x²(2x), xy(y), xz(z)
        j[eq][0] = c[o + 1]
            + 2.0 * c[o + 4] * x
            + c[o + 7] * y
            + c[o + 8] * z;
        // d/dy: terms y(1), y²(2y), xy(x), yz(z)
        j[eq][1] = c[o + 2]
            + 2.0 * c[o + 5] * y
            + c[o + 7] * x
            + c[o + 9] * z;
        // d/dz: terms z(1), z²(2z), xz(x), yz(y)
        j[eq][2] = c[o + 3]
            + 2.0 * c[o + 6] * z
            + c[o + 8] * x
            + c[o + 9] * y;
    }
    j
}

/// Single RK4 step. No allocation.
#[inline(always)]
pub fn rk4_step(
    c: &Coeffs,
    s: &State,
    dt: f64,
) -> State {
    let k1 = rhs(c, s);
    let s2 = [
        s[0] + 0.5 * dt * k1[0],
        s[1] + 0.5 * dt * k1[1],
        s[2] + 0.5 * dt * k1[2],
    ];
    let k2 = rhs(c, &s2);
    let s3 = [
        s[0] + 0.5 * dt * k2[0],
        s[1] + 0.5 * dt * k2[1],
        s[2] + 0.5 * dt * k2[2],
    ];
    let k3 = rhs(c, &s3);
    let s4 = [
        s[0] + dt * k3[0],
        s[1] + dt * k3[1],
        s[2] + dt * k3[2],
    ];
    let k4 = rhs(c, &s4);
    [
        s[0] + dt / 6.0
            * (k1[0] + 2.0*k2[0] + 2.0*k3[0] + k4[0]),
        s[1] + dt / 6.0
            * (k1[1] + 2.0*k2[1] + 2.0*k3[1] + k4[1]),
        s[2] + dt / 6.0
            * (k1[2] + 2.0*k2[2] + 2.0*k3[2] + k4[2]),
    ]
}

/// Integrate for n steps, return final state.
pub fn integrate(
    c: &Coeffs,
    s0: &State,
    dt: f64,
    n: usize,
) -> State {
    let mut s = *s0;
    for _ in 0..n {
        s = rk4_step(c, &s, dt);
    }
    s
}

/// Integrate and collect trajectory.
pub fn integrate_traj(
    c: &Coeffs,
    s0: &State,
    dt: f64,
    n: usize,
) -> Vec<State> {
    let mut traj = Vec::with_capacity(n + 1);
    let mut s = *s0;
    traj.push(s);
    for _ in 0..n {
        s = rk4_step(c, &s, dt);
        traj.push(s);
    }
    traj
}

/// Build Lorenz system coefficients.
/// dx/dt = sigma*(y - x)
/// dy/dt = x*(rho - z) - y
/// dz/dt = x*y - beta*z
pub fn lorenz_coeffs(
    sigma: f64,
    rho: f64,
    beta: f64,
) -> Coeffs {
    let mut c = [0.0; NCOEFFS];
    // eq0: sigma*y - sigma*x
    c[1] = -sigma; // x
    c[2] = sigma;  // y
    // eq1: rho*x - y - x*z
    c[NTERMS + 1] = rho;    // x
    c[NTERMS + 2] = -1.0;   // y
    c[NTERMS + 8] = -1.0;   // xz
    // eq2: x*y - beta*z
    c[2 * NTERMS + 3] = -beta; // z
    c[2 * NTERMS + 7] = 1.0;   // xy
    c
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rk4_harmonic_oscillator() {
        // dx/dt = y, dy/dt = -x (2D in 3D space)
        // z unused. Energy = x² + y² conserved.
        let mut c = [0.0; NCOEFFS];
        c[2] = 1.0;           // eq0: y
        c[NTERMS + 1] = -1.0; // eq1: -x
        let s0: State = [1.0, 0.0, 0.0];
        let dt = 0.001;
        let n = 10000; // 10 units of time
        let sf = integrate(&c, &s0, dt, n);
        let e0 = s0[0] * s0[0] + s0[1] * s0[1];
        let ef = sf[0] * sf[0] + sf[1] * sf[1];
        let err = (ef - e0).abs() / e0;
        assert!(
            err < 1e-6,
            "energy drift {err} too large"
        );
    }

    #[test]
    fn lorenz_bounded() {
        let c = lorenz_coeffs(10.0, 28.0, 8.0 / 3.0);
        let s0: State = [1.0, 1.0, 1.0];
        let sf = integrate(&c, &s0, 0.005, 10000);
        let r = (sf[0]*sf[0]
            + sf[1]*sf[1]
            + sf[2]*sf[2]).sqrt();
        assert!(r < 100.0, "Lorenz diverged: r={r}");
    }

    #[test]
    fn jacobian_check() {
        let c = lorenz_coeffs(10.0, 28.0, 8.0 / 3.0);
        let s: State = [1.0, 2.0, 3.0];
        let j = jacobian(&c, &s);
        // df0/dx = -sigma = -10
        assert!((j[0][0] - (-10.0)).abs() < 1e-10);
        // df0/dy = sigma = 10
        assert!((j[0][1] - 10.0).abs() < 1e-10);
        // df1/dx = rho - z = 28 - 3 = 25
        assert!((j[1][0] - 25.0).abs() < 1e-10);
        // df2/dz = -beta = -8/3
        let expect = -8.0 / 3.0;
        assert!((j[2][2] - expect).abs() < 1e-10);
    }
}
