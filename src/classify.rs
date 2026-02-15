// Classify dynamics of a 3D quadratic ODE system.

use crate::lyapunov::{self, LyapData};
use crate::ode::{self, Coeffs, State};

#[derive(Debug, Clone)]
pub enum Dynamics {
    Divergent,
    FixedPoint,
    Cycle,
    Chaotic(LyapData),
}

/// Thresholds.
const DIV_THRESH: f64 = 1e6;
const FP_VAR: f64 = 1e-4;
const CHAOS_THRESH: f64 = 0.01;
const DT: f64 = 0.005;
const TRANSIENT: usize = 1000;
const CLASSIFY_STEPS: usize = 5000;
const LYAP_STEPS: usize = 30000;
const RENORM: usize = 10;

/// Classify a system from its coefficients.
pub fn classify(c: &Coeffs) -> Dynamics {
    let s0: State = [0.1, 0.1, 0.1];
    let mut s = s0;

    // Transient integration.
    for _ in 0..TRANSIENT {
        s = ode::rk4_step(c, &s, DT);
        let r2 = s[0]*s[0]+s[1]*s[1]+s[2]*s[2];
        if r2 > DIV_THRESH || !r2.is_finite() {
            return Dynamics::Divergent;
        }
    }

    // Collect trajectory stats.
    let mut mean = [0.0f64; 3];
    let mut var = [0.0f64; 3];
    let n = CLASSIFY_STEPS;
    let mut traj_s = s;
    for _ in 0..n {
        traj_s = ode::rk4_step(c, &traj_s, DT);
        let r2 = traj_s[0]*traj_s[0]
            + traj_s[1]*traj_s[1]
            + traj_s[2]*traj_s[2];
        if r2 > DIV_THRESH || !r2.is_finite() {
            return Dynamics::Divergent;
        }
        for i in 0..3 {
            mean[i] += traj_s[i];
        }
    }
    for i in 0..3 {
        mean[i] /= n as f64;
    }

    // Variance pass.
    traj_s = s;
    for _ in 0..n {
        traj_s = ode::rk4_step(c, &traj_s, DT);
        for i in 0..3 {
            let d = traj_s[i] - mean[i];
            var[i] += d * d;
        }
    }
    let total_var = (var[0] + var[1] + var[2])
        / n as f64;
    if total_var < FP_VAR {
        return Dynamics::FixedPoint;
    }

    // Lyapunov spectrum.
    match lyapunov::full_spectrum(
        c, &s0, DT, TRANSIENT, LYAP_STEPS, RENORM,
    ) {
        None => Dynamics::Divergent,
        Some(data) => {
            if data.spectrum[0] > CHAOS_THRESH {
                Dynamics::Chaotic(data)
            } else {
                Dynamics::Cycle
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ode::lorenz_coeffs;

    #[test]
    fn lorenz_is_chaotic() {
        let c = lorenz_coeffs(10.0, 28.0, 8.0/3.0);
        match classify(&c) {
            Dynamics::Chaotic(d) => {
                assert!(d.spectrum[0] > 0.3);
            }
            other => {
                panic!("expected chaotic, got {other:?}");
            }
        }
    }

    #[test]
    fn damped_is_fixed_point() {
        // dx/dt = -x, dy/dt = -y, dz/dt = -z
        let mut c = [0.0; 30];
        c[1] = -1.0;    // eq0: -x
        c[12] = -1.0;   // eq1: -y
        c[23] = -1.0;   // eq2: -z
        match classify(&c) {
            Dynamics::FixedPoint => {}
            other => {
                panic!(
                    "expected fixed point, got {other:?}"
                );
            }
        }
    }
}
