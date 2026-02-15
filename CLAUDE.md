# CLAUDE.md

This file provides guidance to Claude Code
(claude.ai/code) when working with code in
this repository.

## Style

- All lines under 80 characters. No exceptions.
- Short functions, short names.
- Zero fluff. No verbose comments.
- Zero-allocation inner loops.

## Commands

```bash
# build
~/.cargo/bin/cargo build --release

# test (8 tests: ODE, Lyapunov, classify, catalog)
~/.cargo/bin/cargo test

# single test
~/.cargo/bin/cargo test lorenz_is_chaotic

# run: random search
./target/release/attractor-search random --count 10000

# run: evolutionary search
./target/release/attractor-search evolve \
  --generations 500 --pop 200
```

Results go to `results/` as JSON files.

## Architecture

3D quadratic ODE systems defined by 30 f64
coefficients (10 monomials per equation:
`1, x, y, z, x², y², z², xy, xz, yz`).

Core data types:
- `Coeffs = [f64; 30]` — system definition
- `State = [f64; 3]` — phase space point
- `LyapData` — spectrum `[f64; 3]` + KY dim
- `Dynamics` — enum: Divergent/Fixed/Cycle/Chaotic

Pipeline per candidate system:

```
coefficients
 → RK4 integrate (ode.rs)
 → classify (classify.rs)
   → transient: diverge? discard
   → variance: fixed point? discard
   → Lyapunov spectrum (lyapunov.rs)
   → λ₁ > 0.01 + bounded → chaotic
 → catalogue (catalog.rs)
```

**ode.rs** — RK4 integrator + analytic Jacobian.
Hot path. Stack-only, inlined.

**lyapunov.rs** — Full spectrum via tangent vector
co-integration + Gram-Schmidt QR. Also computes
Kaplan-Yorke dimension.

**classify.rs** — Multi-stage pipeline with early
exits. Constants: DT=0.005, TRANSIENT=1000,
CLASSIFY_STEPS=5000, LYAP_STEPS=30000, RENORM=10.

**search.rs** — Random sampling (uniform [-2,2])
and evolutionary (Gaussian mutation, top-half
selection, fitness = bounded positive λ₁).

**catalog.rs** — JSON per entry in results/.
FNV hash of coefficients for dedup/naming.

## Extending

To add search strategies: new fn in search.rs,
call classify::classify(), save chaotic hits.

To extend ODE order: modify NTERMS, NCOEFFS,
basis(), eval_eq() — evaluation auto-scales.

To tune speed vs accuracy: adjust constants in
classify.rs (TRANSIENT, LYAP_STEPS, RENORM).
