# Attractor Search

Find strange attractors that don't exist in the
literature yet. Automatically. Fast.

## Core loop

```
parameterize ODE system
    → integrate forward (RK4, microseconds)
    → classify: diverge / fixed / cycle / chaos
    → if chaos: compute Lyapunov spectrum
    → if novel: catalogue it
    → repeat
```

Thousands of evaluations per second on one core.
The whole thing runs on a laptop.

## Why this is interesting

Strange attractors live on a fractal boundary
between "blows up" and "collapses to nothing."
That boundary has measure zero in parameter space.
Most of it has never been explored.

Sprott found ~19 new chaotic systems in the 90s
by brute-forcing quadratic 3D ODEs with small
integer coefficients. Tiny slice. Simple heuristics.
Nobody has done it systematically since.

## Search space

3D polynomial ODE systems, quadratic to start:

```
dx/dt = Σ a_ij * x^i * y^j * z^k
dy/dt = Σ b_ij * x^i * y^j * z^k
dz/dt = Σ c_ij * x^i * y^j * z^k
```

~30-60 free coefficients per system.
Mostly boring. Chaos is rare. Finding it is
the whole problem.

## Search strategy

Random sampling alone won't work — chaos is
too rare. Need smarter approaches:

- **Boundary tracking**: find one chaotic system,
  perturb params, follow the edge of chaos
- **Evolutionary search**: fitness = positive
  max Lyapunov exponent + bounded trajectory
- **Continuation**: follow bifurcation curves
  through parameter space
- Probably a combination of all three.

## Classification pipeline

For each candidate system:

1. Integrate from random IC, 10k steps
2. Fast reject: did it diverge? → discard
3. Fast reject: did it converge? → discard
4. Compute max Lyapunov exponent
5. If λ₁ > 0 and trajectory bounded → chaotic
6. Full Lyapunov spectrum (3 exponents)
7. Kaplan-Yorke dimension
8. Render the attractor

## Novelty detection

Two systems can look different but be the same
attractor up to smooth deformation (topologically
conjugate). Need invariants:

- Lyapunov spectrum
- Fractal dimension
- Topological template (linking numbers)
- Correlation dimension

Open question: what's sufficient to confirm
genuinely new vs. reparameterization of known?

## Tech

- **Rust** for the core: integration, Lyapunov
  computation, search, classification
- Tight, zero-allocation inner loops
- SIMD where it helps
- Flat-file or SQLite catalogue
- Small visualization tool (TBD)

## Open questions

- Polynomial order: quadratic first, but cubic
  or mixed-order sparse systems might yield
  qualitatively different attractors
- Dashboard vs batch? Live search visualization
  or just dump results?
- How to render attractors? 3D point clouds,
  projected orbits, or something else?
- Parallelism strategy: trivially parallel per
  candidate, but boundary-tracking is sequential

## Success

Run search. Find attractors not in Sprott's
catalogue or anywhere published. Characterize
them. Render them. Ship the catalogue.
