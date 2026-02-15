# attractor-search

Autonomous discovery of strange attractors that
don't exist in the literature. Searches the space
of 3D quadratic polynomial ODE systems, classifies
dynamics via Lyapunov exponents, catalogues novel
chaotic objects.

## Why

Strange attractors live on a fractal boundary
between divergence and collapse. That boundary
has measure zero. Most of it has never been
explored. Sprott found ~19 in the 90s by brute
force. This is the systematic version.

## Run

```bash
cargo build --release

# random search — sample 50k systems
./target/release/attractor-search random \
  --count 50000

# evolutionary — breed for chaos
./target/release/attractor-search evolve \
  --generations 200 --pop 100
```

Discovered attractors saved as JSON in `results/`.

Visualize the top finds:

```bash
python visualize.py
# gallery at gallery/index.html
```

## How it works

Each candidate system has 30 coefficients defining
a 3D quadratic ODE. The pipeline:

1. RK4 integrate forward
2. Early reject: divergent or collapsed
3. Compute full Lyapunov spectrum
4. Positive λ₁ + bounded = chaotic → catalogue it

Random search hits chaos ~0.5% of the time.
Evolutionary search hits ~29% by breeding for
positive Lyapunov exponents.

## Results so far

**2,045 chaotic systems** catalogued from ~60k
evaluations across two search methods:

- **Random**: 246 from 50k evals (0.49%)
- **Evolutionary**: 1,799 from 10k evals (17.9%)

Top finds have Kaplan-Yorke dimension >2.3 and
Lyapunov exponents in the 0.1–0.5 range (moderate
chaos with visible geometric structure).

Each entry is tagged with its search method.
Gallery renders the top 15 from each, ranked by
fractal dimension + moderate λ₁.

## Novelty

These systems use arbitrary real-valued
coefficients in [-2, 2] — a vastly larger space
than Sprott's integer-coefficient search. The
specific coefficient combinations found are almost
certainly not in any published catalogue.

Topological novelty (genuinely new attractor
classes vs. deformations of known types) is an
open question requiring comparison of topological
invariants against known systems.

## Prior art

- Sprott (1994) — found 19 simple chaotic systems
  with quadratic nonlinearities
- Sprott (2000) — systematic search of 3D quadratic
  systems with fewer terms
- Gilmore & Letellier — topological characterization
  of attractors via template theory
