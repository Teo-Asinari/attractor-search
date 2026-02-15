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

# random search — sample 10k systems
./target/release/attractor-search random --count 10000

# evolutionary — breed for chaos
./target/release/attractor-search evolve \
  --generations 500 --pop 200
```

Discovered attractors saved as JSON in `results/`.

## How it works

Each candidate system has 30 coefficients defining
a 3D quadratic ODE. The pipeline:

1. RK4 integrate forward
2. Early reject: divergent or collapsed
3. Compute full Lyapunov spectrum
4. Positive λ₁ + bounded = chaotic → catalogue it

Random search hits chaos ~0.5% of the time.
Evolutionary search converges faster by breeding
for positive Lyapunov exponents.

## First results

10k random evaluations found 48 chaotic systems.
Some with Kaplan-Yorke dimension >2 and Lyapunov
exponents rivaling the Lorenz attractor.
