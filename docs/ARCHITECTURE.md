# Architecture

This document describes the design decisions behind `field-core`.

## Why a Continuous Field?

Most constraint solvers use discrete graphs — nodes for variables, edges for constraints. This works but has fundamental limitations:

1. **No gradients**: You can't differentiate through a SAT solver
2. **Binary satisfaction**: Constraints are satisfied or not — no graceful degradation
3. **Topology blindness**: The graph doesn't tell you about holes or connectedness
4. **Brittle propagation**: Local changes require full re-solve

A continuous constraint field addresses all four:

- **Field values are differentiable** — gradient is a first-class output
- **Confidence is continuous** — constraints degrade gracefully [0, 1]
- **Topology emerges** — Delaunay triangulation gives you Betti numbers for free
- **Propagation is local** — influence falls off as a Gaussian, no global re-solve needed

## Core Types

```
FieldPosition         A single constraint: id, coords, weight, stiffness
FieldQuery            A point in space-time to evaluate: point, time, tolerance
FieldReading          Field value, gradient, confidence, nearby positions
ConstraintField       Container: all positions + time
FieldTopology         Topology state: connectedness, breaches, Betti numbers
Tile                  A handle returned when embedding a new position
FieldStatus           Summary stats for the field
```

### Stiffness (σ)

Stiffness controls the spread of the Gaussian basis. Higher stiffness = wider influence. It maps to the standard deviation of the Gaussian:

```
f_i(x) = w_i · exp(-||x - p_i||² / σ_i²)
```

At distance σ, the influence drops to ~36.8% of max. At 3σ, ~1.1%.

### Time Decay

Confidence decays exponentially with time:

```
decay(t) = exp(-t / 1000)
```

After 1000 time steps, confidence is ~36.8%. This lets the field "cool" so stale constraints fade.

## Topology Detection

### 2D (Delaunay Triangulation)

For 2D fields, we use the `delaunator` crate to compute Delaunay triangulation, then:

1. Build adjacency from triangle edges
2. Find connected components via BFS (β₀)
3. Compute Euler characteristic: β₁ = E - V + C

This gives us the number of holes (β₁) in the constraint arrangement.

### N-Dimensional (Proximity Graph)

For N > 2 dimensions, we fall back to a proximity graph: positions within `2 × (σ₁ + σ₂)` are connected. This is less precise than Delaunay but works in arbitrary dimensions.

## Propagation

Each propagation step:

1. For each position, compute the net "force" from neighbors within 3σ radius
2. Force = weighted attraction toward neighbor, scaled by Gaussian influence
3. Apply force with damping factor (0.1) to avoid oscillation
4. Track displacement and topology changes

This is essentially gradient descent on the field energy, with the positions as free parameters. Over multiple steps, the field converges to a consensus configuration.

## CLI Persistence

The CLI saves field state to `.field-state.json` (configurable via `--state`). This is the full `ConstraintField` serialized as JSON — everything needed to resume a session.

## Trade-offs

| Decision | Rationale | Cost |
|----------|-----------|------|
| Gaussian bases | Smooth, differentiable, local influence | Can't represent hard constraints |
| Delaunay topology | Automatic, no manual labeling | Only works in 2D (falls back to proximity) |
| Force-based propagation | Simple, local, no global solve | May converge slowly |
| JSON serialization | Human-readable, debuggable | Not optimal for speed |
| u64 IDs | No collisions, simple | Memory overhead for millions of positions |
