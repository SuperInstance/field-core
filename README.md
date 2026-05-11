# field-core 🧭

**Continuous constraint field library for topological constraint satisfaction.**

The constraint field is a paradigm shift from discrete constraint graphs to continuous, differentiable fields. Each constraint is a Gaussian basis function embedded in N-dimensional space. The field naturally captures topology (connectedness, holes) and propagates influences between neighboring constraints.

## Core Idea

```
    w·exp(-||x-p||²/σ²)          Σ f_i(x) · decay(t)
    ──────────────────     →     ──────────────────────
    Single Gaussian              Field Value at Point

    positions ──Delaunay──► β₀ = components
                           β₁ = holes (E - V + C)
```

Instead of solving constraint graphs with SAT solvers or SMT, you embed constraints in a continuous field and read the field value, gradient, and topology at any point. Propagation lets constraints adjust toward consensus — like relaxation but without discrete search.

## Quick Start

```bash
# Install
cargo install field-core

# Create a constraint field and embed positions
field embed --position "0.0,0.0" --weight 1.0 --stiffness 200
field embed --position "5.0,5.0" --weight 0.9 --stiffness 150

# Read the field at a point
field read --query "2.5,2.5" --time 0

# Run propagation
field propagate

# Check topology
field topology

# List all nails (positions)
field nails

# Field health
field status
```

**State** is automatically persisted to `.field-state.json` in the current directory.

## Library Usage

```rust
use field_core::{ConstraintField, FieldQuery};

let mut field = ConstraintField::new();

// Embed constraints
field.embed(vec![0.0, 0.0], 1.0, 200.0);
field.embed(vec![5.0, 5.0], 0.9, 150.0);

// Read the field
let reading = field.read(&FieldQuery {
    point: vec![2.5, 2.5],
    time: 0.0,
    tolerance: 3.0,
});
println!("Field value: {}", reading.value);
println!("Gradient: {:?}", reading.gradient);
println!("Confidence: {}", reading.confidence);

// Analyze topology
let topo = field.topology();
println!("Connected: {}", topo.is_connected);
println!("Betti numbers (β₀, β₁): {:?}", topo.betti_numbers);

// Propagate
let result = field.propagate();
println!("Updated: {}, max Δ: {}", result.updated_count, result.max_distance);
```

## The Field Paradigm

| Concept | Constraint Graph | Constraint Field |
|---------|-----------------|-----------------|
| Representation | Nodes + Edges | Gaussian bases |
| Solving | SAT/SMT search | Gradient descent |
| Topology | Manual labeling | Automatic (Delaunay) |
| Change propagation | Event-based | Field diffusion |
| Confidence | Boolean | Continuous [0, 1] |
| Gradient | N/A | Available at any point |

## Crate

- **GitHub**: [SuperInstance/field-core](https://github.com/SuperInstance/field-core)
- **License**: MIT

## Architecture

See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for design decisions.
See [docs/EXAMPLES.md](docs/EXAMPLES.md) for usage scenarios.
