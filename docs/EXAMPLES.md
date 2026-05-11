# Examples

## 1. Simple Gradient Estimation

Embed a single constraint and read the field off-center to get gradient:

```rust
use field_core::{ConstraintField, FieldQuery};

let mut field = ConstraintField::new();
field.embed(vec![0.0, 0.0], 1.0, 2.0);

let reading = field.read(&FieldQuery {
    point: vec![1.0, 0.0],
    time: 0.0,
    tolerance: 3.0,
});

println!("Value: {:.4}", reading.value);   // ~0.7788
println!("Gradient: {:.4}, {:.4}", reading.gradient[0], reading.gradient[1]);
// Gradient at (1,0) points toward origin: ~-0.3894, 0.0
```

## 2. Multi-Constraint Field With Propagation

Embed multiple constraints and let them propagate to consensus:

```rust
use field_core::ConstraintField;

let mut field = ConstraintField::new();
field.embed(vec![0.0, 0.0], 0.9, 5.0);
field.embed(vec![10.0, 0.0], 0.8, 5.0);
field.embed(vec![5.0, 8.0], 0.7, 5.0);

println!("Before: {:?}", field.positions().iter().map(|p| &p.coords).collect::<Vec<_>>());

for step in 0..10 {
    let result = field.propagate();
    println!("Step {}: {} updated, max Δ {:.4}", step, result.updated_count, result.max_distance);
}

println!("After: {:?}", field.positions().iter().map(|p| &p.coords).collect::<Vec<_>>());
// Positions should have converged toward each other
```

## 3. Topology Detection — Circle With a Hole

A ring of 16 points forms a topological hole:

```rust
use field_core::ConstraintField;

let mut field = ConstraintField::new();
let n = 16;
for i in 0..n {
    let angle = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
    field.embed(vec![10.0 * angle.cos(), 10.0 * angle.sin()], 1.0, 5.0);
}

let topo = field.topology();
println!("Connected: {}", topo.is_connected);           // true
println!("Betti numbers: {:?}", topo.betti_numbers);    // [1, 1] — 1 hole
println!("Breaches: {:?}", topo.breaches);
```

## 4. Disconnected Clusters in 3D

Two clusters far apart are detected as disconnected:

```rust
use field_core::ConstraintField;

let mut field = ConstraintField::new();
// Cluster A
field.embed(vec![0.0, 0.0, 0.0], 1.0, 2.0);
field.embed(vec![1.0, 1.0, 0.0], 0.9, 2.0);
// Cluster B (far away)
field.embed(vec![100.0, 100.0, 100.0], 1.0, 2.0);
field.embed(vec![101.0, 101.0, 101.0], 0.9, 2.0);

let topo = field.topology();
println!("Connected: {}", topo.is_connected);  // false
println!("Components: {}", topo.betti_numbers[0]);  // 2
```

## 5. CLI Workflow

```bash
# Start fresh
rm -f .field-state.json

# Embed nails
field embed --position "0,0" --weight 1.0 --stiffness 10
field embed --position "10,0" --weight 0.9 --stiffness 10
field embed --position "5,8" --weight 0.8 --stiffness 10

# Check status
field status

# Read at center
field read --query "5,4" --time 0

# Propagate 5 times
for i in 1 2 3 4 5; do
  field propagate
done

# Check final topology
field topology

# List all nails
field nails
```

## 6. Programmatic Classification

Use field reading confidence as a classification signal:

```rust
use field_core::{ConstraintField, FieldQuery};

let mut field = ConstraintField::new();
// Train: embed known positions
field.embed(vec![2.0, 3.0], 1.0, 1.0); // Class A
field.embed(vec![8.0, 7.0], 1.0, 1.0); // Class B

// Classify a new point
let reading = field.read(&FieldQuery {
    point: vec![2.5, 3.5],
    time: 0.0,
    tolerance: 3.0,
});

if reading.confidence > 0.5 {
    println!("Class A (confidence: {:.2})", reading.confidence);
}
```
