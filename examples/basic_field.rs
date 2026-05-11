/// Basic field example: embed, read, propagate, topology.
use field_core::{ConstraintField, FieldQuery};

fn main() {
    println!("=== field-core: Basic Example ===\n");

    // Create a field
    let mut field = ConstraintField::new();

    // Embed three positions
    field.embed(vec![0.0, 0.0], 1.0, 5.0);
    field.embed(vec![5.0, 0.0], 0.9, 5.0);
    field.embed(vec![2.5, 5.0], 0.8, 5.0);

    println!("Positions: {}", field.len());
    for p in field.positions() {
        println!("  #{} @ {:?} (w={:.2}, σ={:.2})", p.id, p.coords, p.weight, p.stiffness);
    }

    // Read at center
    let reading = field.read(&FieldQuery {
        point: vec![2.5, 2.5],
        time: 0.0,
        tolerance: 3.0,
    });
    println!("\nRead at (2.5, 2.5):");
    println!("  Value: {:.6}", reading.value);
    println!("  Gradient: ({:.4}, {:.4})", reading.gradient[0], reading.gradient[1]);
    println!("  Confidence: {:.4}", reading.confidence);
    println!("  Nearby IDs: {:?}", reading.nearby_positions);

    // Topology
    let topo = field.topology();
    println!("\nTopology:");
    println!("  Connected: {}", topo.is_connected);
    println!("  Betti numbers: {:?}", topo.betti_numbers);
    println!("  Breaches: {} found", topo.breaches.len());

    // Propagate
    println!("\nPropagating...");
    for step in 0..5 {
        let result = field.propagate();
        println!("  Step {}: {} updated, max Δ = {:.6}", step, result.updated_count, result.max_distance);
    }

    // Final positions
    println!("\nFinal positions:");
    for p in field.positions() {
        println!("  #{} @ {:?}", p.id, p.coords);
    }
}
