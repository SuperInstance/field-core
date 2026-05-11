/// Topology demonstration: embedding a circle and detecting the hole.
use field_core::ConstraintField;

fn main() {
    println!("=== Topology Demo: Circle Detection ===\n");

    // Create a ring of points — this creates a topological hole
    let mut field = ConstraintField::new();
    let n = 16;
    for i in 0..n {
        let angle = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
        let x = 10.0 * angle.cos();
        let y = 10.0 * angle.sin();
        field.embed(vec![x, y], 1.0, 3.0);
    }

    println!("Circle with {} points", n);
    let topo = field.topology();
    println!("Connected: {}", topo.is_connected);
    println!("Betti numbers: {:?}", topo.betti_numbers);
    println!("  β₀ (components): {}", topo.betti_numbers[0]);
    println!("  β₁ (holes): {}", topo.betti_numbers.get(1).copied().unwrap_or(0));
    println!("Breaches: {}", topo.breaches.len());

    for breach in &topo.breaches {
        println!("  - {:?}: severity={:.2}, positions={}", breach.breach_type, breach.severity, breach.positions.len());
    }

    // Now add a point in the center — should become connected
    println!("\nAdding center point...");
    field.embed(vec![0.0, 0.0], 1.0, 5.0);

    let topo2 = field.topology();
    println!("Connected: {}", topo2.is_connected);
    println!("Betti numbers: {:?}", topo2.betti_numbers);
    // The center point connects to the ring, potentially filling the hole
    // but the ring itself still encloses a region
}
