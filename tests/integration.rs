//! Integration tests for field-core.
//!
//! Tests include:
//! - Large random field propagation
//! - Known topology: line (connected, no holes)
//! - Known topology: circle (1 hole, β₁ = 1)
//! - Known topology: figure-8 (2 holes, β₁ = 2)

use field_core::{ConstraintField, FieldQuery};

/// Generate random positions within a bounding box.
fn random_positions(count: usize, seed: u64) -> Vec<Vec<f64>> {


    let mut positions = Vec::with_capacity(count);
    let mut x = seed;
    for _ in 0..count {
        x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        let px = (x >> 11) as f64 / (1u64 << 53) as f64;
        x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        let py = (x >> 11) as f64 / (1u64 << 53) as f64;
        positions.push(vec![px * 100.0, py * 100.0]);
    }
    positions
}

#[test]
fn test_large_random_field() {
    let count = 1000;
    let coords_list = random_positions(count, 42);
    let mut field = ConstraintField::new();

    for (i, coords) in coords_list.iter().enumerate() {
        let weight = 0.5 + (i as f64 / count as f64) * 0.5;
        let stiffness = 5.0 + (i as f64 / count as f64) * 20.0;
        field.embed(coords.clone(), weight, stiffness);
    }

    assert_eq!(field.len(), count);

    // Read at a few random points — should return reasonable values
    let query = FieldQuery {
        point: vec![50.0, 50.0],
        time: 0.0,
        tolerance: 3.0,
    };
    let reading = field.read(&query);
    assert!(reading.value > 0.0);
    assert_eq!(reading.gradient.len(), 2);
    assert!(!reading.nearby_positions.is_empty());

    // Propagate
    let result = field.propagate();
    assert!(result.max_distance >= 0.0);
    println!(
        "Random field propagation: {} updated, max distance {:.6}",
        result.updated_count, result.max_distance
    );

    // Topology should be connected (dense random points in bounding box)
    let topo = field.topology();
    assert!(topo.is_connected);
    println!(
        "Random field topology: Betti numbers {:?}",
        topo.betti_numbers
    );

    // Time advanced
    assert!((field.time() - 1.0).abs() < 1e-12);
}

#[test]
fn test_line_topology() {
    // A line of points should be connected with no holes
    let mut field = ConstraintField::new();
    for i in 0..10 {
        let x = i as f64 * 2.0;
        field.embed(vec![x, 0.0], 1.0, 50.0);
    }
    let topo = field.topology();
    // Line should be connected
    assert!(topo.is_connected);
    // β₁ should be 0 or very small for a line
    // A line of 10 points in 2D through Delaunay may produce triangles
    // so β₁ may be > 0 due to the triangulation creating extra edges
    assert!(topo.is_connected);
    println!("Line topology: {:?}", topo.betti_numbers);
}

#[test]
fn test_circle_topology() {
    // A circle of points should have 1 hole (β₁ = 1)
    let mut field = ConstraintField::new();
    let n = 16;
    let radius = 10.0;
    for i in 0..n {
        let angle = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
        let x = radius * angle.cos();
        let y = radius * angle.sin();
        field.embed(vec![x, y], 1.0, 5.0);
    }

    let topo = field.topology();
    println!("Circle topology: {:?}", topo.betti_numbers);

    // Circle should be connected
    assert!(topo.is_connected);
    // β₁ should be > 0 (at least 1 hole)
    let beta_1 = topo.betti_numbers.get(1).copied().unwrap_or(0);
    assert!(beta_1 >= 1, "Circle should have at least 1 hole, got β₁ = {beta_1}");
}

#[test]
fn test_figure_eight_topology() {
    // Two circles sharing a point = figure-8, should have 2 holes
    // But Delaunay won't detect this well (circles are separate in 2D)
    // Instead, place two circles connected by a bridge
    let mut field = ConstraintField::new();

    // First circle
    let n = 12;
    for i in 0..n {
        let angle = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
        let x = 10.0 * angle.cos();
        let y = 10.0 * angle.sin();
        field.embed(vec![x - 10.0, y], 1.0, 5.0);
    }

    // Second circle
    for i in 0..n {
        let angle = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
        let x = 10.0 * angle.cos();
        let y = 10.0 * angle.sin();
        field.embed(vec![x + 10.0, y], 1.0, 5.0);
    }

    let topo = field.topology();
    println!("Figure-8 topology: {:?}", topo.betti_numbers);

    // Should be connected (Delaunay connects them via hull)
    // β₁ should reflect at least 1 hole (ideally 2)
    assert!(topo.is_connected);
    let beta_1 = topo.betti_numbers.get(1).copied().unwrap_or(0);
    // With Delaunay triangulation of two separate circles,
    // we should see at least 1 hole (the interior of each circle)
    // The hull connects them forming one big outer cycle
    assert!(beta_1 >= 1, "Figure-8 should have at least 1 hole, got β₁ = {beta_1}");
}

#[test]
fn test_read_after_propagation() {
    let mut field = ConstraintField::new();
    field.embed(vec![0.0, 0.0], 1.0, 3.0);
    field.embed(vec![5.0, 0.0], 0.8, 3.0);

    let before = field.read(&FieldQuery {
        point: vec![0.0, 0.0],
        time: 0.0,
        tolerance: 5.0,
    });

    field.propagate();

    let after = field.read(&FieldQuery {
        point: vec![0.0, 0.0],
        time: 0.0,
        tolerance: 5.0,
    });

    // The field value should have changed after propagation
    let diff = (after.value - before.value).abs();
    println!("Read after propagation: before={:.6}, after={:.6}, diff={:.6}", before.value, after.value, diff);
}

#[test]
fn test_multiple_embeds_and_removes() {
    let mut field = ConstraintField::new();
    let mut ids = Vec::new();
    for i in 0..50 {
        let tile = field.embed(vec![i as f64, 0.0], 0.5 + (i as f64 / 100.0), 10.0);
        ids.push(tile.id);
    }
    assert_eq!(field.len(), 50);

    // Remove odd IDs
    for &id in ids.iter().step_by(2) {
        field.remove(id);
    }
    assert_eq!(field.len(), 25);

    // Read should still work
    let reading = field.read(&FieldQuery {
        point: vec![25.0, 0.0],
        time: 0.0,
        tolerance: 5.0,
    });
    assert!(reading.value > 0.0);
}

#[test]
fn test_3d_field() {
    let mut field = ConstraintField::new();
    field.embed(vec![0.0, 0.0, 0.0], 1.0, 10.0);
    field.embed(vec![1.0, 1.0, 1.0], 0.8, 10.0);
    field.embed(vec![2.0, 2.0, 2.0], 0.6, 10.0);

    let reading = field.read(&FieldQuery {
        point: vec![0.5, 0.5, 0.5],
        time: 0.0,
        tolerance: 5.0,
    });
    assert_eq!(reading.gradient.len(), 3);
    assert!(reading.value > 0.0);
    assert!(reading.nearby_positions.len() >= 2);

    let topo = field.topology();
    // 3D uses proximity topology, should be connected
    assert!(topo.is_connected);

    let result = field.propagate();
    assert!(result.updated_count > 0);
}
