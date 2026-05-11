//! # field-core
//!
//! A continuous constraint field library for topological constraint satisfaction.
//!
//! ## Overview
//!
//! Each constraint in the field is a Gaussian basis function:
//!
//! ```text
//!   f_i(x) = w_i · exp(-||x - p_i||² / σ_i²)
//! ```
//!
//! where `w_i` is the weight, `p_i` is the position, and `σ_i` is the stiffness
//! (spread) of the constraint. The field at any point is the sum over all
//! basis functions, weighted by confidence decay over time.
//!
//! Topology is detected via Delaunay triangulation of positions, computing
//! Betti numbers (β₀ = components, β₁ = holes/cycles).

mod topology;

use topology::compute_topology;

/// A single constraint position in the continuous field.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct FieldPosition {
    /// Unique identifier
    pub id: u64,
    /// Coordinates in the field (supports N dimensions)
    pub coords: Vec<f64>,
    /// Confidence weight [0, 1]
    pub weight: f64,
    /// Stiffness (spread) of the Gaussian — higher = wider influence
    pub stiffness: f64,
}

/// A query against the field at a given point in space and time.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct FieldQuery {
    /// Query point coordinates
    pub point: Vec<f64>,
    /// Time offset for confidence decay
    pub time: f64,
    /// Tolerance for nearby position matching
    pub tolerance: f64,
}

/// The result of reading the field at a query point.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct FieldReading {
    /// Scalar field value at the query point
    pub value: f64,
    /// Gradient vector (partial derivatives) at the query point
    pub gradient: Vec<f64>,
    /// Confidence of the reading [0, 1]
    pub confidence: f64,
    /// IDs of nearby positions influencing this reading
    pub nearby_positions: Vec<u64>,
}

/// The topological state of the constraint field.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct FieldTopology {
    /// Whether the field forms a single connected component
    pub is_connected: bool,
    /// Topology breaches (disconnected regions, holes)
    pub breaches: Vec<Breach>,
    /// Betti numbers: β₀ = components, β₁ = cycles/holes
    pub betti_numbers: Vec<u64>,
}

/// A topology breach in the constraint field.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Breach {
    /// Type of breach
    pub breach_type: BreachType,
    /// Position IDs involved
    pub positions: Vec<u64>,
    /// Severity of the breach [0, 1]
    pub severity: f64,
}

/// Types of topology breaches.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub enum BreachType {
    /// A disconnected component (β₀ breach)
    Disconnected,
    /// A cycle or hole (β₁ breach)
    Hole,
}

/// Result of a propagation step.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct PropagationResult {
    /// Number of positions that were updated
    pub updated_count: usize,
    /// Maximum distance any position moved
    pub max_distance: f64,
    /// Whether the topology changed as a result
    pub topology_changed: bool,
}

/// The continuous constraint field — the central data structure.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ConstraintField {
    positions: Vec<FieldPosition>,
    next_id: u64,
    time: f64,
}

impl ConstraintField {
    /// Create a new empty constraint field.
    pub fn new() -> Self {
        ConstraintField {
            positions: Vec::new(),
            next_id: 1,
            time: 0.0,
        }
    }

    /// Create a constraint field with pre-existing positions.
    pub fn with_positions(positions: Vec<FieldPosition>) -> Self {
        let max_id = positions.iter().map(|p| p.id).max().unwrap_or(0);
        ConstraintField {
            positions,
            next_id: max_id + 1,
            time: 0.0,
        }
    }

    /// Embed a new constraint position in the field.
    ///
    /// Returns the generated [`Tile`] representing the new position.
    pub fn embed(&mut self, coords: Vec<f64>, weight: f64, stiffness: f64) -> Tile {
        let id = self.next_id;
        self.next_id += 1;
        let position = FieldPosition {
            id,
            coords,
            weight: weight.clamp(0.0, 1.0),
            stiffness: stiffness.max(1e-10),
        };
        self.positions.push(position);
        Tile { id }
    }

    /// Read the field at a query point, returning the field value, gradient, confidence,
    /// and nearby position IDs.
    pub fn read(&self, query: &FieldQuery) -> FieldReading {
        if self.positions.is_empty() {
            return FieldReading {
                value: 0.0,
                gradient: vec![],
                confidence: 0.0,
                nearby_positions: vec![],
            };
        }

        // Determine dimensionality from first position
        let dim = self.positions[0].coords.len();
        let mut total_value = 0.0;
        let mut gradient = vec![0.0; dim];
        let mut total_weight = 0.0;
        let mut nearby = Vec::new();

        for pos in &self.positions {
            if pos.coords.len() != dim {
                continue; // skip mismatched dimensions
            }

            let sq_dist = squared_distance(&query.point, &pos.coords);
            let sigma = pos.stiffness;
            let sigma_sq = sigma * sigma;
            let gaussian = pos.weight * (-sq_dist / sigma_sq).exp();

            // Apply time confidence decay
            let decay = (-query.time / 1000.0).exp();
            let contribution = gaussian * decay;

            total_value += contribution;
            total_weight += pos.weight;

            // Gradient: partial derivatives of Gaussian
            // d/dx_i = -2*(x_i - p_i)/σ² * gaussian * decay
            for i in 0..dim {
                let dx = query.point[i] - pos.coords[i];
                gradient[i] += (-2.0 * dx / sigma_sq) * contribution;
            }

            // Check if this position is "nearby" (within tolerance * stiffness)
            let radius = sigma * query.tolerance;
            if sq_dist.sqrt() <= radius {
                nearby.push(pos.id);
            }
        }

        let confidence = if total_weight > 0.0 {
            let max_contrib = self
                .positions
                .iter()
                .map(|p| {
                    let sq_d = squared_distance(&query.point, &p.coords);
                    p.weight * (-sq_d / (p.stiffness * p.stiffness)).exp()
                })
                .fold(0.0_f64, f64::max);
            max_contrib / total_weight
        } else {
            0.0
        };

        FieldReading {
            value: total_value,
            gradient,
            confidence,
            nearby_positions: nearby,
        }
    }

    /// Propagate changes across the field.
    ///
    /// Each position influences its neighbors within stiffness radius.
    /// Positions with overlapping influence regions adjust toward consensus.
    pub fn propagate(&mut self) -> PropagationResult {
        if self.positions.len() < 2 {
            return PropagationResult {
                updated_count: 0,
                max_distance: 0.0,
                topology_changed: false,
            };
        }

        let dim = self.positions[0].coords.len();
        let positions_copy = self.positions.clone();
        let mut updated_count = 0;
        let mut max_distance: f64 = 0.0;
        let old_topology = self.topology();
        let mut topology_changed = false;

        for pos in &mut self.positions {
            let old_coords = pos.coords.clone();
            let sigma = pos.stiffness;
            let influence_radius = sigma * 3.0; // 3σ captures ~99% of influence

            // Gather neighbors within influence radius
            let mut force = vec![0.0; dim];
            let mut neighbor_count = 0;

            for other in &positions_copy {
                if other.id == pos.id {
                    continue;
                }
                let dist = euclidean_distance(&pos.coords, &other.coords);
                if dist <= influence_radius {
                    // Attractive force toward neighbor, weighted by neighbor's influence
                    let strength = other.weight * (-(dist * dist) / (sigma * sigma)).exp();
                    for i in 0..dim {
                        force[i] += (other.coords[i] - pos.coords[i]) * strength;
                    }
                    neighbor_count += 1;
                }
            }

            if neighbor_count > 0 {
                // Apply force (damped)
                let damping = 0.1;
                for i in 0..dim {
                    pos.coords[i] += force[i] * damping / neighbor_count as f64;
                }

                // Track displacement
                let displacement = euclidean_distance(&old_coords, &pos.coords);
                if displacement > 1e-12 {
                    updated_count += 1;
                    max_distance = max_distance.max(displacement);
                }
            }
        }

        self.time += 1.0;

        // Check topology change
        let new_topology = self.topology();
        if old_topology.is_connected != new_topology.is_connected
            || old_topology.betti_numbers != new_topology.betti_numbers
        {
            topology_changed = true;
        }

        PropagationResult {
            updated_count,
            max_distance,
            topology_changed,
        }
    }

    /// Analyze the topology of the constraint field.
    ///
    /// Uses Delaunay triangulation for 2D positions, then computes
    /// connected components (β₀) and cycle rank (β₁ = E - V + C).
    pub fn topology(&self) -> FieldTopology {
        compute_topology(&self.positions)
    }

    /// Return all positions in the field.
    pub fn positions(&self) -> &[FieldPosition] {
        &self.positions
    }

    /// Return the number of positions.
    pub fn len(&self) -> usize {
        self.positions.len()
    }

    /// Check whether the field has no positions.
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }

    /// Return the current time step.
    pub fn time(&self) -> f64 {
        self.time
    }

    /// Remove a position by ID.
    pub fn remove(&mut self, id: u64) -> Option<FieldPosition> {
        if let Some(idx) = self.positions.iter().position(|p| p.id == id) {
            Some(self.positions.remove(idx))
        } else {
            None
        }
    }
}

impl Default for ConstraintField {
    fn default() -> Self {
        Self::new()
    }
}

/// A tile representing a constraint embedded in the field.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Tile {
    pub id: u64,
}

/// Field status summary.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct FieldStatus {
    pub num_positions: usize,
    pub time: f64,
    pub is_connected: bool,
    pub betti_numbers: Vec<u64>,
    pub average_weight: f64,
    pub average_stiffness: f64,
}

impl From<&ConstraintField> for FieldStatus {
    fn from(field: &ConstraintField) -> Self {
        let n = field.positions.len();
        let topo = field.topology();
        let avg_w = if n > 0 {
            field.positions.iter().map(|p| p.weight).sum::<f64>() / n as f64
        } else {
            0.0
        };
        let avg_s = if n > 0 {
            field.positions.iter().map(|p| p.stiffness).sum::<f64>() / n as f64
        } else {
            0.0
        };
        FieldStatus {
            num_positions: n,
            time: field.time,
            is_connected: topo.is_connected,
            betti_numbers: topo.betti_numbers,
            average_weight: avg_w,
            average_stiffness: avg_s,
        }
    }
}

// --- Utility functions ---

/// Squared Euclidean distance between two vectors.
pub(crate) fn squared_distance(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y) * (x - y))
        .sum()
}

/// Euclidean distance between two vectors.
pub(crate) fn euclidean_distance(a: &[f64], b: &[f64]) -> f64 {
    squared_distance(a, b).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_field_is_empty() {
        let field = ConstraintField::new();
        assert!(field.is_empty());
        assert_eq!(field.len(), 0);
    }

    #[test]
    fn test_embed_increases_count() {
        let mut field = ConstraintField::new();
        let tile = field.embed(vec![0.0, 0.0], 1.0, 1.0);
        assert_eq!(tile.id, 1);
        assert_eq!(field.len(), 1);
    }

    #[test]
    fn test_embed_assigns_incrementing_ids() {
        let mut field = ConstraintField::new();
        let t1 = field.embed(vec![0.0, 0.0], 0.5, 2.0);
        let t2 = field.embed(vec![1.0, 1.0], 0.8, 1.0);
        assert_eq!(t1.id, 1);
        assert_eq!(t2.id, 2);
    }

    #[test]
    fn test_embed_clamps_weight() {
        let mut field = ConstraintField::new();
        field.embed(vec![0.0], 1.5, 1.0); // weight > 1 should clamp
        assert!((field.positions()[0].weight - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_embed_ensures_min_stiffness() {
        let mut field = ConstraintField::new();
        field.embed(vec![0.0], 0.5, 0.0); // stiffness > 0
        assert!(field.positions()[0].stiffness > 0.0);
    }

    #[test]
    fn test_read_empty_field() {
        let field = ConstraintField::new();
        let query = FieldQuery {
            point: vec![0.0, 0.0],
            time: 0.0,
            tolerance: 1.0,
        };
        let reading = field.read(&query);
        assert!((reading.value - 0.0).abs() < 1e-12);
        assert_eq!(reading.gradient.len(), 0);
    }

    #[test]
    fn test_read_single_position() {
        let mut field = ConstraintField::new();
        field.embed(vec![0.0, 0.0], 1.0, 1.0);
        let query = FieldQuery {
            point: vec![0.0, 0.0],
            time: 0.0,
            tolerance: 3.0,
        };
        let reading = field.read(&query);
        // At the center, value = w * exp(0) = 1.0
        assert!((reading.value - 1.0).abs() < 1e-6);
        assert_eq!(reading.gradient.len(), 2);
        // Gradient should be ~0 at center of Gaussian
        assert!(reading.gradient.iter().all(|g| g.abs() < 1e-6));
        // Nearby should include position 1
        assert!(reading.nearby_positions.contains(&1));
    }

    #[test]
    fn test_read_with_time_decay() {
        let mut field = ConstraintField::new();
        field.embed(vec![0.0, 0.0], 1.0, 1.0);
        let query = FieldQuery {
            point: vec![0.0, 0.0],
            time: 1000.0,
            tolerance: 3.0,
        };
        let reading = field.read(&query);
        // After 1000 time units: value = 1.0 * exp(-1) ≈ 0.3679
        assert!((reading.value - 0.367879).abs() < 1e-4);
    }

    #[test]
    fn test_gradient_off_center() {
        let mut field = ConstraintField::new();
        field.embed(vec![0.0, 0.0], 1.0, 2.0);
        let query = FieldQuery {
            point: vec![1.0, 0.0],
            time: 0.0,
            tolerance: 3.0,
        };
        let reading = field.read(&query);
        // Gradient should point from query toward the position center
        // d/dx at (1,0) for Gaussian centered at (0,0) with σ=2:
        // = -2*(1-0)/4 * exp(-1/4) = -0.5 * 0.7788 = -0.3894
        assert!(reading.gradient[0] < 0.0); // should point toward origin
        assert!((reading.gradient[1]).abs() < 1e-6); // y gradient ~ 0
    }

    #[test]
    fn test_propagate_no_positions() {
        let mut field = ConstraintField::new();
        let result = field.propagate();
        assert_eq!(result.updated_count, 0);
        assert!((result.max_distance - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_propagate_single_position_no_op() {
        let mut field = ConstraintField::new();
        field.embed(vec![0.0, 0.0], 1.0, 1.0);
        let result = field.propagate();
        assert_eq!(result.updated_count, 0);
    }

    #[test]
    fn test_propagate_two_positions() {
        let mut field = ConstraintField::new();
        field.embed(vec![0.0, 0.0], 1.0, 5.0);
        field.embed(vec![10.0, 10.0], 1.0, 5.0);
        let result = field.propagate();
        // With large stiffness, they should attract
        assert!(result.updated_count > 0);
    }

    #[test]
    fn test_topology_empty_field() {
        let field = ConstraintField::new();
        let topo = field.topology();
        assert!(topo.is_connected);
        assert!(topo.breaches.is_empty());
        assert_eq!(topo.betti_numbers, vec![1, 0]);
    }

    #[test]
    fn test_topology_single_position() {
        let mut field = ConstraintField::new();
        field.embed(vec![0.0, 0.0], 1.0, 1.0);
        let topo = field.topology();
        assert!(topo.is_connected);
        assert_eq!(topo.betti_numbers, vec![1, 0]);
    }

    #[test]
    fn test_topology_two_known_positions() {
        let mut field = ConstraintField::new();
        field.embed(vec![0.0, 0.0], 1.0, 1.0);
        field.embed(vec![1.0, 1.0], 1.0, 1.0);
        let topo = field.topology();
        // Two positions in 2D should be connected by Delaunay edge
        assert!(topo.is_connected);
    }

    #[test]
    fn test_remove_position() {
        let mut field = ConstraintField::new();
        let tile = field.embed(vec![0.0, 0.0], 1.0, 1.0);
        field.embed(vec![1.0, 1.0], 0.5, 1.0);
        assert_eq!(field.len(), 2);

        let removed = field.remove(tile.id);
        assert!(removed.is_some());
        assert_eq!(removed.unwrap().id, 1);
        assert_eq!(field.len(), 1);
    }

    #[test]
    fn test_remove_nonexistent() {
        let mut field = ConstraintField::new();
        field.embed(vec![0.0, 0.0], 1.0, 1.0);
        assert!(field.remove(999).is_none());
    }

    #[test]
    fn test_field_status() {
        let mut field = ConstraintField::new();
        field.embed(vec![0.0, 0.0], 1.0, 1.0);
        field.embed(vec![1.0, 0.0], 0.5, 2.0);
        let status: FieldStatus = (&field).into();
        assert_eq!(status.num_positions, 2);
        assert!((status.average_weight - 0.75).abs() < 1e-12);
        assert!((status.average_stiffness - 1.5).abs() < 1e-12);
    }
}
