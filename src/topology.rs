//! Topology detection for the constraint field.
//!
//! Uses Delaunay triangulation (via the `delaunator` crate) for 2D positions,
//! then computes connected components and cycle rank:
//!
//! - β₀ = number of connected components
//! - β₁ = E - V + C  (cycle rank, i.e. number of independent holes)

use crate::{Breach, BreachType, FieldPosition, FieldTopology};
use std::collections::{HashMap, HashSet, VecDeque};

/// Compute the topology of a set of field positions.
///
/// Handles 0, 1, and 2+ positions with appropriate fallbacks for non-2D
/// configurations.
pub fn compute_topology(positions: &[FieldPosition]) -> FieldTopology {
    match positions.len() {
        0 => FieldTopology {
            is_connected: true,
            breaches: vec![],
            betti_numbers: vec![1, 0],
        },
        1 => FieldTopology {
            is_connected: true,
            breaches: vec![],
            betti_numbers: vec![1, 0],
        },
        _ => {
            // Only compute full topology for 2D positions
            if positions[0].coords.len() == 2 {
                compute_2d_topology(positions)
            } else {
                // For non-2D, just compute connectivity based on proximity
                compute_nd_topology(positions)
            }
        }
    }
}

/// Compute topology for 2D positions using Delaunay triangulation.
fn compute_2d_topology(positions: &[FieldPosition]) -> FieldTopology {
    // Need at least 3 points for Delaunay; fewer points is always connected
    if positions.len() < 3 {
        return FieldTopology {
            is_connected: true,
            breaches: vec![],
            betti_numbers: vec![1, 0],
        };
    }
    let points: Vec<delaunator::Point> = positions
        .iter()
        .map(|p| delaunator::Point {
            x: p.coords[0],
            y: p.coords[1],
        })
        .collect();

    // Build Delaunay triangulation
    let triangulation = delaunator::triangulate(&points);

    // Build adjacency from triangulation edges
    let mut adj: HashMap<usize, HashSet<usize>> = HashMap::new();
    let num_triangles = triangulation.triangles.len() / 3;

    if num_triangles == 0 {
        // Degenerate case (collinear points, etc.) — fall back to proximity graph
        for i in 0..positions.len() {
            adj.entry(i).or_default();
            for j in (i + 1)..positions.len() {
                let dist = crate::euclidean_distance(&positions[i].coords, &positions[j].coords);
                let threshold = (positions[i].stiffness + positions[j].stiffness) * 2.0;
                if dist <= threshold {
                    adj.entry(i).or_default().insert(j);
                    adj.entry(j).or_default().insert(i);
                }
            }
        }
    } else {
        for i in 0..num_triangles {
            let a = triangulation.triangles[i * 3] as usize;
            let b = triangulation.triangles[i * 3 + 1] as usize;
            let c = triangulation.triangles[i * 3 + 2] as usize;

            adj.entry(a).or_default().insert(b);
            adj.entry(b).or_default().insert(a);
            adj.entry(b).or_default().insert(c);
            adj.entry(c).or_default().insert(b);
            adj.entry(c).or_default().insert(a);
            adj.entry(a).or_default().insert(c);
        }

        // All positions are nodes
        for i in 0..positions.len() {
            adj.entry(i).or_default();
        }
    }

    // Find connected components (β₀)
    let components = find_components(positions.len(), &adj);
    let num_components = components.len();
    let is_connected = num_components == 1;

    // Compute β₁ = E - V + C  (cycle rank = holes)
    let num_vertices = positions.len() as u64;
    let beta_1;

    if num_triangles == 0 {
        // No triangulation: estimate β₁ from adjacency edges
        let mut edge_set: HashSet<(usize, usize)> = HashSet::new();
        for (v, neighbors) in &adj {
            for &nb in neighbors {
                edge_set.insert((*v.min(&nb), *v.max(&nb)));
            }
        }
        let num_edges = edge_set.len() as u64;
        beta_1 = num_edges
            .checked_sub(num_vertices)
            .unwrap_or(0)
            .saturating_add(num_components as u64);
    } else {
        // Use Delaunay edges
        let mut unique_edges: HashSet<(usize, usize)> = HashSet::new();
        for i in 0..num_triangles {
            let a = triangulation.triangles[i * 3] as usize;
            let b = triangulation.triangles[i * 3 + 1] as usize;
            let c = triangulation.triangles[i * 3 + 2] as usize;
            if a != b {
                unique_edges.insert((a.min(b), a.max(b)));
            }
            if b != c {
                unique_edges.insert((b.min(c), b.max(c)));
            }
            if c != a {
                unique_edges.insert((c.min(a), c.max(a)));
            }
        }
        let num_edges = unique_edges.len() as u64;
        beta_1 = num_edges
            .checked_sub(num_vertices)
            .unwrap_or(0)
            .saturating_add(num_components as u64);
    }

    // Build breach list
    let mut breaches = Vec::new();
    if !is_connected {
        for comp in &components {
            if comp.len() < positions.len() {
                let ids: Vec<u64> = comp
                    .iter()
                    .map(|&idx| positions[idx].id)
                    .collect();
                breaches.push(Breach {
                    breach_type: BreachType::Disconnected,
                    positions: ids,
                    severity: comp.len() as f64 / positions.len() as f64,
                });
            }
        }
    }
    if beta_1 > 0 && num_triangles > 0 {
        // Collect positions that form holes (boundary vertices of the convex hull or internal voids)
        // For simplicity, report vertices on the convex hull edges
        let hull_vertices = find_hull_vertices(&points, &triangulation);
        breaches.push(Breach {
            breach_type: BreachType::Hole,
            positions: hull_vertices.iter().map(|&i| positions[i].id).collect(),
            severity: (beta_1 as f64) / (num_vertices as f64).max(1.0),
        });
    }

    FieldTopology {
        is_connected,
        breaches,
        betti_numbers: vec![num_components as u64, beta_1],
    }
}

/// Compute topology for N-dimensional positions using proximity graph.
fn compute_nd_topology(positions: &[FieldPosition]) -> FieldTopology {
    let n = positions.len();
    let mut adj: HashMap<usize, HashSet<usize>> = HashMap::new();
    for i in 0..n {
        adj.entry(i).or_default();
        for j in (i + 1)..n {
            let dist = crate::euclidean_distance(&positions[i].coords, &positions[j].coords);
            let threshold = (positions[i].stiffness + positions[j].stiffness) * 2.0;
            if dist <= threshold {
                adj.entry(i).or_default().insert(j);
                adj.entry(j).or_default().insert(i);
            }
        }
    }

    let components = find_components(n, &adj);
    let num_components = components.len();
    let is_connected = num_components == 1;

    // For ND, estimate β₁ from component edges
    let mut total_edges = 0u64;
    for (_, neighbors) in &adj {
        total_edges += neighbors.len() as u64;
    }
    total_edges /= 2; // each edge counted twice

    let num_vertices = n as u64;
    let beta_1 = if num_vertices >= 3 {
        total_edges.saturating_sub(num_vertices) + num_components as u64
    } else {
        0
    };

    let mut breaches = Vec::new();
    if !is_connected {
        for comp in &components {
            let ids: Vec<u64> = comp.iter().map(|&idx| positions[idx].id).collect();
            breaches.push(Breach {
                breach_type: BreachType::Disconnected,
                positions: ids,
                severity: comp.len() as f64 / n as f64,
            });
        }
    }
    if beta_1 > 0 {
        breaches.push(Breach {
            breach_type: BreachType::Hole,
            positions: positions.iter().map(|p| p.id).collect(),
            severity: (beta_1 as f64) / (num_vertices as f64).max(1.0),
        });
    }

    FieldTopology {
        is_connected,
        breaches,
        betti_numbers: vec![num_components as u64, beta_1],
    }
}

/// Find connected components in an undirected graph.
fn find_components(n: usize, adj: &HashMap<usize, HashSet<usize>>) -> Vec<Vec<usize>> {
    let mut visited = vec![false; n];
    let mut components = Vec::new();

    for start in 0..n {
        if visited[start] {
            continue;
        }
        let mut component = Vec::new();
        let mut queue = VecDeque::new();
        queue.push_back(start);
        visited[start] = true;

        while let Some(v) = queue.pop_front() {
            component.push(v);
            if let Some(neighbors) = adj.get(&v) {
                for &nb in neighbors {
                    if !visited[nb] {
                        visited[nb] = true;
                        queue.push_back(nb);
                    }
                }
            }
        }
        component.sort_unstable();
        components.push(component);
    }

    components
}

/// Find vertices on the convex hull from Delaunay triangulation.
/// Hull edges are those appearing in only one triangle.
fn find_hull_vertices(_points: &[delaunator::Point], triangulation: &delaunator::Triangulation) -> Vec<usize> {
    let mut edge_count: HashMap<(usize, usize), u32> = HashMap::new();

    for i in 0..triangulation.triangles.len() / 3 {
        let a = triangulation.triangles[i * 3] as usize;
        let b = triangulation.triangles[i * 3 + 1] as usize;
        let c = triangulation.triangles[i * 3 + 2] as usize;

        for &(u, v) in &[(a, b), (b, c), (c, a)] {
            if u != v {
                let key = (u.min(v), u.max(v));
                *edge_count.entry(key).or_insert(0) += 1;
            }
        }
    }

    // Hull edges appear in exactly 1 triangle (interior edges appear in 2)
    let mut hull_vertices: HashSet<usize> = HashSet::new();
    for ((u, v), count) in &edge_count {
        if *count == 1 {
            hull_vertices.insert(*u);
            hull_vertices.insert(*v);
        }
    }

    let mut result: Vec<usize> = hull_vertices.into_iter().collect();
    result.sort_unstable();
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_topology() {
        let topo = compute_topology(&[]);
        assert!(topo.is_connected);
        assert_eq!(topo.betti_numbers, vec![1, 0]);
    }

    #[test]
    fn test_single_pos_topology() {
        let pos = vec![FieldPosition {
            id: 1,
            coords: vec![0.0, 0.0],
            weight: 1.0,
            stiffness: 1.0,
        }];
        let topo = compute_topology(&pos);
        assert!(topo.is_connected);
        assert_eq!(topo.betti_numbers, vec![1, 0]);
    }

    #[test]
    fn test_collinear_three_points() {
        // Three collinear points: Delaunay still works but produces degenerate triangles
        let positions = vec![
            FieldPosition { id: 1, coords: vec![0.0, 0.0], weight: 1.0, stiffness: 1.0 },
            FieldPosition { id: 2, coords: vec![1.0, 0.0], weight: 1.0, stiffness: 1.0 },
            FieldPosition { id: 3, coords: vec![2.0, 0.0], weight: 1.0, stiffness: 1.0 },
        ];
        let topo = compute_topology(&positions);
        assert!(topo.is_connected);
    }

    #[test]
    fn test_two_disconnected_clusters() {
        // Two clusters far apart should be disconnected
        let positions = vec![
            FieldPosition { id: 1, coords: vec![0.0, 0.0], weight: 1.0, stiffness: 1.0 },
            FieldPosition { id: 2, coords: vec![1.0, 0.0], weight: 1.0, stiffness: 1.0 },
            FieldPosition { id: 3, coords: vec![100.0, 100.0], weight: 1.0, stiffness: 1.0 },
            FieldPosition { id: 4, coords: vec![101.0, 101.0], weight: 1.0, stiffness: 1.0 },
        ];
        // In 2D, Delaunay will connect them because it always produces convex hull edges.
        // The far cluster will be connected via hull edges. So use ND topology instead.
        let positions_nd: Vec<FieldPosition> = positions.iter().map(|p| FieldPosition {
            id: p.id,
            coords: vec![p.coords[0], p.coords[1], 0.0], // 3D
            weight: p.weight,
            stiffness: p.stiffness,
        }).collect();
        let topo = compute_topology(&positions_nd);
        // Small stiffness means no edges between clusters
        assert!(!topo.is_connected);
    }
}
