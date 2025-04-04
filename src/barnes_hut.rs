use rand::Rng;
use std::{u64::MAX, vec};
use wasm_bindgen::prelude::*;

use web_sys::console;

/// @brief Point is a struct that holds the x and y coordinates of a point.
#[derive(Clone, Debug)]
#[wasm_bindgen]
pub struct Point {
    x: f32,
    y: f32,
}

impl Point {
    fn distance(&self, other: &Point) -> f32 {
        ((self.x - other.x).powi(2) + (self.y - other.y).powi(2)).sqrt()
    }
}

/// @brief PointView is a struct that holds references to the x and y
/// coordinates of a point from the particle buffer.
struct PointView<'a> {
    x: &'a f32,
    y: &'a f32,
}

struct ParticleView<'a> {
    x: &'a f32,
    y: &'a f32,
    mass: &'a f32,
}

impl ParticleView<'_> {
    fn point(&self) -> PointView {
        PointView {
            x: self.x,
            y: self.y,
        }
    }

    fn mass(&self) -> f32 {
        *self.mass
    }
}

impl Point {
    fn view(&self) -> PointView {
        PointView {
            x: &self.x,
            y: &self.y,
        }
    }
}

/// @brief BoundingBox is a struct that holds the top left and bottom right points of a
/// bounding box. This is used to determine the bounds of a quadtree node.
struct BoundingBox {
    top_right: Point,
    bottom_left: Point,
}

impl BoundingBox {
    fn width(&self) -> f32 {
        self.top_right.x - self.bottom_left.x
    }

    fn height(&self) -> f32 {
        self.top_right.y - self.bottom_left.y
    }
}

/// @brief ParticleBuffer is a struct that holds the x and y coordinates of particles.
/// The coordinates are stored in separate vectors in order to optimize iteration.
struct ParticleBuffer {
    x: Vec<f32>,
    y: Vec<f32>,

    velocity_x: Vec<f32>,
    velocity_y: Vec<f32>,

    force_x: Vec<f32>,
    force_y: Vec<f32>,

    mass: Vec<f32>,
}

impl ParticleBuffer {
    fn new(capacity: usize) -> Self {
        ParticleBuffer {
            x: Vec::with_capacity(capacity),
            y: Vec::with_capacity(capacity),
            velocity_x: Vec::with_capacity(capacity),
            velocity_y: Vec::with_capacity(capacity),
            force_x: Vec::with_capacity(capacity),
            force_y: Vec::with_capacity(capacity),
            mass: Vec::with_capacity(capacity),
        }
    }

    /// @brief Adds a particle to the buffer. In theory particles only should added at
    /// initialization time. Consider pushing this method to the constructor.
    fn add_particle(&mut self, x: f32, y: f32, mass: f32) {
        self.x.push(x);
        self.y.push(y);
        self.velocity_x.push(0.0);
        self.velocity_y.push(0.0);
        self.force_x.push(0.0);
        self.force_y.push(0.0);
        self.mass.push(mass);
    }

    /// @brief Returns a PointView for the particle at the given index.
    fn get_position(&self, index: usize) -> PointView {
        PointView {
            x: &self.x[index],
            y: &self.y[index],
        }
    }

    /// @brief Sets the position of the particle at the given index.
    fn set_position(&mut self, index: usize, x: f32, y: f32) {
        self.x[index] = x;
        self.y[index] = y;
    }

    /// @brief Returns a ParticleView for the particle at the given index.
    fn get_particle(&self, index: usize) -> ParticleView {
        ParticleView {
            x: &self.x[index],
            y: &self.y[index],
            mass: &self.mass[index],
        }
    }

    fn apply_force(&mut self, index: usize, force: &Point) {
        self.force_x[index] += force.x;
        self.force_y[index] += force.y;
    }

    /// @brief Computes the min and max of the x coordinates in the particle buffer.
    fn bounding_box(&self) -> BoundingBox {
        let (min_x, max_x) = self
            .x
            .iter()
            .fold((f32::INFINITY, f32::NEG_INFINITY), |(min, max), &val| {
                (min.min(val), max.max(val))
            });
        let (min_y, max_y) = self
            .y
            .iter()
            .fold((f32::INFINITY, f32::NEG_INFINITY), |(min, max), &val| {
                (min.min(val), max.max(val))
            });

        return BoundingBox {
            top_right: Point { x: max_x, y: max_y },
            bottom_left: Point { x: min_x, y: min_y },
        };
    }

    fn get_raw_x_coordinates(&self) -> *const f32 {
        self.x.as_ptr()
    }

    fn get_raw_y_coordinates(&self) -> *const f32 {
        self.y.as_ptr()
    }

    fn size(&self) -> usize {
        self.x.len()
    }
}
struct QuadTreeNode {
    bounding_box: BoundingBox,
    center_of_mass: Point,
    mass: f32,
    /// @brief The index of the particle that is contained in this node.
    /// If this node is a leaf, then this is the index of the particle.
    particle_indices: Vec<usize>,
    children: [Option<Box<QuadTreeNode>>; 4],
}

struct QuadTree {
    root: QuadTreeNode,
    max_depth: usize,
    max_particles_per_node: usize,
}

impl QuadTree {
    fn new(
        particles: &ParticleBuffer,
        max_depth: usize,
        max_particles_per_node: usize,
    ) -> QuadTree {
        let mass = particles.mass.iter().sum::<f32>();
        let center_of_mass = Point {
            x: particles
                .x
                .iter()
                .zip(particles.mass.iter())
                .map(|(&x, &m)| x * m)
                .sum::<f32>()
                / mass,
            y: particles
                .y
                .iter()
                .zip(particles.mass.iter())
                .map(|(&y, &m)| y * m)
                .sum::<f32>()
                / mass,
        };

        let mut tree = QuadTree {
            root: QuadTreeNode {
                bounding_box: particles.bounding_box(),
                particle_indices: (0..particles.x.len()).collect(),
                center_of_mass,
                mass,
                children: [None, None, None, None],
            },
            max_depth,
            max_particles_per_node,
        };

        tree.build(particles);

        tree
    }

    fn get_quad(center: PointView, point: PointView) -> usize {
        // encode the quadrant as a 2 bit number; this simplifies storing and accessing quandrants
        // given a bounding box and point.
        // |---|---|
        // | 2 | 3 |
        // |---|---|
        // | 0 | 1 |
        // |---|---|
        return (*point.x > *center.x) as usize | ((*point.y > *center.y) as usize) << 1;
    }

    fn get_quad_bounding_box(
        center: PointView,
        quad: usize,
        bounding_box: &BoundingBox,
    ) -> BoundingBox {
        let bottom_left = Point {
            x: if quad & 0b01 == 0 {
                bounding_box.bottom_left.x
            } else {
                *center.x
            },
            y: if quad & 0b10 == 0 {
                bounding_box.bottom_left.y
            } else {
                *center.y
            },
        };

        let top_right = Point {
            x: if quad & 0b01 == 0 {
                *center.x
            } else {
                bounding_box.top_right.x
            },
            y: if quad & 0b10 == 0 {
                *center.y
            } else {
                bounding_box.top_right.y
            },
        };
        return BoundingBox {
            top_right,
            bottom_left,
        };
    }

    fn build(&mut self, particles: &ParticleBuffer) {
        QuadTree::split(
            &mut self.root,
            particles,
            0,
            self.max_depth,
            self.max_particles_per_node,
        );
    }

    fn split(
        node: &mut QuadTreeNode,
        particles: &ParticleBuffer,
        depth: usize,
        max_depth: usize,
        max_particles_per_node: usize,
    ) {
        if node.particle_indices.len() <= max_particles_per_node || depth >= max_depth {
            // base case: no need to split further
            return;
        }

        let center = Point {
            x: (node.bounding_box.top_right.x + node.bounding_box.bottom_left.x) / 2.0,
            y: (node.bounding_box.top_right.y + node.bounding_box.bottom_left.y) / 2.0,
        };

        for index in &node.particle_indices {
            let particle = particles.get_particle(*index);
            let quad = QuadTree::get_quad(center.view(), particle.point());
            let child = node.children[quad].get_or_insert(Box::new(QuadTreeNode {
                bounding_box: QuadTree::get_quad_bounding_box(
                    center.view(),
                    quad,
                    &node.bounding_box,
                ),
                particle_indices: Vec::new(),
                children: [None, None, None, None],
                center_of_mass: Point { x: 0.0, y: 0.0 },
                mass: 0.0,
            }));
            child.particle_indices.push(*index);
            child.mass += particle.mass();
            child.center_of_mass.x += particle.point().x * particle.mass();
            child.center_of_mass.y += particle.point().y * particle.mass();
        }

        // clear the particle indices of the current node
        node.particle_indices.clear();

        for child in &mut node.children {
            if let Some(child_node) = child {
                child_node.center_of_mass.x /= child_node.mass;
                child_node.center_of_mass.y /= child_node.mass;

                QuadTree::split(
                    child_node,
                    particles,
                    depth + 1,
                    max_depth,
                    max_particles_per_node,
                );
            }
        }
    }
}

#[wasm_bindgen]
pub struct BarnesHut {
    particles: ParticleBuffer,
    theta: f32,
    max_depth: usize,
    max_particles_per_node: usize,
    gravitational_constant: f32,
}

#[wasm_bindgen]
impl BarnesHut {
    pub fn new(
        num_particles: usize,
        theta: f32,
        max_depth: usize,
        max_particles_per_node: usize,
        gravitational_constant: f32,
    ) -> Self {
        BarnesHut {
            particles: ParticleBuffer::new(num_particles),
            theta,
            max_depth,
            max_particles_per_node,
            gravitational_constant,
        }
    }

    pub fn set_gravitational_constant(&mut self, value: f32) {
        self.gravitational_constant = value;
    }

    pub fn set_theta(&mut self, value: f32) {
        self.theta = value;
    }

    pub fn set_max_particles_per_node(&mut self, value: usize) {
        self.max_particles_per_node = value;
    }

    pub fn reinitialize(
        &mut self,
        num_particles: usize,
        theta: f32,
        max_depth: usize,
        max_particles_per_node: usize,
        gravitational_constant: f32,
    ) {
        self.theta = theta;
        self.max_depth = max_depth;
        self.max_particles_per_node = max_particles_per_node;
        self.gravitational_constant = gravitational_constant;

        self.particles = ParticleBuffer::new(num_particles);
    }

    pub fn step(&mut self, time_step: f32) {
        // TODO (sanjay) : there are too many dynamic allocations occuring in each time step; Create
        //  * Create a DFS iterator that reuses a stack and takes a hint to the max iteration depth.
        //  * Similarly, create a reusable buffer to contain the leaf nodes at each step.
        //  * There are two DFS passes occuring here; one to find the leaf nodes and another to compute
        //    the forces. Try to combine these into a single pass.
        //  * QuadTree should be able to reuse memory between iterations.

        let quad_tree = QuadTree::new(&self.particles, self.max_depth, self.max_particles_per_node);

        // find all the leaf nodes
        let leaf_nodes = self.find_leaf_nodes(&quad_tree.root, self.particles.x.len());
        let mut net_forces = vec![Point { x: 0.0, y: 0.0 }; leaf_nodes.len()];

        // for each leaf node compute the force
        for (i, node) in leaf_nodes.iter().enumerate() {
            // TODO (sanjay) : make efficient
            // At the moment, this is still an O(n^2) algorithm; in order to reduce the complexity
            // we need to compute the reciprocal forces on each body or internal node and then
            // apply the net forces to bodies within the internal node. This way we do not need
            // to revisit leaf nodes.

            self.compute_forces(&mut net_forces[i], &node, &quad_tree.root);

            for p in node.particle_indices.iter() {
                self.particles.apply_force(*p, &net_forces[i]);
            }
        }

        // update positions
        for i in 0..self.particles.x.len() {
            self.particles.velocity_x[i] +=
                self.particles.force_x[i] / self.particles.mass[i] * time_step;
            self.particles.velocity_y[i] +=
                self.particles.force_y[i] / self.particles.mass[i] * time_step;
            self.particles.x[i] += self.particles.velocity_x[i] * time_step
                + 0.5 * self.particles.force_x[i] / self.particles.mass[i] * time_step.powi(2);
            self.particles.y[i] += self.particles.velocity_y[i] * time_step
                + 0.5 * self.particles.force_y[i] / self.particles.mass[i] * time_step.powi(2);

            self.particles.force_x[i] = 0.0;
            self.particles.force_y[i] = 0.0;
        }
    }

    pub fn get_raw_x_coordinates(&self) -> *const f32 {
        self.particles.get_raw_x_coordinates()
    }
    pub fn get_raw_y_coordinates(&self) -> *const f32 {
        self.particles.get_raw_y_coordinates()
    }

    pub fn get_particle(&self, index: usize) -> Point {
        Point {
            x: *self.particles.get_particle(index).x,
            y: *self.particles.get_particle(index).y,
        }
    }

    pub fn init_cube(&mut self, num_particles: usize) {
        let mut rng = rand::thread_rng();

        for i in 0..num_particles {
            self.particles
                .add_particle(rng.gen::<f32>() * 25.0, rng.gen::<f32>() * 25.0, 5.0);
        }
    }

    fn find_leaf_nodes<'a>(
        &self,
        root: &'a QuadTreeNode,
        particle_size_hint: usize,
    ) -> Vec<&'a QuadTreeNode> {
        let mut leaf_nodes = Vec::with_capacity(particle_size_hint);
        let mut stack = vec![root];
        while let Some(current_node) = stack.pop() {
            if !current_node.particle_indices.is_empty() {
                leaf_nodes.push(current_node);
            } else {
                for child in &current_node.children {
                    if let Some(child_node) = child {
                        stack.push(child_node);
                    }
                }
            }
        }
        leaf_nodes
    }

    fn compute_forces(&self, force: &mut Point, node: &QuadTreeNode, root: &QuadTreeNode) {
        let mut stack = vec![root];

        while let Some(current_node) = stack.pop() {
            let distance = node.center_of_mass.distance(&current_node.center_of_mass);

            if distance < 1e-3 {
                continue;
            }

            let width = current_node.bounding_box.width();
            let quotient = width / distance;
            if current_node.particle_indices.len() > 0 || quotient < self.theta {
                // the particle is far enough away, so we can treat the node as a single particle
                // TODO (sanjay) : are the forces reciprocal
                // TODO (sanjay) : compute the forces
                let force_magnitude =
                    self.gravitational_constant * node.mass * current_node.mass / distance.powi(2);

                // TODO (sanjay): double check this math
                force.x -= force_magnitude
                    * (node.center_of_mass.x - current_node.center_of_mass.x)
                    / distance;
                force.y -= force_magnitude
                    * (node.center_of_mass.y - current_node.center_of_mass.y)
                    / distance;
                continue;
            }

            for child in &current_node.children {
                if let Some(child_node) = child {
                    stack.push(child_node);
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Basic tests for the ParticleBuffer struct; testing insertion and access.
    #[test]
    fn test_new() {
        let buffer = ParticleBuffer::new(10);
        assert_eq!(buffer.x.capacity(), 10);
        assert_eq!(buffer.y.capacity(), 10);
        assert_eq!(buffer.mass.capacity(), 10);
    }

    #[test]
    fn test_add_particle() {
        let mut buffer = ParticleBuffer::new(10);
        buffer.add_particle(1.0, 2.0, 3.0);
        assert_eq!(buffer.x.len(), 1);
        assert_eq!(buffer.y.len(), 1);
        assert_eq!(buffer.mass.len(), 1);
        assert_eq!(buffer.x[0], 1.0);
        assert_eq!(buffer.y[0], 2.0);
        assert_eq!(buffer.mass[0], 3.0);
    }

    #[test]
    fn test_get_point() {
        let mut buffer = ParticleBuffer::new(10);
        buffer.add_particle(1.0, 2.0, 3.0);
        let point = buffer.get_position(0);
        assert_eq!(*point.x, 1.0);
        assert_eq!(*point.y, 2.0);
    }

    #[test]
    fn test_get_particle() {
        let mut buffer = ParticleBuffer::new(10);
        buffer.add_particle(1.0, 2.0, 3.0);
        let particle = buffer.get_particle(0);
        assert_eq!(*particle.x, 1.0);
        assert_eq!(*particle.y, 2.0);
        assert_eq!(particle.mass(), 3.0);
    }

    #[test]
    fn test_bounding_box() {
        let mut buffer = ParticleBuffer::new(10);
        buffer.add_particle(1.0, 2.0, 1.0);
        buffer.add_particle(3.0, 4.0, 1.0);
        buffer.add_particle(-1.0, 0.0, 1.0);
        let bounding_box = buffer.bounding_box();
        assert_eq!(bounding_box.top_right.x, 3.0);
        assert_eq!(bounding_box.top_right.y, 4.0);
        assert_eq!(bounding_box.bottom_left.x, -1.0);
        assert_eq!(bounding_box.bottom_left.y, 0.0);
    }

    #[test]
    fn test_get_quad() {
        let center = Point { x: 0.0, y: 0.0 };

        let point = Point { x: 0.75, y: 0.75 };
        assert_eq!(QuadTree::get_quad(center.view(), point.view()), 3);

        let point = Point { x: -0.75, y: 0.75 };
        assert_eq!(QuadTree::get_quad(center.view(), point.view()), 2);

        let point = Point { x: 0.75, y: -0.75 };
        assert_eq!(QuadTree::get_quad(center.view(), point.view()), 1);

        let point = Point { x: -0.75, y: -0.75 };
        assert_eq!(QuadTree::get_quad(center.view(), point.view()), 0);
    }

    #[test]
    fn test_get_quad_bounding_box() {
        let center = Point { x: 0.0, y: 0.0 };
        let bounding_box = BoundingBox {
            top_right: Point { x: 1.0, y: 1.0 },
            bottom_left: Point { x: -1.0, y: -1.0 },
        };

        let quad_bounding_box = QuadTree::get_quad_bounding_box(center.view(), 3, &bounding_box);
        assert_eq!(quad_bounding_box.top_right.x, 1.0);
        assert_eq!(quad_bounding_box.top_right.y, 1.0);
        assert_eq!(quad_bounding_box.bottom_left.x, 0.0);
        assert_eq!(quad_bounding_box.bottom_left.y, 0.0);

        let quad_bounding_box = QuadTree::get_quad_bounding_box(center.view(), 2, &bounding_box);
        assert_eq!(quad_bounding_box.top_right.x, 0.0);
        assert_eq!(quad_bounding_box.top_right.y, 1.0);
        assert_eq!(quad_bounding_box.bottom_left.x, -1.0);
        assert_eq!(quad_bounding_box.bottom_left.y, 0.0);

        let quad_bounding_box = QuadTree::get_quad_bounding_box(center.view(), 1, &bounding_box);
        assert_eq!(quad_bounding_box.top_right.x, 1.0);
        assert_eq!(quad_bounding_box.top_right.y, 0.0);
        assert_eq!(quad_bounding_box.bottom_left.x, 0.0);
        assert_eq!(quad_bounding_box.bottom_left.y, -1.0);

        let quad_bounding_box = QuadTree::get_quad_bounding_box(center.view(), 0, &bounding_box);
        assert_eq!(quad_bounding_box.top_right.x, 0.0);
        assert_eq!(quad_bounding_box.top_right.y, 0.0);
        assert_eq!(quad_bounding_box.bottom_left.x, -1.0);
        assert_eq!(quad_bounding_box.bottom_left.y, -1.0);
    }

    #[test]
    fn test_quad_tree_split() {
        let mut buffer = ParticleBuffer::new(10);
        buffer.add_particle(0.75, 0.75, 1.0);
        buffer.add_particle(-0.75, 0.75, 1.0);
        buffer.add_particle(0.75, -0.75, 1.0);
        buffer.add_particle(-0.75, -0.75, 1.0);

        let quad_tree = QuadTree::new(&buffer, 4, 1);
        assert_eq!(quad_tree.root.particle_indices.len(), 0);
        assert_eq!(quad_tree.root.mass, 4.0);
        assert_eq!(quad_tree.root.center_of_mass.x, 0.0);
        assert_eq!(quad_tree.root.center_of_mass.y, 0.0);
        assert_eq!(quad_tree.root.bounding_box.top_right.x, 0.75);
        assert_eq!(quad_tree.root.bounding_box.top_right.y, 0.75);
        assert_eq!(quad_tree.root.bounding_box.bottom_left.x, -0.75);
        assert_eq!(quad_tree.root.bounding_box.bottom_left.y, -0.75);
        assert_eq!(quad_tree.root.children[0].is_some(), true);
        assert_eq!(quad_tree.root.children[1].is_some(), true);
        assert_eq!(quad_tree.root.children[2].is_some(), true);
        assert_eq!(quad_tree.root.children[3].is_some(), true);

        let quadrant0 = quad_tree.root.children[0].as_ref().unwrap();
        let quadrant1 = quad_tree.root.children[1].as_ref().unwrap();
        let quadrant2 = quad_tree.root.children[2].as_ref().unwrap();
        let quadrant3 = quad_tree.root.children[3].as_ref().unwrap();

        assert_eq!(quadrant0.particle_indices.len(), 1);
        assert_eq!(quadrant1.particle_indices.len(), 1);
        assert_eq!(quadrant2.particle_indices.len(), 1);
        assert_eq!(quadrant3.particle_indices.len(), 1);

        assert_eq!(quad_tree.root.mass, 4.0);
        assert_eq!(quad_tree.root.center_of_mass.x, 0.0);
        assert_eq!(quad_tree.root.center_of_mass.y, 0.0);
        assert_eq!(quad_tree.root.bounding_box.top_right.x, 0.75);
        assert_eq!(quad_tree.root.bounding_box.top_right.y, 0.75);
        assert_eq!(quad_tree.root.bounding_box.bottom_left.x, -0.75);
        assert_eq!(quad_tree.root.bounding_box.bottom_left.y, -0.75);

        assert_eq!(quadrant0.particle_indices[0], 3);
        assert_eq!(quadrant1.particle_indices[0], 2);
        assert_eq!(quadrant2.particle_indices[0], 1);
        assert_eq!(quadrant3.particle_indices[0], 0);
    }

    #[test]
    fn test_quad_tree_split_with_mass() {
        let mut buffer = ParticleBuffer::new(10);
        buffer.add_particle(0.75, 0.75, 1.0);
        buffer.add_particle(-0.75, 0.75, 1.0);
        buffer.add_particle(0.75, -0.75, 1.0);
        buffer.add_particle(-0.75, -0.75, 1.0);
        buffer.add_particle(0.1, 0.3, 1.0);
        buffer.add_particle(0.50, 0.5, 1.0);
        buffer.add_particle(0.50, 0.25, 1.0);

        let quad_tree = QuadTree::new(&buffer, 1, 1);

        assert_eq!(quad_tree.root.mass, 7.0);
        assert_eq!(quad_tree.root.bounding_box.top_right.x, 0.75);
        assert_eq!(quad_tree.root.bounding_box.top_right.y, 0.75);
        assert_eq!(quad_tree.root.bounding_box.bottom_left.x, -0.75);
        assert_eq!(quad_tree.root.bounding_box.bottom_left.y, -0.75);

        assert_eq!(quad_tree.root.particle_indices.len(), 0);
        assert_eq!(quad_tree.root.children[0].is_some(), true);
        assert_eq!(quad_tree.root.children[1].is_some(), true);
        assert_eq!(quad_tree.root.children[2].is_some(), true);
        assert_eq!(quad_tree.root.children[3].is_some(), true);

        let quadrant0 = quad_tree.root.children[0].as_ref().unwrap();
        let quadrant1 = quad_tree.root.children[1].as_ref().unwrap();
        let quadrant2 = quad_tree.root.children[2].as_ref().unwrap();
        let quadrant3 = quad_tree.root.children[3].as_ref().unwrap();

        assert_eq!(quadrant0.particle_indices.len(), 1);
        assert_eq!(quadrant1.particle_indices.len(), 1);
        assert_eq!(quadrant2.particle_indices.len(), 1);
        assert_eq!(quadrant3.particle_indices.len(), 4);

        assert_eq!(quadrant3.mass, 4.0);
        assert_eq!(quadrant3.bounding_box.top_right.x, 0.75);
        assert_eq!(quadrant3.bounding_box.top_right.y, 0.75);
        assert_eq!(quadrant3.bounding_box.bottom_left.x, 0.0);
        assert_eq!(quadrant3.bounding_box.bottom_left.y, 0.0);
        assert_eq!(quadrant3.center_of_mass.x, 0.4625);
        assert_eq!(quadrant3.center_of_mass.y, 0.45);

        let expected_indices = vec![0, 6, 5, 4];
        assert!(quadrant3
            .particle_indices
            .iter()
            .all(|item| expected_indices.contains(item)));
    }

    #[test]
    fn barnes_hut_basic_test() {
        let mut buffer = ParticleBuffer::new(10);
        buffer.add_particle(0.75, 0.75, 1.0);
        buffer.add_particle(-0.75, 0.75, 1.0);
        buffer.add_particle(0.75, -0.75, 1.0);
        buffer.add_particle(-0.75, -0.75, 1.0);

        let mut barnes_hut = BarnesHut::new(buffer, 1.0, 1, 1);
        barnes_hut.step(1.0);

        assert!(true);
        // TODO (sanjay) : write tests
    }
}
