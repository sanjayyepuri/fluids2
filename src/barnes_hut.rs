use rand::Rng;
use std::vec;
use wasm_bindgen::prelude::*;


/// @brief Point is a struct that holds the x and y coordinates of a point.
#[derive(Clone, Debug)]
#[wasm_bindgen]
pub struct Point {
    x: f32,
    y: f32,
    z: f32,
}

impl Point {
    fn distance(&self, other: &Point) -> f32 {
        ((self.x - other.x).powi(2) + (self.y - other.y).powi(2) + (self.z - other.z).powi(2)).sqrt()
    }
}

/// @brief PointView is a struct that holds references to the x, y and z
/// coordinates of a point from the particle buffer.
struct PointView<'a> {
    x: &'a f32,
    y: &'a f32,
    z: &'a f32,
}

struct ParticleView<'a> {
    x: &'a f32,
    y: &'a f32,
    z: &'a f32,
    mass: &'a f32,
}

impl ParticleView<'_> {
    fn point(&self) -> PointView {
        PointView {
            x: self.x,
            y: self.y,
            z: self.z,
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
            z: &self.z,
        }
    }
}

/// @brief BoundingBox is a struct that holds the min and max points of a
/// bounding box. This is used to determine the bounds of a quadtree node.
struct BoundingBox {
    min: Point,
    max: Point,
}

impl BoundingBox {
    fn width(&self) -> f32 {
        self.max.x - self.min.x
    }

    fn height(&self) -> f32 {
        self.max.y - self.min.y
    }

    fn depth(&self) -> f32 {
        self.max.z - self.min.z
    }
}

/// @brief ParticleBuffer is a struct that holds the x, y and z coordinates of particles.
/// The coordinates are stored in separate vectors in order to optimize iteration.
struct ParticleBuffer {
    x: Vec<f32>,
    y: Vec<f32>,
    z: Vec<f32>,

    velocity_x: Vec<f32>,
    velocity_y: Vec<f32>,
    velocity_z: Vec<f32>,

    force_x: Vec<f32>,
    force_y: Vec<f32>,
    force_z: Vec<f32>,

    mass: Vec<f32>,
}

impl ParticleBuffer {
    fn new(capacity: usize) -> Self {
        ParticleBuffer {
            x: Vec::with_capacity(capacity),
            y: Vec::with_capacity(capacity),
            z: Vec::with_capacity(capacity),
            velocity_x: Vec::with_capacity(capacity),
            velocity_y: Vec::with_capacity(capacity),
            velocity_z: Vec::with_capacity(capacity),
            force_x: Vec::with_capacity(capacity),
            force_y: Vec::with_capacity(capacity),
            force_z: Vec::with_capacity(capacity),
            mass: Vec::with_capacity(capacity),
        }
    }

    /// @brief Adds a particle to the buffer. In theory particles only should added at
    /// initialization time. Consider pushing this method to the constructor.
    fn add_particle(&mut self, x: f32, y: f32, z: f32, mass: f32) {
        self.x.push(x);
        self.y.push(y);
        self.z.push(z);
        self.velocity_x.push(0.0);
        self.velocity_y.push(0.0);
        self.velocity_z.push(0.0);
        self.force_x.push(0.0);
        self.force_y.push(0.0);
        self.force_z.push(0.0);
        self.mass.push(mass);
    }

    /// @brief Returns a PointView for the particle at the given index.
    fn get_position(&self, index: usize) -> PointView {
        PointView {
            x: &self.x[index],
            y: &self.y[index],
            z: &self.z[index],
        }
    }

    /// @brief Sets the position of the particle at the given index.
    fn set_position(&mut self, index: usize, x: f32, y: f32, z: f32) {
        self.x[index] = x;
        self.y[index] = y;
        self.z[index] = z;
    }

    /// @brief Returns a ParticleView for the particle at the given index.
    fn get_particle(&self, index: usize) -> ParticleView {
        ParticleView {
            x: &self.x[index],
            y: &self.y[index],
            z: &self.z[index],
            mass: &self.mass[index],
        }
    }

    fn apply_force(&mut self, index: usize, force: &Point) {
        self.force_x[index] += force.x;
        self.force_y[index] += force.y;
        self.force_z[index] += force.z;
    }

    /// @brief Computes the min and max of the x, y and z coordinates in the particle buffer.
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
        let (min_z, max_z) = self
            .z
            .iter()
            .fold((f32::INFINITY, f32::NEG_INFINITY), |(min, max), &val| {
                (min.min(val), max.max(val))
            });

        return BoundingBox {
            max: Point { x: max_x, y: max_y, z: max_z },
            min: Point { x: min_x, y: min_y, z: min_z },
        };
    }

    fn get_raw_x_coordinates(&self) -> *const f32 {
        self.x.as_ptr()
    }

    fn get_raw_y_coordinates(&self) -> *const f32 {
        self.y.as_ptr()
    }

    fn get_raw_z_coordinates(&self) -> *const f32 {
        self.z.as_ptr()
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
    children: [Option<Box<QuadTreeNode>>; 8],
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
            z: particles
                .z
                .iter()
                .zip(particles.mass.iter())
                .map(|(&z, &m)| z * m)
                .sum::<f32>()
                / mass,
        };

        let mut tree = QuadTree {
            root: QuadTreeNode {
                bounding_box: particles.bounding_box(),
                particle_indices: (0..particles.x.len()).collect(),
                center_of_mass,
                mass,
                children: [None, None, None, None, None, None, None, None],
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
        // Encode the octant as a 3-bit number.
        // Bit 0: x > center.x
        // Bit 1: y > center.y
        // Bit 2: z > center.z
        ((*point.x > *center.x) as usize)
            | (((*point.y > *center.y) as usize) << 1)
            | (((*point.z > *center.z) as usize) << 2)
    }

    fn get_octant_bounding_box(
        center: PointView,
        octant: usize,
        bounding_box: &BoundingBox,
    ) -> BoundingBox {
        let min = Point {
            x: if octant & 0b001 == 0 {
                bounding_box.min.x
            } else {
                *center.x
            },
            y: if octant & 0b010 == 0 {
                bounding_box.min.y
            } else {
                *center.y
            },
            z: if octant & 0b100 == 0 {
                bounding_box.min.z
            } else {
                *center.z
            },
        };

        let max = Point {
            x: if octant & 0b001 == 0 {
                *center.x
            } else {
                bounding_box.max.x
            },
            y: if octant & 0b010 == 0 {
                *center.y
            } else {
                bounding_box.max.y
            },
            z: if octant & 0b100 == 0 {
                *center.z
            } else {
                bounding_box.max.z
            },
        };
        return BoundingBox {
            min,
            max,
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
            x: (node.bounding_box.max.x + node.bounding_box.min.x) / 2.0,
            y: (node.bounding_box.max.y + node.bounding_box.min.y) / 2.0,
            z: (node.bounding_box.max.z + node.bounding_box.min.z) / 2.0,
        };

        for index in &node.particle_indices {
            let particle = particles.get_particle(*index);
            let octant = QuadTree::get_quad(center.view(), particle.point());
            let child = node.children[octant].get_or_insert(Box::new(QuadTreeNode {
                bounding_box: QuadTree::get_octant_bounding_box(
                    center.view(),
                    octant,
                    &node.bounding_box,
                ),
                particle_indices: Vec::new(),
                children: [None, None, None, None, None, None, None, None],
                center_of_mass: Point { x: 0.0, y: 0.0, z: 0.0 },
                mass: 0.0,
            }));
            child.particle_indices.push(*index);
            child.mass += particle.mass();
            child.center_of_mass.x += particle.point().x * particle.mass();
            child.center_of_mass.y += particle.point().y * particle.mass();
            child.center_of_mass.z += particle.point().z * particle.mass();
        }

        // clear the particle indices of the current node
        node.particle_indices.clear();

        for child in &mut node.children {
            if let Some(child_node) = child {
                child_node.center_of_mass.x /= child_node.mass;
                child_node.center_of_mass.y /= child_node.mass;
                child_node.center_of_mass.z /= child_node.mass;

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
        let mut net_forces = vec![Point { x: 0.0, y: 0.0, z: 0.0 }; leaf_nodes.len()];

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
            self.particles.velocity_z[i] +=
                self.particles.force_z[i] / self.particles.mass[i] * time_step;
            self.particles.x[i] += self.particles.velocity_x[i] * time_step
                + 0.5 * self.particles.force_x[i] / self.particles.mass[i] * time_step.powi(2);
            self.particles.y[i] += self.particles.velocity_y[i] * time_step
                + 0.5 * self.particles.force_y[i] / self.particles.mass[i] * time_step.powi(2);
            self.particles.z[i] += self.particles.velocity_z[i] * time_step
                + 0.5 * self.particles.force_z[i] / self.particles.mass[i] * time_step.powi(2);

            self.particles.force_x[i] = 0.0;
            self.particles.force_y[i] = 0.0;
            self.particles.force_z[i] = 0.0;
        }
    }

    pub fn get_raw_x_coordinates(&self) -> *const f32 {
        self.particles.get_raw_x_coordinates()
    }
    pub fn get_raw_y_coordinates(&self) -> *const f32 {
        self.particles.get_raw_y_coordinates()
    }
    #[wasm_bindgen]
    pub fn get_raw_z_coordinates(&self) -> *const f32 {
        self.particles.get_raw_z_coordinates()
    }

    pub fn get_particle(&self, index: usize) -> Point {
        Point {
            x: *self.particles.get_particle(index).x,
            y: *self.particles.get_particle(index).y,
            z: *self.particles.get_particle(index).z,
        }
    }

    #[wasm_bindgen]
    pub fn get_system_energy(&self) -> f32 {
        let mut kinetic_energy = 0.0;
        for i in 0..self.particles.size() {
            let vx = self.particles.velocity_x[i];
            let vy = self.particles.velocity_y[i];
            let vz = self.particles.velocity_z[i];
            let speed_sq = vx.powi(2) + vy.powi(2) + vz.powi(2);
            kinetic_energy += 0.5 * self.particles.mass[i] * speed_sq;
        }

        let mut potential_energy = 0.0;
        let num_particles = self.particles.size();
        for i in 0..num_particles {
            let p1 = Point {
                x: self.particles.x[i],
                y: self.particles.y[i],
                z: self.particles.z[i],
            };
            for j in (i + 1)..num_particles {
                let p2 = Point {
                    x: self.particles.x[j],
                    y: self.particles.y[j],
                    z: self.particles.z[j],
                };
                let distance = p1.distance(&p2);
                // Add a small epsilon to avoid division by zero or extremely large energies.
                let effective_distance = distance.max(1e-3);
                potential_energy -= self.gravitational_constant * self.particles.mass[i] * self.particles.mass[j] / effective_distance;
            }
        }
        kinetic_energy + potential_energy
    }

    pub fn init_cube(&mut self, num_particles: usize) {
        let mut rng = rand::thread_rng();

        for i in 0..num_particles {
            self.particles
                .add_particle(rng.gen::<f32>() * 25.0, rng.gen::<f32>() * 25.0, rng.gen::<f32>() * 25.0, 5.0);
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

            let bb = &current_node.bounding_box;
            let max_dim = bb.width().max(bb.height()).max(bb.depth());
            let quotient = max_dim / distance;

            if current_node.particle_indices.len() > 0 || quotient < self.theta {
                // The node is far enough away, or it's a leaf node.
                // Treat the current_node as a single particle.
                let force_magnitude =
                    self.gravitational_constant * node.mass * current_node.mass / distance.powi(2);

                force.x -= force_magnitude
                    * (node.center_of_mass.x - current_node.center_of_mass.x)
                    / distance;
                force.y -= force_magnitude
                    * (node.center_of_mass.y - current_node.center_of_mass.y)
                    / distance;
                force.z -= force_magnitude
                    * (node.center_of_mass.z - current_node.center_of_mass.z)
                    / distance;
                continue;
            }
            // Node is too close, recursively explore its children.

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
        assert_eq!(buffer.z.capacity(), 10);
        assert_eq!(buffer.mass.capacity(), 10);
    }

    #[test]
    fn test_point_distance_3d() {
        let p1 = Point { x: 1.0, y: 2.0, z: 3.0 };
        let p2 = Point { x: 4.0, y: 6.0, z: 8.0 }; // dx=3, dy=4, dz=5
        assert_eq!(p1.distance(&p2), ((3.0_f32.powi(2) + 4.0_f32.powi(2) + 5.0_f32.powi(2)).sqrt())); // sqrt(9+16+25) = sqrt(50)

        let p3 = Point { x: -1.0, y: -2.0, z: -3.0 };
        let p4 = Point { x: -4.0, y: -6.0, z: -8.0 }; // dx=-3, dy=-4, dz=-5
        assert_eq!(p3.distance(&p4), (( (-3.0_f32).powi(2) + (-4.0_f32).powi(2) + (-5.0_f32).powi(2)).sqrt()));

        let p5 = Point { x: 0.0, y: 0.0, z: 0.0 };
        let p6 = Point { x: 1.0, y: 1.0, z: 1.0 };
        assert_eq!(p5.distance(&p6), (3.0_f32).sqrt());

        let p7 = Point { x: 1.0, y: 1.0, z: 1.0 };
        assert_eq!(p7.distance(&p7), 0.0);
    }

    #[test]
    fn test_bounding_box_dimensions_3d() {
        let bb = BoundingBox {
            min: Point { x: -1.0, y: -2.0, z: -3.0 },
            max: Point { x: 1.0, y: 2.0, z: 3.0 },
        };
        assert_eq!(bb.width(), 2.0);
        assert_eq!(bb.height(), 4.0);
        assert_eq!(bb.depth(), 6.0);

        let bb_zero = BoundingBox {
            min: Point { x: 0.0, y: 0.0, z: 0.0 },
            max: Point { x: 0.0, y: 0.0, z: 0.0 },
        };
        assert_eq!(bb_zero.width(), 0.0);
        assert_eq!(bb_zero.height(), 0.0);
        assert_eq!(bb_zero.depth(), 0.0);
    }

    #[test]
    fn test_add_particle() {
        let mut buffer = ParticleBuffer::new(10);
        buffer.add_particle(1.0, 2.0, 3.0, 4.0);
        assert_eq!(buffer.x.len(), 1);
        assert_eq!(buffer.y.len(), 1);
        assert_eq!(buffer.z.len(), 1);
        assert_eq!(buffer.mass.len(), 1);
        assert_eq!(buffer.x[0], 1.0);
        assert_eq!(buffer.y[0], 2.0);
        assert_eq!(buffer.z[0], 3.0);
        assert_eq!(buffer.mass[0], 4.0);
    }

    #[test]
    fn test_get_point() {
        let mut buffer = ParticleBuffer::new(10);
        buffer.add_particle(1.0, 2.0, 3.0, 4.0);
        let point = buffer.get_position(0);
        assert_eq!(*point.x, 1.0);
        assert_eq!(*point.y, 2.0);
        assert_eq!(*point.z, 3.0);
    }

    #[test]
    fn test_get_particle() {
        let mut buffer = ParticleBuffer::new(10);
        buffer.add_particle(1.0, 2.0, 3.0, 4.0);
        let particle = buffer.get_particle(0);
        assert_eq!(*particle.x, 1.0);
        assert_eq!(*particle.y, 2.0);
        assert_eq!(*particle.z, 3.0);
        assert_eq!(particle.mass(), 4.0);
    }

    #[test]
    fn test_bounding_box() {
        let mut buffer = ParticleBuffer::new(10);
        buffer.add_particle(1.0, 2.0, 3.0, 1.0);
        buffer.add_particle(3.0, 4.0, 5.0, 1.0);
        buffer.add_particle(-1.0, 0.0, -2.0, 1.0);
        let bounding_box = buffer.bounding_box();
        assert_eq!(bounding_box.max.x, 3.0);
        assert_eq!(bounding_box.max.y, 4.0);
        assert_eq!(bounding_box.max.z, 5.0);
        assert_eq!(bounding_box.min.x, -1.0);
        assert_eq!(bounding_box.min.y, 0.0);
        assert_eq!(bounding_box.min.z, -2.0);
    }

    #[test]
    fn test_get_octant() {
        let center = Point { x: 0.0, y: 0.0, z: 0.0 };

        // Test all 8 octants
        let points = [
            Point { x: -0.75, y: -0.75, z: -0.75 }, // Octant 0 (---)
            Point { x: 0.75, y: -0.75, z: -0.75 },  // Octant 1 (+--)
            Point { x: -0.75, y: 0.75, z: -0.75 },  // Octant 2 (-+-)
            Point { x: 0.75, y: 0.75, z: -0.75 },   // Octant 3 (++-)
            Point { x: -0.75, y: -0.75, z: 0.75 },  // Octant 4 (--+)
            Point { x: 0.75, y: -0.75, z: 0.75 },   // Octant 5 (+-+)
            Point { x: -0.75, y: 0.75, z: 0.75 },   // Octant 6 (-++)
            Point { x: 0.75, y: 0.75, z: 0.75 },    // Octant 7 (+++)
        ];

        for (i, point) in points.iter().enumerate() {
            assert_eq!(QuadTree::get_quad(center.view(), point.view()), i);
        }
    }

    #[test]
    fn test_get_octant_bounding_box() {
        let center = Point { x: 0.0, y: 0.0, z: 0.0 };
        let bounding_box = BoundingBox {
            max: Point { x: 1.0, y: 1.0, z: 1.0 },
            min: Point { x: -1.0, y: -1.0, z: -1.0 },
        };

        // Octant 0 (---)
        let octant_bb = QuadTree::get_octant_bounding_box(center.view(), 0, &bounding_box);
        assert_eq!(octant_bb.min.x, -1.0);
        assert_eq!(octant_bb.min.y, -1.0);
        assert_eq!(octant_bb.min.z, -1.0);
        assert_eq!(octant_bb.max.x, 0.0);
        assert_eq!(octant_bb.max.y, 0.0);
        assert_eq!(octant_bb.max.z, 0.0);

        // Octant 1 (+--)
        let octant_bb = QuadTree::get_octant_bounding_box(center.view(), 1, &bounding_box);
        assert_eq!(octant_bb.min.x, 0.0);
        assert_eq!(octant_bb.min.y, -1.0);
        assert_eq!(octant_bb.min.z, -1.0);
        assert_eq!(octant_bb.max.x, 1.0);
        assert_eq!(octant_bb.max.y, 0.0);
        assert_eq!(octant_bb.max.z, 0.0);
        
        // Octant 6 (-++)
        let octant_bb = QuadTree::get_octant_bounding_box(center.view(), 6, &bounding_box);
        assert_eq!(octant_bb.min.x, -1.0);
        assert_eq!(octant_bb.min.y, 0.0);
        assert_eq!(octant_bb.min.z, 0.0);
        assert_eq!(octant_bb.max.x, 0.0);
        assert_eq!(octant_bb.max.y, 1.0);
        assert_eq!(octant_bb.max.z, 1.0);

        // Octant 7 (+++)
        let octant_bb = QuadTree::get_octant_bounding_box(center.view(), 7, &bounding_box);
        assert_eq!(octant_bb.min.x, 0.0);
        assert_eq!(octant_bb.min.y, 0.0);
        assert_eq!(octant_bb.min.z, 0.0);
        assert_eq!(octant_bb.max.x, 1.0);
        assert_eq!(octant_bb.max.y, 1.0);
        assert_eq!(octant_bb.max.z, 1.0);
    }

    #[test]
    fn test_octree_split() {
        let mut buffer = ParticleBuffer::new(10);
        // Add 8 particles, one for each octant
        buffer.add_particle(0.75, 0.75, 0.75, 1.0);    // Octant 7
        buffer.add_particle(-0.75, 0.75, 0.75, 1.0);   // Octant 6
        buffer.add_particle(0.75, -0.75, 0.75, 1.0);   // Octant 5
        buffer.add_particle(-0.75, -0.75, 0.75, 1.0);  // Octant 4
        buffer.add_particle(0.75, 0.75, -0.75, 1.0);   // Octant 3
        buffer.add_particle(-0.75, 0.75, -0.75, 1.0);  // Octant 2
        buffer.add_particle(0.75, -0.75, -0.75, 1.0);  // Octant 1
        buffer.add_particle(-0.75, -0.75, -0.75, 1.0); // Octant 0


        let octree = QuadTree::new(&buffer, 4, 1); // Max depth 4, max 1 particle per node
        assert_eq!(octree.root.particle_indices.len(), 0); // Root node should be split
        assert_eq!(octree.root.mass, 8.0); // Total mass
        assert_eq!(octree.root.center_of_mass.x, 0.0);
        assert_eq!(octree.root.center_of_mass.y, 0.0);
        assert_eq!(octree.root.center_of_mass.z, 0.0); // Adjusted for symmetry

        for i in 0..8 {
            assert!(octree.root.children[i].is_some(), "Child {} should exist", i);
            let child_node = octree.root.children[i].as_ref().unwrap();
            assert_eq!(child_node.particle_indices.len(), 1, "Child {} should have 1 particle", i);
            // The exact particle index depends on insertion order and octant calculation,
            // so we verify properties of the particle's position relative to the octant.
            let particle_view = buffer.get_particle(child_node.particle_indices[0]);
            let octant = QuadTree::get_quad(octree.root.center_of_mass.view(), particle_view.point());
            assert_eq!(octant, i, "Particle in child {} is not in the correct octant", i);
        }
    }

    #[test]
    fn test_octree_split_with_mass() {
        let mut buffer = ParticleBuffer::new(20); // Increased capacity
        // Particles in Octant 7 (+++)
        buffer.add_particle(0.75, 0.75, 0.75, 1.0);
        buffer.add_particle(0.50, 0.5, 0.5, 1.0);
        buffer.add_particle(0.25, 0.25, 0.25, 1.0);
        // Particle in Octant 0 (---)
        buffer.add_particle(-0.75, -0.75, -0.75, 3.0); // Heavier particle

        // Other octants to ensure splitting
        buffer.add_particle(0.1, -0.1, -0.1, 1.0); // Octant 1
        buffer.add_particle(-0.1, 0.1, -0.1, 1.0); // Octant 2
        buffer.add_particle(0.1, 0.1, -0.1, 1.0);  // Octant 3
        buffer.add_particle(-0.1, -0.1, 0.1, 1.0); // Octant 4
        buffer.add_particle(0.1, -0.1, 0.1, 1.0);  // Octant 5
        buffer.add_particle(-0.1, 0.1, 0.1, 1.0);  // Octant 6


        let octree = QuadTree::new(&buffer, 1, 1); // Max depth 1, max 1 particle per node

        assert_eq!(octree.root.mass, 10.0); // Total mass
        // Center of mass calculation would be more complex here, skip direct check for now
        assert_eq!(octree.root.particle_indices.len(), 0); // Root should split

        // Check Octant 7
        let octant7_node = octree.root.children[7].as_ref().unwrap();
        assert_eq!(octant7_node.particle_indices.len(), 3); // All 3 particles in this octant
        assert_eq!(octant7_node.mass, 3.0);
        // Verify center of mass for Octant 7
        let expected_com_x = (0.75 * 1.0 + 0.50 * 1.0 + 0.25 * 1.0) / 3.0;
        let expected_com_y = (0.75 * 1.0 + 0.50 * 1.0 + 0.25 * 1.0) / 3.0;
        let expected_com_z = (0.75 * 1.0 + 0.50 * 1.0 + 0.25 * 1.0) / 3.0;
        assert!((octant7_node.center_of_mass.x - expected_com_x).abs() < 1e-6);
        assert!((octant7_node.center_of_mass.y - expected_com_y).abs() < 1e-6);
        assert!((octant7_node.center_of_mass.z - expected_com_z).abs() < 1e-6);


        // Check Octant 0
        let octant0_node = octree.root.children[0].as_ref().unwrap();
        assert_eq!(octant0_node.particle_indices.len(), 1);
        assert_eq!(octant0_node.mass, 3.0);
        assert_eq!(octant0_node.center_of_mass.x, -0.75);
        assert_eq!(octant0_node.center_of_mass.y, -0.75);
        assert_eq!(octant0_node.center_of_mass.z, -0.75);

        let expected_indices_octant7 = vec![0, 1, 2]; // Assuming these were the first three added
        assert!(octant7_node
            .particle_indices
            .iter()
            .all(|item| expected_indices_octant7.contains(item)));
    }

    #[test]
    fn barnes_hut_basic_test() {
        let mut buffer = ParticleBuffer::new(10);
        buffer.add_particle(0.75, 0.75, 0.75, 1.0);
        buffer.add_particle(-0.75, 0.75, 0.75, 1.0);
        buffer.add_particle(0.75, -0.75, 0.75, 1.0);
        buffer.add_particle(-0.75, -0.75, 0.75, 1.0);

        // TODO(sanjay): The following fields are not available in ParticleBuffer
        // let mut barnes_hut = BarnesHut::new(buffer, 1.0, 1, 1);
        // barnes_hut.step(1.0);

        assert!(true);
        // TODO (sanjay) : write tests
    }

    #[test]
    fn test_system_energy_simple() {
        let mut sim = BarnesHut::new(3, 0.5, 5, 1, 1.0); // G = 1.0

        // Particle 1: m=1, p=(0,0,0), v=(0,0,0)
        sim.particles.add_particle(0.0, 0.0, 0.0, 1.0);
        // Velocities are initialized to 0, so no need to set them for this particle

        // Particle 2: m=2, p=(1,0,0), v=(1,0,0)
        sim.particles.add_particle(1.0, 0.0, 0.0, 2.0);
        sim.particles.velocity_x[1] = 1.0;

        // Particle 3: m=3, p=(0,1,0), v=(0,1,1)
        sim.particles.add_particle(0.0, 1.0, 0.0, 3.0);
        sim.particles.velocity_y[2] = 1.0;
        sim.particles.velocity_z[2] = 1.0;
        
        // Expected Kinetic Energy:
        // P1_ke = 0.5 * 1 * (0^2) = 0
        // P2_ke = 0.5 * 2 * (1^2 + 0^2 + 0^2) = 0.5 * 2 * 1 = 1.0
        // P3_ke = 0.5 * 3 * (0^2 + 1^2 + 1^2) = 0.5 * 3 * 2 = 3.0
        // Total KE = 0 + 1.0 + 3.0 = 4.0
        let expected_ke = 4.0;

        // Expected Potential Energy (G=1):
        // PE_12 = -1 * 1 * 2 / distance((0,0,0), (1,0,0)) = -2 / 1 = -2.0
        // PE_13 = -1 * 1 * 3 / distance((0,0,0), (0,1,0)) = -3 / 1 = -3.0
        // PE_23 = -1 * 2 * 3 / distance((1,0,0), (0,1,0)) = -6 / sqrt(2) = -6 / 1.41421356 approx
        let dist_23 = ((1.0_f32 - 0.0_f32).powi(2) + (0.0_f32 - 1.0_f32).powi(2) + (0.0_f32 - 0.0_f32).powi(2)).sqrt(); // sqrt(1+1+0) = sqrt(2)
        let expected_pe = -2.0 / 1.0 - 3.0 / 1.0 - (2.0 * 3.0) / dist_23; // -2.0 - 3.0 - 6.0/sqrt(2) = -5.0 - 4.242640687
                                                                      // = -9.242640687

        let expected_total_energy = expected_ke + expected_pe;
        let calculated_total_energy = sim.get_system_energy();
        assert!((calculated_total_energy - expected_total_energy).abs() < 1e-5, "Energy mismatch: expected {}, got {}", expected_total_energy, calculated_total_energy);

        // Test case: Zero velocity (only potential energy)
        let mut sim_zero_vel = BarnesHut::new(2, 0.5, 5, 1, 1.0);
        sim_zero_vel.particles.add_particle(0.0, 0.0, 0.0, 1.0); // m1=1, p1=(0,0,0)
        sim_zero_vel.particles.add_particle(1.0, 0.0, 0.0, 2.0); // m2=2, p2=(1,0,0)
        // Velocities are 0 by default.
        let expected_pe_zero_vel = -1.0 * 1.0 * 2.0 / 1.0; // -2.0
        let calculated_energy_zero_vel = sim_zero_vel.get_system_energy();
        assert!((calculated_energy_zero_vel - expected_pe_zero_vel).abs() < 1e-5, "Zero velocity energy: expected {}, got {}", expected_pe_zero_vel, calculated_energy_zero_vel);
    
        // Test case: Particles at the same location (distance clamping)
        let mut sim_same_loc = BarnesHut::new(2, 0.5, 5, 1, 1.0);
        sim_same_loc.particles.add_particle(0.0, 0.0, 0.0, 1.0);
        sim_same_loc.particles.add_particle(0.0, 0.0, 0.0, 2.0);
        // Expected PE with clamping at 1e-3
        let expected_pe_clamped = -1.0 * 1.0 * 2.0 / 1e-3; // -2000.0
        let calculated_energy_clamped = sim_same_loc.get_system_energy();
         assert!((calculated_energy_clamped - expected_pe_clamped).abs() < 1e-5, "Clamped energy: expected {}, got {}", expected_pe_clamped, calculated_energy_clamped);
    }
}
