/// @brief Point is a struct that holds the x and y coordinates of a point.
struct Point {
    x: f32,
    y: f32,
}

/// @brief PointView is a struct that holds references to the x and y
/// coordinates of a point from the particle buffer.
struct PointView<'a> {
    x: &'a f32,
    y: &'a f32,
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

/// @brief ParticleBuffer is a struct that holds the x and y coordinates of particles.
/// The coordinates are stored in separate vectors in order to optimize iteration.
struct ParticleBuffer {
    x: Vec<f32>,
    y: Vec<f32>,
}
impl ParticleBuffer {
    fn new(capacity: usize) -> Self {
        ParticleBuffer {
            x: Vec::with_capacity(capacity),
            y: Vec::with_capacity(capacity),
        }
    }

    /// @brief Adds a particle to the buffer. In theory particles only should added at
    /// initialization time. Consider pushing this method to the constructor.
    fn add_particle(&mut self, x: f32, y: f32) {
        self.x.push(x);
        self.y.push(y);
    }

    /// @brief Returns a PointView for the particle at the given index.
    fn get_point(&self, index: usize) -> PointView {
        PointView {
            x: &self.x[index],
            y: &self.y[index],
        }
    }

    /// @brief Computes the min and max of the x coordinates in the particle buffer.
    fn bounding_box(&self) -> BoundingBox {
        let (min_x, max_x) = self.x.iter().fold((f32::INFINITY, f32::NEG_INFINITY), |(min, max), &val| {
            (min.min(val), max.max(val))
        });
        let (min_y, max_y) = self.y.iter().fold((f32::INFINITY, f32::NEG_INFINITY), |(min, max), &val| {
            (min.min(val), max.max(val))
        });

        return BoundingBox {
            top_right: Point {
                x: max_x,
                y: max_y,
            },
            bottom_left: Point {
                x: min_x,
                y: min_y,
            },
        };
    }
}
struct QuadTreeNode {
    bounding_box: BoundingBox,
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
    fn new (particles: &ParticleBuffer, max_depth: usize, max_particles_per_node: usize) -> QuadTree {
        let mut tree = QuadTree{
            root: QuadTreeNode {
                bounding_box: particles.bounding_box(),
                particle_indices: (0..particles.x.len()).collect(),
                children: [None, None, None, None],
            },
            max_depth,
            max_particles_per_node
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
        let bottom_left = Point { x: if quad & 0b01 == 0 { bounding_box.bottom_left.x } else { *center.x },
                                 y: if quad & 0b10 == 0 { bounding_box.bottom_left.y } else { *center.y } };

        let top_right = Point { x: if quad & 0b01 == 0 { *center.x } else { bounding_box.top_right.x },
                                 y: if quad & 0b10 == 0 { *center.y } else { bounding_box.top_right.y } };
        return BoundingBox {
            top_right,
            bottom_left,
        };
    }

    fn build(&mut self, particles: &ParticleBuffer) {
        QuadTree::split(&mut self.root, particles, 0, self.max_depth, self.max_particles_per_node);
    }

    fn split(node: &mut QuadTreeNode, particles: &ParticleBuffer, depth: usize, max_depth: usize, max_particles_per_node: usize) {
        if node.particle_indices.len() <= max_particles_per_node || depth > max_depth {
            // base case: no need to split further
            return;
        }

        let center = Point {
            x: (node.bounding_box.top_right.x + node.bounding_box.bottom_left.x) / 2.0,
            y: (node.bounding_box.top_right.y + node.bounding_box.bottom_left.y) / 2.0,
        };

        for index in &node.particle_indices {
            let point = particles.get_point(*index);
            let quad = QuadTree::get_quad(center.view(), point);
            let child = node.children[quad].get_or_insert(Box::new(QuadTreeNode {
                bounding_box: QuadTree::get_quad_bounding_box(center.view(), quad, &node.bounding_box),
                particle_indices: Vec::new(),
                children: [None, None, None, None],
            }));

            child.particle_indices.push(*index);
        }

        // clear the particle indices of the current node
        node.particle_indices.clear();

        for child in &mut node.children {
            if let Some(child_node) = child {
                QuadTree::split(child_node, particles, depth+1, max_depth, max_particles_per_node);
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
    }

    #[test]
    fn test_add_particle() {
        let mut buffer = ParticleBuffer::new(10);
        buffer.add_particle(1.0, 2.0);
        assert_eq!(buffer.x.len(), 1);
        assert_eq!(buffer.y.len(), 1);
        assert_eq!(buffer.x[0], 1.0);
        assert_eq!(buffer.y[0], 2.0);
    }

    #[test]
    fn test_get_point() {
        let mut buffer = ParticleBuffer::new(10);
        buffer.add_particle(1.0, 2.0);
        let point = buffer.get_point(0);
        assert_eq!(*point.x, 1.0);
        assert_eq!(*point.y, 2.0);
    }

    #[test]
    fn test_bounding_box() {
        let mut buffer = ParticleBuffer::new(10);
        buffer.add_particle(1.0, 2.0);
        buffer.add_particle(3.0, 4.0);
        buffer.add_particle(-1.0, 0.0);
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
        buffer.add_particle(0.75, 0.75);
        buffer.add_particle(-0.75, 0.75);
        buffer.add_particle(0.75, -0.75);
        buffer.add_particle(-0.75, -0.75);

        let quad_tree = QuadTree::new(&buffer, 4, 1);
        assert_eq!(quad_tree.root.particle_indices.len(), 0);
        assert_eq!(quad_tree.root.children[0].is_some(), true);
        assert_eq!(quad_tree.root.children[1].is_some(), true);
        assert_eq!(quad_tree.root.children[2].is_some(), true);
        assert_eq!(quad_tree.root.children[3].is_some(), true);

        assert_eq!(quad_tree.root.children[0].as_ref().unwrap().particle_indices.len(), 1);
        assert_eq!(quad_tree.root.children[1].as_ref().unwrap().particle_indices.len(), 1);
        assert_eq!(quad_tree.root.children[2].as_ref().unwrap().particle_indices.len(), 1);
        assert_eq!(quad_tree.root.children[3].as_ref().unwrap().particle_indices.len(), 1);
    }
}
