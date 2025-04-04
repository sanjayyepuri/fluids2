use rand::Rng;
use wasm_bindgen::prelude::*;

pub const BOUNDING_WIDTH: f32 = 50.0;

pub struct BoundingBox {
    pub x0: f32,
    pub x1: f32,
    pub y0: f32,
    pub y1: f32,
    pub z0: f32,
    pub z1: f32,
}

pub struct ParticleProperties {
    pub radius_: f32,
    pub elasticity_: f32,
    pub mass_: f32,
}


pub struct ParticleBuffer {
    position_buffer_: Vec<f32>,
    velocity_buffer_: Vec<f32>,
    velocity_normalized_buffer_: Vec<f32>,
    // force_buffer_: Vec<f32>, // TODO (sanjay) add force
    net_kinetic_energy: f32,
    net_potential_energy: f32,
    num_particles_: usize,
}


impl ParticleBuffer {
    pub fn new(num_particles: usize) -> ParticleBuffer {
        ParticleBuffer {
            position_buffer_: vec![0.0; num_particles * 3],
            velocity_buffer_: vec![0.0; num_particles * 3],
            velocity_normalized_buffer_: vec![0.0; num_particles],
            // force_buffer_: vec![0.0; num_particles * 3],
            net_kinetic_energy: 0.0,
            net_potential_energy: 0.0,
            num_particles_: num_particles,
        }
    }

    pub fn position_buffer(&self) -> *const f32 {
        return self.position_buffer_.as_ptr();
    }

    pub fn position_mut(&mut self, idx: usize) -> &mut [f32]{
        let i = idx * 3;
        return &mut self.position_buffer_[i..i+3];
    }

    pub fn velocity_mut(&mut self, idx: usize) -> &mut [f32]{
        let i = idx * 3;
        return &mut self.velocity_buffer_[i..i+3];
    }

    pub fn position(&self, idx: usize) -> &[f32]{
        let i = idx * 3;
        return &self.position_buffer_[i..i+3];
    }

    pub fn velocity(&self, idx: usize) -> &[f32]{
        let i = idx * 3;
        return &self.velocity_buffer_[i..i+3];
    }

    pub fn velocity_buffer(&self) -> *const f32 {
        return self.velocity_buffer_.as_ptr();
    }

    pub fn velocity_normalized_buffer(&self) -> *const f32 {
        return self.velocity_normalized_buffer_.as_ptr();
    }

    pub fn net_particle_energy(&self) -> f32 {
        return self.net_particle_kinetic_energy() + self.net_particle_potential_energy();
    }

    pub fn net_particle_kinetic_energy(&self) -> f32 {
        return self.net_kinetic_energy;
    }

    pub fn net_particle_potential_energy(&self) -> f32 {
        return self.net_potential_energy;
    }

    pub fn compute_system_attributes(&mut self, gravity: f32) {
        // TODO (sanjay) Likely want to split computing these metrics into separate functions.
        self.net_kinetic_energy = 0.0;
        self.net_potential_energy = 0.0;

        for i in 0..self.size() {
            let v = &self.velocity_buffer_[i*3..i*3+3];
            let p = &self.position_buffer_[i*3..i*3+3];
            let v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
            self.net_kinetic_energy += 0.5 * 1.0 * v2;
            self.net_potential_energy += 1.0 * gravity * p[1];
            self.velocity_normalized_buffer_[i] = v2.sqrt();
            self.velocity_normalized_buffer_[i] = self.velocity_normalized_buffer_[i] / (self.velocity_normalized_buffer_[i] + 5.0);
        }
    }

    pub fn set_position(&mut self, idx: usize, x: f32, y: f32, z: f32) {
        let i = idx * 3;
        self.position_buffer_[i] = x;
        self.position_buffer_[i + 1] = y;
        self.position_buffer_[i + 2] = z;
    }

    pub fn set_velocity(&mut self, idx: usize, u: f32, v: f32, w: f32) {
        let i = idx * 3;
        self.velocity_buffer_[i] = u;
        self.velocity_buffer_[i + 1] = v;
        self.velocity_buffer_[i + 2] = w;
    }

    pub fn size(&self) -> usize {
        self.num_particles_
    }

    pub fn advect_particles(&mut self, t: f32) {
        for i in 0..self.size() * 3 {
            self.position_buffer_[i] += t * self.velocity_buffer_[i];
        }
    }

    pub fn apply_centroid_gravity(&mut self, gravity: f32, t: f32) {
        for i in 0..self.size() {
            let p = &self.position_buffer_[i*3..i*3+3];
            let v = &mut self.velocity_buffer_[i*3..i*3+3];
            let rs = p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
            let factor = t * gravity / (rs + 1e-10);
            v[0] += -p[0] * factor * 1.1;
            v[1] += -p[1] * factor * 1.1;
            v[2] += -p[2] * factor * 1.1;
        }
    }

    pub fn apply_uniform_gravity(&mut self, gravity: f32, t: f32) {
        for i in 0..self.size() {
            let v = &mut self.velocity_buffer_[i*3..i*3+3];
            v[1] += -gravity * t;
        }
    }
}

#[wasm_bindgen]
pub struct FluidSimulation {
    particle_buffer_: ParticleBuffer,
    particle_properties_: ParticleProperties,
    bounding_box_: BoundingBox,
    gravity_: f32,
    boundary_elasticity_: f32,
}

#[wasm_bindgen]
impl FluidSimulation {
    pub fn new(num_particles: usize) -> FluidSimulation {
        FluidSimulation {
            particle_buffer_: ParticleBuffer::new(num_particles),
            particle_properties_: ParticleProperties {
                radius_: 5.0,
                elasticity_: 0.05,
                mass_: 1.0,
            },
            bounding_box_: BoundingBox {
                x0: 0.0,
                x1: BOUNDING_WIDTH,
                y0: 0.0,
                y1: BOUNDING_WIDTH,
                z0: 0.0,
                z1: BOUNDING_WIDTH,
            },
            gravity_: 0.0,
            boundary_elasticity_: 0.05,
        }
    }

    pub fn update(&mut self, t: f32) {
        self.particle_buffer_
            .apply_uniform_gravity(self.gravity_, t);
        self.compute_particle_interactions(t);
        self.particle_buffer_.advect_particles(t);
        self.compute_boundary_condition();
        self.particle_buffer_
            .compute_system_attributes(self.gravity_);
    }

    pub fn position_buffer(&self) -> *const f32 {
        return self.particle_buffer_.position_buffer();
    }

    pub fn velocity_buffer(&self) -> *const f32 {
        return self.particle_buffer_.velocity_buffer();
    }

    pub fn velocity_normalized_buffer(&self) -> *const f32 {
        return self.particle_buffer_.velocity_normalized_buffer();
    }

    pub fn net_particle_energy(&self) -> f32 {
        return self.particle_buffer_.net_particle_energy();
    }

    pub fn set_gravity(&mut self, gravity: f32) {
        self.gravity_ = gravity;
    }

    pub fn set_boundary_elasticity(&mut self, elasticity: f32) {
        self.boundary_elasticity_ = elasticity;
    }

    pub fn set_particle_elasiticity(&mut self, elasticity: f32) {
        self.particle_properties_.elasticity_ = elasticity;
    }

    pub fn set_particle_radius(&mut self, radius: f32) {
        self.particle_properties_.radius_ = radius;
    }

    pub fn compute_boundary_condition(&mut self) {
        for i in 0..self.particle_buffer_.size() {
            let position = self.particle_buffer_.position_mut(i);

            let mut x_b = false;
            let mut y_b = false;
            let mut z_b = false;

            if position[0] < self.bounding_box_.x0 {
                position[0] = self.bounding_box_.x0;
                x_b = true;
            }
            if position[0] > self.bounding_box_.x1 {
                position[0] = self.bounding_box_.x1;
                x_b = true;
            }
            if position[1] < self.bounding_box_.y0 {
                position[1] = self.bounding_box_.y0;
                y_b = true;
            }
            if position[1] > self.bounding_box_.y1 {
                position[1] = self.bounding_box_.y1;
                y_b = true;
            }
            if position[2] < self.bounding_box_.z0 {
                position[2] = self.bounding_box_.z0;
                z_b = true;
            }
            if position[2] > self.bounding_box_.z1 {
                position[2] = self.bounding_box_.z1;
                z_b = true;
            }

            if x_b {
                self.particle_buffer_.velocity_mut(i)[0] *= -1.0 * self.boundary_elasticity_;
            }

            if y_b {
                self.particle_buffer_.velocity_mut(i)[1] *= -1.0 * self.boundary_elasticity_;
            }

            if z_b {
                self.particle_buffer_.velocity_mut(i)[2] *= -1.0 * self.boundary_elasticity_;
            }
        }
    }

    pub fn compute_particle_interactions(&mut self, t: f32) {
        for p1 in 0..self.particle_buffer_.size() {
            for p2 in (p1 + 1)..self.particle_buffer_.size() {
                let position1 = self.particle_buffer_.position(p1);
                let position2 = self.particle_buffer_.position(p2);

                let dx = position2[0] - position1[0];
                let dy = position2[1] - position1[1];
                let dz = position2[2] - position1[2];

                let distance = (dx * dx + dy * dy + dz * dz).sqrt() + 1e-8;

                let n = [dx / distance, dy / distance, dz / distance];

                if distance > self.particle_properties_.radius_ * 2.0 {
                    continue;
                }

                let u1 = self.particle_buffer_.velocity(p1).to_owned();
                let u2 = self.particle_buffer_.velocity(p2).to_owned();

                let c1 = ((u1[0] - u2[0]) * dx) + ((u1[1] - u2[1]) * dy) + ((u1[2] - u2[2]) * dz);
                let c1 = c1 / distance;
                for axis in 0..3 {
                    self.particle_buffer_.velocity_mut(p1)[axis] =
                        (u1[axis] - c1 * n[axis]) * self.particle_properties_.elasticity_;
                    self.particle_buffer_.velocity_mut(p2)[axis] =
                        (u2[axis] + c1 * n[axis]) * self.particle_properties_.elasticity_;
                }
            }
        }
    }

    pub fn init_cube(&mut self) {
        let mut rng = rand::thread_rng();

        for i in 0..self.particle_buffer_.size() {
            self.particle_buffer_.set_position(
                i,
                rng.gen::<f32>() * BOUNDING_WIDTH,
                rng.gen::<f32>() * BOUNDING_WIDTH,
                rng.gen::<f32>() * BOUNDING_WIDTH,
            );
        }
    }

    pub fn init_random_velocity(&mut self, max_velocity: f32) {
        let mut rng = rand::thread_rng();
        for i in 0..self.particle_buffer_.size() {
            self.particle_buffer_.set_velocity(
                i,
                rng.gen::<f32>() * max_velocity - max_velocity / 2.0,
                rng.gen::<f32>() * max_velocity - max_velocity / 2.0,
                rng.gen::<f32>() * max_velocity - max_velocity / 2.0,
            )
        }
    }

    pub fn init_uniform_velocity(&mut self, max_velocity: f32) {
        for i in 0..self.particle_buffer_.size() {
            self.particle_buffer_
                .set_velocity(i, max_velocity, max_velocity, max_velocity);
        }
    }
}
