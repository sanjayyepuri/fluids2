use rand::Rng;
use wasm_bindgen::prelude::*;

use crate::particle_buffer::{BoundingBox, ParticleBuffer, ParticleProperties, BOUNDING_WIDTH};

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
