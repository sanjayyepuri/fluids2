use rand::Rng;
use wasm_bindgen::prelude::*;

use crate::particle_buffer::ParticleBuffer;

#[wasm_bindgen]
pub struct FluidSimulation {
    particle_buffer_: ParticleBuffer,
    gravity_: f32
}

#[wasm_bindgen]
impl FluidSimulation {
    pub fn new(num_particles: usize) -> FluidSimulation {
        FluidSimulation {
            particle_buffer_: ParticleBuffer::new(num_particles),
            gravity_: 0.0,
        }
    }

    pub fn update(&mut self, t: f32) {
        self.particle_buffer_.apply_gravity(self.gravity_, t);
        self.particle_buffer_.advect_particles(t);
    }

    pub fn position_buffer(&self) -> *const f32 {
        return self.particle_buffer_.position_buffer();
    }

    pub fn set_gravity(&mut self, gravity: f32) {
        self.gravity_ = gravity;
    }

    pub fn init_cube(&mut self, size: f32) {
        let mut rng = rand::thread_rng();

        for i in 0..self.particle_buffer_.size() {
            self.particle_buffer_.set_position(i,
                rng.gen::<f32>() * size - size / 2.0,
                rng.gen::<f32>() * size - size / 2.0,
                rng.gen::<f32>() * size - size / 2.0,
            );
        }
    }

    pub fn init_random_velocity(&mut self, max_velocity: f32) {
        let mut rng = rand::thread_rng();
        for i in 0..self.particle_buffer_.size() {
            self.particle_buffer_.set_velocity(i,
                rng.gen::<f32>() * max_velocity - max_velocity / 2.0,
                rng.gen::<f32>() * max_velocity - max_velocity / 2.0,
                rng.gen::<f32>() * max_velocity - max_velocity / 2.0,
            )
        }
    }
}
