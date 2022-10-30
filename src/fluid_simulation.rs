use rand::Rng;
use wasm_bindgen::prelude::*;

use crate::particle_buffer::ParticleBuffer;

#[wasm_bindgen]
pub struct FluidSimulation {
    particle_buffer_: ParticleBuffer,
}

#[wasm_bindgen]
impl FluidSimulation {
    pub fn new(num_particles: usize) -> FluidSimulation {
        FluidSimulation {
            particle_buffer_: ParticleBuffer::new(num_particles),
        }
    }

    pub fn particle_buffer(&self) -> *const f32 {
        return self.particle_buffer_.buffer();
    }
    
    pub fn init_cube(&mut self, size: f32) {
        let mut rng = rand::thread_rng();

        for i in 0..self.particle_buffer_.size() {
            self.particle_buffer_.set_point(i, 
                rng.gen::<f32>() * size - size/2.0,
                rng.gen::<f32>() * size - size/2.0,
                rng.gen::<f32>() * size - size/2.0,
            );
        }
    }    
}
