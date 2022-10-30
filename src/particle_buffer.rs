use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub struct ParticleBuffer {
    buffer_: Vec<f32>,
    num_particles_: usize,
}

#[wasm_bindgen]
impl ParticleBuffer { 
    pub fn new(num_particles: usize) -> ParticleBuffer {
        ParticleBuffer {
            buffer_:  vec![0.0; num_particles * 3],
            num_particles_: num_particles,
        }
    }
    
    pub fn buffer(&self) -> *const f32 {
        return self.buffer_.as_ptr();
    }

    pub fn set_point(&mut self, idx: usize, x: f32, y: f32, z: f32) {
        let i = idx * 3;
        self.buffer_[i] = x;
        self.buffer_[i+1] = y;
        self.buffer_[i+2] = z;
    }

    pub fn size(&self) -> usize {
        self.num_particles_
    }
}
