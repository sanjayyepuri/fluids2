use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub struct ParticleBuffer {
    position_buffer_: Vec<f32>,
    velocity_buffer_: Vec<f32>,
    num_particles_: usize,
}

#[wasm_bindgen]
impl ParticleBuffer {
    pub fn new(num_particles: usize) -> ParticleBuffer {
        ParticleBuffer {
            position_buffer_: vec![0.0; num_particles * 3],
            velocity_buffer_: vec![0.0; num_particles * 3],
            num_particles_: num_particles,
        }
    }

    pub fn position_buffer(&self) -> *const f32 {
        return self.position_buffer_.as_ptr();
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

    pub fn advect_particles(&mut self) {
        for i in 0..self.size() * 3 {
            self.position_buffer_[i] += self.velocity_buffer_[i];
        }
    }

    pub fn apply_gravity(&mut self, gravity: f32) {
        for i in 0..self.size() {
            let p = &self.position_buffer_[i*3..i*3+3];
            let v = &mut self.velocity_buffer_[i*3..i*3+3];
            let factor = gravity / (p[0] * p[0] + p[1] * p[1] + p[2] * p[2] + 1e-10);
            v[0] += -p[0] * factor;
            v[1] += -p[1] * factor;
            v[2] += -p[2] * factor;
        }
    }
}
