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
