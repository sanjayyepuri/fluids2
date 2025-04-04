import * as THREE from "three";
import { FluidSimulation } from "fluids2";
import { memory } from "fluids2/fluids2_bg.wasm";

import GUI from "lil-gui";
import { Particle, ParticleSimulation } from "./particleVisualization";

export class SimulationConfig {
  numParticles: number;
  gravity: number;
  boundaryElasticity: number;
  particleElasticity: number;
  particleRadius: number;
  initialVelocity: number;
}

class SimulationBuffers {
  position_: Float32Array;
  velocity_: Float32Array;
  normalizedVelocity_: Float32Array;

  constructor(
    position: Float32Array,
    velocity: Float32Array,
    normalizedVelocity: Float32Array,
  ) {
    this.position_ = position;
    this.velocity_ = velocity;
    this.normalizedVelocity_ = normalizedVelocity;
  }

  get position() {
    return this.position_;
  }

  get velocity() {
    return this.velocity_;
  }

  get normalizedVelocity() {
    return this.normalizedVelocity_;
  }
}

export class SimulationCradle implements ParticleSimulation {
  constructor(config: SimulationConfig) {
    this.config_ = config;

    this.color1 = new THREE.Color("#ff0000");
    this.color2 = new THREE.Color("#049ef4");

    this.initialize();
  }

  initialize(): void {
    this.simulation_ = FluidSimulation.new(this.config_.numParticles);
    this.buffers_ = new SimulationBuffers(
      new Float32Array(
        memory.buffer,
        this.simulation_.position_buffer(),
        this.config_.numParticles * 3,
      ),
      new Float32Array(
        memory.buffer,
        this.simulation_.velocity_buffer(),
        this.config_.numParticles * 3,
      ),
      new Float32Array(
        memory.buffer,
        this.simulation_.velocity_normalized_buffer(),
        this.config_.numParticles,
      ),
    );
    this.simulation_.init_cube();
    this.simulation_.init_uniform_velocity(this.config_.initialVelocity);
  }

  restart(): void {
    this.initialize();
  }

  updateParameters(): void {
    this.simulation_.set_gravity(this.config_.gravity);
    this.simulation_.set_boundary_elasticity(this.config_.boundaryElasticity);
    this.simulation_.set_particle_elasiticity(this.config_.particleElasticity);
    this.simulation_.set_particle_radius(this.config_.particleRadius);
  }

  step(time: number): void {
    this.simulation_.update(time);
  }

  bind(gui: GUI): void {
    gui.add(this.config_, "gravity", 0, 200);
    gui.add(this.config_, "boundaryElasticity", 0, 1);
    gui.add(this.config_, "particleElasticity", 0, 1);
    gui.add(this.config_, "particleRadius", 1, 10);

    const simulationConfig = gui.addFolder("");
    simulationConfig
      .add(this.config_, "numParticles", 0, 1000)
      .onFinishChange(() => {
        this.initialize();
      });
    simulationConfig
      .add(this.config_, "initialVelocity", 0, 500)
      .onFinishChange(() => {
        this.initialize();
      });
    simulationConfig.add(this, "restart");
  }

  reinitialize(): void {
    this.initialize();
  }

  getParticle(index: number): Particle {
    return {
      x: this.buffers_.position[index * 3] - 25,
      y: this.buffers_.position[index * 3 + 1] - 25,
      z: this.buffers_.position[index * 3 + 2] - 25,
    };
  }

  particleColor(index: number): number {
    let c = new THREE.Color(this.color2);
    c.lerp(this.color1, this.buffers_.normalizedVelocity[index]);
    return c.getHex();
  }

  numParticles(): number {
    return this.config_.numParticles;
  }

  simulation_: FluidSimulation;
  config_: SimulationConfig;
  buffers_: SimulationBuffers;

  color1: THREE.Color;
  color2: THREE.Color;

  onInitialize_: (config: SimulationConfig) => void;
}
