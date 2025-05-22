import { BarnesHut } from "fluids2";
import { memory } from "fluids2/fluids2_bg.wasm";
import { Particle, ParticleSimulation } from "./particleVisualization";
import GUI from "lil-gui";

export class SimulationConfig {
  numParticles: number;
  theta: number;
  maxDepth: number;
  maxParticlesPerNode: number;
  gravitationalConstant: number;
}

class SimulationBuffer {
  positionX: Float32Array;
  positionY: Float32Array;
  positionZ: Float32Array;

  constructor(positionX: Float32Array, positionY: Float32Array, positionZ: Float32Array) {
    this.positionX = positionX;
    this.positionY = positionY;
    this.positionZ = positionZ;
  }
}

export class BarnesHutCradle implements ParticleSimulation {
  constructor(config: SimulationConfig) {
    this.config_ = config;
    this.initialize();
  }

  initializeBuffers(): void {
    this.buffers_ = new SimulationBuffer(
      new Float32Array(
        memory.buffer,
        this.simulation_.get_raw_x_coordinates(),
        this.config_.numParticles,
      ),
      new Float32Array(
        memory.buffer,
        this.simulation_.get_raw_y_coordinates(),
        this.config_.numParticles,
      ),
      new Float32Array(
        memory.buffer,
        this.simulation_.get_raw_z_coordinates(),
        this.config_.numParticles,
      ),
    );
  }

  initialize(): void {
    this.simulation_ = BarnesHut.new(
      this.config_.numParticles,
      this.config_.theta,
      this.config_.maxDepth,
      this.config_.maxParticlesPerNode,
      this.config_.gravitationalConstant,
    );

    this.simulation_.init_cube(this.config_.numParticles);
    this.initializeBuffers();
  }

  bind(gui: GUI): void {
    gui.add(this.config_, "gravitationalConstant", 0.1, 50);
    gui.add(this.config_, "theta", 0, 1);
    gui.add(this.config_, "maxParticlesPerNode", 1, 5);

    const simulationGui = gui.addFolder("Simulation Parameters");
    simulationGui.add(this.config_, "numParticles");
    simulationGui.add(this.config_, "maxDepth");
    simulationGui.add(this, "initializeBuffers");
  }

  updateParameters(): void {
    this.simulation_.set_gravitational_constant(
      this.config_.gravitationalConstant,
    );
    this.simulation_.set_theta(this.config_.theta);
    this.simulation_.set_max_particles_per_node(
      this.config_.maxParticlesPerNode,
    );
  }

  latestEnergy: number = 0;
  latestMomentum: { x: number; y: number; z: number } | null = null;

  step(time: number): void {
    this.simulation_.step(time);
    // TODO (sanjay) for some reason for large numbers of particles the buffers die.
    this.initializeBuffers();

    this.latestEnergy = this.simulation_.get_system_energy();
    const momentumPoint = this.simulation_.get_system_momentum();
    this.latestMomentum = { x: momentumPoint.x, y: momentumPoint.y, z: momentumPoint.z };
    // Ensure to free the Point manually if it's not done by wasm-bindgen automatically
    // and if Point is not a plain JS object returned by value.
    // Assuming Point from wasm_bindgen is a JS object that can be GC'd, or if it's on Rust heap,
    // that it's freed by its destructor when the JS wrapper is GC'd.
    // If `get_system_momentum` returns a pointer or an object that needs explicit freeing,
    // and if `momentumPoint` is such an object, then `momentumPoint.free()` might be needed.
    // For now, assuming it's a JS object or handled by wasm-bindgen's wrappers.
  }

  reinitialize(): void {
    this.simulation_.reinitialize(
      this.config_.numParticles,
      this.config_.theta,
      this.config_.maxDepth,
      this.config_.numParticles,
      this.config_.gravitationalConstant,
    );
    this.initializeBuffers();
    // TODO (sanjay) why does init cube take the number of particles as a configuration
    this.simulation_.init_cube(this.config_.numParticles);
  }

  getParticle(index: number): Particle {
    return {
      x: this.buffers_.positionX[index],
      y: this.buffers_.positionY[index],
      z: this.buffers_.positionZ[index],
    };
  }

  getSystemEnergy(): number {
    return this.latestEnergy;
  }

  getSystemMomentum(): { x: number; y: number; z: number } {
    return this.latestMomentum ?? { x: 0, y: 0, z: 0 };
  }

  particleColor(index: number): number {
    return 0xff0000;
  }

  numParticles(): number {
    return this.config_.numParticles;
  }

  config_: SimulationConfig;
  simulation_: BarnesHut;
  buffers_: SimulationBuffer;
  latestEnergy: number = 0;
  latestMomentum: { x: number; y: number; z: number } | null = null;
}
