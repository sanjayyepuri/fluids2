import { BarnesHut } from "fluids2";
import { memory } from "fluids2/fluids2_bg.wasm";
import { Particle, ParticleSimulation } from "./particleVisualization";

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

  constructor(positionX: Float32Array, positionY: Float32Array) {
    this.positionX = positionX;
    this.positionY = positionY;
  }
}

export class BarnesHutCradle implements ParticleSimulation {
  constructor(config: SimulationConfig) {
    this.config_ = config;
    this.initialize();
  }

  initialize(): void {
    this.simulation_ = BarnesHut.new(
      this.config_.numParticles,
      this.config_.theta,
      this.config_.maxDepth,
      this.config_.maxParticlesPerNode,
    );

    this.simulation_.init_cube(this.config_.numParticles);

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
    );
  }

  step(time: number): void {
    this.simulation_.step(time);
  }

  getParticle(index: number): Particle {
    return {
      x: this.buffers_.positionX[index],
      y: this.buffers_.positionY[index],
      z: 0,
    };
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
}
