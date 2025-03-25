import { BarnesHut } from "fluids2";
import { memory } from "fluids2/fluids2_bg.wasm";


export class SimulationConfig {
    numParticles: number;
    theta: number;
    maxDepth: number;
    maxParticlesPerNode: number;
    gravitionalConstant: number;
}

class SimulationBuffer {
    positionX: Float32Array;
    positionY: Float32Array;

    constructor(positionX: Float32Array, positionY: Float32Array) {
        this.positionX = positionX;
        this.positionY = positionY;
    }

    getPosition(index: number): [number, number] {
        return [this.positionX[index], this.positionY[index]];
    }
}

export class BarnesHutCradle {
    constructor(config: SimulationConfig) {
        this.config_ = config;
    }

    initialize(): void {
        this.simulation_ = BarnesHut.new(
            this.config_.numParticles,
            this.config_.theta,
            this.config_.maxDepth,
            this.config_.maxParticlesPerNode);

        this.simulation_.init_cube(this.config_.numParticles);

        this.buffers_ = new SimulationBuffer(
            new Float32Array(memory.buffer, this.simulation_.get_raw_x_coordinates(), this.config_.numParticles),
            new Float32Array(memory.buffer, this.simulation_.get_raw_y_coordinates(), this.config_.numParticles)
        );

        // log all positions
        for (let i = 0; i < this.config_.numParticles; i++) {
            const [x, y] = this.buffers_.getPosition(i);
            console.log(`Particle ${i}: (${x}, ${y})`);
        }
    }

    step(time: number): void {
        this.simulation_.step(time);
    }

    config_: SimulationConfig;
    simulation_: BarnesHut;
    buffers_: SimulationBuffer;
}

