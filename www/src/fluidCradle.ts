
import { FluidSimulation } from "fluids2";
import { memory } from "fluids2/fluids2_bg.wasm";

import GUI from "lil-gui";

export class SimulationConfig {
    numParticles: number;
    gravity: number;
    boundaryElasticity: number;
    particleElasticity: number;
    particleRadius: number;
    boundingBoxSize: number;
    initialVelocity: number;
}

class SimulationBuffers {
    position_: Float32Array;
    velocity_: Float32Array;
    normalizedVelocity_: Float32Array;

    constructor(position: Float32Array, velocity: Float32Array, normalizedVelocity: Float32Array) {
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

export class SimulationCradle {
    constructor(config: SimulationConfig) {
        this.config_ = config;
    };

    initialize(): void {
        this.simulation_ = FluidSimulation.new(this.config_.numParticles);
        this.buffers_ = new SimulationBuffers(
            new Float32Array(memory.buffer, this.simulation_.position_buffer(), this.config_.numParticles * 3),
            new Float32Array(memory.buffer, this.simulation_.velocity_buffer(), this.config_.numParticles * 3),
            new Float32Array(memory.buffer, this.simulation_.velocity_normalized_buffer(), this.config_.numParticles)
        );
        this.simulation_.init_cube(this.config_.boundingBoxSize);
        this.simulation_.init_uniform_velocity(this.config_.initialVelocity);
        this.onInitialize_(this.config_);
    }

    restart(): void {
        this.initialize();
    }

    step(timeStep: number): void {
        this.simulation_.set_gravity(this.config_.gravity);
        this.simulation_.set_boundary_elasticity(this.config_.boundaryElasticity);
        this.simulation_.set_particle_elasiticity(this.config_.particleElasticity);
        this.simulation_.set_particle_radius(this.config_.particleRadius);
        this.simulation_.update(timeStep);
    }

    bind(gui: GUI): void {
        gui.add(this.config_, "gravity", 0, 200);
        gui.add(this.config_, "boundaryElasticity", 0, 1);
        gui.add(this.config_, "particleElasticity", 0, 1);
        gui.add(this.config_, "particleRadius", 1, 10);

        const simulationConfig = gui.addFolder("");
        simulationConfig.add(this.config_, "numParticles", 0, 1000).onFinishChange(() => { this.initialize() });
        simulationConfig.add(this.config_, "initialVelocity", 0, 500).onFinishChange(() => { this.initialize() });
        simulationConfig.add(this, "restart");
    }

    addInitializeListener(func: (config: SimulationConfig) => void) {
        this.onInitialize_ = func;
    }

    get buffers() {
        return this.buffers_;
    }

    get config() {
        return this.config_;
    }

    set config(config: SimulationConfig) {
        this.config_ = config;
    }


    simulation_: FluidSimulation;
    config_: SimulationConfig;
    buffers_: SimulationBuffers;

    onInitialize_: (config: SimulationConfig) => void;
}