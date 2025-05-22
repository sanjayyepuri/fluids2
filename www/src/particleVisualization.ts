import GUI from "lil-gui";
import * as THREE from "three";
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls";
import Stats from "three/examples/jsm/libs/stats.module";

export type Particle = {
  x: number;
  y: number;
  z: number;
};

export interface ParticleSimulation {
  // updates any dynamic values in the simulation prior to the step
  // TODO (sanjay) should this be part of the step?
  updateParameters: () => void;
  // steps the simulation forward in time
  step: (dt: number) => void;
  // returns the color to render the particle as
  particleColor: (index: number) => number;
  // returns the number of particles
  numParticles: () => number;
  // returns the particle at the given index
  getParticle: (index: number) => Particle;
  // binds the configuration object for the simulation and manages updates.
  bind: (gui: GUI) => void;
  // returns the current system energy
  getSystemEnergy: () => number;
  // returns the current system momentum
  getSystemMomentum: () => { x: number; y: number; z: number };
  // Resets the simulation to a known state
  reinitialize: () => void;
}

export class ParticleVisualization {
  constructor(simulation: ParticleSimulation) {
    this.simulation = simulation;

    this.renderer = new THREE.WebGLRenderer();
    this.renderer.setSize(window.innerWidth, window.innerHeight);

    // Setup Camera
    this.camera = new THREE.PerspectiveCamera(
      75,
      window.innerWidth / window.innerHeight,
      0.1,
      3500,
    );
    this.camera.position.z = 100;

    // Setup scene geometry
    this.scene = new THREE.Scene();
    this.initializeParticles();

    // Add scene lights
    const directionalLight = new THREE.DirectionalLight(0xffffff, 0.5);
    const ambientLight = new THREE.AmbientLight();
    this.scene.add(directionalLight);
    this.scene.add(ambientLight);

    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.clock = new THREE.Clock();

    this.stats = Stats();
    // Default panels: 0: fps, 1: ms, 2: mb
    this.energyPanel = new Stats.Panel('Energy', '#ff8', '#221'); // Name, fg color, bg color
    this.stats.addPanel(this.energyPanel);
    this.momentumXPanel = new Stats.Panel('Px', '#f8f', '#212');
    this.stats.addPanel(this.momentumXPanel);
    this.momentumYPanel = new Stats.Panel('Py', '#f8f', '#212');
    this.stats.addPanel(this.momentumYPanel);
    this.momentumZPanel = new Stats.Panel('Pz', '#f8f', '#212');
    this.stats.addPanel(this.momentumZPanel);
    this.stats.showPanel(0); // Start with FPS, user can cycle through (Energy is 3, Px 4, Py 5, Pz 6)

    this.initializeControlPanel();
  }

  updateEnergyStat(energy: number) {
    if (this.energyPanel) {
      // For stats.js, the first argument to update is the value,
      // and the second is the maxValue for the graph.
      // We might need to normalize or find a dynamic maxValue later.
      // For now, let's pass a reasonably high number for maxValue or just the value itself
      // if the panel text display is sufficient.
      // Let's assume a large max value for the graph, or just update the text part.
      // The text part will show the 'energy' value regardless of maxValue for graph.
      this.energyPanel.update(energy, 2000); // Assuming energy values might go up to 2000 for graph scale
    }
  }

  updateMomentumStats(momentum: { x: number; y: number; z: number }) {
    if (this.momentumXPanel) {
      // Assuming momentum values could be e.g. +/- 500 for graph scale
      this.momentumXPanel.update(momentum.x, 500);
    }
    if (this.momentumYPanel) {
      this.momentumYPanel.update(momentum.y, 500);
    }
    if (this.momentumZPanel) {
      this.momentumZPanel.update(momentum.z, 500);
    }
  }

  initializeParticles() {
    if (this.particles != null) {
      this.scene.remove(this.particles);
    }
    const sphere = new THREE.SphereGeometry(1, 16, 8);
    const material = new THREE.MeshStandardMaterial();
    this.particles = new THREE.InstancedMesh(
      sphere,
      material,
      this.simulation.numParticles(),
    );
    this.scene.add(this.particles);
  }

  initializeControlPanel() {
    this.configPanel = new GUI();
    this.configPanel.add(this, "pause");
    this.configPanel.add(this, "resume");
    this.configPanel.add(this, "step");
    this.configPanel.add(this, "stepSize");
    this.configPanel.add(this, "reset");

    const simulationPanel = this.configPanel.addFolder("Simulation");
    this.simulation.bind(simulationPanel);
  }


  initializeDom() {
    document.body.appendChild(this.renderer.domElement);
    document.body.appendChild(this.stats.dom);
  }

  pause() {
    this.paused = true;
  }

  resume() {
    this.paused = false;
  }

  step() {
    if (this.paused) {
      this.simulation.updateParameters();
      this.simulation.step(this.stepSize);
    }
  }

  render() {
    const delta = this.clock.getDelta();

    if (!this.paused) {
      this.simulation.updateParameters();
      this.simulation.step(delta);
    }

    let matrix = new THREE.Matrix4();
    const baseScale = 0.5; // Base size of the particle
    const focalLength = 500; // Affects perspective effect; adjust as needed

    for (let i = 0; i < this.simulation.numParticles(); i++) {
      const particle = this.simulation.getParticle(i);

      // Simple perspective scaling: scale = focalLength / (focalLength + z)
      // Add a small epsilon to particle.z if it can be zero and focalLength is also small, to avoid division by zero.
      // However, with typical z values and focalLength, (focalLength + particle.z) should be positive.
      // We also clamp the scale to avoid extremely large/small particles.
      let scale = focalLength / (focalLength + particle.z);
      scale = Math.max(0.1, Math.min(scale, 5.0)) * baseScale; // Clamp scale and apply base scale

      // Create a matrix that includes translation and scaling
      matrix.compose(
        new THREE.Vector3(particle.x, particle.y, particle.z), // position
        new THREE.Quaternion(), // rotation (identity)
        new THREE.Vector3(scale, scale, scale), // scale
      );
      this.particles.setMatrixAt(i, matrix);

      this.particles.setColorAt(
        i,
        new THREE.Color(this.simulation.particleColor(i)),
      );
    }

    this.particles.instanceMatrix.needsUpdate = true;
    this.particles.instanceColor.needsUpdate = true;

    this.controls.update();
    this.renderer.render(this.scene, this.camera);

    const energy = this.simulation.getSystemEnergy();
    this.updateEnergyStat(energy);

    const momentum = this.simulation.getSystemMomentum();
    this.updateMomentumStats(momentum);

    this.stats.update(); // Update all stats panels
  }

  resize() {
    this.camera.aspect = window.innerWidth / window.innerHeight;
    this.camera.updateProjectionMatrix();
    this.renderer.setSize(window.innerWidth, window.innerHeight);
  }

  reset() {
    // invariant: the simulation should reinitialize with the latest values of the configuration.
    this.initializeParticles();
    this.simulation.reinitialize();
  }

  simulation: ParticleSimulation;

  renderer: THREE.WebGLRenderer;
  scene: THREE.Scene;
  camera: THREE.PerspectiveCamera;
  particles: THREE.InstancedMesh;
  controls: OrbitControls;
  clock: THREE.Clock;
  stats: Stats;
  energyPanel: Stats.Panel;
  momentumXPanel: Stats.Panel;
  momentumYPanel: Stats.Panel;
  momentumZPanel: Stats.Panel;
  configPanel: GUI;

  paused: boolean = false;
  stepSize: number = 0.016;
}
