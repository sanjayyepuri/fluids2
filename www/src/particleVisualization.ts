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
    this.stats.showPanel(0);

    this.initializeControlPanel();
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

    let translation = new THREE.Matrix4();

    console.log(this.simulation.getParticle(0));

    for (let i = 0; i < this.simulation.numParticles(); i++) {
      const particle = this.simulation.getParticle(i);
      this.particles.setMatrixAt(
        i,
        translation.makeTranslation(particle.x, particle.y, particle.z),
      );
      this.particles.setColorAt(
        i,
        new THREE.Color(this.simulation.particleColor(i)),
      );
    }

    this.particles.scale.setScalar(1);
    this.particles.instanceMatrix.needsUpdate = true;
    this.particles.instanceColor.needsUpdate = true;

    this.controls.update();
    this.renderer.render(this.scene, this.camera);
    this.stats.update();
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
  configPanel: GUI;

  paused: boolean = false;
  stepSize: number = 0.016;
}
