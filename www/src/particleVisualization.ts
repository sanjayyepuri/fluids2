import * as THREE from "three";
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls";
import Stats from "three/examples/jsm/libs/stats.module"

export type Particle = {
  x: number;
  y: number;
  z: number;
};

export interface ParticleSimulation {
  // steps the simulation forward in time
  step: (dt: number) => void;

  // returns the particle at the given index
  getParticle: (index: number) => Particle;

  particleColor: (index: number) => number;

  numParticles: () => number;
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

    const sphere = new THREE.SphereGeometry(1, 16, 8);
    const material = new THREE.MeshStandardMaterial();
    this.particles = new THREE.InstancedMesh(
      sphere,
      material,
      this.simulation.numParticles(),
    );
    this.scene.add(this.particles);

    // Add scene lights
    const directionalLight = new THREE.DirectionalLight(0xffffff, 0.5);
    const ambientLight = new THREE.AmbientLight();
    this.scene.add(directionalLight);
    this.scene.add(ambientLight);

    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.clock = new THREE.Clock();

    this.stats = Stats();
    this.stats.showPanel(0);
  }

  initializeDom() {
    document.body.appendChild(this.renderer.domElement);
    document.body.appendChild(this.stats.dom);
  }

  render() {
    this.simulation.step(this.clock.getDelta());

    let translation = new THREE.Matrix4();

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

  simulation: ParticleSimulation;

  renderer: THREE.WebGLRenderer;
  scene: THREE.Scene;
  camera: THREE.PerspectiveCamera;
  particles: THREE.InstancedMesh;
  controls: OrbitControls;
  clock: THREE.Clock;
  stats: Stats;
}
