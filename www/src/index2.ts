import * as THREE from "three";
import Stats from "three/examples/jsm/libs/stats.module"
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls";
import GUI from "lil-gui";

import { SimulationConfig, SimulationCradle } from "./fluidCradle";

const config = new SimulationConfig();
config.gravity = 30.0;
config.numParticles = 100;
config.boundaryElasticity = 0.5;
config.particleElasticity = 0.95;
config.particleRadius = 1.0;
config.initialVelocity = 0.0;

// Setup Scene
const renderer = new THREE.WebGLRenderer();
renderer.setSize(window.innerWidth, window.innerHeight);

const color1 = new THREE.Color("#ff0000");
const color2 = new THREE.Color("#049ef4");

const geometry = new THREE.SphereGeometry(config.particleRadius, 16, 8);
const material = new THREE.MeshStandardMaterial();

const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(
    75, window.innerWidth / window.innerHeight, 0.1, 3500
);
camera.position.z = 100;
const directionalLight = new THREE.DirectionalLight(0xffffff, 0.5);
const ambientLight = new THREE.AmbientLight();
scene.add(directionalLight);
scene.add(ambientLight);
const controls = new OrbitControls(camera, renderer.domElement);

const gui = new GUI();
const stats = Stats();

var particleEnergyPanel = stats.addPanel(Stats.Panel('energy', '#ff8', '#221'));

stats.showPanel(3);

document.body.appendChild(renderer.domElement);
document.body.appendChild(stats.dom);

// Setup Simulation
const simulationCradle = new SimulationCradle(config);
simulationCradle.bind(gui);

let points: THREE.InstancedMesh = null;
simulationCradle.addInitializeListener(config => {
    if (points != null) {
        scene.remove(points);
    }
    points = new THREE.InstancedMesh(geometry, material, config.numParticles);
    points.setColorAt(0, new THREE.Color());
    scene.add(points);
});

simulationCradle.initialize();

const clock = new THREE.Clock();
const dummy = new THREE.Object3D();
function animate() {
    requestAnimationFrame(animate);

    simulationCradle.step(Math.min(clock.getDelta(), 1));
    particleEnergyPanel.update(simulationCradle.simulation_.net_particle_energy(), 100000);

    for (let i = 0; i < simulationCradle.config.numParticles; ++i) {
        dummy.position.x = simulationCradle.buffers.position[i * 3] - 25;
        dummy.position.y = simulationCradle.buffers.position[i * 3 + 1] - 25;
        dummy.position.z = simulationCradle.buffers.position[i * 3 + 2] - 25;
        dummy.updateMatrix()

        points.setMatrixAt(i, dummy.matrix);
        let c = new THREE.Color(color2);
        points.setColorAt(i, c.lerp(color1, simulationCradle.buffers.normalizedVelocity[i]));
        points.scale.setScalar(simulationCradle.config.particleRadius);
        points.instanceMatrix.needsUpdate = true;
        points.instanceColor.needsUpdate = true;
    }

    controls.update()
    renderer.render(scene, camera);

    stats.update();
}

function onWindowResize() {
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(window.innerWidth, window.innerHeight);
}

window.addEventListener('resize', onWindowResize);
animate();