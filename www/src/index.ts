import * as THREE from "three";
import Stats from "three/examples/jsm/libs/stats.module"
import GUI from "lil-gui";

import { FluidSimulation } from "fluids2";
import { memory } from "fluids2/fluids2_bg.wasm";


const renderer = new THREE.WebGLRenderer();
renderer.setSize(window.innerWidth, window.innerHeight);
document.body.appendChild(renderer.domElement);

const gui = new GUI();
const stats = Stats()
document.body.appendChild(stats.dom);

const config = { speed: 0.01, };
gui.add(config, "speed", 0, 0.10, 0.001);

const numParticles = 5000;
const n = 500;

const fluidSimulation = FluidSimulation.new(numParticles);
fluidSimulation.init_cube(n);

const positions = new Float32Array(memory.buffer, fluidSimulation.particle_buffer(), numParticles * 3);
const colors = [];
for ( let i = 0; i < positions.length; i ++ ) {
    colors.push((positions[i] / n) + 0.5);
}

const geometry = new THREE.BufferGeometry();
geometry.setAttribute("position", new THREE.Float32BufferAttribute(positions, 3));
geometry.setAttribute("color", new THREE.Float32BufferAttribute(colors, 3));

const material = new THREE.PointsMaterial( { size: 15, vertexColors: true } );
const points = new THREE.Points(geometry, material);

const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(
    75, window.innerWidth / window.innerHeight, 0.1, 3500 
);
scene.add(points);
camera.position.z = 2000;

function animate() {
    requestAnimationFrame(animate);
    points.rotation.x += config.speed;
    points.rotation.y += config.speed;
    renderer.render(scene, camera);
    stats.update();
}

animate();