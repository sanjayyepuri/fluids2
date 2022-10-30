import * as THREE from "three";
import * as fluids from "fluids2";
import { memory } from "fluids2/fluids2_bg.wasm";
import GUI from "lil-gui";

const gui = new GUI();

const config = {
    speed: 0.01,
};

gui.add(config, "speed", 0, 0.10, 0.001);

const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(
    75, window.innerWidth / window.innerHeight, 0.1, 3500 
);
const renderer = new THREE.WebGLRenderer();
renderer.setSize(window.innerWidth, window.innerHeight);
document.body.appendChild(renderer.domElement);

const geometry = new THREE.BufferGeometry();

const colors = [];

const numParts = 5000;
const n = 500;

const fluidSimulation = fluids.FluidSimulation.new(numParts);
fluidSimulation.init_cube(n);

const positions = new Float32Array(memory.buffer, fluidSimulation.particle_buffer(), numParts * 3);
for ( let i = 0; i < positions.length; i ++ ) {
    colors.push((positions[i] / n) + 0.5);
}

console.log(positions);

geometry.setAttribute("position", new THREE.Float32BufferAttribute(positions, 3));
geometry.setAttribute("color", new THREE.Float32BufferAttribute(colors, 3));

const material = new THREE.PointsMaterial( { size: 15, vertexColors: true } );

const points = new THREE.Points(geometry, material);

scene.add(points);

camera.position.z = 2000;

function animate() {
    requestAnimationFrame(animate);
    points.rotation.x += config.speed;
    points.rotation.y += config.speed;
    renderer.render(scene, camera);
}

animate();