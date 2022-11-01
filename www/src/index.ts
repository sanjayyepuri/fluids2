import * as THREE from "three";
import Stats from "three/examples/jsm/libs/stats.module"
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls";
import GUI from "lil-gui";

import { FluidSimulation } from "fluids2";
import { memory } from "fluids2/fluids2_bg.wasm";


const renderer = new THREE.WebGLRenderer();
renderer.setSize(window.innerWidth, window.innerHeight);
document.body.appendChild(renderer.domElement);

const gui = new GUI();
const stats = Stats()
document.body.appendChild(stats.dom);

const config = { gravity: 10.0 };
gui.add(config, "gravity", 0, 100);

const numParticles = 100;
const n = 500;

const fluidSimulation = FluidSimulation.new(numParticles);
fluidSimulation.init_cube(n);
fluidSimulation.init_random_velocity(5);

const positions = new Float32Array(memory.buffer, fluidSimulation.position_buffer(), numParticles * 3);

const geometry = new THREE.SphereGeometry( 10, 16, 8 );
const material = new THREE.MeshStandardMaterial( { color: 0x049ef4 } );

const points = new THREE.InstancedMesh(geometry, material, numParticles);


const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(
    75, window.innerWidth / window.innerHeight, 0.1, 3500
);

const directionalLight = new THREE.DirectionalLight( 0xffffff, 1 );
const ambientLight = new THREE.AmbientLight();
scene.add( directionalLight );
scene.add(ambientLight);

const controls = new OrbitControls( camera, renderer.domElement );


scene.add(points);

camera.position.z = 2000;


const dummy = new THREE.Object3D();

function animate() {
    requestAnimationFrame(animate);

    fluidSimulation.set_gravity(config.gravity);
    fluidSimulation.update()

    for (let i = 0; i < numParticles; ++i) {
        dummy.position.x = positions[i * 3];
        dummy.position.y = positions[i*3 + 1];
        dummy.position.z = positions[i*3 + 2];
        dummy.updateMatrix()

        points.setMatrixAt(i, dummy.matrix);
    }

    points.instanceMatrix.needsUpdate = true;

    controls.update()
    renderer.render(scene, camera);

    stats.update();
}

function onWindowResize() {
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();

    renderer.setSize( window.innerWidth, window.innerHeight );
}

window.addEventListener( 'resize', onWindowResize );
animate();