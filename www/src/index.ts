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

const config = { gravity: 30.0, numParticles: 100, boundary_elasticity: 0.5, particle_elasticity: 0.95, particle_radius: 1.0 };
gui.add(config, "gravity", 0, 200);
gui.add(config, "numParticles", 0, 5000);
gui.add(config, "boundary_elasticity", 0, 1);
gui.add(config, "particle_elasticity", 0, 1);
gui.add(config, "particle_radius", 1, 10);

const n = 50;

const fluidSimulation = FluidSimulation.new(config.numParticles);
fluidSimulation.init_cube(n);
fluidSimulation.init_uniform_velocity(0);

const positions = new Float32Array(memory.buffer, fluidSimulation.position_buffer(), config.numParticles * 3);
const velocities = new Float32Array(memory.buffer, fluidSimulation.velocity_buffer(), config.numParticles * 3);
const velocityNormalized = new Float32Array(memory.buffer, fluidSimulation.velocity_normalized_buffer(), config.numParticles);


const color1 = new THREE.Color("#ff0000");
const color2 = new THREE.Color("#049ef4");

const geometry = new THREE.SphereGeometry(config.particle_radius, 16, 8);
const material = new THREE.MeshStandardMaterial();

const points = new THREE.InstancedMesh(geometry, material, config.numParticles);
points.setColorAt(0, new THREE.Color());


const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(
    75, window.innerWidth / window.innerHeight, 0.1, 3500
);
camera.position.z = 100;

const directionalLight = new THREE.DirectionalLight(0xffffff, 0.5);
const ambientLight = new THREE.AmbientLight();

scene.add(directionalLight);
scene.add(ambientLight);
scene.add(points);

const controls = new OrbitControls( camera, renderer.domElement );

const clock = new THREE.Clock();



const dummy = new THREE.Object3D();
function animate() {
    requestAnimationFrame(animate);

    fluidSimulation.set_gravity(config.gravity);
    fluidSimulation.set_boundary_elasticity(config.boundary_elasticity);
    fluidSimulation.set_particle_elasiticity(config.particle_elasticity);
    fluidSimulation.set_particle_radius(config.particle_radius);
    fluidSimulation.update(clock.getDelta());


    for (let i = 0; i < config.numParticles; ++i) {
        dummy.position.x = positions[i * 3];
        dummy.position.y = positions[i*3 + 1];
        dummy.position.z = positions[i*3 + 2];
        dummy.updateMatrix()

        points.setMatrixAt(i, dummy.matrix);
        let c = new THREE.Color(color2);
        // points.setColorAt(i, color1.setHex( 0xffffff * Math.random() ));
        points.setColorAt(i, c.lerp(color1, velocityNormalized[i]));
        points.scale.setScalar(config.particle_radius);
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

    renderer.setSize( window.innerWidth, window.innerHeight );
}

window.addEventListener( 'resize', onWindowResize );
animate();