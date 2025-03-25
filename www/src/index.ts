import * as THREE from "three";
import Stats from "three/examples/jsm/libs/stats.module";
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls";
import GUI from "lil-gui";


import { SimulationConfig, BarnesHutCradle } from "./barnesHut";

const config = new SimulationConfig();
config.numParticles = 10;
config.theta = 0.5;
config.maxDepth = 4;
config.maxParticlesPerNode = 4;
config.gravitionalConstant = 0.1;

// Setup Scene
const renderer = new THREE.WebGLRenderer();
renderer.setSize(window.innerWidth, window.innerHeight);
const color1 = new THREE.Color("#ff0000");
const color2 = new THREE.Color("#049ef4");
const geometry = new THREE.SphereGeometry(1, 16, 8);
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


stats.showPanel(3);

document.body.appendChild(renderer.domElement);
document.body.appendChild(stats.dom);

const barnesHut = new BarnesHutCradle(config);
barnesHut.initialize();

let points: THREE.InstancedMesh = null;
points = new THREE.InstancedMesh(geometry, material, config.numParticles);
points.setColorAt(0, new THREE.Color());
scene.add(points);

const clock = new THREE.Clock();
const dummy = new THREE.Object3D();

function animate() {
    requestAnimationFrame(animate);

    barnesHut.step(clock.getDelta());

    for (let i = 0; i < config.numParticles; i++) {
        const [x, y] = barnesHut.buffers_.getPosition(i);
        dummy.position.x = x - 25;
        dummy.position.y = y - 25;
        dummy.position.z = 0;
        dummy.updateMatrix();

        points.setMatrixAt(i, dummy.matrix);
        points.setColorAt(i, new THREE.Color(color1.r, color1.g, color1.b));
        console.log(`Particle ${i}: (${x}, ${y})`);

        points.scale.setScalar(1);
        points.instanceMatrix.needsUpdate = true;
        points.instanceColor.needsUpdate = true;
    }

    controls.update();
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