import * as THREE from "three";
import * as fluids from "fluids2";
import GUI from "lil-gui";

const gui = new GUI();

const config = {
    speed: 0.01,
    radius: 1,
    tube: 0.4,
    tubularSegments: 64,
    radialSegments:  8,
};

gui.add(config, "speed", 0, 0.10, 0.001);

const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(
    75, window.innerWidth / window.innerHeight, 0.1, 1000
);
const renderer = new THREE.WebGLRenderer();
renderer.setSize(window.innerWidth, window.innerHeight);
document.body.appendChild(renderer.domElement);

const light = new THREE.DirectionalLight( 0xffffff );
				light.position.set( 0, 0, 1 );
				scene.add( light );

const geometry = new THREE.TorusKnotGeometry(config.radius, 
    config.tube, config.tubularSegments, config.radialSegments);
const material = new THREE.MeshToonMaterial(
    {color: 0x90C2E7}
);
const cube = new THREE.Mesh(geometry, material);
scene.add(cube);
scene.add(light);
camera.position.z = 5;

function animate() {
    requestAnimationFrame(animate);
    cube.rotation.x = fluids.add(cube.rotation.x, config.speed);
    cube.rotation.y = fluids.add(cube.rotation.y, config.speed);
    renderer.render(scene, camera);
}

animate();