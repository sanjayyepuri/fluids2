import * as THREE from "three";
import * as fluids from "fluids2";
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

const positions = [];
const colors = [];

const particles = 5000;

const color = new THREE.Color();
const n = 500, n2 = n / 2; // particles spread in the cube

for ( let i = 0; i < particles; i ++ ) {

    // positions

    const x = Math.random() * n - n2;
    const y = Math.random() * n - n2;
    const z = Math.random() * n - n2;

    positions.push( x, y, z );

    const vx = ( x / n ) + 0.5;
    const vy = ( y / n ) + 0.5;
    const vz = ( z / n ) + 0.5;

    color.setRGB( vx, vy, vz );

    colors.push( color.r, color.g, color.b );
}


geometry.setAttribute("position", new THREE.Float32BufferAttribute(positions, 3));
geometry.setAttribute("color", new THREE.Float32BufferAttribute(colors, 3));

const material = new THREE.PointsMaterial( { size: 15, vertexColors: true } );

const points = new THREE.Points(geometry, material);

scene.add(points);

camera.position.z = 2000;

function animate() {
    requestAnimationFrame(animate);
    points.rotation.x = fluids.add(points.rotation.x, config.speed);
    points.rotation.y = fluids.add(points.rotation.y, config.speed);
    renderer.render(scene, camera);
}

animate();