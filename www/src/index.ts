import { SimulationConfig as BhSimulationConfig, BarnesHutCradle } from "./barnesHut";
import { SimulationConfig as FSimulationConfig, FluidCradle } from "./fluidSimulation";
import { ParticleVisualization } from "./particleVisualization";

const bhConfig = new BhSimulationConfig();
bhConfig.numParticles = 10;
bhConfig.theta = 0.5;
bhConfig.maxDepth = 4;
bhConfig.maxParticlesPerNode = 4;
bhConfig.gravitationalConstant = 0.1;


const fConfig = new FSimulationConfig();
fConfig.gravity = 30.0;
fConfig.numParticles = 100;
fConfig.boundaryElasticity = 0.5;
fConfig.particleElasticity = 0.95;
fConfig.particleRadius = 1.0;
fConfig.initialVelocity = 0.0;

const barnesHut = new BarnesHutCradle(bhConfig);
const fluidSimulation = new FluidCradle(fConfig);

const particleVisualization = new ParticleVisualization(barnesHut);

particleVisualization.initializeDom()

window.addEventListener('resize', () => { particleVisualization.resize() });

const animate = () => {
  requestAnimationFrame(animate);
  particleVisualization.render();
};

animate();
