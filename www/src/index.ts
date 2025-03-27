import { SimulationConfig, BarnesHutCradle } from "./barnesHut";
import { ParticleVisualization } from "./particleVisualization";

const config = new SimulationConfig();
config.numParticles = 10;
config.theta = 0.5;
config.maxDepth = 4;
config.maxParticlesPerNode = 4;
config.gravitationalConstant = 0.1;

const barnesHut = new BarnesHutCradle(config);
const particleVisualization = new ParticleVisualization(barnesHut);

particleVisualization.initializeDom()

window.addEventListener('resize', () => { particleVisualization.resize() });

const animate = () => {
  requestAnimationFrame(animate);
  particleVisualization.render();
};

animate();
