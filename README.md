# fluids2


1. cargo install  wasm-pack
2. cd www && npm install 




# Metrics to look at

Objectives are:
1. assessing if the fluid / gas is in equilibrium
2. assesing properties of the system



### Equilibrium

Disable external forces (gravity), and maybe reduce to a 2d slice for simplicity.
Slide the area into a 10x10 grid (e.g.) and 
- look at the distribution of the particle count in each bin. Expect this to be stationary,
- look at the distribution of velocities within each bin, should also be similar
  - another way of looking at this is plotting a chart of velocity as a fucntion of (e.g.) x. expect this to be stationary


### Energy conservation
- energy of the system $E_m = E_k + E_p$
  - $E_k = \frac{1}{2} m v^2$
  - $E_p = m g h$

### Fluid properties

- size of particle relative to size of container
- volume of all particles relative to volume of container
- time between consecutive collisions for any given particle
- Knudsen number $Kn = \frac{\lambda}{L}$ where L is a cacarteristic of the size of the space, lambda the mean free path

$$Kn = \frac{k_b T}{\sqrt{2 \pi} d^2 p L}$$

where
- $k_{\text{B}}$ is the Boltzmann constant (1.380649 × 10−23 J/K in SI units) [M1 L2 T−2 Θ−1]
- $T$ is the thermodynamic temperature 
- $d$ is the particle hard-shell diameter 
- $p$ is the static pressure
- $R_{s}$ is the specific gas constant (287.05 J/(kg K) for air),
- $\rho$ is the density