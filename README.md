# `pmsimulator`
2D gravitational N body code using a particle-mesh method

# Introduction
`pmsimulator` is an object-oriented code for running gravitational N body simulations on a 2D mesh.  In each time step:
* The density field is calculated on the mesh
* The gravitional potential field is calculated from the density field by solving the Poisson equation in Fourier space
* The acceleration vector for each particle is calculated from the gravitional potential field on the mesh
* The particle positions and velocities are updated

The following features cannot be controlled by the user:
* The simulation grid is a square
* The simulation grid is periodic

The dimensions and resolution of the simuation grid, as well as the time step, can be set by the user.  Adaptive time stepping may be used as well.

## Units

This code uses Planck units, such that $G = c = 1$.  Thus, Newton's gravitational force law is written
$$\frac{F}{F_\mathrm{P}} = \frac{\left(\frac{m_1}{m_\mathrm{P}}\right)\left(\frac{m_2}{m_\mathrm{P}}\right)}{\left(\frac{r}{l_\mathrm{P}}\right)^2},$$
where the subscript P denotes the Planck unit.

We solve the Poisson equation:
$$\nabla^2 \phi(\vec{x}) = 4\pi\rho(\vec{x})$$
which in Fourier space becomes
$$\phi(\vec{k}) = -\frac{4\pi\rho(\vec{k})}{k^2} = -\frac{4\pi\rho(\vec{k})}{k_x^2 + k_y^2}$$
where $k_i = 2\pi/\lambda_i$ are the wavenumbers in the $i$ direction.
The gravitational acceleration vector at each grid cell is then given by
$$\vec{a} = -\vec{\nabla}\phi.$$

In the above, $\rho$ is in units of $m_\mathrm{P} l_\mathrm{P}^{-3} = c^5 \hbar^{-1} G^{-2}$, $a$ is in units of $l_\mathrm{P} t_\mathrm{P}^{-2} = \sqrt{c^7 \hbar^{-1} G^{-1}}$, and $\phi$ is in units of $l_\mathrm{P}^2 t_\mathrm{P}^{-2} = c^2$.

## Constants

The `pmsimulator.constants` class includes some numerical constants that may be of use working in Planck units.
* `pmsimulator.constants.msun = 9.138e37` ($1 M_\odot = 9.138\times10^{37} m_\mathrm{P}$)
* `pmsimulator.constants.kpc = 1.9092e54` (1 kpc $= 1.9092\times10^{54} l_\mathrm{P}$)
* `pmsimulator.constants.year = 5.849e50` (1 year $= 5.849\times10^{50} t_\mathrm{P}$)
* `pmsimulator.constants.hub = 1.747e-61` (100 $\frac{\mathrm{km}}{\mathrm{Mpc}\ \mathrm{s}}= 1.747\times10^{-61} t_\mathrm{P}^{-1}$)

# Structure

The primary class is the `Simulator` class found in `pmsimulator.simulator`.  An instantiation of this class will be an object that holds:
* the `simulation_settings` attribute, which is an instantiation of `pmsimulator.simulator.SimulationSettings`.  This stores the global properties of the simulation: `domain_size`, `grid_size`, and `dt`.  (Note that if adaptive time stepping is used, then `dt` is overridden).
* the `grid` attribute, , which is an instantiation of `pmsimulator.grid.Grid`.  This is automatically generated upon calling `pmsimulator.simulator.Simulator` based on the the `simulation_settings` attributes.  The simulator has only one `grid` in which all particle species live.
* any number of particle species, which are stored in the list `species_list`.  Each in an instantiation of the `pmsimulator.particles.Particles` class.  Upon construction, no particle species are part of the simuation.  The `create_particles` method can be called to add a particle species to the `species_list`.
* a set of attributes to keep track of the simulation:
  - `integration_time`: The total integration time that the simulator has evolved through
  - `n_species`: The total number of particle species in the simulator
  - `particles_per_species`: A Python list that keeps track of the number of particles in each species
