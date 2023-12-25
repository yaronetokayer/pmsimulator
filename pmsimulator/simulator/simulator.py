import numpy as np
import matplotlib.pyplot as plt

from pmsimulator.particles.particles import Particles
from pmsimulator.grid.grid import Grid

class SimulationSettings:
    def __init__(self, domain_size=1.0, grid_size=0.01, dt=0.01):
        self.domain_size = domain_size
        self.grid_size = grid_size
        self.dt = dt # Time resolution of the simulation

class PMSimulator:
    def __init__(self, **kwargs):
        '''
        Instatiate a simulator object
        **kwarg - optional keyword arguments for the simulation_settings attribute
        '''

        self.simulation_settings = SimulationSettings(**kwargs)

        # Keep track of simulation properties
        self.integration_time = 0.0 # Keeps track of how much time the simulator has evolved through
        self.n_species = 0 # The total number of particle species in the simulator
        self.particles_per_species = []
        self.species_list = []

        # Store energy calculations
        self.potential_energy_array = None
        self.kinetic_energy_array = None
        self.time_array = None
        
        self.grid = Grid(self.simulation_settings) # Create the simulation grid

    def reset_integration_time(self):
        self.integration_time = 0.0

    def create_particles(self, num_particles, species_name=None):
        '''
        Create a particle species with 'num_particles' particles
        The species is made a simulator attribute, which is appended to species_list

        Inputs:
        num_particles - the number of particles for the species
        species_name - (optional) give the species a name for easy access of the simulator attribute 
        '''
        if species_name is None:
            # Assign generic attribute name if none is passed by user
            species_name = f'particles_{len([attr for attr in self.__dict__ if "particles_" in attr])}'

        particles_instance = Particles(num_particles, self.simulation_settings, species_name)
            
        setattr(self, species_name, particles_instance)
        
        self.n_species += 1
        self.particles_per_species.append(num_particles)
        self.species_list.append(getattr(self, species_name))

    def plot_snapshot(self, save_path=None, **kwargs):
        '''
        Method to plot current particle positions
        '''
        fig, ax = plt.subplots(1, 1, figsize=(8,8))

        for particle_pop in self.species_list:
            scatter_kwargs = {
                'marker': kwargs.get('marker', 'o'),
                's': kwargs.get('s', 2),
                'alpha': kwargs.get('alpha', 0.5),
                'label': particle_pop.species_name
            }
            ax.scatter(
                particle_pop.positions.x.data, particle_pop.positions.y.data, **scatter_kwargs
            )
        ax.legend(loc='lower right')

        # Cosmetics
        ax.set_title(f't={self.integration_time:.4f}')
        ax.set_xlabel('X Position')
        ax.set_ylabel('Y Position')
        ax.set_xlim(0, self.simulation_settings.domain_size)
        ax.set_ylim(0, self.simulation_settings.domain_size)

        if save_path is not None:
            fig.savefig(save_path, dpi=300)
            plt.close('all')

    def advance_one(self, dt=None, density_method='tsc', accel_method='ngp'):
        '''
        Advance the entire simulation by a single time step

        Parameters:
        dt (float) - Time step length.
        density_method (str) - Density assignment method.
        accel_method (str) - Acceleration calculation method.
        '''
        if dt is None:
            dt = self.simulation_settings.dt
            
        self.grid.assign_density(self.species_list, method=density_method)
        self.grid.compute_potential()
        for particle_pop in self.species_list:
            particle_pop.calculate_accels(self.grid, method=accel_method)
            particle_pop.update(dt)

        self.integration_time += dt

    def advance(
        self, int_time, 
        density_method='tsc', accel_method='ngp', adaptive_time=False, energy=False,
        ):
        '''
        Advance the entire simulation by 'int_time' time.

        Parameters:
        int_time (float) - Time to advance the simulation.
        density_method (str) - Density assignment method.
        accel_method (str) - Acceleration calculation method.
        adaptive_time (bool) - If True, uses adaptive time steps.
        energy (bool) - If True, includes energy calculations and returns lists of values.
        '''
        time_elapsed = 0.0
        first = True

        if energy:
            time_array = [0.0]
            potential_energy_array = []
            kinetic_energy_array = []
            potential_energy_array.append(0.0)
            kinetic_energy_array.append(0.0)
            for particle_pop in self.species_list:
                particle_pop.calculate_potential_energy(self.grid)
                particle_pop.calculate_kinetic_energy()
                potential_energy_array[-1] += particle_pop.potential_energy
                kinetic_energy_array[-1] += particle_pop.kinetic_energy

        while time_elapsed < int_time:
            if adaptive_time:
                dt = self.compute_adaptive_dt(first=first)
                first = False
            else:
                dt = self.simulation_settings.dt 

            self.advance_one(dt, density_method=density_method, accel_method=accel_method)

            time_elapsed += dt

            if energy:
                potential_energy_array.append(0.0)
                kinetic_energy_array.append(0.0)
                time_array.append(self.integration_time)
                for particle_pop in self.species_list:
                    particle_pop.calculate_potential_energy(self.grid)
                    particle_pop.calculate_kinetic_energy()
                    potential_energy_array[-1] += particle_pop.potential_energy
                    kinetic_energy_array[-1] += particle_pop.kinetic_energy

        if energy:
            self.potential_energy_array = np.array(potential_energy_array)
            self.kinetic_energy_array = np.array(kinetic_energy_array)
            self.time_array = np.array(time_array)

    def compute_adaptive_dt(self, first=False):
        '''
        Compute adaptive time step for the simulation given the current configuration

        Inputs:
        first - Flag indicating if this is the first time step
        '''
        if first: # At first iteration, require a step size of at most dt
            max_v = self.simulation_settings.grid_size / self.simulation_settings.dt 
        else:
            max_v = 1e-5 # Avoid division by zero

        for particle_pop in self.species_list:
            max_v = max_v = np.max([
                np.max(np.abs(particle_pop.velocities.x.data)), 
                np.max(np.abs(particle_pop.velocities.y.data)), 
                max_v
                ])

        return self.simulation_settings.grid_size / max_v

    def export_snapshot(self):
        '''
        Function to export a current snapshot of the simulation as a Python dictionary.  Includes all particle positions and velocities.
        '''
        pass