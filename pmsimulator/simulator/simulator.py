import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from pmsimulator.particles.particles import *
from pmsimulator.grid.grid import *
from pmsimulator.simulator.animate import *

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

    def plot_snapshot(self, **kwargs):
        '''
        Method to plot current particle positions
        '''
        plt.figure(figsize=(8, 8))

        for particle_pop in self.species_list:
            scatter_kwargs = {
                'marker': kwargs.get('marker', 'o'),
                's': kwargs.get('s', 2),
                'alpha': kwargs.get('alpha', 0.5),
                'label': particle_pop.species_name
            }
            plt.scatter(
                particle_pop.x_positions.data, particle_pop.y_positions.data, **scatter_kwargs
            )
        plt.legend(loc='lower right')

        # Cosmetics
        plt.title(f't={self.integration_time:.4f}')
        plt.xlabel('X Position')
        plt.ylabel('Y Position')
        plt.xlim(0, self.simulation_settings.domain_size)
        plt.ylim(0, self.simulation_settings.domain_size)
    
        plt.show()

    def advance_one(self, dt, density_method='tsc', accel_method='ngp'):
       '''
        Advance the entire simulation by a single time step

        Parameters:
        dt (float) - Time step length.
        density_method (str) - Density assignment method.
        accel_method (str) - Acceleration calculation method.
        '''
        self.grid.assign_density(self.species_list, method=density_method)
        self.grid.compute_potential()
        for particle_pop in self.species_list:
            particle_pop.calculate_accels(self.grid, method=accel_method)
            particle_pop.update(dt)

        self.integration_time += dt

    def advance(self, int_time, density_method='tsc', accel_method='ngp', adaptive_time=False, energy=False):
        '''
        Advance the entire simulation by 'int_time' time.

        Parameters:
        int_time (float) - Time to advance the simulation.
        density_method (str) - Density assignment method.
        accel_method (str) - Acceleration calculation method.
        adaptive_time (bool) - If True, uses adaptive time steps.
        energy (bool) - If True, includes energy calculations and returns lists of values.
        animate (bool) - If True, creates an animation. Default is False.
        save_path (str) - If provided, saves the animation to a file.
        fps (int) - Frames per second for the animation. Default is 30.
        '''
        time_elapsed = 0.0
        counter = 0

        while time_elapsed < int_time:
            if adaptive_time:
                dt = self.compute_adaptive_dt(first=not counter)
            else:
                dt = self.simulation_settings.dt 

            self.advance_one(self, dt, density_method=density_method, accel_method=accel_method)

            time_elapsed += dt
            counter += 1

            if energy:
                '''TO BE FILLED IN'''
                pass

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
                np.max(np.abs(particle_pop.x_velocities.data)), 
                np.max(np.abs(particle_pop.y_velocities.data)), 
                max_v
                ])

        return self.simulation_settings.grid_size / max_v

