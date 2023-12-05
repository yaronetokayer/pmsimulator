import numpy as np

from pmsimulator.particles.distributions import *
from pmsimulator.particles.accel_schemes import *

class Particles:
    
    def __init__(self, num_particles, simulation_settings, species_name):
        
        if not isinstance(num_particles, int):
            raise TypeError("num_particles must be an integer.")
        
        self.species_name = species_name
        self.num_particles = num_particles
        self.positions = Vector(num_particles)
        self.velocities = Vector(num_particles)
        # self.x_positions = ParticleArray(num_particles)
        # self.y_positions = ParticleArray(num_particles)
        # self.x_velocities = ParticleArray(num_particles)
        # self.y_velocities = ParticleArray(num_particles)
        self.masses = ParticleArray(num_particles)

        self.accels = Vector(num_particles)

        self.domain_size = simulation_settings.domain_size

    def apply_periodic_bc(self):
        '''
        Apply periodic boundary conditions for position arrays
        '''
        for comp in self.positions.comps:
            comp.data = comp.data % self.domain_size
        # self.positions.x.data = self.positions.x.data % self.domain_size
        # self.positions.y.data = self.positions.y.data % self.domain_size

    def kinetic_energy(self):
        '''
        Calculate and return the kinetic energy of particle species in current state
        '''
        v_squared = self.velocities[0].data**2 + self.velocities[1].data**2

        return 0.5 * np.sum( self.masses.data * v_squared )

    def potential_energy(self, grid):
        '''
        Calculate and return the potential energy of particle species in current state
        PROBLEM: THIS DOESN'T WORK WITH THE WAY WE RESCALE THE POTENTIAL FIELD AT EACH TIME STEP
        '''
        pe = np.zeros(self.num_particles)

        for i in range(self.num_particles):
            x_particle = self.positions.x.data[i]
            y_particle = self.positions.y.data[i]
            mass_particle = self.masses.data[i]

            # Find the corresponding grid cell for the particle
            grid_x = int(x_particle / grid.simulation_settings.grid_size)
            grid_y = int(y_particle / grid.simulation_settings.grid_size)

            pe += mass_particle * grid.potential_field[grid_y, grid_x]

    def calculate_accels(self, grid, method='tsc'):
        '''
        Calculate the acceleration of each particle given a force field from the grid object

        Inputs:
        grid - A Grid object with a force field already calculated
        '''
        if method not in ACCEL_METHODS:
            raise ValueError(f"Invalid force assignment method. Allowed methods: {', '.join(ACCEL_METHODS)}")

        accels_x = np.empty(self.num_particles)
        accels_y = np.empty(self.num_particles)

        for i in range(self.num_particles):
            x_particle = self.positions.x.data[i]
            y_particle = self.positions.y.data[i]

            # Convert particle positions to grid indices
            grid_x = int(x_particle / grid.simulation_settings.grid_size)
            grid_y = int(y_particle / grid.simulation_settings.grid_size)

            # Calculate the force from the gravitational potential in the grid cell
            if method == 'ngp':
                accels_x[i], accels_y[i] = calculate_accel_ngp(grid_x, grid_y, grid.force_field)

        self.accels.x.data, self.accels.y.data = accels_x, accels_y

    def update(self, time_step):
        '''
        Update the positions and velocities of the species with current values through a time step 'time_step' in length, using leapfrog integration

        Inputs:
        time_step - time step size
        '''
        for pos, vel, accel in zip(self.positions.comps, self.velocities.comps, self.accels.comps):
            pos.data += vel.data * time_step
            vel.data += accel.data * time_step
        # self.x_positions.data += self.x_velocities.data * time_step
        # self.y_positions.data += self.y_velocities.data * time_step
        self.apply_periodic_bc()

        # self.x_velocities.data += self.accels[0] * time_step
        # self.y_velocities.data += self.accels[1] * time_step

    def polar_coords(self):
        '''
        Export the particle positions and velocities in polar coordinates
        '''
        pass

class Vector:
    '''
    Class for vector quantities: position, velocity, and acceleration
    '''

    def __init__(self, num_particles):
        if not isinstance(num_particles, int):
            raise TypeError("num_particles must be an integer.")

        self.x = ParticleArray(num_particles)
        self.y = ParticleArray(num_particles)
        self.comps = (self.x, self.y)

class ParticleArray:
        
    def __init__(self, num_particles):
        if not isinstance(num_particles, int):
            raise TypeError("num_particles must be an integer.")
        
        self.num_particles = num_particles
        self.data = np.zeros(num_particles, dtype='float64') # Initialize to zero

    def initialize(self, distribution='uniform', **kwargs):
        '''
        Method to initialize an array according to some specified distribution

        Inputs:
        distribution - the distribution according which to intialize the array.  Must be from ALLOWED_DISTRIBUTIONS
            'uniform': a uniform distribution, range specified by the tuple 'bounds'
            'normal': a Gaussian distribution, mean and std can be specified by the user
            'allones': set all entries to a single value, specified by the float 'scale'
        
        Returns:
        None
        '''
        if distribution not in ALLOWED_DISTRIBUTIONS:
            raise ValueError(f"Invalid distribution type. Allowed distributions: {', '.join(ALLOWED_DISTRIBUTIONS)}")

        # Initialize according to chosen distribution
        if distribution == 'uniform':
            self.data = initialize_uniform(self.num_particles, **kwargs)
        elif distribution == 'normal':
            self.data = initialize_normal(self.num_particles, **kwargs)
        elif distribution == 'allones':
            self.data = initialize_all_ones(self.num_particles, **kwargs)