import numpy as np

class Particles:
    
    def __init__(self, num_particles, species_name=None):
        
        if not isinstance(num_particles, int):
            raise TypeError("num_particles must be an integer.")
            
        self.num_particles = num_particles
        self.x_positions = ParticleArray(num_particles)
        self.y_positions = ParticleArray(num_particles)
        self.x_velocities = ParticleArray(num_particles)
        self.y_velocities = ParticleArray(num_particles)
        self.masses = ParticleArray(num_particles)

    def apply_periodic_bc(self):
        '''
        Apply periodic boundary conditions for positions arrays
        '''
        self.x_positions.data = self.x_positions.data % SIMULATION_PARAMETERS['DOMAIN_SIZE']
        self.y_positions.data = self.y_positions.data % SIMULATION_PARAMETERS['DOMAIN_SIZE']

    def kinetic_energy(self):
        '''
        Calculate and return the kinetic energy of particle species in current state
        '''
        v_squared = self.x_velocities.data**2 + self.y_velocities.data**2

        return 0.5 * np.sum( self.masses.data * v_squared )

class ParticleArray:

    ALLOWED_DISTRIBUTIONS = ['uniform', 'normal', 'allones']
        
    def __init__(self, num_particles):
        if not isinstance(num_particles, int):
            raise TypeError("num_particles must be an integer.")
        
        self.num_particles = num_particles
        self.data = np.zeros(num_particles, dtype='float64')

    def initialize_uniform(self, bounds=(0, 1)):
        self.data = np.random.uniform(low=bounds[0], high=bounds[1], size=self.num_particles)

    def initialize_normal(self, mean=0, std=1):
        self.data = np.random.normal(mean, std, self.num_particles)

    def initialize_all_ones(self, scale=1.0):
        self.data = scale * np.ones(self.num_particles)

    def initialize(self, distribution='uniform', **kwargs):
        '''
        Method to initialize an array according to some specifies distribution

        Inputs:
        distribution - the distribution according which to intialize the array.  Must be from ALLOWED_DISTRIBUTIONS
            'uniform': a uniform distribution, range specified by the tuple 'bounds'
            'normal': a Gaussian distribution, mean and std can be specified by the user
            'allones': set all entries to a single value, specified by the float 'scale'
        
        Returns:
        None
        '''
        if distribution not in self.ALLOWED_DISTRIBUTIONS:
            raise ValueError(f"Invalid distribution type. Allowed distributions: {', '.join(self.ALLOWED_DISTRIBUTIONS)}")

        # Initialize according to chosen distribution
        if distribution == 'uniform':
            self.initialize_uniform(**kwargs)
        elif distribution == 'normal':
            self.initialize_normal(**kwargs)
        elif distribution == 'allones':
            self.initialize_all_ones(**kwargs)