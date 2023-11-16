import numpy as np
import matplotlib.pyplot as plt

class SimulationSettings:
    def __init__(self, domain_size=1.0, grid_size=0.01, dt=dt):
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
        self.time_steps = 0 # Counter for how many time steps the simulator has evolved through
        self.n_species = 0 # The total number of particle species in the simulator
        self.particles_per_species = []
        self.species_list = []
        
        self.grid = Grid(self.simulation_settings) # Create the simulation grid

    def create_particles(self, num_particles, self.simulation_settings, species_name=None):
        particles_instance = Particles(num_particles, self.simulation_settings)
        if species_name is None:
            # Assign generic attribute name if none is passed by user
            species_name = f'particles_{len([attr for attr in self.__dict__ if "particles_" in attr])}'
            
        setattr(self, species_name, particles_instance)
        
        self.n_species += 1
        self.particles_per_species.append(num_particles)
        self.species_list.append(getattr(self, species_name))

    def reset_all(self):
        '''
        Delete all particle attributes and reset the simulation grid according to current values of SIMULATION_PARAMETERS
        Update PARTICLE_SPECIES accordingly
        '''
        # Delete particle attributes
        for attr in list(self.__dict__):
            if attr != 'grid':
                delattr(self, attr)
    
        # Reset the simulation grid
        self.grid = Grid(SIMULATION_PARAMETERS['DOMAIN_SIZE'], SIMULATION_PARAMETERS['GRID_SIZE'])

        # Reset particle species information
        PARTICLE_SPECIES['N'] = 0
        PARTICLE_SPECIES['PARTICLES_PER_SPECIES'] = []
        PARTICLE_SPECIES['SPECIES_LIST'] = []

class Particles:
    
    def __init__(self, num_particles, simulation_settings, species_name=None):
        
        if not isinstance(num_particles, int):
            raise TypeError("num_particles must be an integer.")
            
        self.num_particles = num_particles
        self.x_positions = ParticleArray(num_particles)
        self.y_positions = ParticleArray(num_particles)
        self.x_velocities = ParticleArray(num_particles)
        self.y_velocities = ParticleArray(num_particles)
        self.masses = ParticleArray(num_particles)

        self.domain_size = simulation_settings.domain_size

    def apply_periodic_bc(self):
        '''
        Apply periodic boundary conditions for positions arrays
        '''
        self.x_positions.data = self.x_positions.data % self.domain_size
        self.y_positions.data = self.y_positions.data % self.domain_size

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

class Grid:

    DENSITY_METHODS = ['tsc', 'cic', 'ngp']
    
    def __init__(self, simulation_settings):
        self.simulation_settings = simulation_settings
        
        # Number of grid cells along each axis
        self.grid_cells_x = int(simulation_settings.domain_size / simulation_settings.grid_size)
        self.grid_cells_y = int(simulation_settings.domain_size / simulation_settings.grid_size)

        # Initialize the density field to zero
        self.density_field = np.zeros((self.grid_cells_x, self.grid_cells_y))

    def assign_density(self, particle_populations, method='tsc'):
        '''
        Assign a grid density based on the particles on the grid

        Inputs:
        particle_populations - an instantiation of Particles, or a list of instantiations, representing the particles on the grid
        method - the density assignment method.  Default is TSC

        Returns:
        None
        '''
        
        if method not in self.DENSITY_METHODS:
            raise ValueError(f"Invalid density assignment method. Allowed methods: {', '.join(self.DENSITY_METHODS)}")

        # Ensure particle_populations is a list
        if not isinstance(particle_populations, list):
            particle_populations = [particle_populations]

        # Assign density to each grid cell
        for particles in particle_populations:    
            for x, y, mass in zip(particles.x_positions.data, 
                                  particles.y_positions.data, 
                                  particles.masses.data):
                if method == 'tsc':
                    self.assign_density_tsc(x, y, mass)
                elif method == 'cic':
                    self.assign_density_cic(x, y, mass)
                elif method == 'ngp':
                    self.assign_density_ngp(x, y, mass)
    
    def assign_density_tsc(self, x, y, mass):
        '''
        Triangular Shaped Cloud (TSC) scheme for applying density to grid cells

        Inputs:
        x - x particle positions
        y - y particle positions
        mass - mass of each particle

        Returns:
        None
        '''
        # Convert particle positions to grid indices
        grid_x = int(x / self.simulation_settings.grid_size)
        grid_y = int(y / self.simulation_settings.grid_size)
        
        # Loop over the 3x3 grid around the particle's grid cell
        for i in range(grid_x - 1, grid_x + 2):
            for j in range(grid_y - 1, grid_y + 2):
                # Apply periodic boundary conditions to grid indices
                i_periodic = i % self.grid_cells_x
                j_periodic = j % self.grid_cells_y

                # Distance between particle and grid cell center
                dx = (i_periodic + 0.5) * self.simulation_settings.grid_size - x
                dy = (j_periodic + 0.5) * self.simulation_settings.grid_size - y
                    
                # TSC weights along x and y directions
                weight_x = max(1 - abs(2 * dx / self.simulation_settings.grid_size), 0)
                weight_y = max(1 - abs(2 * dy / self.simulation_settings.grid_size), 0)

                weight = weight_x * weight_y

                self.density_field[j_periodic, i_periodic] += mass * weight

    def assign_density_cic(self, grid_x, grid_y, x, y, mass):
        '''
        Cloud in Cell (CIC) scheme for applying density to grid cells
        '''
        
        # CIC density assignment logic
        # ... (implementation needed)
        
    def assign_density_ngp(self, grid_x, grid_y, x, y, mass):
        '''
        Nearest Grid Point (NGP) scheme for applying density to grid cells
        '''
        # NGP density assignment logic
        # ... (implementation needed)

    def plot_density_heatmap(self, particles=None, **kwargs):
        '''
        Method to plot the current grid density as a heatmap.
        Note that the only kwargs that will get noticed are:
        -cmap for the heatmap
        -c, marker, s, and alpha for the scatter
    
        Parameters:
        particles - particles to overlay on the heatmap. Can be a list. Default is None.
    
        Returns:
        None
        '''
        plt.figure(figsize=(8, 8))
    
        # Plot the grid density as a heatmap
        heatmap_kwargs = {
            'cmap': kwargs.get('cmap', 'viridis'),
            'origin': 'lower',
            'extent': [0, self.simulation_settings.domain_size, 0, self.simulation_settings.domain_size]
        }
        heatmap = plt.imshow(self.denisty_field, **heatmap_kwargs)
    
        # Overlay particle positions if provided
        if particles is not None:
            # Ensure particles is a list
            if not isinstance(particles, list):
                particles = [particles]

            for i, particle_pop in enumerate(particles):
                scatter_kwargs = {
                    'marker': kwargs.get('marker', 'o'),
                    's': kwargs.get('s', 2),
                    'alpha': kwargs.get('alpha', 0.5),
                    'label': f'Particles {i + 1}'
                }
                plt.scatter(
                    particle_pop.x_positions.data, particle_pop.y_positions.data, **scatter_kwargs
                )
            plt.legend()
    
        colorbar = plt.colorbar(heatmap, label='Density')

        # Cosmetics
        plt.title('Grid Density Heatmap')
        plt.xlabel('X Position')
        plt.ylabel('Y Position')
        plt.xlim(0, self.simulation_settings.domain_size)
        plt.ylim(0, self.simulation_settings.domain_size)
    
        plt.show()
