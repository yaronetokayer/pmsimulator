from pmsimulator.particles import *
from pmsimulator.grid import *

class SimulationSettings:
    def __init__(self, domain_size=1.0, grid_size=0.01, dt=0.01H):
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

    def create_particles(self, num_particles, species_name=None):
        particles_instance = Particles(num_particles, self.simulation_settings)
        if species_name is None:
            # Assign generic attribute name if none is passed by user
            species_name = f'particles_{len([attr for attr in self.__dict__ if "particles_" in attr])}'
            
        setattr(self, species_name, particles_instance)
        
        self.n_species += 1
        self.particles_per_species.append(num_particles)
        self.species_list.append(getattr(self, species_name))
