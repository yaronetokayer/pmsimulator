from .particles import *
from .grid import *

SIMULATION_PARAMETERS = {
    'DOMAIN_SIZE': 1.0,
    'GRID_SIZE': 0.1,
    'TIME_STEPS': 0 # How many time steps has the simulator evolved through
}

PARTICLE_SPECIES = {
    'N': 0, # The total number of particle species in the simulator
    'PARTICLES_PER_SPECIES': [], # How many particles are in each species
    'SPECIES_LIST': [], # Pointers to the Particle objects
}

class PMSimulator:
    def __init__(self):
        
        self.grid = Grid(SIMULATION_PARAMETERS['DOMAIN_SIZE'], SIMULATION_PARAMETERS['GRID_SIZE']) # Create the simulation grid

    def create_particles(self, num_particles, species_name=None):
        particles_instance = Particles(num_particles)
        if species_name is None:
            # Assign generic attribute name if none is passed by user
            species_name = f'particles_{len([attr for attr in self.__dict__ if "particles_" in attr])}'
            
        setattr(self, species_name, particles_instance)
        
        PARTICLE_SPECIES['N'] += 1
        PARTICLE_SPECIES['PARTICLES_PER_SPECIES'].append(num_particles)
        PARTICLE_SPECIES['SPECIES_LIST'].append(getattr(self, species_name))

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
