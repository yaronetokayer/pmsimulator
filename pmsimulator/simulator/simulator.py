from pmsimulator.particles.particles import *
from pmsimulator.grid.grid import *

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
        plt.title('Snapshot t=' + str(self.integration_time))
        plt.xlabel('X Position')
        plt.ylabel('Y Position')
        plt.xlim(0, self.simulation_settings.domain_size)
        plt.ylim(0, self.simulation_settings.domain_size)
    
        plt.show()

    def advance(self, int_time, density_method='tsc', accel_method='ngp', adaptive_time=False, energy=False):
        '''
        Advance the entire simulation by 'int_time' time.
        Option to include energy calculations - KE, PE, and total E for each particle species at each time step
        '''
        time_elapsed = 0.0
        counter = 0

        while time_elapsed < int_time:
            if adaptive_time:
                # Find maximum velocity component
                if counter == 0: # At first iteration, require a step size of at most dt
                    max_v = self.simulation_settings.grid_size / self.simulation_settings.dt 
                else:
                    max_v = 1e-5 # Avoid division by zero
                for particle_pop in self.species_list:
                    max_v = max_v = np.max([
                        np.max(np.abs(particle_pop.x_velocities.data)), 
                        np.max(np.abs(particle_pop.y_velocities.data)), 
                        max_v
                        ])
                dt = self.simulation_settings.grid_size / max_v
            else:
                dt = self.simulation_settings.dt 

            self.grid.assign_density(self.species_list, method=density_method)
            self.grid.compute_potential()
            for particle_pop in self.species_list:
                particle_pop.calculate_accels(self.grid, method=accel_method)
                particle_pop.update(dt)
            time_elapsed += dt
            self.integration_time += dt
            counter += 1