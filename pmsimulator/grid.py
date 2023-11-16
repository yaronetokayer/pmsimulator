import numpy as np
import matplotlib.pyplot as plt

class Grid:

    DENSITY_METHODS = ['tsc', 'cic', 'ngp']
    
    def __init__(self, domain_size, grid_size):
        
        # Number of grid cells along each axis
        self.grid_cells_x = int(domain_size / grid_size)
        self.grid_cells_y = int(domain_size / grid_size)

        # Initialize the density field to zero
        self.denisty_field = np.zeros((self.grid_cells_x, self.grid_cells_y))

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
        grid_x = int(x / SIMULATION_PARAMETERS['GRID_SIZE'])
        grid_y = int(y / SIMULATION_PARAMETERS['GRID_SIZE'])
        
        # Loop over the 3x3 grid around the particle's grid cell
        for i in range(grid_x - 1, grid_x + 2):
            for j in range(grid_y - 1, grid_y + 2):
                # Apply periodic boundary conditions to grid indices
                i_periodic = i % self.grid_cells_x
                j_periodic = j % self.grid_cells_y

                # Distance between particle and grid cell center
                dx = (i_periodic + 0.5) * SIMULATION_PARAMETERS['GRID_SIZE'] - x
                dy = (j_periodic + 0.5) * SIMULATION_PARAMETERS['GRID_SIZE'] - y
                    
                # TSC weights along x and y directions
                weight_x = max(1 - abs(2 * dx / SIMULATION_PARAMETERS['GRID_SIZE']), 0)
                weight_y = max(1 - abs(2 * dy / SIMULATION_PARAMETERS['GRID_SIZE']), 0)

                weight = weight_x * weight_y

                self.denisty_field[j_periodic, i_periodic] += mass * weight

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
        -c, maker, s, and alpha for the scatter
    
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
            'extent': [0, SIMULATION_PARAMETERS['DOMAIN_SIZE'], 0, SIMULATION_PARAMETERS['DOMAIN_SIZE']]
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
        plt.xlim(0, SIMULATION_PARAMETERS['DOMAIN_SIZE'])
        plt.ylim(0, SIMULATION_PARAMETERS['DOMAIN_SIZE'])
    
        plt.show()