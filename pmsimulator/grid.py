from pmsimulator.particles import *

import numpy as np
import matplotlib.pyplot as plt

class Grid:

    DENSITY_METHODS = ['tsc', 'cic', 'ngp']
    
    def __init__(self, simulation_settings):
        self.simulation_settings = simulation_settings
        
        # Number of grid cells along each axis
        self.grid_cells_x = int(simulation_settings.domain_size / simulation_settings.grid_size)
        self.grid_cells_y = int(simulation_settings.domain_size / simulation_settings.grid_size)

        # Initialize the density field
        self.density_field = np.zeros((self.grid_cells_x, self.grid_cells_y))

    def assign_density(self, particle_populations, method='tsc'):
        '''
        Assign a grid density based on the particles on the grid

        Inputs:
        particle_populations - an instantiation of Particles, or a list of instantiations, 
                               representing the particles on the grid
        method - the density assignment method.  Default is TSC

        Returns:
        None
        '''
        # Reset the density field
        self.density_field = np.zeros((self.grid_cells_x, self.grid_cells_y))
        
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

    def assign_density_cic(self, x, y, mass):
        '''
        Cloud in Cell (CIC) scheme for applying density to grid cells
        '''
        
        # CIC density assignment logic
        # ... (implementation needed)
        
    def assign_density_ngp(self, x, y, mass):
        '''
        Nearest Grid Point (NGP) scheme for applying density to grid cells
        '''
        # NGP density assignment logic
        # ... (implementation needed)

    def compute_potential(self):
        '''
        Compute the gravitational potential field from the density field using FFT
        We solve the 2D Poisson eq in Fourier space: -k^2 \\phi(k) = -2*\\pi G \\rho(k)
        '''

        # Perform FFT on the density field
        density_fft = np.fft.fft2(self.density_field)

        # Define wave numbers along x and y axes
        kx = 2 * np.pi * np.fft.fftfreq(self.grid_cells_x, d=self.simulation_settings.grid_size)
        ky = 2 * np.pi * np.fft.fftfreq(self.grid_cells_y, d=self.simulation_settings.grid_size)

        # Create 2D arrays of wave numbers
        kx, ky = np.meshgrid(kx, ky, indexing='ij')

        # Avoid division by zero by adding a small constant to the denominators
        epsilon = 1e-10
        denominator = kx**2 + ky**2 + epsilon

        # Solve Poisson equation in Fourier space
        phi_fft = 2 * np.pi * density_fft / denominator

        # Perform inverse FFT to obtain the potential in real space
        self.potential_field = np.fft.ifft2(phi_fft).real  * (2 * np.pi)**2

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
        heatmap = plt.imshow(self.density_field, **heatmap_kwargs)
    
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