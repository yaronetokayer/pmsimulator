# from pmsimulator.particles.particles import *
from pmsimulator.grid.density_schemes import DENSITY_METHODS, assign_density_tsc, assign_density_cic, assign_density_ngp

import numpy as np
import matplotlib.pyplot as plt

class Grid:
    
    def __init__(self, simulation_settings):
        self.simulation_settings = simulation_settings
        
        # Number of grid cells along each axis
        self.grid_cells_x = int(simulation_settings.domain_size / simulation_settings.grid_size)
        self.grid_cells_y = int(simulation_settings.domain_size / simulation_settings.grid_size)

        # Initialize the fields
        self.density_field = np.zeros((self.grid_cells_x, self.grid_cells_y))
        self.potential_field = np.zeros((self.grid_cells_x, self.grid_cells_y))
        self.accel_field = (
            np.zeros((self.grid_cells_x, self.grid_cells_y)), 
            np.zeros((self.grid_cells_x, self.grid_cells_y))
            )

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

        # Ensure particle_populations is a list
        if not isinstance(particle_populations, list):
            particle_populations = [particle_populations]

        # Calculate mean density 
        total_particles = 0
        for particle_pop in particle_populations:
            total_particles += particle_pop.num_particles
        mean_density = total_particles / self.simulation_settings.domain_size**2

        # Reset the density field
        self.density_field = np.zeros((self.grid_cells_x, self.grid_cells_y))
        
        if method not in DENSITY_METHODS:
            raise ValueError(f"Invalid density assignment method. Allowed methods: {', '.join(DENSITY_METHODS)}")

        # Assign density to each grid cell
        for particles in particle_populations:    
            for x, y, mass in zip(particles.positions.x.data, 
                                  particles.positions.y.data, 
                                  particles.masses.data):
                if method == 'tsc':
                    assign_density_tsc(self, x, y, mass)
                elif method == 'cic':
                    assign_density_cic(self, x, y, mass)
                elif method == 'ngp':
                    assign_density_ngp(self, x, y, mass)

        # Scale by cell size and subtract off the mean
        self.density_field /= (self.simulation_settings.grid_size**2 / self.simulation_settings.domain_size**2)
        self.density_field -= mean_density

    def compute_potential(self):
        '''
        Compute the gravitational potential field from the density field using FFT
        We solve the 2D Poisson eq in Fourier space: -k^2 \\phi(k) = -2*\\pi G \\rho(k)

        The potential field is stored in self.potential_field
        The negative gradient of the potential field is stored in self.accel_field
        '''

        # Perform FFT on the density field
        density_fft = np.fft.fft2(self.density_field)

        # Define wave numbers along x and y axes
        kx = 2 * np.pi * np.fft.fftfreq(self.grid_cells_x, d=self.simulation_settings.grid_size)
        ky = 2 * np.pi * np.fft.fftfreq(self.grid_cells_y, d=self.simulation_settings.grid_size)

        # Create 2D arrays of wave numbers
        kxkx, kyky = np.meshgrid(kx, ky, indexing='xy')

        # Avoid division by zero by adding a small constant to the denominators
        epsilon = 1e-10
        denominator = kxkx**2 + kyky**2 + epsilon

        # Solve Poisson equation in Fourier space
        phi_fft = -4 * np.pi * density_fft / denominator

        # Set the (0,0) element to 0 to properly normalize the density
        phi_fft[0,0] = 0.0

        # # Compute the acceleration field by taking the gradient in Fourier space
        # grad_x_fft = 1j * kxkx * phi_fft
        # grad_y_fft = 1j * kyky * phi_fft

        # Perform inverse FFT to obtain the potential in real space
        phi_real = np.fft.ifft2(phi_fft).real
        self.potential_field = phi_real

        # Compute the gradient by taking the gradient in coordinate space
        # Create a ghost layer to account for periodic BCs in the gradient calculation
        potential_with_ghost = np.zeros((self.grid_cells_x + 2, self.grid_cells_y + 2))
        potential_with_ghost[1:-1, 1:-1] = self.potential_field
        potential_with_ghost[0, :] = potential_with_ghost[-2, :]
        potential_with_ghost[-1, :] = potential_with_ghost[1, :]
        potential_with_ghost[:, 0] = potential_with_ghost[:, -2]
        potential_with_ghost[:, -1] = potential_with_ghost[:, 1]
        dx = self.simulation_settings.grid_size
        grad_x = np.gradient(potential_with_ghost, axis=1) / dx
        grad_y = np.gradient(potential_with_ghost, axis=0) / dx

        self.accel_field = (-grad_x[1:-1, 1:-1], -grad_y[1:-1, 1:-1])

        # self.accel_field = (-np.fft.ifft2(grad_x_fft).real, -np.fft.ifft2(grad_y_fft).real)

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

            for particle_pop in particles:
                scatter_kwargs = {
                    'marker': kwargs.get('marker', 'o'),
                    's': kwargs.get('s', 2),
                    'alpha': kwargs.get('alpha', 0.5),
                    'label': particle_pop.species_name
                }
                plt.scatter(
                    particle_pop.positions.x.data, particle_pop.positions.y.data, **scatter_kwargs
                )
            plt.legend()
    
        colorbar = plt.colorbar(heatmap, label='Density - Mean Density')

        # Cosmetics
        plt.title('Grid Density Heatmap')
        plt.xlabel('X Position')
        plt.ylabel('Y Position')
        plt.xlim(0, self.simulation_settings.domain_size)
        plt.ylim(0, self.simulation_settings.domain_size)
    
        plt.show()

    def plot_potential_heatmap(self, force=True, particles=None, **kwargs):
        '''
        Method to plot the current potential field as a heatmap.
        Note that the only kwargs that will get noticed are:
        -cmap for the heatmap
        -c, marker, s, and alpha for the scatter
    
        Parameters:
        force - option to include force vectors. Deault is True.
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
        heatmap = plt.imshow(self.potential_field, **heatmap_kwargs)

        if force:
            x = np.linspace(0, self.simulation_settings.domain_size, self.grid_cells_x)
            y = np.linspace(0, self.simulation_settings.domain_size, self.grid_cells_y)
            X, Y = np.meshgrid(x, y)
            plt.quiver(
                X, Y, self.accel_field[0], self.accel_field[1],
                headlength=3, headaxislength=3, 
                scale=None, color='white', width=0.005
                )
    
        # Overlay particle positions if provided
        if particles is not None:
            # Ensure particles is a list
            if not isinstance(particles, list):
                particles = [particles]

            for particle_pop in particles:
                scatter_kwargs = {
                    'marker': kwargs.get('marker', 'o'),
                    's': kwargs.get('s', 2),
                    'alpha': kwargs.get('alpha', 0.5),
                    'label': particle_pop.species_name
                }
                plt.scatter(
                    particle_pop.positions.x.data, particle_pop.positions.y.data, **scatter_kwargs
                )
            plt.legend(loc='lower right')
    
        colorbar = plt.colorbar(heatmap, label='Potential')

        # Cosmetics
        plt.title('Potential Field Heatmap')
        plt.xlabel('X Position')
        plt.ylabel('Y Position')
        plt.xlim(0, self.simulation_settings.domain_size)
        plt.ylim(0, self.simulation_settings.domain_size)
    
        plt.show()