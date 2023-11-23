DENSITY_METHODS = ['tsc', 'cic', 'ngp']

def assign_density_tsc(grid, x, y, mass):
    '''
    Triangular Shaped Cloud (TSC) scheme for applying density to grid cells

    Inputs:
    grid - grid object
    x - array of x particle positions
    y - array of y particle positions
    mass - mass of each particle

    Returns:
    None
    '''
    # Convert particle positions to grid indices
    grid_x = int(x / grid.simulation_settings.grid_size)
    grid_y = int(y / grid.simulation_settings.grid_size)
    
    # Loop over the 3x3 grid around the particle's grid cell
    for i in range(grid_x - 1, grid_x + 2):
        for j in range(grid_y - 1, grid_y + 2):
            # Apply periodic boundary conditions to grid indices
            i_periodic = i % grid.grid_cells_x
            j_periodic = j % grid.grid_cells_y

            # Distance between particle and grid cell center
            dx = (i_periodic + 0.5) * grid.simulation_settings.grid_size - x
            dy = (j_periodic + 0.5) * grid.simulation_settings.grid_size - y
                
            # TSC weights along x and y directions
            weight_x = max(1 - abs(2 * dx / grid.simulation_settings.grid_size), 0)
            weight_y = max(1 - abs(2 * dy / grid.simulation_settings.grid_size), 0)

            weight = weight_x * weight_y

            grid.density_field[j_periodic, i_periodic] += mass * weight

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