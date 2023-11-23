ACCEL_METHODS = ['tsc', 'cic', 'ngp']

def calculate_accel_ngp(x, y, force_field):
    '''
    Nearest Grid Point (NGP) scheme for calculating acceleration of each particle

    Inputs:
    x - the x grid index of the particle
    x - the y grid index of the particle
    force_field - force_field attribute of Grid object

    Returns:
    accel_x, accel_y
    '''
    
    return force_field[0][y, x], force_field[1][y, x]